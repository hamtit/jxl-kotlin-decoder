package io.github.hamtit.jxl

/** デコード中に壊れた・対応していないデータを見つけたときの例外 */
public class JxlException(message: String) : Exception(message)

internal fun fail(message: String): Nothing = throw JxlException(message)

/**
 * LSB 先頭のビットリーダー。[start] から [end] までのバイト列を読む。
 * 終端を越えた読み出しは 0 を返し、[checkOverrun] で検出する。
 */
internal class BitReader(private val data: ByteArray, private val start: Int, private val end: Int) {
    private var pos = start
    private var buf = 0L
    private var bits = 0

    constructor(data: ByteArray) : this(data, 0, data.size)

    /** これまでに消費したビット数 */
    val bitsRead: Long get() = (pos - start).toLong() * 8 - bits

    val totalBits: Long get() = (end - start).toLong() * 8

    private val view: java.nio.ByteBuffer = java.nio.ByteBuffer.wrap(data).order(java.nio.ByteOrder.LITTLE_ENDIAN)

    private fun refill() {
        val n = (63 - bits) ushr 3
        if (pos + 8 <= end) {
            // 8 バイトまとめて読む。buf の有効なビットより上には次のバイトのビットが入るが、
            // 次の補充で同じ位置に同じビットを OR するので壊れない
            buf = buf or (view.getLong(pos) shl bits)
            pos += n
            bits += n shl 3
            return
        }
        if (pos + n <= end) {
            val d = data
            var p = pos
            var b = bits
            var v = buf
            for (i in 0 until n) {
                v = v or ((d[p].toLong() and 0xff) shl b)
                p++
                b += 8
            }
            buf = v
            pos = p
            bits = b
            return
        }
        while (bits <= 56) {
            val b = if (pos < end) data[pos].toLong() and 0xff else 0L
            pos++
            buf = buf or (b shl bits)
            bits += 8
        }
    }

    /** n (0..32) ビットを先読みする */
    fun peek(n: Int): Int {
        if (bits < n) refill()
        return (buf and ((1L shl n) - 1)).toInt()
    }

    fun consume(n: Int) {
        if (bits < n) refill()
        buf = buf ushr n
        bits -= n
    }

    fun readBits(n: Int): Int {
        if (n == 0) return 0
        if (bits < n) refill()
        val v = (buf and ((1L shl n) - 1)).toInt()
        buf = buf ushr n
        bits -= n
        return v
    }

    fun readBool(): Boolean = readBits(1) == 1

    /** U32(d0, d1, d2, d3)。各分布は (offset, bits) で表す */
    fun readU32(o0: Int, b0: Int, o1: Int, b1: Int, o2: Int, b2: Int, o3: Int, b3: Int): Int =
        when (readBits(2)) {
            0 -> o0 + readBits(b0)
            1 -> o1 + readBits(b1)
            2 -> o2 + readBits(b2)
            else -> o3 + readBits(b3)
        }

    fun readU64(): Long = when (readBits(2)) {
        0 -> 0L
        1 -> 1L + readBits(4)
        2 -> 17L + readBits(8)
        else -> {
            var value = readBits(12).toLong()
            var shift = 12
            while (readBits(1) == 1) {
                if (shift == 60) {
                    value = value or (readBits(4).toLong() shl shift)
                    break
                }
                value = value or (readBits(8).toLong() shl shift)
                shift += 8
            }
            value
        }
    }

    fun readF16(): Float {
        val v = readBits(16)
        val sign = (v and 0x8000) shl 16
        if (v and 0x7fff == 0) return Float.fromBits(sign)
        val mantissa = v and 0x3ff
        val exponent = (v shr 10) and 0x1f
        if (exponent == 0x1f) fail("invalid f16 value")
        if (exponent == 0) {
            val value = (1.0f / 16384.0f) * (mantissa / 1024.0f)
            return if (sign != 0) -value else value
        }
        return Float.fromBits((mantissa shl 13) or ((exponent + 112) shl 23) or sign)
    }

    /** Enum: U32(0, 1, 2 + u(4), 18 + u(6)) */
    fun readEnum(): Int = readU32(0, 0, 1, 0, 2, 4, 18, 6)

    fun zeroPadToByte() {
        val rem = (bitsRead % 8).toInt()
        if (rem != 0) {
            if (readBits(8 - rem) != 0) fail("non-zero padding")
        }
    }

    /** n ビット読み飛ばす。残りより長ければ (壊れたデータ) 読まずに失敗する */
    fun skipBits(n: Long) {
        if (n < 0 || n > totalBits - bitsRead) fail("unexpected end of bitstream")
        if (n <= bits) {
            consume(n.toInt())
            return
        }
        // バッファにあるビットを捨て、残りはバイト単位で位置を進める
        val left = n - bits
        buf = 0L
        bits = 0
        pos += (left / 8).toInt()
        readBits((left % 8).toInt())
    }

    fun checkOverrun() {
        if (bitsRead > totalBits) fail("unexpected end of bitstream")
    }

    /** 現在位置 (バイト境界であること) の data 内での絶対オフセット */
    fun bytePosition(): Int {
        val b = bitsRead
        if (b % 8 != 0L) fail("not byte aligned")
        return start + (b / 8).toInt()
    }

    fun extensions() {
        val ext = readU64()
        var sum = 0L
        var bitsLeft = ext
        var i = 0
        while (bitsLeft != 0L && i < 64) {
            if (bitsLeft and 1L != 0L) {
                // 長さは入力の残りを越えられない (越えるなら壊れている)。合計の桁あふれもここで防ぐ
                val len = readU64()
                if (len < 0 || len > totalBits) fail("invalid extension length")
                sum += len
            }
            bitsLeft = bitsLeft ushr 1
            i++
        }
        skipBits(sum)
    }
}

internal fun unpackSigned(v: Int): Int = (v ushr 1) xor -(v and 1)
