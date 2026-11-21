package io.github.hamtit.jxl

import java.util.concurrent.ExecutorService
import java.util.concurrent.Executors
import java.util.concurrent.ThreadFactory
import java.util.concurrent.atomic.AtomicInteger

/** デコードの設定 */
public class JxlDecodeOptions(
    /** 並列処理に使う ExecutorService。null なら [threads] 本の共有スレッドプールを使う */
    public val executor: ExecutorService? = null,
    /** 1 なら呼び出したスレッドだけで処理する */
    public val threads: Int = Runtime.getRuntime().availableProcessors(),
    /** EXIF 相当の向き (orientation) を出力に反映するか */
    public val applyOrientation: Boolean = true,
    /** これより画素数の多い画像は復号しない (壊れたデータで巨大な領域を確保しないため) */
    public val maxPixels: Long = 1L shl 28,
)

/** 画像の基本情報 */
public class JxlInfo internal constructor(
    /** 向きを反映した後の幅 */
    public val width: Int,
    /** 向きを反映した後の高さ */
    public val height: Int,
    public val hasAlpha: Boolean,
    public val isAnimated: Boolean,
    /** アニメーションのループ回数 (0 は無限) */
    public val loopCount: Int,
    public val bitsPerSample: Int,
    public val orientation: Int,
)

/** 表示する 1 フレーム。pixels は ARGB (アルファは非乗算) */
public class JxlFrame internal constructor(
    public val width: Int,
    public val height: Int,
    public val pixels: IntArray,
    /** 表示時間 (ミリ秒)。静止画では 0 */
    public val durationMillis: Int,
)

/**
 * JPEG XL のデコーダー。フレームを先頭から順に [nextFrame] で取り出す。
 * スレッドセーフではない。
 */
public class JxlDecoder(data: ByteArray, private val options: JxlDecodeOptions = JxlDecodeOptions()) {
    private val codestream: ByteArray = Container.extractCodestream(data)
    private lateinit var header: ImageHeader
    private var firstFramePos: Int = 0
    private var pos: Int = 0
    private var refs = arrayOfNulls<ReferenceFrame>(4)
    private var lfFrames = arrayOfNulls<FrameImage>(4)
    private var finished = false
    private var visibleFrames = 0
    private var invisibleSinceVisible = 0
    /** 次に読むフレームの番号 (表示しないフレームも数える) */
    private var frameIndex = 0
    private val executor: ExecutorService? = options.executor ?: if (options.threads <= 1) null else sharedPool(options.threads)

    public lateinit var info: JxlInfo
        private set

    init {
        try {
            init0()
        } catch (e: JxlException) {
            throw e
        } catch (e: RuntimeException) {
            throw JxlException("invalid JPEG XL data (${e.javaClass.simpleName})")
        }
    }

    private fun init0() {
        val br = BitReader(codestream)
        if (br.readBits(16) != 0x0aff) fail("not a JPEG XL codestream")
        header = ImageHeader.parse(br)
        if (header.width.toLong() * header.height > options.maxPixels) fail("image too large")
        if (header.colourEncoding.wantIcc) Icc.skip(br)
        br.zeroPadToByte()
        if (header.previewWidth > 0) {
            // プレビューのフレームは読み飛ばす
            val ph = header.withSize(header.previewWidth, header.previewHeight)
            val fh = FrameHeader.parse(br, ph)
            val toc = Toc.parse(br, fh)
            val end = toc.end
            if (end > codestream.size) fail("truncated preview")
            firstFramePos = end.toInt()
        } else {
            firstFramePos = br.bytePosition()
        }
        pos = firstFramePos
        val swap = header.orientation in 5..8 && options.applyOrientation
        info = JxlInfo(
            if (swap) header.height else header.width,
            if (swap) header.width else header.height,
            header.alphaChannel >= 0,
            header.animation != null,
            header.animation?.numLoops ?: 0,
            header.bitDepth.bitsPerSample,
            header.orientation,
        )
    }

    /** 最初のフレームから読み直す */
    public fun rewind() {
        pos = firstFramePos
        frameIndex = 0
        refs = arrayOfNulls(4)
        lfFrames = arrayOfNulls(4)
        finished = false
        visibleFrames = 0
        invisibleSinceVisible = 0
    }

    /** 次に表示するフレーム。終わりなら null */
    public fun nextFrame(): JxlFrame? {
        try {
            while (!finished) {
                val frame = decodeOneFrame()
                if (frame != null) return frame
            }
            return null
        } catch (e: JxlException) {
            finished = true
            throw e
        } catch (e: RuntimeException) {
            // 壊れたデータで配列の範囲外などになった場合
            finished = true
            throw JxlException("invalid JPEG XL data (${e.javaClass.simpleName})")
        }
    }

    private fun decodeOneFrame(): JxlFrame? {
        if (pos >= codestream.size) fail("unexpected end of codestream")
        val br = BitReader(codestream, pos, codestream.size)
        val fh = FrameHeader.parse(br, header)
        // 切り出したフレームは画像と別に大きさを持つので、それも上限を確かめる
        if (fh.width.toLong() * fh.height > minOf(options.maxPixels, Int.MAX_VALUE.toLong())) fail("frame too large")
        val toc = Toc.parse(br, fh)
        val end = toc.end
        if (end > codestream.size) fail("truncated frame")
        // libjxl と同じく、表示フレームの番号は 1 から、非表示フレームの番号は直前の表示フレームから数える
        val seed = if (fh.isKeyframe) (visibleFrames + 1).toLong() shl 32
        else (visibleFrames.toLong() shl 32) + invisibleSinceVisible + 1
        val index = frameIndex
        // 保存しても後のフレームが使わない参照は持たない
        val keepRef = fh.canReference && referenceUsedLater(index, fh.saveAsReference)
        // 表示するフレームで、参照として持つ必要がなく、キャンバス全体を置き換えるなら逐次描画できる
        // (全体の浮動小数点の平面を持たない)。静止画の最後の 1 枚と、前のフレームに重ねないアニメーションのフレーム
        val streamable = fh.isKeyframe && !keepRef && fh.x0 == 0 && fh.y0 == 0 &&
            fh.width == header.width && fh.height == header.height &&
            fh.blending.mode == BlendMode.REPLACE && fh.ecBlending.all { it.mode == BlendMode.REPLACE }
        val stream = if (streamable) FrameDecoder.StreamOutput() else null
        val decoded = FrameDecoder(header, fh, codestream, toc, refs, lfFrames, executor, seed).decode(stream)
        pos = end.toInt()
        frameIndex++
        if (fh.isLast) finished = true
        val streamed = stream?.pixels
        // 使われない参照は、このフレームの復号 (パッチ・ブレンド) が終わったら手放す
        val dropRef = fh.canReference && !keepRef
        if (streamed != null) {
            if (dropRef) refs[fh.saveAsReference] = null
            visibleFrames++
            invisibleSinceVisible = 0
            return orient(streamed, header.width, header.height, frameDuration(fh))
        }
        val img = decoded!!

        if (fh.frameType == FrameType.LF) {
            lfFrames[fh.lfLevel - 1] = img
            invisibleSinceVisible++
            return null
        }

        if (fh.saveBeforeCt && fh.canReference) {
            if (keepRef) refs[fh.saveAsReference] = ReferenceFrame(img)
            if (!fh.isKeyframe) {
                invisibleSinceVisible++
                return null
            }
        }
        if (fh.frameType == FrameType.REFERENCE_ONLY) {
            invisibleSinceVisible++
            return null
        }

        // 色変換 (ブレンドは出力の色空間で行う)
        var t0 = System.nanoTime()
        val converted = if (fh.saveBeforeCt && keepRef) copyImage(img) else img
        convertColor(converted)
        if (FrameDecoder.TIMING) System.err.println("  %-14s %7.2f ms".format("color", (System.nanoTime() - t0) / 1e6))
        t0 = System.nanoTime()
        val canvas = blend(fh, converted)
        if (keepRef && !fh.saveBeforeCt) refs[fh.saveAsReference] = ReferenceFrame(canvas)
        if (dropRef) refs[fh.saveAsReference] = null

        if (!fh.isKeyframe) {
            invisibleSinceVisible++
            return null
        }
        visibleFrames++
        invisibleSinceVisible = 0
        val out = toOutput(canvas, frameDuration(fh))
        if (FrameDecoder.TIMING) System.err.println("  %-14s %7.2f ms".format("blend+output", (System.nanoTime() - t0) / 1e6))
        return out
    }

    // 参照を使うか決めるために先読みしたフレームのヘッダー (先頭から順に)。必要になった分だけ読む
    private val scannedHeaders = ArrayList<FrameHeader>()
    private var scanPos = -1
    private var scanDone = false

    /** j 番目のフレームのヘッダー。読めない (壊れている・終わりを越える) なら null */
    private fun headerAt(j: Int): FrameHeader? {
        if (scanPos < 0) scanPos = firstFramePos
        while (scannedHeaders.size <= j && !scanDone) {
            try {
                if (scanPos >= codestream.size) {
                    scanDone = true
                    break
                }
                val br = BitReader(codestream, scanPos, codestream.size)
                val fh = FrameHeader.parse(br, header)
                val end = Toc.parse(br, fh).end
                if (end > codestream.size) {
                    scanDone = true
                    break
                }
                scannedHeaders.add(fh)
                if (fh.isLast) scanDone = true
                scanPos = end.toInt()
            } catch (e: JxlException) {
                // 壊れたところは、そのフレームを復号するときにエラーにする
                scanDone = true
            } catch (e: RuntimeException) {
                scanDone = true
            }
        }
        return scannedHeaders.getOrNull(j)
    }

    /**
     * index 番目のフレームが参照 [slot] に保存する内容を、後のフレームが使うか。
     * 同じ場所に次の参照が保存されるフレームまでを見る。ヘッダーを読めなかったフレームやパッチを使うフレームがあれば、使うものとする
     */
    private fun referenceUsedLater(index: Int, slot: Int): Boolean {
        var j = index + 1
        while (true) {
            val f = headerAt(j++) ?: return true
            if (f.hasPatches) return true
            if (f.isNormalFrame) {
                val fullReplace = f.x0 == 0 && f.y0 == 0 && f.width == header.width && f.height == header.height &&
                    f.blending.mode == BlendMode.REPLACE && f.ecBlending.all { it.mode == BlendMode.REPLACE }
                if (!fullReplace && (f.blending.source == slot || f.ecBlending.any { it.source == slot })) return true
            }
            if (f.canReference && f.saveAsReference == slot) return false
            if (f.isLast) return false
        }
    }

    private fun frameDuration(fh: FrameHeader): Int {
        val anim = header.animation ?: return 0
        if (anim.tpsNumerator <= 0) return 0
        return (fh.duration.toLong() * 1000L * anim.tpsDenominator / anim.tpsNumerator).toInt()
    }

    private fun copyImage(img: FrameImage) = FrameImage(
        img.header, img.width, img.height,
        Array(img.color.size) { img.color[it].copyOf() },
        Array(img.extra.size) { img.extra[it].copyOf() },
        img.colorTransformPending,
    )

    private fun convertColor(img: FrameImage) {
        if (header.xybEncoded) {
            if (img.colorTransformPending) {
                ColorConvert.xybToSrgb(img, header, executor)
                if (header.isGrey) img.color = arrayOf(img.color[1])
            }
        } else {
            if (img.header.doYCbCr) ColorConvert.ycbcrToRgb(img)
            ColorConvert.toSrgb(img, header.colourEncoding)
        }
        img.colorTransformPending = false
    }

    /** 参照フレームの上にフレームを重ねて、画像全体の大きさのキャンバスを作る */
    private fun blend(fh: FrameHeader, img: FrameImage): FrameImage {
        val w = header.width
        val h = header.height
        val colorCount = img.color.size
        val numExtra = img.extra.size
        if (fh.x0 == 0 && fh.y0 == 0 && img.width == w && img.height == h &&
            fh.blending.mode == BlendMode.REPLACE && fh.ecBlending.all { it.mode == BlendMode.REPLACE }
        ) {
            return img
        }
        val outColor = arrayOfNulls<FloatArray>(colorCount)
        val outExtra = arrayOfNulls<FloatArray>(numExtra)
        // 元のアルファ (ブレンド中に書き換わる前の値)
        val baseAlphaCache = HashMap<Int, FloatArray?>()
        fun baseChannel(source: Int, ch: Int): FloatArray? {
            val r = refs[source]?.image ?: return null
            val arr = if (ch < colorCount) r.color.getOrNull(ch) else r.extra.getOrNull(ch - colorCount)
            if (arr == null || r.width != w || r.height != h) {
                if (arr == null) return null
                // 大きさが違う参照 (切り出し) はキャンバスに置き直す
                val out = FloatArray(w * h)
                for (y in 0 until minOf(h, r.height)) System.arraycopy(arr, y * r.width, out, y * w, minOf(w, r.width))
                return out
            }
            return arr
        }
        val fx0 = maxOf(0, fh.x0)
        val fy0 = maxOf(0, fh.y0)
        val fx1 = minOf(w, fh.x0 + img.width)
        val fy1 = minOf(h, fh.y0 + img.height)
        val bw = fx1 - fx0
        val bh = fy1 - fy0
        for (ch in 0 until colorCount + numExtra) {
            val bi = if (ch < colorCount) fh.blending else fh.ecBlending[ch - colorCount]
            val src = baseChannel(bi.source, ch)
            val base = src?.copyOf() ?: FloatArray(w * h)
            if (bw > 0 && bh > 0) {
                val newData = if (ch < colorCount) img.color[ch] else img.extra[ch - colorCount]
                var mode = when (bi.mode) {
                    BlendMode.REPLACE -> Features.BLEND_REPLACE
                    BlendMode.ADD -> Features.BLEND_ADD
                    BlendMode.BLEND -> Features.BLEND_ALPHA
                    BlendMode.MUL_ADD -> Features.BLEND_MULADD
                    else -> Features.BLEND_MUL
                }
                var baseAlpha: FloatArray? = null
                var newAlpha: FloatArray? = null
                var premul = false
                val useAlpha = (bi.mode == BlendMode.BLEND || bi.mode == BlendMode.MUL_ADD) && numExtra > 0
                if (useAlpha) {
                    val a = bi.alphaChannel
                    if (ch == colorCount + a) {
                        mode = if (bi.mode == BlendMode.BLEND) Features.BLEND_MIX_ALPHA else Features.BLEND_SKIP
                    } else {
                        baseAlpha = baseAlphaCache.getOrPut(bi.source * 1000 + a) { baseChannel(bi.source, colorCount + a) }
                        newAlpha = img.extra[a]
                        premul = header.extraChannels[a].alphaAssociated
                    }
                }
                val nx = fx0 - fh.x0
                val ny = fy0 - fh.y0
                if (newAlpha != null && (baseAlpha == null)) baseAlpha = null
                Features.blend(
                    mode, bi.clamp, false, premul,
                    base, w, fx0, fy0, baseAlpha,
                    newData, img.width, nx, ny, newAlpha,
                    bw, bh,
                )
            }
            if (ch < colorCount) outColor[ch] = base else outExtra[ch - colorCount] = base
        }
        return FrameImage(fh, w, h, Array(colorCount) { outColor[it]!! }, Array(numExtra) { outExtra[it]!! }, false)
    }

    private fun toOutput(img: FrameImage, duration: Int): JxlFrame {
        val w = img.width
        val h = img.height
        val n = w * h
        val pixels = IntArray(n)
        val r = img.color[0]
        val g = if (img.color.size >= 3) img.color[1] else r
        val b = if (img.color.size >= 3) img.color[2] else r
        val alphaIdx = header.alphaChannel
        val alpha = if (alphaIdx >= 0) img.extra[alphaIdx] else null
        val premul = alphaIdx >= 0 && header.extraChannels[alphaIdx].alphaAssociated
        parallelRows(executor, h) { y0, y1 ->
            for (i in y0 * w until y1 * w) {
                var a = 1f
                if (alpha != null) a = alpha[i].coerceIn(0f, 1f)
                var rv = r[i]
                var gv = g[i]
                var bv = b[i]
                if (premul && a > 0f) {
                    rv /= a; gv /= a; bv /= a
                }
                pixels[i] = (to8(a) shl 24) or (to8(rv) shl 16) or (to8(gv) shl 8) or to8(bv)
            }
        }
        return orient(pixels, w, h, duration)
    }

    /** 向き (orientation) を反映する */
    private fun orient(pixels: IntArray, w: Int, h: Int, duration: Int): JxlFrame {
        val n = w * h
        val o = if (options.applyOrientation) header.orientation else 1
        if (o == 1) return JxlFrame(w, h, pixels, duration)
        val swap = o >= 5
        val ow = if (swap) h else w
        val oh = if (swap) w else h
        val out = IntArray(n)
        for (y in 0 until h) for (x in 0 until w) {
            val (dx, dy) = when (o) {
                2 -> Pair(w - 1 - x, y)
                3 -> Pair(w - 1 - x, h - 1 - y)
                4 -> Pair(x, h - 1 - y)
                5 -> Pair(y, x)
                6 -> Pair(h - 1 - y, x)
                7 -> Pair(h - 1 - y, w - 1 - x)
                else -> Pair(y, w - 1 - x)
            }
            out[dy * ow + dx] = pixels[y * w + x]
        }
        return JxlFrame(ow, oh, out, duration)
    }

    public companion object {
        private fun to8(v: Float): Int {
            val f = v * 255f + 0.5f
            return if (f <= 0f) 0 else if (f >= 255f) 255 else f.toInt()
        }

        @Volatile
        private var pool: ExecutorService? = null
        private var poolThreads = 0

        @Synchronized
        private fun sharedPool(threads: Int): ExecutorService {
            val p = pool
            if (p != null && poolThreads == threads) return p
            val counter = AtomicInteger()
            val factory = ThreadFactory { r ->
                Thread(r, "jxl-decoder-${counter.incrementAndGet()}").apply { isDaemon = true }
            }
            val np = Executors.newFixedThreadPool(threads, factory)
            pool = np
            poolThreads = threads
            return np
        }

        /** 先頭のバイト列が JPEG XL (コードストリームまたはコンテナ) か */
        @JvmStatic
        public fun isJxl(head: ByteArray): Boolean {
            if (head.size >= 2 && head[0] == 0xff.toByte() && head[1] == 0x0a.toByte()) return true
            return head.size >= 12 && Container.SIGNATURE.indices.all { head[it] == Container.SIGNATURE[it] }
        }

        /** 最初のフレームだけを復号する (静止画用) */
        @JvmStatic
        @JvmOverloads
        public fun decode(data: ByteArray, options: JxlDecodeOptions = JxlDecodeOptions()): JxlFrame =
            JxlDecoder(data, options).nextFrame() ?: throw JxlException("no frames")
    }
}

internal fun ImageHeader.withSize(w: Int, h: Int) = ImageHeader(
    w, h, orientation, animation, bitDepth, modular16BitBuffers, extraChannels, xybEncoded, colourEncoding,
    toneMapping, opsinInverseMatrix, opsinBias, quantBias, quantBiasNumerator, up2Weights, up4Weights, up8Weights, 0, 0,
)

/** ISOBMFF 形式のコンテナ */
internal object Container {
    val SIGNATURE = byteArrayOf(0, 0, 0, 0x0c, 0x4a, 0x58, 0x4c, 0x20, 0x0d, 0x0a, 0x87.toByte(), 0x0a)

    fun extractCodestream(data: ByteArray): ByteArray {
        if (data.size >= 2 && data[0] == 0xff.toByte() && data[1] == 0x0a.toByte()) return data
        if (data.size < 12 || !SIGNATURE.indices.all { data[it] == SIGNATURE[it] }) fail("not a JPEG XL file")
        var p = 0
        val parts = ArrayList<Pair<Int, Int>>()
        var whole: Pair<Int, Int>? = null
        while (p + 8 <= data.size) {
            var size = readBe32(data, p).toLong() and 0xffffffffL
            val type = String(data, p + 4, 4, Charsets.ISO_8859_1)
            var header = 8
            if (size == 1L) {
                if (p + 16 > data.size) fail("truncated box")
                size = (readBe32(data, p + 8).toLong() shl 32) or (readBe32(data, p + 12).toLong() and 0xffffffffL)
                header = 16
            } else if (size == 0L) {
                size = (data.size - p).toLong()
            }
            if (size < header || p + size > data.size) fail("invalid box size")
            val start = p + header
            val end = (p + size).toInt()
            when (type) {
                "jxlc" -> whole = Pair(start, end)
                "jxlp" -> if (end - start >= 4) parts.add(Pair(start + 4, end))
            }
            p = end
        }
        whole?.let { return data.copyOfRange(it.first, it.second) }
        if (parts.isEmpty()) fail("no codestream in container")
        val total = parts.sumOf { it.second - it.first }
        val out = ByteArray(total)
        var o = 0
        for ((s, e) in parts) {
            System.arraycopy(data, s, out, o, e - s)
            o += e - s
        }
        return out
    }

    private fun readBe32(d: ByteArray, p: Int): Int =
        ((d[p].toInt() and 0xff) shl 24) or ((d[p + 1].toInt() and 0xff) shl 16) or
            ((d[p + 2].toInt() and 0xff) shl 8) or (d[p + 3].toInt() and 0xff)
}

/** ICC プロファイルのストリームを読み飛ばす (色変換には使わない) */
internal object Icc {
    private fun kind1(b: Int): Int = when {
        b in 'a'.code..'z'.code || b in 'A'.code..'Z'.code -> 0
        b in '0'.code..'9'.code || b == '.'.code || b == ','.code -> 1
        b == 0 -> 2
        b == 1 -> 3
        b < 16 -> 4
        b == 255 -> 6
        b > 240 -> 5
        else -> 7
    }

    private fun kind2(b: Int): Int = when {
        b in 'a'.code..'z'.code || b in 'A'.code..'Z'.code -> 0
        b in '0'.code..'9'.code || b == '.'.code || b == ','.code -> 1
        b < 16 -> 2
        b > 240 -> 3
        else -> 4
    }

    fun skip(br: BitReader) {
        val size = br.readU64()
        if (size > (1L shl 28)) fail("ICC profile too large")
        val code = EntropyCode.parse(br, 41)
        val r = code.newReader()
        r.begin(br)
        var b1 = 0
        var b2 = 0
        for (i in 0 until size) {
            val ctx = if (i <= 128) 0 else 1 + kind1(b1) + kind2(b2) * 8
            val v = r.read(br, ctx) and 0xff
            b2 = b1
            b1 = v
            if ((i and 0xffff) == 0L) br.checkOverrun()
        }
        r.finish()
    }
}
