package io.github.hamtit.jxl

import java.io.File
import javax.imageio.ImageIO
import kotlin.math.abs
import kotlin.test.Test
import kotlin.test.assertEquals
import kotlin.test.assertFailsWith
import kotlin.test.assertTrue
import kotlin.test.fail

/**
 * fixtures の各 .jxl を復号し、libjxl (djxl) で 8 ビットに復号した fixtures/<名前>.<フレーム番号>.png と比べる。
 * 非可逆は丸めの違いで ±1 まで、可逆 (lossless_*, anim_lossless, *_q100) は完全一致を求める。
 */
class DecoderTest {
    private val dir = File(javaClass.getResource("/fixtures")!!.toURI())

    private fun fixtures(): List<File> = dir.listFiles { f -> f.name.endsWith(".jxl") }!!.sortedBy { it.name }

    private fun compare(name: String, threads: Int) {
        val data = File(dir, "$name.jxl").readBytes()
        val decoder = JxlDecoder(data, JxlDecodeOptions(threads = threads))
        val lossless = name.startsWith("lossless") || name == "anim_lossless" || name.endsWith("q100")
        var index = 0
        while (true) {
            val frame = decoder.nextFrame() ?: break
            val refFile = File(dir, "$name.$index.png")
            if (!refFile.exists()) fail("$name: extra frame $index")
            val ref = ImageIO.read(refFile)
            assertEquals(ref.width, frame.width, "$name width")
            assertEquals(ref.height, frame.height, "$name height")
            var maxDiff = 0
            for (y in 0 until ref.height) for (x in 0 until ref.width) {
                val a = frame.pixels[y * frame.width + x]
                val b = ref.getRGB(x, y)
                for (c in 0 until 4) {
                    val d = abs(((a ushr (8 * c)) and 0xff) - ((b ushr (8 * c)) and 0xff))
                    if (d > maxDiff) maxDiff = d
                }
            }
            assertTrue(maxDiff <= (if (lossless) 0 else 1), "$name frame $index: max difference $maxDiff")
            index++
        }
        assertTrue(!File(dir, "$name.$index.png").exists(), "$name: missing frame $index")
        assertTrue(index > 0)
    }

    @Test
    fun allFixturesMatchLibjxl() {
        val names = fixtures().map { it.name.removeSuffix(".jxl") }
        assertTrue(names.size >= 20)
        for (n in names) compare(n, threads = 1)
    }

    @Test
    fun multiThreadedOutputIsTheSame() {
        for (f in fixtures()) {
            val data = f.readBytes()
            val a = JxlDecoder(data, JxlDecodeOptions(threads = 1))
            val b = JxlDecoder(data, JxlDecodeOptions(threads = 4))
            while (true) {
                val fa = a.nextFrame()
                val fb = b.nextFrame()
                if (fa == null || fb == null) {
                    assertTrue(fa == null && fb == null, f.name)
                    break
                }
                assertTrue(fa.pixels.contentEquals(fb.pixels), f.name)
            }
        }
    }

    @Test
    fun animationInfoAndDurations() {
        val decoder = JxlDecoder(File(dir, "anim_lossy.jxl").readBytes())
        assertTrue(decoder.info.isAnimated)
        assertTrue(decoder.info.hasAlpha)
        assertEquals(40, decoder.info.width)
        assertEquals(30, decoder.info.height)
        val durations = generateSequence { decoder.nextFrame() }.map { it.durationMillis }.toList()
        assertEquals(listOf(40, 50, 60, 70), durations)
        // 巻き戻すと最初から同じフレームが出る
        val first = JxlDecoder(File(dir, "anim_lossy.jxl").readBytes()).nextFrame()!!
        decoder.rewind()
        assertTrue(first.pixels.contentEquals(decoder.nextFrame()!!.pixels))
    }

    @Test
    fun returnedFramesStayUnchangedAndRewindRepeatsThem() {
        // 置き換えで重ねるアニメーションはキャンバスをその場で書き換えるので、
        // 返したフレームが後のフレームの復号で変わらないこと、巻き戻すと同じフレームが出ることを確かめる
        for (name in listOf("anim_lossy", "anim_lossless")) {
            val decoder = JxlDecoder(File(dir, "$name.jxl").readBytes(), JxlDecodeOptions(threads = 1))
            // 返された時点の写しを取っておく
            val returned = generateSequence { decoder.nextFrame() }.map { it to it.pixels.copyOf() }.toList()
            val frames = returned.map { it.first }
            val snapshots = returned.map { it.second }
            decoder.rewind()
            val again = generateSequence { decoder.nextFrame() }.toList()
            assertEquals(frames.size, again.size, name)
            for (i in frames.indices) {
                assertTrue(frames[i].pixels.contentEquals(snapshots[i]), "$name frame $i changed")
                assertTrue(again[i].pixels.contentEquals(snapshots[i]), "$name frame $i after rewind")
            }
            // 出力配列を呼び出し元が変更しても、後続フレームの参照は書き換わらない
            decoder.rewind()
            for ((i, expected) in snapshots.withIndex()) {
                val pixels = decoder.nextFrame()!!.pixels
                assertTrue(pixels.contentEquals(expected), "$name frame $i after editing output")
                pixels.fill(0x12345678)
            }
        }
    }

    @Test
    fun detectsSignatures() {
        assertTrue(JxlDecoder.isJxl(byteArrayOf(0xff.toByte(), 0x0a)))
        assertTrue(JxlDecoder.isJxl(File(dir, "vardct_container.jxl").readBytes().copyOf(12)))
        assertTrue(!JxlDecoder.isJxl(byteArrayOf(0x89.toByte(), 0x50, 0x4e, 0x47)))
    }

    @Test
    fun brokenDataThrowsJxlException() {
        val data = File(dir, "vardct_d1.jxl").readBytes()
        // 途中で切れたデータ
        for (len in listOf(1, 2, 10, data.size / 3, data.size - 1)) {
            assertFailsWith<JxlException>("truncated to $len") {
                val d = JxlDecoder(data.copyOf(len))
                while (d.nextFrame() != null) {
                }
            }
        }
        // 中身を壊したデータは、例外になるか何らかの画像になる (固まったり別の例外を投げたりしない)
        val rnd = java.util.Random(7)
        repeat(200) {
            val broken = data.copyOf()
            repeat(3) { broken[16 + rnd.nextInt(broken.size - 16)] = rnd.nextInt(256).toByte() }
            try {
                val d = JxlDecoder(broken)
                while (d.nextFrame() != null) {
                }
            } catch (e: JxlException) {
                // 想定どおり
            }
        }
    }

    @Test
    fun skipBitsAdvancesAndRejectsSkipsPastTheEnd() {
        val data = ByteArray(64) { (it * 37 + 11).toByte() }
        for (n in listOf(0L, 3L, 13L, 63L, 64L, 100L, 301L)) {
            val skipped = BitReader(data).apply { readBits(5); skipBits(n) }
            val read = BitReader(data).apply { readBits(5) }
            var left = n
            while (left > 0) {
                val k = minOf(left, 16L).toInt()
                read.readBits(k)
                left -= k
            }
            assertEquals(read.readBits(16), skipped.readBits(16), "skip $n")
        }
        // 入力より長い読み飛ばし (壊れた拡張データ) は、ループで回さずにすぐ失敗する
        assertFailsWith<JxlException> { BitReader(data).skipBits(1L shl 60) }
        assertFailsWith<JxlException> { BitReader(data).skipBits(64L * 8 + 1) }
    }

    @Test
    fun floatSamplesIncludingSubnormals() {
        val f16 = BitDepth(16, 5)
        assertEquals(1.0f, f16.parseFloatSample(0x3c00))
        assertEquals(-2.0f, f16.parseFloatSample(0xc000))
        assertEquals(Math.scalb(1f, -24), f16.parseFloatSample(0x0001))
        assertEquals(Math.scalb(1023f, -24), f16.parseFloatSample(0x03ff))
        assertEquals(-0.0f, f16.parseFloatSample(0x8000))
        val f32 = BitDepth(32, 8)
        for (v in listOf(1.5f, -3.25f, Float.MIN_VALUE, 1e-40f, -1e-39f)) {
            assertEquals(v, f32.parseFloatSample(v.toRawBits()))
        }
    }
}
