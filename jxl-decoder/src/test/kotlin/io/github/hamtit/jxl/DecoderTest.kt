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
 * 非可逆は丸めの違いで ±1 まで、可逆 (lossless_*, *_q100) は完全一致を求める。
 */
class DecoderTest {
    private val dir = File(javaClass.getResource("/fixtures")!!.toURI())

    private fun fixtures(): List<File> = dir.listFiles { f -> f.name.endsWith(".jxl") }!!.sortedBy { it.name }

    private fun compare(name: String, threads: Int) {
        val data = File(dir, "$name.jxl").readBytes()
        val decoder = JxlDecoder(data, JxlDecodeOptions(threads = threads))
        val lossless = name.startsWith("lossless") || name.endsWith("q100")
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
}
