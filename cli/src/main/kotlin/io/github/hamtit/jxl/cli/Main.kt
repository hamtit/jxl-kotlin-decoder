package io.github.hamtit.jxl.cli

import io.github.hamtit.jxl.JxlDecodeOptions
import io.github.hamtit.jxl.JxlDecoder
import io.github.hamtit.jxl.JxlException
import io.github.hamtit.jxl.JxlFrame
import java.awt.image.BufferedImage
import java.awt.image.ColorModel
import java.awt.image.DataBufferByte
import java.awt.image.DataBufferInt
import java.awt.image.DirectColorModel
import java.awt.image.Raster
import java.io.File
import javax.imageio.IIOImage
import javax.imageio.ImageIO
import javax.imageio.ImageWriteParam
import kotlin.system.exitProcess

private fun usage(): Nothing {
    System.err.println(
        """
        使い方:
          decode <in.jxl> <out.png> [--threads N]     最初のフレームを PNG に書き出す
          frames <in.jxl> <outdir> [--threads N]      全フレームを PNG に書き出す
          bench <in.jxl>... [--iterations N] [--threads N]  デコード時間を測る
          compare <in.jxl> <ref.png> [--threads N]    参照画像との PSNR・最大誤差を表示する
        """.trimIndent(),
    )
    exitProcess(2)
}

private fun options(args: List<String>): JxlDecodeOptions =
    JxlDecodeOptions(threads = intOpt(args, "--threads", Runtime.getRuntime().availableProcessors()))

/** 数値のオプション。値が無い・数でない・1 未満なら使い方を表示して終わる */
private fun intOpt(args: List<String>, name: String, def: Int): Int {
    val i = args.indexOf(name)
    if (i < 0) return def
    val v = args.getOrNull(i + 1)?.toIntOrNull() ?: usage()
    if (v < 1) usage()
    return v
}

/** 端末に出す前に制御文字 (ESC・改行・C1 制御文字など) を \\uXXXX にする */
private fun escapeControl(s: String): String = buildString {
    for (c in s) if (Character.isISOControl(c)) append("\\u%04x".format(c.code)) else append(c)
}

private fun positional(args: List<String>): List<String> {
    val out = ArrayList<String>()
    var i = 0
    while (i < args.size) {
        if (args[i].startsWith("--")) i += 2 else out.add(args[i++])
    }
    return out
}

private fun toImage(f: JxlFrame, hasAlpha: Boolean): BufferedImage {
    val px = f.pixels
    var allOpaque = !hasAlpha
    var allGray = true
    if (hasAlpha) {
        allOpaque = true
        for (p in px) {
            if ((p ushr 24) != 0xff) {
                allOpaque = false
                allGray = false
                break
            }
            val b = p and 0xff
            if (((p ushr 8) and 0xff) != b || ((p ushr 16) and 0xff) != b) {
                allGray = false
            }
        }
    } else {
        for (p in px) {
            val b = p and 0xff
            if (((p ushr 8) and 0xff) != b || ((p ushr 16) and 0xff) != b) {
                allGray = false
                break
            }
        }
    }
    if (allOpaque && allGray) {
        val img = BufferedImage(f.width, f.height, BufferedImage.TYPE_BYTE_GRAY)
        val dst = (img.raster.dataBuffer as DataBufferByte).data
        for (i in px.indices) dst[i] = (px[i] and 0xff).toByte()
        return img
    }
    val db = DataBufferInt(px, px.size)
    if (allOpaque) {
        val masks = intArrayOf(0x00ff0000, 0x0000ff00, 0x000000ff)
        val raster = Raster.createPackedRaster(db, f.width, f.height, f.width, masks, null)
        val cm = DirectColorModel(24, 0x00ff0000, 0x0000ff00, 0x000000ff)
        return BufferedImage(cm, raster, false, null)
    }
    val masks = intArrayOf(0x00ff0000, 0x0000ff00, 0x000000ff, 0xff000000.toInt())
    val raster = Raster.createPackedRaster(db, f.width, f.height, f.width, masks, null)
    return BufferedImage(ColorModel.getRGBdefault(), raster, false, null)
}

private fun writePng(img: BufferedImage, file: File) {
    val writers = ImageIO.getImageWritersByFormatName("png")
    if (!writers.hasNext()) {
        ImageIO.write(img, "png", file)
        return
    }
    val writer = writers.next()
    try {
        ImageIO.createImageOutputStream(file).use { ios ->
            writer.output = ios
            val param = writer.defaultWriteParam
            if (param.canWriteCompressed()) {
                param.compressionMode = ImageWriteParam.MODE_EXPLICIT
                param.compressionQuality = 0.0f
            }
            writer.write(null, IIOImage(img, null, null), param)
        }
    } finally {
        writer.dispose()
    }
}

public fun main(argv: Array<String>) {
    if (argv.isEmpty()) usage()
    val args = argv.drop(1)
    val pos = positional(args)
    when (argv[0]) {
        "decode" -> {
            if (pos.size < 2) usage()
            val t0 = System.nanoTime()
            val dec = JxlDecoder(File(pos[0]).readBytes(), options(args))
            val frame = dec.nextFrame() ?: throw JxlException("no frames")
            val t1 = System.nanoTime()
            writePng(toImage(frame, dec.info.hasAlpha), File(pos[1]))
            println("${frame.width}x${frame.height} ${"%.1f".format((t1 - t0) / 1e6)} ms")
        }
        "frames" -> {
            if (pos.size < 2) usage()
            val dir = File(pos[1]).apply { mkdirs() }
            val dec = JxlDecoder(File(pos[0]).readBytes(), options(args))
            var i = 0
            while (true) {
                val f = dec.nextFrame() ?: break
                writePng(toImage(f, dec.info.hasAlpha), File(dir, "frame%04d.png".format(i)))
                println("frame $i ${f.width}x${f.height} ${f.durationMillis} ms")
                i++
            }
        }
        "bench" -> {
            if (pos.isEmpty()) usage()
            val iterations = intOpt(args, "--iterations", 10)
            val opts = options(args)
            for (path in pos) {
                val data = File(path).readBytes()
                // ウォームアップ
                repeat(3) { decodeAll(data, opts) }
                val times = DoubleArray(iterations)
                for (k in 0 until iterations) {
                    val t0 = System.nanoTime()
                    decodeAll(data, opts)
                    times[k] = (System.nanoTime() - t0) / 1e6
                }
                times.sort()
                println("%s median %.2f ms min %.2f ms".format(File(path).name, times[iterations / 2], times[0]))
            }
        }
        "compare" -> {
            if (pos.size < 2) usage()
            val frame = JxlDecoder.decode(File(pos[0]).readBytes(), options(args))
            val ref = ImageIO.read(File(pos[1])) ?: run {
                System.err.println("読めない画像形式: ${escapeControl(pos[1])}")
                exitProcess(1)
            }
            if (ref.width != frame.width || ref.height != frame.height) {
                println("size mismatch: ${frame.width}x${frame.height} vs ${ref.width}x${ref.height}")
                exitProcess(1)
            }
            var se = 0.0
            var maxErr = 0
            var n = 0L
            val hasAlpha = ref.colorModel.hasAlpha()
            for (y in 0 until ref.height) for (x in 0 until ref.width) {
                val a = frame.pixels[y * frame.width + x]
                val b = ref.getRGB(x, y)
                val channels = if (hasAlpha) 4 else 3
                for (c in 0 until channels) {
                    val va = (a ushr (c * 8)) and 0xff
                    val vb = (b ushr (c * 8)) and 0xff
                    val d = va - vb
                    se += (d * d).toDouble()
                    if (kotlin.math.abs(d) > maxErr) maxErr = kotlin.math.abs(d)
                    n++
                }
            }
            val mse = se / n
            val psnr = if (mse == 0.0) Double.POSITIVE_INFINITY else 10 * Math.log10(255.0 * 255.0 / mse)
            println("PSNR %.2f dB, max error %d".format(psnr, maxErr))
        }
        else -> usage()
    }
}

private fun decodeAll(data: ByteArray, opts: JxlDecodeOptions) {
    val dec = JxlDecoder(data, opts)
    while (dec.nextFrame() != null) {
    }
}
