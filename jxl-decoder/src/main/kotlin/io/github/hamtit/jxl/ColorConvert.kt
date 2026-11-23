package io.github.hamtit.jxl

import java.util.concurrent.ExecutorService

/** 色変換 (XYB / YCbCr / 各種伝達関数 → sRGB) */
internal object ColorConvert {

    fun srgbEncode(v: Float): Float {
        if (v <= 0.0031308f) return if (v <= 0f) 0f else v * 12.92f
        if (v >= 1f) return 1f
        return (1.055 * Math.pow(v.toDouble(), 1.0 / 2.4) - 0.055).toFloat()
    }

    /** sRGB の符号化 (libjxl と同じ有理多項式近似、誤差 5e-7 程度)。0..1 に切り詰める */
    fun srgbEncodeFast(v: Float): Float {
        if (!(v > 0.0031308f)) return if (v > 0f) v * 12.92f else 0f
        if (v >= 1f) return 1f
        val s = Math.sqrt(v.toDouble()).toFloat()
        val p = -5.135152395e-04f + s * (5.287254571e-03f + s * (3.903842876e-01f + s * (1.474205315e+00f + s * 7.352629620e-01f)))
        val q = 1.004519624e-02f + s * (3.036675394e-01f + s * (1.340816930e+00f + s * (9.258482155e-01f + s * 2.424867759e-02f)))
        return p / q
    }

    /*
     * 線形の値から 8 ビットの sRGB への変換表。2^-16 .. 1 の値を、浮動小数点のビット列の指数と仮数の上位 9 ビット
     * (8192 区間 = 8 KB) で引き、各区間の中央値を四捨五入した 8 ビット値を持つ。
     * 99.2% 以上は厳密な四捨五入と一致し、残りの境界付近でも差は最大 ±1。
     * 2^-16 未満は 0 になる (2^-16 * 12.92 * 255 < 0.5)。
     */
    private val SRGB_MIN_BITS = java.lang.Float.floatToRawIntBits(1f / 65536f)
    private val SRGB_ONE_BITS = java.lang.Float.floatToRawIntBits(1f)
    private val SRGB_TABLE = ByteArray((SRGB_ONE_BITS - SRGB_MIN_BITS) ushr 14)

    private fun srgbExact8(v: Double): Int {
        val e = if (v <= 0.0031308) v * 12.92 else 1.055 * Math.pow(v, 1.0 / 2.4) - 0.055
        return Math.round(e * 255.0).toInt().coerceIn(0, 255)
    }

    init {
        for (i in SRGB_TABLE.indices) {
            val midBits = SRGB_MIN_BITS + (i shl 14) + (1 shl 13)
            SRGB_TABLE[i] = srgbExact8(Float.fromBits(midBits).toDouble()).toByte()
        }
    }

    /** 線形の値を 8 ビットの sRGB にする (sRGB の伝達関数をかけて四捨五入。0..1 の外は切り詰める) */
    fun linearToSrgb8(v: Float): Int {
        val bits = java.lang.Float.floatToRawIntBits(v)
        val idx = (bits - SRGB_MIN_BITS) ushr 14
        // 2^-16 .. 1 未満なら idx は 0 until SRGB_TABLE.size (負数・極小値・1 以上・NaN はすべて範囲外になる)
        if (idx < SRGB_TABLE.size) return SRGB_TABLE[idx].toInt() and 0xff
        return if (bits >= SRGB_ONE_BITS && bits <= 0x7f800000) 255 else 0
    }

    /** XYB を sRGB (伝達関数適用済み、0..1 付近) に変換する */
    fun xybToSrgb(img: FrameImage, image: ImageHeader, executor: ExecutorService?) {
        val x = img.color[0]
        val y = img.color[1]
        val b = img.color[2]
        val ob = image.opsinBias
        val cb0 = Math.cbrt(ob[0].toDouble()).toFloat()
        val cb1 = Math.cbrt(ob[1].toDouble()).toFloat()
        val cb2 = Math.cbrt(ob[2].toDouble()).toFloat()
        val its = 255f / image.toneMapping.intensityTarget
        val m = image.opsinInverseMatrix
        val target = image.colourEncoding
        val grey = !target.wantIcc && target.colourSpace == ColourSpace.GREY
        val n = x.size
        parallelRows(executor, img.height) { y0, y1 ->
            for (i in y0 * img.width until minOf(n, y1 * img.width)) {
                val gl = y[i] + x[i] - cb0
                val gm = y[i] - x[i] - cb1
                val gs = b[i] - cb2
                val ml = (gl * gl * gl + ob[0]) * its
                val mm = (gm * gm * gm + ob[1]) * its
                val ms = (gs * gs * gs + ob[2]) * its
                var r = m[0] * ml + m[1] * mm + m[2] * ms
                var g = m[3] * ml + m[4] * mm + m[5] * ms
                var bb = m[6] * ml + m[7] * mm + m[8] * ms
                if (grey) {
                    val l = 0.2126f * r + 0.7152f * g + 0.0722f * bb
                    r = l; g = l; bb = l
                }
                // 出力は常に 8 ビット sRGB なので、画像の色空間が線形でも sRGB の伝達関数をかける
                x[i] = srgbEncodeFast(r)
                y[i] = srgbEncodeFast(g)
                b[i] = srgbEncodeFast(bb)
            }
        }
    }

    fun ycbcrToRgb(img: FrameImage) {
        val cb = img.color[0]
        val yy = img.color[1]
        val cr = img.color[2]
        for (i in cb.indices) {
            val c1 = cb[i]
            val y = yy[i] + 128f / 255f
            val c2 = cr[i]
            cb[i] = c2 * 1.402f + y
            yy[i] = c1 * (-0.114f * 1.772f / 0.587f) + (c2 * (-0.299f * 1.402f / 0.587f) + y)
            cr[i] = c1 * 1.772f + y
        }
    }

    // ---- XYB 以外の画像を sRGB に変換する ----

    private fun toLinear(v: Float, enc: ColourEncoding): Float {
        val s = if (v < 0) -1f else 1f
        val a = kotlin.math.abs(v).toDouble()
        val r = when (enc.transferFunction) {
            TransferFunction.LINEAR -> a
            TransferFunction.SRGB -> if (a <= 0.04045) a / 12.92 else Math.pow((a + 0.055) / 1.055, 2.4)
            TransferFunction.BT709 -> if (a < 0.081) a / 4.5 else Math.pow((a + 0.099) / 1.099, 1.0 / 0.45)
            TransferFunction.GAMMA -> Math.pow(a, 1.0 / enc.gamma)
            TransferFunction.DCI -> Math.pow(a, 2.6)
            TransferFunction.PQ -> {
                val m1 = 2610.0 / 16384
                val m2 = 2523.0 / 4096 * 128
                val c1 = 3424.0 / 4096
                val c2 = 2413.0 / 4096 * 32
                val c3 = 2392.0 / 4096 * 32
                val p = Math.pow(a, 1 / m2)
                // 10000 nits を 1 とし、SDR の白 (203 nits) を 1 に合わせる
                Math.pow(maxOf(p - c1, 0.0) / (c2 - c3 * p), 1 / m1) * 10000.0 / 203.0
            }
            TransferFunction.HLG -> {
                val aa = 0.17883277
                val bb = 1 - 4 * aa
                val cc = 0.5 - aa * Math.log(4 * aa)
                if (a <= 0.5) a * a / 3 else (Math.exp((a - cc) / aa) + bb) / 12
            }
            else -> if (a <= 0.04045) a / 12.92 else Math.pow((a + 0.055) / 1.055, 2.4)
        }
        return (s * r).toFloat()
    }

    private fun primariesToXyz(p: FloatArray, white: FloatArray): DoubleArray {
        val xr = p[0].toDouble(); val yr = p[1].toDouble()
        val xg = p[2].toDouble(); val yg = p[3].toDouble()
        val xb = p[4].toDouble(); val yb = p[5].toDouble()
        val m = doubleArrayOf(
            xr / yr, xg / yg, xb / yb,
            1.0, 1.0, 1.0,
            (1 - xr - yr) / yr, (1 - xg - yg) / yg, (1 - xb - yb) / yb,
        )
        val wx = white[0] / white[1].toDouble()
        val wz = (1 - white[0] - white[1]) / white[1].toDouble()
        val inv = invert3(m)
        val s = mul3v(inv, doubleArrayOf(wx, 1.0, wz))
        return doubleArrayOf(
            m[0] * s[0], m[1] * s[1], m[2] * s[2],
            m[3] * s[0], m[4] * s[1], m[5] * s[2],
            m[6] * s[0], m[7] * s[1], m[8] * s[2],
        )
    }

    private fun bradford(srcWhite: FloatArray, dstWhite: FloatArray): DoubleArray {
        val b = doubleArrayOf(0.8951, 0.2664, -0.1614, -0.7502, 1.7135, 0.0367, 0.0389, -0.0685, 1.0296)
        fun xyz(w: FloatArray) = doubleArrayOf(w[0] / w[1].toDouble(), 1.0, (1 - w[0] - w[1]) / w[1].toDouble())
        val s = mul3v(b, xyz(srcWhite))
        val d = mul3v(b, xyz(dstWhite))
        val diag = doubleArrayOf(d[0] / s[0], 0.0, 0.0, 0.0, d[1] / s[1], 0.0, 0.0, 0.0, d[2] / s[2])
        return mul3(invert3(b), mul3(diag, b))
    }

    private fun mul3(a: DoubleArray, b: DoubleArray) = DoubleArray(9) { i ->
        val r = i / 3
        val c = i % 3
        a[r * 3] * b[c] + a[r * 3 + 1] * b[3 + c] + a[r * 3 + 2] * b[6 + c]
    }

    private fun mul3v(a: DoubleArray, v: DoubleArray) = DoubleArray(3) { r -> a[r * 3] * v[0] + a[r * 3 + 1] * v[1] + a[r * 3 + 2] * v[2] }

    private fun invert3(m: DoubleArray): DoubleArray {
        val a = m[0]; val b = m[1]; val c = m[2]
        val d = m[3]; val e = m[4]; val f = m[5]
        val g = m[6]; val h = m[7]; val i = m[8]
        val det = a * (e * i - f * h) - b * (d * i - f * g) + c * (d * h - e * g)
        return doubleArrayOf(
            (e * i - f * h) / det, (c * h - b * i) / det, (b * f - c * e) / det,
            (f * g - d * i) / det, (a * i - c * g) / det, (c * d - a * f) / det,
            (d * h - e * g) / det, (b * g - a * h) / det, (a * e - b * d) / det,
        )
    }

    /** 列挙型の色空間で表された RGB / グレーを sRGB に変換する (sRGB ならそのまま) */
    fun toSrgb(img: FrameImage, enc: ColourEncoding) {
        if (enc.wantIcc) return
        val needTf = enc.transferFunction != TransferFunction.SRGB
        val needPrim = enc.colourSpace == ColourSpace.RGB && (!enc.isSrgbPrimaries || !enc.isD65)
        if (!needTf && !needPrim) return
        val matrix: DoubleArray? = if (needPrim) {
            val toXyz = primariesToXyz(enc.primaries, enc.whitePoint)
            val adapt = bradford(enc.whitePoint, ColourEncoding.D65)
            val srgbToXyz = primariesToXyz(ColourEncoding.SRGB_PRIMARIES, ColourEncoding.D65)
            mul3(invert3(srgbToXyz), mul3(adapt, toXyz))
        } else null
        val n = img.color[0].size
        if (img.color.size == 3) {
            val r = img.color[0]
            val g = img.color[1]
            val b = img.color[2]
            for (i in 0 until n) {
                var lr = toLinear(r[i], enc)
                var lg = toLinear(g[i], enc)
                var lb = toLinear(b[i], enc)
                if (matrix != null) {
                    val nr = (matrix[0] * lr + matrix[1] * lg + matrix[2] * lb).toFloat()
                    val ng = (matrix[3] * lr + matrix[4] * lg + matrix[5] * lb).toFloat()
                    val nb = (matrix[6] * lr + matrix[7] * lg + matrix[8] * lb).toFloat()
                    lr = nr; lg = ng; lb = nb
                }
                r[i] = srgbEncodeFast(lr)
                g[i] = srgbEncodeFast(lg)
                b[i] = srgbEncodeFast(lb)
            }
        } else {
            val gch = img.color[0]
            for (i in 0 until n) gch[i] = srgbEncodeFast(toLinear(gch[i], enc))
        }
    }
}
