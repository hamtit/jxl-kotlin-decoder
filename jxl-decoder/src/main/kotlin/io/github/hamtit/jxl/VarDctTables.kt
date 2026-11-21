package io.github.hamtit.jxl

/** 変換の種類 (DctSelect) ごとの定数表 */
internal object Tt {
    const val DCT8 = 0
    const val HORNUSS = 1
    const val DCT2 = 2
    const val DCT4 = 3
    const val DCT16 = 4
    const val DCT32 = 5
    const val DCT16X8 = 6
    const val DCT8X16 = 7
    const val DCT32X8 = 8
    const val DCT8X32 = 9
    const val DCT32X16 = 10
    const val DCT16X32 = 11
    const val DCT4X8 = 12
    const val DCT8X4 = 13
    const val AFV0 = 14
    const val AFV3 = 17
    const val COUNT = 27

    /** 8x8 ブロック単位の幅 */
    val BW = intArrayOf(1, 1, 1, 1, 2, 4, 1, 2, 1, 4, 2, 4, 1, 1, 1, 1, 1, 1, 8, 4, 8, 16, 8, 16, 32, 16, 32)
    /** 8x8 ブロック単位の高さ */
    val BH = intArrayOf(1, 1, 1, 1, 2, 4, 2, 1, 4, 1, 4, 2, 1, 1, 1, 1, 1, 1, 8, 8, 4, 16, 16, 8, 32, 32, 16)

    /** 量子化行列の番号 */
    val MATRIX = intArrayOf(0, 1, 2, 3, 4, 5, 6, 6, 7, 7, 8, 8, 9, 9, 10, 10, 10, 10, 11, 12, 12, 13, 14, 14, 15, 16, 16)

    /** 係数順序の番号 */
    val ORDER = intArrayOf(0, 1, 1, 1, 2, 3, 4, 4, 5, 5, 6, 6, 1, 1, 1, 1, 1, 1, 7, 8, 8, 9, 10, 10, 11, 12, 12)

    /** 量子化行列の (幅, 高さ) */
    val MATRIX_W = intArrayOf(8, 8, 8, 8, 16, 32, 16, 32, 32, 8, 8, 64, 64, 128, 128, 256, 256)
    val MATRIX_H = intArrayOf(8, 8, 8, 8, 16, 32, 8, 8, 16, 8, 8, 64, 32, 128, 64, 256, 128)

    fun isSmall(t: Int): Boolean = t == HORNUSS || t == DCT2 || t == DCT4 || t == DCT4X8 || t == DCT8X4 || (t in AFV0..AFV3)

    /** 係数を縦横入れ替えて置くか */
    fun needTranspose(t: Int): Boolean = !isSmall(t) && BH[t] >= BW[t]

    /** 係数順序の (幅, 高さ)。幅 >= 高さ */
    val ORDER_W = intArrayOf(8, 8, 16, 32, 16, 32, 32, 64, 64, 128, 128, 256, 256)
    val ORDER_H = intArrayOf(8, 8, 16, 32, 8, 8, 16, 64, 32, 128, 64, 256, 128)

    private val naturalOrders = java.util.concurrent.atomic.AtomicReferenceArray<IntArray>(13)

    /** 自然順序 (x + y * ORDER_W) の並び */
    fun naturalOrder(orderId: Int): IntArray {
        naturalOrders.get(orderId)?.let { return it }
        run {
            val bw = ORDER_W[orderId]
            val bh = ORDER_H[orderId]
            val out = IntArray(bw * bh)
            val yScale = bw / bh
            val lbw = bw / 8
            val lbh = bh / 8
            var idx = 0
            while (idx < lbw * lbh) {
                val x = idx % lbw
                val y = idx / lbw
                out[idx] = x + y * bw
                idx++
            }
            for (dist in 1 until 2 * bw) {
                val margin = maxOf(0, dist - bw)
                for (order in margin until dist - margin) {
                    val x: Int
                    val y: Int
                    if (dist % 2 == 1) {
                        x = order
                        y = dist - 1 - order
                    } else {
                        x = dist - 1 - order
                        y = order
                    }
                    if (x < lbw && y < lbw) continue
                    if (y % yScale != 0) continue
                    out[idx++] = x + (y / yScale) * bw
                }
            }
            naturalOrders.compareAndSet(orderId, null, out)
            return naturalOrders.get(orderId)
        }
    }
}

/** 逆量子化行列 */
internal class QuantMatrices(
    /** [matrix][channel] = 行列 (MATRIX_W x MATRIX_H、行優先) */
    private val matrices: Array<Array<FloatArray>>,
) {
    private val transposed = java.util.concurrent.atomic.AtomicReferenceArray<FloatArray>(17 * 3)

    private fun transposedMatrix(m: Int, c: Int): FloatArray {
        transposed.get(m * 3 + c)?.let { return it }
        val src = matrices[m][c]
        val w = Tt.MATRIX_W[m]
        val h = Tt.MATRIX_H[m]
        val out = FloatArray(src.size)
        // 出力は幅 h、高さ w
        for (y in 0 until w) for (x in 0 until h) out[y * h + x] = src[x * w + y]
        transposed.compareAndSet(m * 3 + c, null, out)
        return transposed.get(m * 3 + c)
    }

    /** 変換 [t] の係数配置 (幅 BW*8) に合わせた行列 */
    fun get(t: Int, c: Int): FloatArray {
        val m = Tt.MATRIX[t]
        return if (Tt.needTranspose(t)) transposedMatrix(m, c) else matrices[m][c]
    }

    companion object {
        private val SEQ_A = floatArrayOf(-1.025f, -0.78f, -0.65012f, -0.19041574f, -0.20819396f, -0.421064f, -0.32733846f)
        private val SEQ_B = floatArrayOf(-0.30419582f, -0.36330363f, -0.3566038f, -0.34430745f, -0.33699593f, -0.30180866f, -0.27321684f)
        private val SEQ_C = floatArrayOf(-1.2f, -1.2f, -0.8f, -0.7f, -0.7f, -0.4f, -0.5f)
        private val DCT4X8_PARAMS = arrayOf(
            floatArrayOf(2198.0505f, -0.96269625f, -0.7619425f, -0.65511405f),
            floatArrayOf(764.36554f, -0.926302f, -0.967523f, -0.2784529f),
            floatArrayOf(527.10754f, -1.4594386f, -1.4500821f, -1.5843723f),
        )
        private val DCT4_PARAMS = arrayOf(
            floatArrayOf(2200.0f, 0.0f, 0.0f, 0.0f),
            floatArrayOf(392.0f, 0.0f, 0.0f, 0.0f),
            floatArrayOf(112.0f, -0.25f, -0.25f, -0.5f),
        )

        private fun seq(a: Float, b: Float, c: Float) = arrayOf(
            floatArrayOf(a) + SEQ_A, floatArrayOf(b) + SEQ_B, floatArrayOf(c) + SEQ_C,
        )

        /** 既定のパラメータ */
        private fun defaultParams(m: Int): Params = when (m) {
            0 -> Params.Dct(arrayOf(
                floatArrayOf(3150.0f, 0.0f, -0.4f, -0.4f, -0.4f, -2.0f),
                floatArrayOf(560.0f, 0.0f, -0.3f, -0.3f, -0.3f, -0.3f),
                floatArrayOf(512.0f, -2.0f, -1.0f, 0.0f, -1.0f, -2.0f),
            ))
            1 -> Params.Hornuss(arrayOf(
                floatArrayOf(280.0f, 3160.0f, 3160.0f),
                floatArrayOf(60.0f, 864.0f, 864.0f),
                floatArrayOf(18.0f, 200.0f, 200.0f),
            ))
            2 -> Params.Dct2(arrayOf(
                floatArrayOf(3840.0f, 2560.0f, 1280.0f, 640.0f, 480.0f, 300.0f),
                floatArrayOf(960.0f, 640.0f, 320.0f, 180.0f, 140.0f, 120.0f),
                floatArrayOf(640.0f, 320.0f, 128.0f, 64.0f, 32.0f, 16.0f),
            ))
            3 -> Params.Dct4(arrayOf(floatArrayOf(1f, 1f), floatArrayOf(1f, 1f), floatArrayOf(1f, 1f)), DCT4_PARAMS)
            4 -> Params.Dct(arrayOf(
                floatArrayOf(8996.873f, -1.3000778f, -0.4942453f, -0.43909377f, -0.6350102f, -0.9017726f, -1.6162099f),
                floatArrayOf(3191.4836f, -0.67424583f, -0.80745816f, -0.4492584f, -0.3586544f, -0.3132239f, -0.37615025f),
                floatArrayOf(1157.504f, -2.0531423f, -1.4f, -0.5068713f, -0.4270873f, -1.4856834f, -4.920914f),
            ))
            5 -> Params.Dct(arrayOf(
                floatArrayOf(15718.408f, -1.025f, -0.98f, -0.9012f, -0.4f, -0.48819396f, -0.421064f, -0.27f),
                floatArrayOf(7305.7637f, -0.8041958f, -0.76330364f, -0.5566038f, -0.49785304f, -0.43699592f, -0.40180868f, -0.27321684f),
                floatArrayOf(3803.5317f, -3.0607336f, -2.041327f, -2.023565f, -0.54953897f, -0.4f, -0.4f, -0.3f),
            ))
            6 -> Params.Dct(arrayOf(
                floatArrayOf(7240.7734f, -0.7f, -0.7f, -0.2f, -0.2f, -0.2f, -0.5f),
                floatArrayOf(1448.1547f, -0.5f, -0.5f, -0.5f, -0.2f, -0.2f, -0.2f),
                floatArrayOf(506.85413f, -1.4f, -0.2f, -0.5f, -0.5f, -1.5f, -3.6f),
            ))
            7 -> Params.Dct(arrayOf(
                floatArrayOf(16283.249f, -1.7812846f, -1.6309059f, -1.0382179f, -0.85f, -0.7f, -0.9f, -1.2360638f),
                floatArrayOf(5089.1577f, -0.3200494f, -0.3536285f, -0.3034f, -0.61f, -0.5f, -0.5f, -0.6f),
                floatArrayOf(3397.7761f, -0.32132736f, -0.3450762f, -0.7034f, -0.9f, -1.0f, -1.0f, -1.1754606f),
            ))
            8 -> Params.Dct(arrayOf(
                floatArrayOf(13844.971f, -0.971138f, -0.658f, -0.42026f, -0.22712f, -0.2206f, -0.226f, -0.6f),
                floatArrayOf(4798.964f, -0.6112531f, -0.8377079f, -0.7901486f, -0.26927274f, -0.38272768f, -0.22924222f, -0.20719099f),
                floatArrayOf(1807.2369f, -1.2f, -1.2f, -0.7f, -0.7f, -0.7f, -0.4f, -0.5f),
            ))
            9 -> Params.Dct4x8(arrayOf(floatArrayOf(1f), floatArrayOf(1f), floatArrayOf(1f)), DCT4X8_PARAMS)
            10 -> Params.Afv(
                arrayOf(
                    floatArrayOf(3072.0f, 3072.0f, 256.0f, 256.0f, 256.0f, 414.0f, 0.0f, 0.0f, 0.0f),
                    floatArrayOf(1024.0f, 1024.0f, 50.0f, 50.0f, 50.0f, 58.0f, 0.0f, 0.0f, 0.0f),
                    floatArrayOf(384.0f, 384.0f, 12.0f, 12.0f, 12.0f, 22.0f, -0.25f, -0.25f, -0.25f),
                ),
                DCT4X8_PARAMS, DCT4_PARAMS,
            )
            11 -> Params.Dct(seq(23966.166f, 8380.191f, 4493.024f))
            12 -> Params.Dct(seq(15358.898f, 5597.3604f, 2919.9617f))
            13 -> Params.Dct(seq(47932.332f, 16760.383f, 8986.048f))
            14 -> Params.Dct(seq(30717.797f, 11194.721f, 5839.9233f))
            15 -> Params.Dct(seq(95864.664f, 33520.766f, 17972.096f))
            else -> Params.Dct(seq(61435.594f, 24209.441f, 12979.847f))
        }

        sealed class Params {
            class Dct(val p: Array<FloatArray>) : Params()
            class Hornuss(val p: Array<FloatArray>) : Params()
            class Dct2(val p: Array<FloatArray>) : Params()
            class Dct4(val p: Array<FloatArray>, val dct: Array<FloatArray>) : Params()
            class Dct4x8(val p: Array<FloatArray>, val dct: Array<FloatArray>) : Params()
            class Afv(val p: Array<FloatArray>, val dct: Array<FloatArray>, val dct4: Array<FloatArray>) : Params()
            class Raw(val matrix: Array<FloatArray>) : Params()
        }

        private fun interpolate(pos: Float, max: Float, bands: FloatArray): Float {
            if (bands.size == 1) return bands[0]
            val scaled = pos * (bands.size - 1) / max
            val idx = scaled.toInt()
            val frac = scaled - idx
            val a = bands[idx]
            val b = bands[idx + 1]
            return a * Math.pow((b / a).toDouble(), frac.toDouble()).toFloat()
        }

        private fun mult(x: Float): Float = if (x > 0f) 1f + x else 1f / (1f - x)

        private fun dctWeights(params: FloatArray, width: Int, height: Int): FloatArray {
            val bands = FloatArray(params.size)
            var last = params[0]
            bands[0] = last
            for (i in 1 until params.size) {
                val b = last * mult(params[i])
                if (b <= 0f) fail("invalid dequant matrix")
                bands[i] = b
                last = b
            }
            val out = FloatArray(width * height)
            for (y in 0 until height) for (x in 0 until width) {
                val dx = x.toFloat() / (width - 1)
                val dy = y.toFloat() / (height - 1)
                val dist = Math.sqrt((dx * dx + dy * dy).toDouble()).toFloat()
                out[y * width + x] = interpolate(dist, 1.4142135f + 1e-6f, bands)
            }
            return out
        }

        private val AFV_FREQS = floatArrayOf(
            0.0f, 0.0f, 0.8517779f, 5.3777843f, 0.0f, 0.0f, 4.734748f, 5.4492455f, 1.659827f, 4.0f,
            7.275749f, 10.423227f, 2.6629324f, 7.6306577f, 8.962389f, 12.971662f,
        )

        private fun toMatrix(m: Int, params: Params): Array<FloatArray> {
            val w = Tt.MATRIX_W[m]
            val h = Tt.MATRIX_H[m]
            if (params is Params.Raw) return params.matrix
            val weights: Array<FloatArray> = when (params) {
                is Params.Dct -> Array(3) { dctWeights(params.p[it], w, h) }
                is Params.Hornuss -> Array(3) { c ->
                    val p = params.p[c]
                    val r = FloatArray(64) { p[0] }
                    r[0] = 1f
                    r[1] = p[1]
                    r[8] = p[1]
                    r[9] = p[2]
                    r
                }
                is Params.Dct2 -> Array(3) { c ->
                    val p = params.p[c]
                    val r = FloatArray(64)
                    r[0] = 1f
                    for (idx in 0 until 6) {
                        val dim = 1 shl (idx / 2)
                        val v = p[idx]
                        if (idx % 2 == 0) {
                            for (y in 0 until dim) for (x in dim until dim * 2) {
                                r[y * 8 + x] = v
                                r[x * 8 + y] = v
                            }
                        } else {
                            for (y in dim until dim * 2) for (x in dim until dim * 2) r[y * 8 + x] = v
                        }
                    }
                    r
                }
                is Params.Dct4 -> Array(3) { c ->
                    val mat = dctWeights(params.dct[c], 4, 4)
                    val r = FloatArray(64)
                    for (y in 0 until 4) for (x in 0 until 4) {
                        val v = mat[y * 4 + x]
                        r[y * 16 + x * 2] = v
                        r[y * 16 + x * 2 + 1] = v
                        r[(y * 2 + 1) * 8 + x * 2] = v
                        r[(y * 2 + 1) * 8 + x * 2 + 1] = v
                    }
                    r[1] /= params.p[c][0]
                    r[8] /= params.p[c][0]
                    r[9] /= params.p[c][1]
                    r
                }
                is Params.Dct4x8 -> Array(3) { c ->
                    val mat = dctWeights(params.dct[c], 8, 4)
                    val r = FloatArray(64)
                    for (y in 0 until 4) for (x in 0 until 8) {
                        r[(2 * y) * 8 + x] = mat[y * 8 + x]
                        r[(2 * y + 1) * 8 + x] = mat[y * 8 + x]
                    }
                    r[8] /= params.p[c][0]
                    r
                }
                is Params.Afv -> Array(3) { c ->
                    val p = params.p[c]
                    val w48 = dctWeights(params.dct[c], 8, 4)
                    val w44 = dctWeights(params.dct4[c], 4, 4)
                    val bands = FloatArray(4)
                    bands[0] = p[5]
                    for (i in 1 until 4) bands[i] = bands[i - 1] * mult(p[5 + i])
                    val r = FloatArray(64)
                    val lo = AFV_FREQS[2]
                    val hi = AFV_FREQS[15]
                    for (y in 0 until 4) for (x in 0 until 4) {
                        r[16 * y + 2 * x] = when {
                            x == 0 && y == 0 -> 1f
                            x == 0 && y == 1 -> p[2]
                            x == 1 && y == 0 -> p[3]
                            x == 1 && y == 1 -> p[4]
                            else -> interpolate(AFV_FREQS[y * 4 + x] - lo, hi - lo + 1e-6f, bands)
                        }
                    }
                    for (y in 0 until 4) {
                        for (x in 0 until 8) {
                            r[16 * y + 8 + x] = if (y == 0 && x == 0) p[0] else w48[y * 8 + x]
                        }
                        for (x in 0 until 4) {
                            r[16 * y + 2 * x + 1] = if (y == 0 && x == 0) p[1] else w44[y * 4 + x]
                        }
                    }
                    r
                }
            }
            for (wa in weights) {
                for (i in wa.indices) {
                    val v = 1f / wa[i]
                    if (!(v > 0f) || v >= 1e8f) fail("invalid dequant matrix")
                    wa[i] = v
                }
            }
            return weights
        }

        /** 既定の行列 (重いので一度だけ作る) */
        private val defaultMatrices = java.util.concurrent.atomic.AtomicReferenceArray<Array<FloatArray>>(17)

        private fun defaultMatrix(m: Int): Array<FloatArray> {
            defaultMatrices.get(m)?.let { return it }
            defaultMatrices.compareAndSet(m, null, toMatrix(m, defaultParams(m)))
            return defaultMatrices.get(m)
        }

        @Volatile
        private var defaultSet: QuantMatrices? = null

        fun parse(br: BitReader, numLfGroups: Int, globalTree: MaTree?, bitDepth: Int): QuantMatrices {
            val all = br.readBool()
            if (all) {
                defaultSet?.let { return it }
                return QuantMatrices(Array(17) { defaultMatrix(it) }).also { defaultSet = it }
            }
            val list = Array(17) { m ->
                val params = readParams(br, m, numLfGroups, globalTree, bitDepth)
                if (params == null) defaultMatrix(m) else toMatrix(m, params)
            }
            return QuantMatrices(list)
        }

        private fun readFixed(br: BitReader, n: Int) = Array(3) { FloatArray(n) { br.readF16() } }

        private fun readDctParams(br: BitReader): Array<FloatArray> {
            val num = br.readBits(4) + 1
            val p = Array(3) { FloatArray(num) { br.readF16() } }
            for (c in 0 until 3) p[c][0] *= 64f
            return p
        }

        private fun readParams(br: BitReader, m: Int, numLfGroups: Int, globalTree: MaTree?, bitDepth: Int): Params? {
            val mode = br.readBits(3)
            if (mode in 1..5 && m !in intArrayOf(0, 1, 2, 3, 9, 10)) fail("invalid dequant encoding mode")
            return when (mode) {
                0 -> null
                1 -> Params.Hornuss(readFixed(br, 3))
                2 -> Params.Dct2(readFixed(br, 6))
                3 -> Params.Dct4(readFixed(br, 2), readDctParams(br))
                4 -> Params.Dct4x8(readFixed(br, 1), readDctParams(br))
                5 -> {
                    val p = readFixed(br, 9)
                    for (c in 0 until 3) for (i in 0 until 6) p[c][i] *= 64f
                    Params.Afv(p, readDctParams(br), readDctParams(br))
                }
                6 -> Params.Dct(readDctParams(br))
                7 -> {
                    val w = Tt.MATRIX_W[m]
                    val h = Tt.MATRIX_H[m]
                    val den = br.readF16()
                    val chans = List(3) { Channel(w, h, 0, 0) }
                    val img = ModularImage.read(br, chans, globalTree, bitDepth)
                    val out = img.decodeAll(br, 1 + numLfGroups * 3 + m)
                    Params.Raw(Array(3) { c ->
                        val d = out[c].data
                        FloatArray(w * h) { d[it] * den }
                    })
                }
                else -> fail("invalid dequant encoding mode")
            }
        }
    }
}
