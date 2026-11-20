package io.github.hamtit.jxl

/** Modular チャンネル 1 枚分のサンプルを MA 木に従って復号する */
internal object ChannelDecoder {

    fun predict(p: Int, w: Int, n: Int, nw: Int, ne: Int, nn: Int, ww: Int, nee: Int, wp: WeightedPredictor?): Int =
        when (p) {
            Predictor.ZERO -> 0
            Predictor.WEST -> w
            Predictor.NORTH -> n
            Predictor.AVG_W_N -> ((w.toLong() + n) / 2).toInt()
            Predictor.SELECT -> {
                val dn = kotlin.math.abs(n.toLong() - nw)
                val dw = kotlin.math.abs(w.toLong() - nw)
                if (dn < dw) w else n
            }
            Predictor.GRADIENT -> gradient(w, n, nw)
            Predictor.WEIGHTED -> ((wp!!.prediction + 3) shr 3).toInt()
            Predictor.NORTH_EAST -> ne
            Predictor.NORTH_WEST -> nw
            Predictor.WEST_WEST -> ww
            Predictor.AVG_W_NW -> ((w.toLong() + nw) / 2).toInt()
            Predictor.AVG_N_NW -> ((n.toLong() + nw) / 2).toInt()
            Predictor.AVG_N_NE -> ((n.toLong() + ne) / 2).toInt()
            else -> ((6L * n - 2L * nn + 7L * w + ww + nee + 3L * ne + 8) / 16).toInt()
        }

    fun gradient(w: Int, n: Int, nw: Int): Int {
        val g = w.toLong() + n - nw
        val lo = minOf(w, n).toLong()
        val hi = maxOf(w, n).toLong()
        return (if (g < lo) lo else if (g > hi) hi else g).toInt()
    }

    fun decode(
        br: BitReader,
        reader: EntropyReader,
        tree: MaTree,
        ch: Channel,
        prev: List<Channel>,
        wpParams: WpParams,
        distMultiplier: Int,
    ) {
        val width = ch.width
        val height = ch.height
        val out = ch.data
        val prop = tree.property

        var usesWp = false
        var maxPrev = 0
        for (i in 0 until tree.size) {
            val p = prop[i]
            if (p == 15) usesWp = true
            if (p >= 16) maxPrev = maxOf(maxPrev, (p - 16) / 4 + 1)
            if (p < 0 && tree.predictor[i] == Predictor.WEIGHTED) usesWp = true
        }

        if (tree.size == 1 && !usesWp) {
            decodeSingleLeaf(br, reader, tree, out, width, height, distMultiplier)
            return
        }

        val wp = if (usesWp) WeightedPredictor(width, wpParams) else null
        val prevData = Array(minOf(maxPrev, prev.size)) { prev[it].data }
        val code = tree.code
        val treeSize = tree.size
        val split = tree.splitValue
        val left = tree.left
        val right = tree.right
        val predictor = tree.predictor
        val offset = tree.offset
        val multiplier = tree.multiplier
        val cluster = tree.cluster
        val lz77 = code.lz77Enabled

        for (y in 0 until height) {
            val row = y * width
            var prevGrad = 0
            for (x in 0 until width) {
                val i = row + x
                val w: Int
                val n: Int
                val nw: Int
                if (y == 0) {
                    w = if (x > 0) out[i - 1] else 0
                    n = w
                    nw = w
                } else {
                    n = out[i - width]
                    w = if (x > 0) out[i - 1] else n
                    nw = if (x > 0) out[i - width - 1] else n
                }
                val ne = if (y > 0 && x + 1 < width) out[i - width + 1] else n
                val nn = if (y > 1) out[i - 2 * width] else n
                val ww = if (x > 1) out[i - 2] else w
                wp?.predict(n, nw, ne, w, nn)

                // 木をたどる
                var node = 0
                while (prop[node] >= 0) {
                    val v: Int = when (val p = prop[node]) {
                        2 -> y
                        3 -> x
                        4 -> kotlin.math.abs(n)
                        5 -> kotlin.math.abs(w)
                        6 -> n
                        7 -> w
                        8 -> w - prevGrad
                        9 -> w + n - nw
                        10 -> w - nw
                        11 -> nw - n
                        12 -> n - ne
                        13 -> n - nn
                        14 -> w - ww
                        15 -> wp!!.maxError
                        else -> prevProperty(p - 16, prevData, x, y, width, i)
                    }
                    node = if (v > split[node]) left[node] else right[node]
                }

                val token = if (lz77) reader.readClustered(br, cluster[node], distMultiplier)
                else reader.readClustered(br, cluster[node], 0)
                val diff = unpackSigned(token) * multiplier[node] + offset[node]
                val pr = predictor[node]
                val pred = when (pr) {
                    Predictor.ZERO -> 0
                    Predictor.WEST -> w
                    Predictor.NORTH -> n
                    Predictor.GRADIENT -> gradient(w, n, nw)
                    Predictor.WEIGHTED -> ((wp!!.prediction + 3) shr 3).toInt()
                    else -> {
                        val nee = if (y > 0 && x + 2 < width) out[i - width + 2] else ne
                        predict(pr, w, n, nw, ne, nn, ww, nee, wp)
                    }
                }
                val value = diff + pred
                out[i] = value
                wp?.record(value)
                prevGrad = w + n - nw
            }
        }
        if (treeSize < 0) error("unreachable")
    }

    private fun prevProperty(pe: Int, prevData: Array<IntArray>, x: Int, y: Int, width: Int, i: Int): Int {
        val idx = pe / 4
        if (idx >= prevData.size) return 0
        val d = prevData[idx]
        val c = d[i]
        val k = pe % 4
        if (k == 0) return kotlin.math.abs(c)
        if (k == 1) return c
        val g = if (x == 0 && y == 0) 0
        else if (x == 0) d[i - width]
        else if (y == 0) d[i - 1]
        else gradient(d[i - 1], d[i - width], d[i - width - 1])
        return if (k == 2) kotlin.math.abs(c - g) else c - g
    }

    /** 分岐のない木 (葉 1 枚) の高速版 */
    private fun decodeSingleLeaf(
        br: BitReader,
        reader: EntropyReader,
        tree: MaTree,
        out: IntArray,
        width: Int,
        height: Int,
        distMultiplier: Int,
    ) {
        val pr = tree.predictor[0]
        val off = tree.offset[0]
        val mul = tree.multiplier[0]
        val cl = tree.cluster[0]
        val code = tree.code
        val single = code.singleToken(cl)
        if (single >= 0 && pr == Predictor.ZERO) {
            java.util.Arrays.fill(out, 0, width * height, unpackSigned(single) * mul + off)
            return
        }
        val dm = if (code.lz77Enabled) distMultiplier else 0
        for (y in 0 until height) {
            val row = y * width
            for (x in 0 until width) {
                val i = row + x
                val diff = unpackSigned(reader.readClustered(br, cl, dm)) * mul + off
                val pred: Int = when (pr) {
                    Predictor.ZERO -> 0
                    Predictor.WEST -> if (x > 0) out[i - 1] else if (y > 0) out[i - width] else 0
                    Predictor.NORTH -> if (y > 0) out[i - width] else if (x > 0) out[i - 1] else 0
                    Predictor.GRADIENT -> {
                        if (y == 0) {
                            if (x > 0) out[i - 1] else 0
                        } else if (x == 0) {
                            out[i - width]
                        } else {
                            gradient(out[i - 1], out[i - width], out[i - width - 1])
                        }
                    }
                    else -> {
                        val w: Int
                        val n: Int
                        val nw: Int
                        if (y == 0) {
                            w = if (x > 0) out[i - 1] else 0
                            n = w
                            nw = w
                        } else {
                            n = out[i - width]
                            w = if (x > 0) out[i - 1] else n
                            nw = if (x > 0) out[i - width - 1] else n
                        }
                        val ne = if (y > 0 && x + 1 < width) out[i - width + 1] else n
                        val nn = if (y > 1) out[i - 2 * width] else n
                        val ww = if (x > 1) out[i - 2] else w
                        val nee = if (y > 0 && x + 2 < width) out[i - width + 2] else ne
                        predict(pr, w, n, nw, ne, nn, ww, nee, null)
                    }
                }
                out[i] = diff + pred
            }
        }
    }
}
