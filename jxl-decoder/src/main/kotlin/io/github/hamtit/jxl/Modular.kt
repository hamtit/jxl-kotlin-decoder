package io.github.hamtit.jxl

/** Modular 画像の 1 チャンネル (整数サンプル) */
internal class Channel(var width: Int, var height: Int, var hshift: Int, var vshift: Int) {
    var data: IntArray = IntArray(0)

    fun allocate() {
        val n = width.toLong() * height
        if (n > Int.MAX_VALUE) fail("channel too large")
        data = IntArray(n.toInt())
    }

    fun copyMeta(): Channel = Channel(width, height, hshift, vshift)
}

internal class WpParams(
    val p1: Int, val p2: Int, val p3a: Int, val p3b: Int, val p3c: Int, val p3d: Int, val p3e: Int,
    val w0: Int, val w1: Int, val w2: Int, val w3: Int,
) {
    companion object {
        val DEFAULT = WpParams(16, 10, 7, 7, 7, 0, 0, 13, 12, 12, 12)
        fun parse(br: BitReader): WpParams {
            if (br.readBool()) return DEFAULT
            return WpParams(
                br.readBits(5), br.readBits(5), br.readBits(5), br.readBits(5), br.readBits(5),
                br.readBits(5), br.readBits(5), br.readBits(4), br.readBits(4), br.readBits(4), br.readBits(4),
            )
        }
    }
}

internal object Predictor {
    const val ZERO = 0
    const val WEST = 1
    const val NORTH = 2
    const val AVG_W_N = 3
    const val SELECT = 4
    const val GRADIENT = 5
    const val WEIGHTED = 6
    const val NORTH_EAST = 7
    const val NORTH_WEST = 8
    const val WEST_WEST = 9
    const val AVG_W_NW = 10
    const val AVG_N_NW = 11
    const val AVG_N_NE = 12
    const val AVG_ALL = 13
}

/**
 * MA 木。ノード i は property[i] >= 0 なら分岐 (値 > splitValue なら left)、-1 なら葉。
 */
internal class MaTree(
    val property: IntArray,
    val splitValue: IntArray,
    val left: IntArray,
    val right: IntArray,
    val predictor: IntArray,
    val offset: IntArray,
    val multiplier: IntArray,
    val cluster: IntArray,
    val code: EntropyCode,
) {
    val size: Int get() = property.size

    companion object {
        fun parse(br: BitReader, nodeLimit: Int): MaTree {
            val treeCode = EntropyCode.parse(br, 6)
            val reader = treeCode.newReader()
            reader.begin(br)
            val prop = IntArrayList()
            val split = IntArrayList()
            val left = IntArrayList()
            val right = IntArrayList()
            val pred = IntArrayList()
            val off = IntArrayList()
            val mul = IntArrayList()
            val ctx = IntArrayList()
            var toDecode = 1
            var numCtx = 0
            while (toDecode > 0) {
                if (prop.size > nodeLimit) fail("MA tree too large")
                toDecode--
                val p = reader.read(br, 1)
                if (p > 0) {
                    val v = unpackSigned(reader.read(br, 0))
                    val l = prop.size + toDecode + 1
                    prop.add(p - 1); split.add(v); left.add(l); right.add(l + 1)
                    pred.add(0); off.add(0); mul.add(0); ctx.add(0)
                    toDecode += 2
                } else {
                    val predictor = reader.read(br, 2)
                    if (predictor > 13) fail("invalid predictor")
                    val offset = unpackSigned(reader.read(br, 3))
                    val mulLog = reader.read(br, 4)
                    if (mulLog > 30) fail("invalid MA tree")
                    val mulBits = reader.read(br, 5)
                    if (mulBits.toLong() > (1L shl (31 - mulLog)) - 2) fail("invalid MA tree")
                    prop.add(-1); split.add(0); left.add(0); right.add(0)
                    pred.add(predictor); off.add(offset); mul.add((mulBits + 1) shl mulLog); ctx.add(numCtx++)
                }
            }
            reader.finish()
            val code = EntropyCode.parse(br, numCtx)
            val clusters = IntArray(prop.size) { if (prop[it] < 0) code.clusters[ctx[it]] else 0 }
            return MaTree(
                prop.toArray(), split.toArray(), left.toArray(), right.toArray(), pred.toArray(), off.toArray(),
                mul.toArray(), clusters, code,
            )
        }
    }

    /** チャンネル番号とストリーム番号で決まる分岐を畳み込んだ木を作る */
    fun specialize(channel: Int, stream: Int): MaTree {
        val map = IntArray(size) { -1 }
        val order = IntArrayList()
        fun resolve(start: Int): Int {
            var n = start
            while (property[n] == 0 || property[n] == 1) {
                val v = if (property[n] == 0) channel else stream
                n = if (v > splitValue[n]) left[n] else right[n]
            }
            return n
        }
        val root = resolve(0)
        // 幅優先で再番号付け
        val queue = IntArrayList()
        queue.add(root)
        var qi = 0
        while (qi < queue.size) {
            val n = queue[qi++]
            if (map[n] >= 0) continue
            map[n] = order.size
            order.add(n)
            if (property[n] >= 0) {
                queue.add(resolve(left[n]))
                queue.add(resolve(right[n]))
            }
        }
        val m = order.size
        val p = IntArray(m)
        val s = IntArray(m)
        val l = IntArray(m)
        val r = IntArray(m)
        val pr = IntArray(m)
        val o = IntArray(m)
        val mu = IntArray(m)
        val c = IntArray(m)
        for (i in 0 until m) {
            val n = order[i]
            p[i] = property[n]
            s[i] = splitValue[n]
            if (property[n] >= 0) {
                l[i] = map[resolve(left[n])]
                r[i] = map[resolve(right[n])]
            }
            pr[i] = predictor[n]
            o[i] = offset[n]
            mu[i] = multiplier[n]
            c[i] = cluster[n]
        }
        return MaTree(p, s, l, r, pr, o, mu, c, code)
    }
}

internal class IntArrayList(capacity: Int = 16) {
    private var a = IntArray(capacity)
    var size = 0
        private set

    fun add(v: Int) {
        if (size == a.size) a = a.copyOf(a.size * 2)
        a[size++] = v
    }

    operator fun get(i: Int): Int = a[i]
    fun toArray(): IntArray = a.copyOf(size)
}

/** 自己補正予測器 (weighted predictor) */
internal class WeightedPredictor(private val width: Int, private val wp: WpParams) {
    private val trueErrRow = IntArray(width)
    private val subErrRow = IntArray(width * 4)
    private var x = 0
    private var y = 0
    private var teW = 0
    private var teNW = 0
    private var teN = 0
    private var teNE = 0
    private val seNWWW = IntArray(4)
    private val seNW = IntArray(4)
    private val seNE = IntArray(4)
    private val subpred = LongArray(4)
    private val weight = IntArray(4)
    var prediction = 0L
        private set
    var maxError = 0
        private set

    fun predict(n: Int, nw: Int, ne: Int, w: Int, nn: Int) {
        val tw = teW.toLong()
        val tnw = teNW.toLong()
        val tn = teN.toLong()
        val tne = teNE.toLong()
        val n3 = n.toLong() shl 3
        val nw3 = nw.toLong() shl 3
        val ne3 = ne.toLong() shl 3
        val w3 = w.toLong() shl 3
        val nn3 = nn.toLong() shl 3
        subpred[0] = w3 + ne3 - n3
        subpred[1] = n3 - (((tw + tn + tne) * wp.p1) shr 5)
        subpred[2] = w3 - (((tw + tn + tnw) * wp.p2) shr 5)
        subpred[3] = n3 - ((tnw * wp.p3a + tn * wp.p3b + tne * wp.p3c + (nn3 - n3) * wp.p3d + (nw3 - w3) * wp.p3e) shr 5)

        var sumW = 0L
        for (i in 0 until 4) {
            val errSum = (seNWWW[i].toLong() and 0xffffffffL) + (seNW[i].toLong() and 0xffffffffL) +
                (seNE[i].toLong() and 0xffffffffL)
            val e = errSum and 0xffffffffL
            val t = (e + 1) shr 5
            val shift = if (t == 0L) 0 else 63 - java.lang.Long.numberOfLeadingZeros(t)
            val maxW = when (i) { 0 -> wp.w0; 1 -> wp.w1; 2 -> wp.w2; else -> wp.w3 }
            val wv = 4 + ((maxW.toLong() * DIV_LOOKUP[((e shr shift) + 1).toInt()]) shr shift)
            weight[i] = wv.toInt()
            sumW += wv
        }
        val logWeight = 63 - java.lang.Long.numberOfLeadingZeros(sumW shr 4)
        var sumW2 = 0
        for (i in 0 until 4) {
            weight[i] = weight[i] ushr logWeight
            sumW2 += weight[i]
        }
        var s = (sumW2.toLong() shr 1) - 1
        for (i in 0 until 4) s += subpred[i] * weight[i]
        var pred = (s * DIV_LOOKUP[sumW2]) shr 24
        if (((tn xor tw) or (tn xor tnw)) <= 0) {
            val mn = minOf(n3, minOf(w3, ne3))
            val mx = maxOf(n3, maxOf(w3, ne3))
            pred = pred.coerceIn(mn, mx)
        }
        prediction = pred
        var me = teW
        if (kotlin.math.abs(teN.toLong()) > kotlin.math.abs(me.toLong())) me = teN
        if (kotlin.math.abs(teNW.toLong()) > kotlin.math.abs(me.toLong())) me = teNW
        if (kotlin.math.abs(teNE.toLong()) > kotlin.math.abs(me.toLong())) me = teNE
        maxError = me
    }

    fun record(sample: Int) {
        val s8 = sample.toLong() shl 3
        val trueErr = (prediction - s8).toInt()
        val base = x * 4
        for (i in 0 until 4) {
            subErrRow[base + i] = ((kotlin.math.abs(subpred[i] - s8) + 3) shr 3).toInt()
        }
        trueErrRow[x] = trueErr
        x++
        if (x >= width) {
            y++
            x = 0
            teW = 0
            teN = trueErrRow[0]
            teNW = teN
            for (i in 0 until 4) seNW[i] = subErrRow[i]
            for (i in 0 until 4) seNWWW[i] = seNW[i]
            if (width <= 1) {
                teNE = teN
                for (i in 0 until 4) seNE[i] = seNW[i]
            } else {
                teNE = trueErrRow[1]
                for (i in 0 until 4) seNE[i] = subErrRow[4 + i]
            }
        } else {
            teW = trueErr
            teNW = teN
            teN = teNE
            for (i in 0 until 4) {
                seNWWW[i] = seNW[i]
                seNW[i] = seNE[i] + subErrRow[base + i]
            }
            if (x + 1 >= width) {
                teNE = teN
                for (i in 0 until 4) seNE[i] = seNW[i]
            } else if (y != 0) {
                teNE = trueErrRow[x + 1]
                for (i in 0 until 4) seNE[i] = subErrRow[(x + 1) * 4 + i]
            }
        }
    }

    companion object {
        val DIV_LOOKUP = LongArray(65).also { for (i in 1..64) it[i] = ((1 shl 24) / i).toLong() }
    }
}

internal sealed class Transform {
    class Rct(val beginC: Int, val type: Int) : Transform()
    class Palette(val beginC: Int, val numC: Int, val nbColours: Int, val nbDeltas: Int, val dPred: Int) : Transform()
    class Squeeze(val params: MutableList<SqueezeParam>) : Transform()
}

internal class SqueezeParam(val horizontal: Boolean, val inPlace: Boolean, val beginC: Int, val numC: Int)

internal class ModularHeader(
    val useGlobalTree: Boolean,
    val wp: WpParams,
    val transforms: List<Transform>,
) {
    companion object {
        fun parse(br: BitReader): ModularHeader {
            val useGlobal = br.readBool()
            val wp = WpParams.parse(br)
            val nb = br.readU32(0, 0, 1, 0, 2, 4, 18, 8)
            if (nb > 512) fail("too many transforms")
            val list = ArrayList<Transform>(nb)
            for (i in 0 until nb) {
                when (br.readBits(2)) {
                    0 -> {
                        val begin = br.readU32(0, 3, 8, 6, 72, 10, 1096, 13)
                        val type = br.readU32(6, 0, 0, 2, 2, 4, 10, 6)
                        if (type >= 42) fail("invalid RCT type")
                        list.add(Transform.Rct(begin, type))
                    }
                    1 -> {
                        val begin = br.readU32(0, 3, 8, 6, 72, 10, 1096, 13)
                        val numC = br.readU32(1, 0, 3, 0, 4, 0, 1, 13)
                        val nbColours = br.readU32(0, 8, 256, 10, 1280, 12, 5376, 16)
                        val nbDeltas = br.readU32(0, 0, 1, 8, 257, 10, 1281, 16)
                        val dPred = br.readBits(4)
                        if (dPred > 13) fail("invalid palette predictor")
                        list.add(Transform.Palette(begin, numC, nbColours, nbDeltas, dPred))
                    }
                    2 -> {
                        val num = br.readU32(0, 0, 1, 4, 9, 6, 41, 8)
                        val sp = ArrayList<SqueezeParam>(num)
                        for (j in 0 until num) {
                            val h = br.readBool()
                            val inPlace = br.readBool()
                            val begin = br.readU32(0, 3, 8, 6, 72, 10, 1096, 13)
                            val numC = br.readU32(1, 0, 2, 0, 3, 0, 4, 4)
                            sp.add(SqueezeParam(h, inPlace, begin, numC))
                        }
                        list.add(Transform.Squeeze(sp))
                    }
                    else -> fail("invalid transform")
                }
            }
            return ModularHeader(useGlobal, wp, list)
        }
    }
}

/**
 * 1 つの Modular サブストリーム。チャンネル情報に変換を適用した後のチャンネル列を復号し、
 * [finish] で逆変換して元のチャンネルに戻す。
 */
internal class ModularImage(
    val header: ModularHeader,
    val tree: MaTree,
    /** 変換前のチャンネル */
    val original: List<Channel>,
    /** 変換後のチャンネル (復号対象) */
    val channels: MutableList<Channel>,
    val nbMetaChannels: Int,
    val bitDepth: Int,
) {
    companion object {
        /**
         * ヘッダと (必要なら) 木を読み、変換後のチャンネルを用意する。
         * [channelsIn] は幅・高さ・シフトだけを使う。
         */
        fun read(br: BitReader, channelsIn: List<Channel>, globalTree: MaTree?, bitDepth: Int): ModularImage {
            val header = ModularHeader.parse(br)
            val chans = channelsIn.map { it.copyMeta() }.toMutableList()
            var nbMeta = 0
            for (tr in header.transforms) {
                nbMeta = applyMeta(tr, chans, nbMeta)
            }
            if (chans.size > (1 shl 16)) fail("too many channels")
            val tree = if (header.useGlobalTree) {
                globalTree ?: fail("global tree not available")
            } else {
                var samples = 0L
                for (c in chans) samples += c.width.toLong() * c.height
                MaTree.parse(br, minOf(1024L + samples, 1L shl 20).toInt())
            }
            for (c in chans) c.allocate()
            return ModularImage(header, tree, channelsIn, chans, nbMeta, bitDepth)
        }

        private fun applyMeta(tr: Transform, chans: MutableList<Channel>, nbMetaIn: Int): Int {
            var nbMeta = nbMetaIn
            when (tr) {
                is Transform.Rct -> {
                    if (tr.beginC + 3 > chans.size) fail("invalid RCT")
                    val c0 = chans[tr.beginC]
                    for (i in 1 until 3) {
                        val c = chans[tr.beginC + i]
                        if (c.width != c0.width || c.height != c0.height) fail("invalid RCT")
                    }
                }
                is Transform.Palette -> {
                    val end = tr.beginC + tr.numC
                    if (end > chans.size || tr.numC == 0) fail("invalid palette")
                    if (tr.beginC < nbMeta) {
                        if (end > nbMeta) fail("invalid palette")
                        nbMeta = nbMeta + 2 - tr.numC
                    } else {
                        nbMeta += 1
                    }
                    val c0 = chans[tr.beginC]
                    for (i in tr.beginC + 1 until end) {
                        val c = chans[i]
                        if (c.width != c0.width || c.height != c0.height) fail("invalid palette")
                    }
                    for (i in end - 1 downTo tr.beginC + 1) chans.removeAt(i)
                    chans.add(0, Channel(tr.nbColours, tr.numC, -1, -1))
                }
                is Transform.Squeeze -> {
                    if (tr.params.isEmpty()) defaultSqueeze(tr.params, chans, nbMeta)
                    for (sp in tr.params) {
                        val begin = sp.beginC
                        val end = begin + sp.numC
                        if (end > chans.size) fail("invalid squeeze")
                        if (begin < nbMeta) {
                            if (!sp.inPlace || end > nbMeta) fail("invalid squeeze")
                            nbMeta += sp.numC
                        }
                        val residuals = ArrayList<Channel>()
                        for (i in begin until end) {
                            val ch = chans[i]
                            if (ch.width == 0 || ch.height == 0) fail("invalid squeeze")
                            if (ch.hshift > 30 || ch.vshift > 30) fail("invalid squeeze")
                            val r = ch.copyMeta()
                            if (sp.horizontal) {
                                val w = ch.width
                                ch.width = (w + 1) / 2
                                r.width = w / 2
                                if (ch.hshift >= 0) { ch.hshift++; r.hshift++ }
                            } else {
                                val h = ch.height
                                ch.height = (h + 1) / 2
                                r.height = h / 2
                                if (ch.vshift >= 0) { ch.vshift++; r.vshift++ }
                            }
                            residuals.add(r)
                        }
                        if (sp.inPlace) {
                            chans.addAll(end, residuals)
                        } else {
                            chans.addAll(residuals)
                        }
                    }
                }
            }
            return nbMeta
        }

        private fun defaultSqueeze(sp: MutableList<SqueezeParam>, chans: List<Channel>, nbMeta: Int) {
            val first = nbMeta
            if (first >= chans.size) fail("invalid squeeze")
            var w = chans[first].width
            var h = chans[first].height
            if (chans.size - first >= 3) {
                val next = chans[first + 1]
                if (next.width == w && next.height == h) {
                    sp.add(SqueezeParam(true, false, first + 1, 2))
                    sp.add(SqueezeParam(false, false, first + 1, 2))
                }
            }
            val numC = chans.size - first
            if (h >= w && h > 8) {
                sp.add(SqueezeParam(false, true, first, numC))
                h = (h + 1) / 2
            }
            while (w > 8 || h > 8) {
                if (w > 8) {
                    sp.add(SqueezeParam(true, true, first, numC))
                    w = (w + 1) / 2
                }
                if (h > 8) {
                    sp.add(SqueezeParam(false, true, first, numC))
                    h = (h + 1) / 2
                }
            }
        }
    }

    /** チャンネル [from] から [to] (含まない) を復号する */
    fun decodeChannels(br: BitReader, streamIdx: Int, from: Int, to: Int) {
        var distMultiplier = 0
        for (i in from until to) distMultiplier = maxOf(distMultiplier, channels[i].width)
        val reader = tree.code.newReader()
        reader.begin(br)
        for (i in from until to) {
            val ch = channels[i]
            if (ch.width == 0 || ch.height == 0) continue
            val t = tree.specialize(i, streamIdx)
            // 同じ大きさ・シフトの前のチャンネル (近い順)
            val prev = ArrayList<Channel>()
            for (j in i - 1 downTo from) {
                val c = channels[j]
                if (c.width == ch.width && c.height == ch.height && c.hshift == ch.hshift && c.vshift == ch.vshift) {
                    prev.add(c)
                }
            }
            ChannelDecoder.decode(br, reader, t, ch, prev, header.wp, distMultiplier)
        }
        reader.finish()
    }

    /** 全チャンネルを復号して逆変換する (グループに分けないサブストリーム用) */
    fun decodeAll(br: BitReader, streamIdx: Int): List<Channel> {
        decodeChannels(br, streamIdx, 0, channels.size)
        return finish()
    }

    /** 逆変換して、変換前のチャンネル配置で返す */
    fun finish(): List<Channel> {
        val chans = channels
        for (ti in header.transforms.indices.reversed()) {
            when (val tr = header.transforms[ti]) {
                is Transform.Rct -> inverseRct(tr, chans)
                is Transform.Palette -> inversePalette(tr, chans, bitDepth, header.wp)
                is Transform.Squeeze -> {
                    for (si in tr.params.indices.reversed()) {
                        val sp = tr.params[si]
                        val begin = sp.beginC
                        val end = begin + sp.numC
                        val residuals = ArrayList<Channel>()
                        if (sp.inPlace) {
                            for (k in 0 until sp.numC) residuals.add(chans.removeAt(end))
                        } else {
                            val start = chans.size - sp.numC
                            if (start < end) fail("invalid squeeze")
                            for (k in 0 until sp.numC) residuals.add(chans.removeAt(start))
                        }
                        for (k in 0 until sp.numC) {
                            val merged = if (sp.horizontal) Squeeze.inverseH(chans[begin + k], residuals[k])
                            else Squeeze.inverseV(chans[begin + k], residuals[k])
                            chans[begin + k] = merged
                        }
                    }
                }
            }
        }
        if (chans.size != original.size) fail("modular channel count mismatch")
        for (i in chans.indices) {
            val c = chans[i]
            val o = original[i]
            if (c.width != o.width || c.height != o.height) fail("modular channel size mismatch")
            c.hshift = o.hshift
            c.vshift = o.vshift
        }
        return chans
    }
}

internal fun inverseRct(tr: Transform.Rct, chans: MutableList<Channel>) {
    val a = chans[tr.beginC].data
    val b = chans[tr.beginC + 1].data
    val c = chans[tr.beginC + 2].data
    val perm = tr.type / 7
    val type = tr.type % 7
    val n = a.size
    for (i in 0 until n) {
        val va = a[i]
        val vb = b[i]
        val vc = c[i]
        var d: Int
        var e: Int
        var f: Int
        if (type == 6) {
            val tmp = va - (vc shr 1)
            e = vc + tmp
            f = tmp - (vb shr 1)
            d = f + vb
        } else {
            d = va
            f = if (type and 1 != 0) vc + va else vc
            e = when (type shr 1) {
                1 -> vb + va
                2 -> vb + ((va + f) shr 1)
                else -> vb
            }
        }
        a[i] = d
        b[i] = e
        c[i] = f
    }
    // チャンネルの並べ替え
    val ca = chans[tr.beginC]
    val cb = chans[tr.beginC + 1]
    val cc = chans[tr.beginC + 2]
    when (perm) {
        1 -> { swapData(ca, cb); swapData(ca, cc) }
        2 -> { swapData(ca, cb); swapData(cb, cc) }
        3 -> swapData(cb, cc)
        4 -> swapData(ca, cb)
        5 -> swapData(ca, cc)
    }
}

private fun swapData(a: Channel, b: Channel) {
    val t = a.data
    a.data = b.data
    b.data = t
}

private val DELTA_PALETTE = intArrayOf(
    0, 0, 0, 4, 4, 4, 11, 0, 0, 0, 0, -13, 0, -12, 0, -10, -10, -10,
    -18, -18, -18, -27, -27, -27, -18, -18, 0, 0, 0, -32, -32, 0, 0, -37, -37, -37,
    0, -32, -32, 24, 24, 45, 50, 50, 50, -45, -24, -24, -24, -45, -45, 0, -24, -24,
    -34, -34, 0, -24, 0, -24, -45, -45, -24, 64, 64, 64, -32, 0, -32, 0, -32, 0,
    -32, 0, 32, -24, -45, -24, 45, 24, 45, 24, -24, -45, -45, -24, 24, 80, 80, 80,
    64, 0, 0, 0, 0, -64, 0, -64, -64, -24, -24, 45, 96, 96, 96, 64, 64, 0,
    45, -24, -24, 34, -34, 0, 112, 112, 112, 24, -45, -45, 45, 45, -24, 0, -32, 32,
    24, -24, 45, 0, 96, 96, 45, -24, 24, 24, -45, -24, -24, -45, 24, 0, -64, 0,
    96, 0, 0, 128, 128, 128, 64, 0, 64, 144, 144, 144, 96, 96, 0, -36, -36, 36,
    45, -24, -45, 45, -45, -24, 0, 0, -96, 0, 128, 128, 0, 96, 0, 45, 24, -45,
    -128, 0, 0, 24, -45, 24, -45, 24, -45, 64, 0, -64, 64, -64, -64, 96, 0, 96,
    45, -45, 24, 24, 45, -45, 64, 64, -64, 128, 128, 0, 0, 0, -128, -24, 45, -45,
)

internal fun inversePalette(tr: Transform.Palette, chans: MutableList<Channel>, bitDepth: Int, wp: WpParams) {
    val palette = chans.removeAt(0)
    val index = chans[tr.beginC]
    val w = index.width
    val h = index.height
    val numC = tr.numC
    val nbColours = tr.nbColours
    val nbDeltas = tr.nbDeltas
    val outs = Array(numC) { if (it == 0) Channel(w, h, index.hshift, index.vshift) else Channel(w, h, index.hshift, index.vshift) }
    for (o in outs) o.allocate()
    val idx = index.data
    var needDelta = false
    val pal = palette.data
    val palW = palette.width
    for (i in 0 until w * h) {
        val v = idx[i]
        if (v < nbDeltas) needDelta = true
        if (v in 0 until nbColours) {
            for (c in 0 until numC) outs[c].data[i] = pal[c * palW + v]
        } else if (v >= nbColours) {
            val vi = v - nbColours
            if (vi < 64) {
                for (c in 0 until numC) {
                    outs[c].data[i] = ((vi shr (2 * c)) % 4) * ((1 shl bitDepth) - 1) / 4 + (1 shl maxOf(0, bitDepth - 3))
                }
            } else {
                var k = vi - 64
                for (c in 0 until numC) {
                    outs[c].data[i] = (k % 5) * ((1 shl bitDepth) - 1) / 4
                    k /= 5
                }
            }
        } else {
            for (c in 0 until numC) {
                if (c >= 3) {
                    outs[c].data[i] = 0
                    continue
                }
                val k = (-(v + 1)) % 143
                var s = DELTA_PALETTE[((k + 1) shr 1) * 3 + c]
                if (k and 1 == 0) s = -s
                if (bitDepth > 8) s = s shl (minOf(bitDepth, 24) - 8)
                outs[c].data[i] = s
            }
        }
    }
    if (needDelta) {
        // nbDeltas 未満の索引は予測値との差分
        for (c in 0 until numC) {
            val d = outs[c].data
            val wpred = if (tr.dPred == Predictor.WEIGHTED) WeightedPredictor(w, wp) else null
            for (y in 0 until h) {
                val row = y * w
                for (x in 0 until w) {
                    val i = row + x
                    val west = if (x > 0) d[i - 1] else if (y > 0) d[i - w] else 0
                    val north = if (y > 0) d[i - w] else west
                    val nw = if (x > 0 && y > 0) d[i - w - 1] else west
                    val ne = if (y > 0 && x + 1 < w) d[i - w + 1] else north
                    val nn = if (y > 1) d[i - 2 * w] else north
                    val ww = if (x > 1) d[i - 2] else west
                    val nee = if (y > 0 && x + 2 < w) d[i - w + 2] else ne
                    wpred?.predict(north, nw, ne, west, nn)
                    if (idx[i] < nbDeltas) {
                        val p = ChannelDecoder.predict(tr.dPred, west, north, nw, ne, nn, ww, nee, wpred)
                        d[i] = d[i] + p
                    }
                    wpred?.record(d[i])
                }
            }
        }
    }
    chans.removeAt(tr.beginC)
    for (c in numC - 1 downTo 0) chans.add(tr.beginC, outs[c])
}

internal object Squeeze {
    private fun tendency(a: Int, b: Int, c: Int): Int {
        if (a >= b && b >= c) {
            var x = (4 * a - 3 * c - b + 6) / 12
            if (x - (x and 1) > 2 * (a - b)) x = 2 * (a - b) + 1
            if (x + (x and 1) > 2 * (b - c)) x = 2 * (b - c)
            return x
        } else if (a <= b && b <= c) {
            var x = (4 * a - 3 * c - b - 6) / 12
            if (x + (x and 1) < 2 * (a - b)) x = 2 * (a - b) - 1
            if (x - (x and 1) < 2 * (b - c)) x = 2 * (b - c)
            return x
        }
        return 0
    }

    fun inverseH(avg: Channel, res: Channel): Channel {
        val w = avg.width + res.width
        val h = avg.height
        if (res.height != h || (avg.width != res.width && avg.width != res.width + 1)) fail("invalid squeeze")
        val out = Channel(w, h, avg.hshift, avg.vshift)
        out.allocate()
        val aw = avg.width
        val rw = res.width
        val a = avg.data
        val r = res.data
        val o = out.data
        for (y in 0 until h) {
            val ar = y * aw
            val rr = y * rw
            val orow = y * w
            if (aw == 0) continue
            var av = a[ar]
            var left = av
            for (x in 0 until rw) {
                val next = if (x + 1 < aw) a[ar + x + 1] else av
                val diff = r[rr + x] + tendency(left, av, next)
                val first = av + diff / 2
                val second = first - diff
                o[orow + 2 * x] = first
                o[orow + 2 * x + 1] = second
                av = next
                left = second
            }
            if (aw > rw) o[orow + w - 1] = a[ar + aw - 1]
        }
        return out
    }

    fun inverseV(avg: Channel, res: Channel): Channel {
        val w = avg.width
        val h = avg.height + res.height
        if (res.width != w || (avg.height != res.height && avg.height != res.height + 1)) fail("invalid squeeze")
        val out = Channel(w, h, avg.hshift, avg.vshift)
        out.allocate()
        val ah = avg.height
        val rh = res.height
        val a = avg.data
        val r = res.data
        val o = out.data
        if (ah == 0) return out
        // 行単位で処理する
        val avgRow = IntArray(w)
        val top = IntArray(w)
        System.arraycopy(a, 0, avgRow, 0, w)
        System.arraycopy(a, 0, top, 0, w)
        for (y in 0 until rh) {
            val nextBase = if (y + 1 < ah) (y + 1) * w else -1
            for (x in 0 until w) {
                val av = avgRow[x]
                val next = if (nextBase >= 0) a[nextBase + x] else av
                val diff = r[y * w + x] + tendency(top[x], av, next)
                val first = av + diff / 2
                val second = first - diff
                o[(2 * y) * w + x] = first
                o[(2 * y + 1) * w + x] = second
                avgRow[x] = next
                top[x] = second
            }
        }
        if (ah > rh) System.arraycopy(a, (ah - 1) * w, o, (h - 1) * w, w)
        return out
    }
}
