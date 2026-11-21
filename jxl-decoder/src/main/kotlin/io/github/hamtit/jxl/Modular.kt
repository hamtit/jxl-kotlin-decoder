package io.github.hamtit.jxl

/** Modular 画像の 1 チャンネル (整数サンプル) */
internal class Channel(var width: Int, var height: Int, var hshift: Int, var vshift: Int) {
    var data: IntArray = IntArray(0)

    /** (x, y) のサンプルは data[offset + y * stride + x]。allocate した配列は offset 0・stride = width */
    var offset = 0
    var stride = 0

    fun allocate() {
        val n = width.toLong() * height
        if (n > Int.MAX_VALUE) fail("channel too large")
        data = IntArray(n.toInt())
        offset = 0
        stride = width
    }

    /** 別の (大きな) 配列の矩形を直接使う */
    fun bind(data: IntArray, offset: Int, stride: Int) {
        this.data = data
        this.offset = offset
        this.stride = stride
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

    /** 自己補正予測器を使うか (プロパティ 15 か、葉の予測器) */
    val usesWp: Boolean

    /** 参照する前のチャンネルの数 (プロパティ 16 以上) */
    val maxPrev: Int

    init {
        var wp = false
        var mp = 0
        for (i in property.indices) {
            val p = property[i]
            if (p == 15) wp = true
            if (p >= 16) mp = maxOf(mp, (p - 16) / 4 + 1)
            if (p < 0 && predictor[i] == Predictor.WEIGHTED) wp = true
        }
        usesWp = wp
        maxPrev = mp
    }

    /** 全ての分岐が同じプロパティ [property] を見る木の、値 → 葉の表。値は [min]..[max] に切り詰めてから引く */
    class SinglePropertyLut(val property: Int, val min: Int, val max: Int, val leaf: IntArray)

    /** 分岐が 1 種類のプロパティ (2〜15) だけで、値の範囲が狭い木なら表を作る ([specialize] した木で使う) */
    val singlePropertyLut: SinglePropertyLut? by lazy(LazyThreadSafetyMode.PUBLICATION) { buildLut() }

    private fun buildLut(): SinglePropertyLut? {
        if (size == 1) return null
        var p = -1
        var lo = Int.MAX_VALUE
        var hi = Int.MIN_VALUE
        for (i in 0 until size) {
            val q = property[i]
            if (q < 0) continue
            if (q < 2 || q > 15) return null
            if (p < 0) p = q else if (p != q) return null
            lo = minOf(lo, splitValue[i])
            hi = maxOf(hi, splitValue[i])
        }
        // 最大の分割値より大きい値は全て同じ葉になる
        val max = hi.toLong() + 1
        if (max - lo >= LUT_LIMIT) return null
        val n = (max - lo + 1).toInt()
        val leaf = IntArray(n)
        for (k in 0 until n) {
            val v = lo + k
            var node = 0
            while (property[node] >= 0) node = if (v > splitValue[node]) left[node] else right[node]
            leaf[k] = node
        }
        return SinglePropertyLut(p, lo, max.toInt(), leaf)
    }

    /** チャンネル番号・ストリーム番号で分岐するノード */
    private val chanStreamNodes: IntArray by lazy(LazyThreadSafetyMode.PUBLICATION) {
        val l = IntArrayList()
        for (i in 0 until size) if (property[i] == 0 || property[i] == 1) l.add(i)
        l.toArray()
    }

    /** [specialize] の結果。チャンネル・ストリームでの分岐の結果が同じなら同じ木になる */
    private val specialized = java.util.concurrent.ConcurrentHashMap<BranchKey, MaTree>()

    private class BranchKey(val bits: LongArray) {
        private val hash = bits.contentHashCode()
        override fun hashCode() = hash
        override fun equals(other: Any?) = other is BranchKey && other.bits.contentEquals(bits)
    }

    companion object {
        private const val LUT_LIMIT = 1L shl 14
        private const val SPECIALIZED_CACHE_LIMIT = 1024

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

    /** チャンネル番号とストリーム番号で決まる分岐を畳み込んだ木を作る (分岐の結果が同じなら作った木を使い回す) */
    fun specialize(channel: Int, stream: Int): MaTree {
        val nodes = chanStreamNodes
        val bits = LongArray((nodes.size + 63) ushr 6)
        for (k in nodes.indices) {
            val n = nodes[k]
            val v = if (property[n] == 0) channel else stream
            if (v > splitValue[n]) bits[k ushr 6] = bits[k ushr 6] or (1L shl k)
        }
        val key = BranchKey(bits)
        specialized[key]?.let { return it }
        val t = specializeUncached(channel, stream)
        if (specialized.size < SPECIALIZED_CACHE_LIMIT) specialized.putIfAbsent(key, t)
        return t
    }

    private fun specializeUncached(channel: Int, stream: Int): MaTree {
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

/**
 * 自己補正予測器 (weighted predictor)。4 つの部分予測は配列にせず個別の変数で持つ (画素ごとに呼ばれるため)
 */
internal class WeightedPredictor(private val width: Int, private val wp: WpParams) {
    private val trueErrRow = IntArray(width)
    private val subErrRow = IntArray(width * 4)
    private var x = 0
    private var y = 0
    private var teW = 0
    private var teNW = 0
    private var teN = 0
    private var teNE = 0
    // 部分予測ごとの誤差 (NWWW: 2 つ前の NW、NW、NE)
    private var nwww0 = 0
    private var nwww1 = 0
    private var nwww2 = 0
    private var nwww3 = 0
    private var nw0 = 0
    private var nw1 = 0
    private var nw2 = 0
    private var nw3 = 0
    private var ne0 = 0
    private var ne1 = 0
    private var ne2 = 0
    private var ne3 = 0
    private var sp0 = 0L
    private var sp1 = 0L
    private var sp2 = 0L
    private var sp3 = 0L
    private val maxW0 = wp.w0.toLong()
    private val maxW1 = wp.w1.toLong()
    private val maxW2 = wp.w2.toLong()
    private val maxW3 = wp.w3.toLong()
    var prediction = 0L
        private set
    var maxError = 0
        private set

    private fun weightOf(a: Int, b: Int, c: Int, maxW: Long): Long {
        val e = ((a.toLong() and 0xffffffffL) + (b.toLong() and 0xffffffffL) + (c.toLong() and 0xffffffffL)) and 0xffffffffL
        val t = (e + 1) shr 5
        val shift = if (t == 0L) 0 else 63 - java.lang.Long.numberOfLeadingZeros(t)
        return 4 + ((maxW * DIV_LOOKUP[((e shr shift) + 1).toInt()]) shr shift)
    }

    fun predict(n: Int, nw: Int, ne: Int, w: Int, nn: Int) {
        val tw = teW.toLong()
        val tnw = teNW.toLong()
        val tn = teN.toLong()
        val tne = teNE.toLong()
        val nS = n.toLong() shl 3
        val nwS = nw.toLong() shl 3
        val neS = ne.toLong() shl 3
        val wS = w.toLong() shl 3
        val nnS = nn.toLong() shl 3
        val s0 = wS + neS - nS
        val s1 = nS - (((tw + tn + tne) * wp.p1) shr 5)
        val s2 = wS - (((tw + tn + tnw) * wp.p2) shr 5)
        val s3 = nS - ((tnw * wp.p3a + tn * wp.p3b + tne * wp.p3c + (nnS - nS) * wp.p3d + (nwS - wS) * wp.p3e) shr 5)
        sp0 = s0
        sp1 = s1
        sp2 = s2
        sp3 = s3

        val w0 = weightOf(nwww0, nw0, ne0, maxW0)
        val w1 = weightOf(nwww1, nw1, ne1, maxW1)
        val w2 = weightOf(nwww2, nw2, ne2, maxW2)
        val w3 = weightOf(nwww3, nw3, ne3, maxW3)
        val logWeight = 63 - java.lang.Long.numberOfLeadingZeros((w0 + w1 + w2 + w3) shr 4)
        // 重みは Int に切り詰めてから (元の実装と同じ) 右シフトする
        val v0 = w0.toInt() ushr logWeight
        val v1 = w1.toInt() ushr logWeight
        val v2 = w2.toInt() ushr logWeight
        val v3 = w3.toInt() ushr logWeight
        val sumW2 = v0 + v1 + v2 + v3
        val sum = (sumW2.toLong() shr 1) - 1 + s0 * v0 + s1 * v1 + s2 * v2 + s3 * v3
        var pred = (sum * DIV_LOOKUP[sumW2]) shr 24
        if (((tn xor tw) or (tn xor tnw)) <= 0) {
            val mn = minOf(nS, minOf(wS, neS))
            val mx = maxOf(nS, maxOf(wS, neS))
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
        val row = subErrRow
        val e0 = ((kotlin.math.abs(sp0 - s8) + 3) shr 3).toInt()
        val e1 = ((kotlin.math.abs(sp1 - s8) + 3) shr 3).toInt()
        val e2 = ((kotlin.math.abs(sp2 - s8) + 3) shr 3).toInt()
        val e3 = ((kotlin.math.abs(sp3 - s8) + 3) shr 3).toInt()
        row[base] = e0
        row[base + 1] = e1
        row[base + 2] = e2
        row[base + 3] = e3
        trueErrRow[x] = trueErr
        x++
        if (x >= width) {
            y++
            x = 0
            teW = 0
            teN = trueErrRow[0]
            teNW = teN
            nw0 = row[0]; nw1 = row[1]; nw2 = row[2]; nw3 = row[3]
            nwww0 = nw0; nwww1 = nw1; nwww2 = nw2; nwww3 = nw3
            if (width <= 1) {
                teNE = teN
                ne0 = nw0; ne1 = nw1; ne2 = nw2; ne3 = nw3
            } else {
                teNE = trueErrRow[1]
                ne0 = row[4]; ne1 = row[5]; ne2 = row[6]; ne3 = row[7]
            }
        } else {
            teW = trueErr
            teNW = teN
            teN = teNE
            nwww0 = nw0; nwww1 = nw1; nwww2 = nw2; nwww3 = nw3
            nw0 = ne0 + e0; nw1 = ne1 + e1; nw2 = ne2 + e2; nw3 = ne3 + e3
            if (x + 1 >= width) {
                teNE = teN
                ne0 = nw0; ne1 = nw1; ne2 = nw2; ne3 = nw3
            } else if (y != 0) {
                teNE = trueErrRow[x + 1]
                val nb = (x + 1) * 4
                ne0 = row[nb]; ne1 = row[nb + 1]; ne2 = row[nb + 2]; ne3 = row[nb + 3]
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
    /** 書き込み先に直接復号しているか (finish の結果は書き込み先の配列を指す) */
    val direct: Boolean = false,
) {
    companion object {
        /**
         * ヘッダと (必要なら) 木を読み、変換後のチャンネルを用意する。
         * [channelsIn] は幅・高さ・シフトだけを使う。
         */
        fun read(
            br: BitReader,
            channelsIn: List<Channel>,
            globalTree: MaTree?,
            bitDepth: Int,
            /**
             * channelsIn と同じ並びの書き込み先 (offset・stride 付き)。変換が RCT だけ (チャンネルの並びと大きさが
             * 変わらない) なら、配列を確保せずここに直接復号する。[ModularImage.direct] でどちらになったか分かる
             */
            targets: List<Channel>? = null,
        ): ModularImage {
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
            val direct = targets != null && header.transforms.all { it is Transform.Rct }
            if (direct) {
                for (i in chans.indices) {
                    val t = targets[i]
                    chans[i].bind(t.data, t.offset, t.stride)
                }
            } else {
                for (c in chans) c.allocate()
            }
            return ModularImage(header, tree, channelsIn, chans, nbMeta, bitDepth, direct)
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

/** RCT の並べ替え: 出力チャンネル k に入るのは、逆変換した 3 つの値 (d, e, f) のうち RCT_PERM[perm][k] 番目 */
private val RCT_PERM = arrayOf(
    intArrayOf(0, 1, 2), intArrayOf(2, 0, 1), intArrayOf(1, 2, 0),
    intArrayOf(0, 2, 1), intArrayOf(1, 0, 2), intArrayOf(2, 1, 0),
)

internal fun inverseRct(tr: Transform.Rct, chans: MutableList<Channel>) {
    val ca = chans[tr.beginC]
    val cb = chans[tr.beginC + 1]
    val cc = chans[tr.beginC + 2]
    val a = ca.data
    val b = cb.data
    val c = cc.data
    val pm = RCT_PERM[tr.type / 7]
    val pa = pm[0]
    val pb = pm[1]
    val type = tr.type % 7
    val w = ca.width
    // チャンネルはグループの矩形 (offset・stride) のこともあるので行ごとに処理する
    for (y in 0 until ca.height) {
        val ia = ca.offset + y * ca.stride
        val ib = cb.offset + y * cb.stride
        val ic = cc.offset + y * cc.stride
        for (x in 0 until w) {
            val va = a[ia + x]
            val vb = b[ib + x]
            val vc = c[ic + x]
            val d: Int
            val e: Int
            val f: Int
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
            // 並べ替え (pa, pb が決まれば残りは 3 - pa - pb)
            a[ia + x] = if (pa == 0) d else if (pa == 1) e else f
            b[ib + x] = if (pb == 0) d else if (pb == 1) e else f
            val pc = 3 - pa - pb
            c[ic + x] = if (pc == 0) d else if (pc == 1) e else f
        }
    }
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
