package io.github.hamtit.jxl

import java.util.concurrent.ExecutorService

/** 浮動小数点の画像平面 */
internal class Plane(val width: Int, val height: Int) {
    val data: FloatArray = FloatArray(checkedSize(width, height))

    companion object {
        fun checkedSize(w: Int, h: Int): Int {
            val n = w.toLong() * h
            if (n > Int.MAX_VALUE - 16) fail("image too large")
            return n.toInt()
        }
    }
}

/** 復号したフレーム (アップサンプリング後、フレーム座標) */
internal class FrameImage(
    val header: FrameHeader,
    val width: Int,
    val height: Int,
    /** 色チャンネル (1 または 3) */
    var color: Array<FloatArray>,
    val extra: Array<FloatArray>,
    /** 色がまだ XYB (または YCbCr) のままか */
    var colorTransformPending: Boolean,
)

internal class LfGlobalData(
    val patches: Patches?,
    val splines: Splines?,
    val noise: FloatArray?,
    val mXLf: Float,
    val mYLf: Float,
    val mBLf: Float,
    val globalScale: Int,
    val quantLf: Int,
    val blockCtx: BlockContextMap?,
    val colourFactor: Int,
    val baseCorrX: Float,
    val baseCorrB: Float,
    val xFactorLf: Int,
    val bFactorLf: Int,
    val globalTree: MaTree?,
    val gmodular: ModularImage?,
    val gmodularDecoded: Int,
)

internal class BlockContextMap(
    val qfThresholds: IntArray,
    val lfThresholds: Array<IntArray>,
    val ctxMap: IntArray,
    val numClusters: Int,
) {
    companion object {
        private val DEFAULT_MAP = intArrayOf(
            0, 1, 2, 2, 3, 3, 4, 5, 6, 6, 6, 6, 6, 7, 8, 9, 9, 10, 11, 12, 13, 14, 14, 14,
            14, 14, 7, 8, 9, 9, 10, 11, 12, 13, 14, 14, 14, 14, 14,
        )

        fun parse(br: BitReader): BlockContextMap {
            if (br.readBool()) return BlockContextMap(IntArray(0), Array(3) { IntArray(0) }, DEFAULT_MAP, 15)
            var bsize = 1
            val lf = Array(3) {
                val n = br.readBits(4)
                bsize *= n + 1
                IntArray(n) { unpackSigned(br.readU32(0, 4, 16, 8, 272, 16, 65808, 32)) }
            }
            val nq = br.readBits(4)
            bsize *= nq + 1
            val qf = IntArray(nq) { 1 + br.readU32(0, 2, 4, 3, 12, 5, 44, 8) }
            if (bsize > 64) fail("too many block contexts")
            val map = EntropyCode.readClusters(br, bsize * 39)
            var n = 0
            for (c in map) n = maxOf(n, c + 1)
            if (n > 16) fail("too many block clusters")
            return BlockContextMap(qf, lf, map, n)
        }
    }
}

internal class HfPass(val orders: Array<Array<IntArray?>>, val code: EntropyCode)

/** 参照フレーム (パッチ・ブレンドの元) */
internal class ReferenceFrame(val image: FrameImage)

internal class FrameDecoder(
    private val image: ImageHeader,
    private val header: FrameHeader,
    private val data: ByteArray,
    private val toc: Toc,
    private val refs: Array<ReferenceFrame?>,
    private val lfFrames: Array<FrameImage?>,
    private val executor: ExecutorService?,
    private val noiseSeed: Long,
) {
    private val numLf = header.numLfGroups
    private val numGroups = header.numGroups
    private var sharedReader: BitReader? = null

    private fun section(idx: Int): BitReader {
        if (toc.singleEntry) {
            return sharedReader ?: BitReader(data, toc.offsets[0], toc.offsets[0] + toc.sizes[0]).also { sharedReader = it }
        }
        return BitReader(data, toc.offsets[idx], toc.offsets[idx] + toc.sizes[idx])
    }

    private fun lfGroupSection(i: Int) = section(1 + i)
    private fun hfGlobalSection() = section(1 + numLf)
    private fun passGroupSection(p: Int, g: Int) = section(2 + numLf + p * numGroups + g)

    // ---- 共通の情報 ----
    private val colorW = header.colorWidth
    private val colorH = header.colorHeight
    private val gd = header.groupDim
    private val hShift = IntArray(3)
    private val vShift = IntArray(3)
    private val maxHs: Int
    private val maxVs: Int

    init {
        val up = header.jpegUpsampling
        val hscale = up.any { it == 1 || it == 2 }
        val vscale = up.any { it == 1 || it == 3 }
        for (c in 0 until 3) {
            when (up[c]) {
                0 -> { hShift[c] = if (hscale) 1 else 0; vShift[c] = if (vscale) 1 else 0 }
                1 -> { hShift[c] = 0; vShift[c] = 0 }
                2 -> { hShift[c] = 0; vShift[c] = if (vscale) 1 else 0 }
                else -> { hShift[c] = if (hscale) 1 else 0; vShift[c] = 0 }
            }
        }
        maxHs = if (hscale) 1 else 0
        maxVs = if (vscale) 1 else 0
    }

    /** JPEG のクロマサブサンプリングを考慮したチャンネルの大きさ */
    private fun jpegShiftSize(w: Int, h: Int, c: Int): IntArray {
        val ww = if (maxHs != 0) { val s = (w + 1) / 2; if (hShift[c] != 0) s else s * 2 } else w
        val hh = if (maxVs != 0) { val s = (h + 1) / 2; if (vShift[c] != 0) s else s * 2 } else h
        return intArrayOf(ww, hh)
    }

    lateinit var lfGlobal: LfGlobalData

    // VarDCT の状態
    private var bw = 0
    private var bh = 0
    private lateinit var lfPlanes: Array<Plane>
    /**
     * ブロックごとの、LF の量子化値から決まる HF の文脈の番号 (BlockContextMap の LF のしきい値で分けたもの)。
     * しきい値がない (既定の文脈の表など) か LF フレームを使うなら null (常に 0)
     */
    private var lfCtx: ByteArray? = null
    private lateinit var dctSelect: ByteArray
    private lateinit var hfMul: IntArray
    private lateinit var epfSigma: FloatArray
    private lateinit var xFromY: IntArray
    private lateinit var bFromY: IntArray
    private var tilesW = 0
    private var colorPlanes: Array<Plane>? = null
    private var quantMatrices: QuantMatrices? = null
    private var numHfPresets = 1
    private lateinit var hfPasses: Array<HfPass>

    private var tMark = System.nanoTime()
    private fun phase(name: String) {
        if (TIMING) {
            val now = System.nanoTime()
            System.err.println("  %-14s %7.2f ms".format(name, (now - tMark) / 1e6))
            tMark = now
        }
    }

    /**
     * 逐次描画したときの出力 (ARGB)。[target] が null ならフレームの大きさの配列を作る。
     * そうでなければ、幅 [targetW]・高さ [targetH] の [target] の (ox, oy) にフレームを描く (はみ出す部分は捨てる)
     */
    class StreamOutput(
        val target: IntArray? = null,
        val targetW: Int = 0,
        val targetH: Int = 0,
        val ox: Int = 0,
        val oy: Int = 0,
    ) {
        var pixels: IntArray? = null

        /** フレームの大きさの ARGB を出力にする */
        fun put(frame: IntArray, w: Int, h: Int) {
            val t = target
            if (t == null) {
                pixels = frame
                return
            }
            val xa = maxOf(0, -ox)
            val xb = minOf(w, targetW - ox)
            if (xb > xa) {
                for (y in maxOf(0, -oy) until minOf(h, targetH - oy)) {
                    System.arraycopy(frame, y * w + xa, t, (y + oy) * targetW + ox + xa, xb - xa)
                }
            }
            pixels = t
        }
    }

    /**
     * フレームを復号する。[stream] を渡し、逐次描画できるフレームなら、全体の浮動小数点の平面を作らずに
     * ARGB を [StreamOutput.pixels] に書いて null を返す。
     */
    fun decode(stream: StreamOutput? = null): FrameImage? {
        tMark = System.nanoTime()
        val br0 = section(0)
        lfGlobal = readLfGlobal(br0)
        if (!toc.singleEntry) br0.checkOverrun()

        val vardct = !header.isModular
        if (vardct) initVarDct()
        if (stream != null && canStreamModular()) {
            phase("lfGlobal")
            stream.pixels = decodeModularStreaming(stream)
            phase("modularGroups")
            return null
        }
        allocateModular()
        phase("lfGlobal")

        // HfGlobal (量子化行列・係数の順序・エントロピー符号) は LF グループと関係なく読めるので、
        // 別のスレッドで LF グループ (1 つだと 1 スレッドでしか読めない) と並行して読む
        val hfGlobalTask = if (vardct && !toc.singleEntry && executor != null) {
            executor.submit {
                val br = hfGlobalSection()
                readHfGlobal(br)
                br.checkOverrun()
            }
        } else {
            null
        }
        try {
            // LF グループ
            if (toc.singleEntry) {
                decodeLfGroup(0, section(1))
            } else {
                parallelFor(executor, numLf) { i ->
                    val br = lfGroupSection(i)
                    decodeLfGroup(i, br)
                    br.checkOverrun()
                }
            }
            phase("lfGroups")
            if (vardct) {
                finishLf()
                phase("lfSmoothing")
            }
        } catch (e: Throwable) {
            // 失敗したら HfGlobal の終わりを待ってから投げる (このフレームの状態を後から書き換えさせない)
            try {
                hfGlobalTask?.get()
            } catch (_: Exception) {
            }
            throw e
        }
        if (hfGlobalTask != null) {
            try {
                hfGlobalTask.get()
            } catch (e: java.util.concurrent.ExecutionException) {
                throw e.cause ?: e
            }
            phase("hfGlobal")
        } else if (vardct) {
            val br = hfGlobalSection()
            readHfGlobal(br)
            if (!toc.singleEntry) br.checkOverrun()
            phase("hfGlobal")
        }

        if (stream != null && canStream()) {
            stream.pixels = decodeStreaming(stream)
            phase("streaming")
            return null
        }
        colorPlanes = Array(3) { c -> Plane((bw shr hShift[c]) * 8, (bh shr vShift[c]) * 8) }

        // パスグループ
        if (toc.singleEntry) {
            decodeGroup(0, null, 0, 0, 0)
        } else {
            parallelFor(executor, numGroups) { g -> decodeGroup(g, null, 0, 0, 0) }
        }
        if (toc.singleEntry) sharedReader?.checkOverrun()

        phase("groups")
        val modularChannels = lfGlobal.gmodular?.finish()
        phase("modularFinish")
        if (stream != null && modularChannels != null && canDirectModular()) {
            stream.put(modularToArgb(modularChannels), colorW, colorH)
            phase("modularToArgb")
            return null
        }
        return render(modularChannels)
    }

    // ---------------- LfGlobal ----------------

    private fun readLfGlobal(br: BitReader): LfGlobalData {
        val patches = if (header.hasPatches) Patches.parse(br, image, header) else null
        val splines = if (header.hasSplines) Splines.parse(br, header) else null
        val noise = if (header.hasNoise) FloatArray(8) { br.readBits(10) / 1024f } else null
        var mx = 1f / 32f
        var my = 1f / 4f
        var mb = 1f / 2f
        if (!br.readBool()) {
            mx = br.readF16()
            my = br.readF16()
            mb = br.readF16()
        }
        var globalScale = 0
        var quantLf = 0
        var blockCtx: BlockContextMap? = null
        var colourFactor = 84
        var baseX = 0f
        var baseB = 1f
        var xF = 128
        var bF = 128
        if (!header.isModular) {
            globalScale = br.readU32(1, 11, 2049, 11, 4097, 12, 8193, 16)
            quantLf = br.readU32(16, 0, 1, 5, 1, 8, 1, 16)
            blockCtx = BlockContextMap.parse(br)
            if (!br.readBool()) {
                colourFactor = br.readU32(84, 0, 256, 0, 2, 8, 258, 16)
                baseX = br.readF16()
                baseB = br.readF16()
                xF = br.readBits(8)
                bF = br.readBits(8)
            }
        }
        // GlobalModular
        var tree: MaTree? = null
        if (br.readBool()) {
            val nch = header.encodedColorChannels + image.extraChannels.size
            val limit = minOf(1L shl 22, 1024L + header.width.toLong() * header.height * nch / 16).toInt()
            tree = MaTree.parse(br, limit)
        }
        val chans = ArrayList<Channel>()
        if (header.isModular) {
            if (header.doYCbCr) {
                for (c in 0 until 3) {
                    val s = jpegShiftSize(colorW, colorH, c)
                    chans.add(Channel(s[0], s[1], hShift[c], vShift[c]))
                }
            } else {
                for (c in 0 until header.encodedColorChannels) chans.add(Channel(colorW, colorH, 0, 0))
            }
        }
        val colorShift = Integer.numberOfTrailingZeros(header.upsampling)
        for ((i, ec) in image.extraChannels.withIndex()) {
            val s = Integer.numberOfTrailingZeros(header.ecUpsampling[i]) + ec.dimShift - colorShift
            if (s < 0) fail("invalid extra channel upsampling")
            val w = ((colorW.toLong() + (1L shl s) - 1) shr s).toInt()
            val h = ((colorH.toLong() + (1L shl s) - 1) shr s).toInt()
            chans.add(Channel(w, h, s, s))
        }
        var gmod: ModularImage? = null
        var decoded = 0
        if (chans.isNotEmpty()) {
            // グループで復号するチャンネルの配列は、グループごとに描画しない場合だけ後で確保する (allocateModular)
            gmod = ModularImage.read(br, chans, tree, image.bitDepth.bitsPerSample, allocate = false)
            val list = gmod.channels
            var end = 0
            while (end < list.size) {
                val c = list[end]
                if (end >= gmod.nbMetaChannels && (c.width > gd || c.height > gd)) break
                end++
            }
            for (i in 0 until end) list[i].allocate()
            gmod.decodeChannels(br, 0, 0, end)
            decoded = end
        }
        return LfGlobalData(
            patches, splines, noise, mx, my, mb, globalScale, quantLf, blockCtx, colourFactor,
            baseX, baseB, xF, bF, tree, gmod, decoded,
        )
    }

    // ---------------- Modular のグループ分割 ----------------

    private class GroupChannel(val index: Int, val x0: Int, val y0: Int, val ch: Channel)

    /** グローバル Modular 画像のうち、矩形 (グループ) に属するチャンネルの部分 */
    private fun groupChannels(lf: Boolean, groupX: Int, groupY: Int, minShift: Int, maxShift: Int): List<GroupChannel> {
        val g = lfGlobal.gmodular ?: return emptyList()
        val out = ArrayList<GroupChannel>()
        val list = g.channels
        for (i in lfGlobal.gmodularDecoded until list.size) {
            val c = list[i]
            if (c.hshift < 0 || c.vshift < 0) continue
            val isLf = c.hshift >= 3 && c.vshift >= 3
            if (lf != isLf) continue
            if (!lf) {
                val s = minOf(c.hshift, c.vshift)
                if (s < minShift || s > maxShift) continue
            }
            val dim = if (lf) gd * 8 else gd
            val gw = dim shr c.hshift
            val gh = dim shr c.vshift
            if (gw == 0 || gh == 0) fail("channel shift too large")
            val x0 = groupX * gw
            val y0 = groupY * gh
            if (x0 >= c.width || y0 >= c.height) continue
            val w = minOf(gw, c.width - x0)
            val h = minOf(gh, c.height - y0)
            out.add(GroupChannel(i, x0, y0, Channel(w, h, c.hshift, c.vshift)))
        }
        return out
    }

    private fun decodeModularGroup(br: BitReader, parts: List<GroupChannel>, stream: Int) {
        if (parts.isEmpty()) return
        val g = lfGlobal.gmodular!!.channels
        // グローバルのチャンネルの矩形をそのまま書き込み先にする (変換が RCT だけなら一時配列とコピーがいらない)
        val targets = parts.map { p ->
            val dst = g[p.index]
            Channel(p.ch.width, p.ch.height, p.ch.hshift, p.ch.vshift).also {
                it.bind(dst.data, p.y0 * dst.width + p.x0, dst.width)
            }
        }
        val img = ModularImage.read(br, parts.map { it.ch }, lfGlobal.globalTree, image.bitDepth.bitsPerSample, targets)
        val result = img.decodeAll(br, stream)
        if (img.direct) return
        for ((k, p) in parts.withIndex()) {
            val dst = g[p.index]
            val src = result[k]
            for (y in 0 until src.height) {
                System.arraycopy(src.data, y * src.width, dst.data, (p.y0 + y) * dst.width + p.x0, src.width)
            }
        }
    }

    // ---------------- LF グループ ----------------

    private fun initVarDct() {
        var bwl = (colorW + 7) / 8
        var bhl = (colorH + 7) / 8
        if (maxHs != 0) bwl = (bwl + 1) / 2 * 2
        if (maxVs != 0) bhl = (bhl + 1) / 2 * 2
        bw = bwl
        bh = bhl
        val n = Plane.checkedSize(bw, bh)
        dctSelect = ByteArray(n) { -2 }
        hfMul = IntArray(n)
        epfSigma = FloatArray(n)
        tilesW = (colorW + 63) / 64
        val tilesH = (colorH + 63) / 64
        xFromY = IntArray(tilesW * tilesH)
        bFromY = IntArray(tilesW * tilesH)
        lfPlanes = Array(3) { c -> Plane(bw shr hShift[c], bh shr vShift[c]) }
        if (!header.useLfFrame && lfGlobal.blockCtx!!.lfThresholds.any { it.isNotEmpty() }) lfCtx = ByteArray(n)
    }

    private fun decodeLfGroup(idx: Int, br: BitReader) {
        val lgx = idx % header.lfGroupsPerRow
        val lgy = idx / header.lfGroupsPerRow
        val lfgd = gd * 8
        val x0 = lgx * lfgd
        val y0 = lgy * lfgd
        val lfw = minOf(lfgd, colorW - x0)
        val lfh = minOf(lfgd, colorH - y0)
        val vardct = !header.isModular
        if (vardct && !header.useLfFrame) {
            val extraPrecision = br.readBits(2)
            val w8 = (lfw + 7) / 8
            val h8 = (lfh + 7) / 8
            val order = intArrayOf(1, 0, 2)
            val chans = order.map { c ->
                val s = jpegShiftSize(w8, h8, c)
                Channel(s[0], s[1], hShift[c], vShift[c])
            }
            val img = ModularImage.read(br, chans, lfGlobal.globalTree, image.bitDepth.bitsPerSample)
            val out = img.decodeAll(br, 1 + idx)
            val m = floatArrayOf(lfGlobal.mXLf, lfGlobal.mYLf, lfGlobal.mBLf)
            val precision = (1 shl (9 - extraPrecision)).toDouble()
            val scaleInv = lfGlobal.globalScale.toDouble() * lfGlobal.quantLf
            for (k in 0 until 3) {
                val c = order[k]
                val ch = out[k]
                val scale = (m[c] * precision / scaleInv).toFloat()
                val plane = lfPlanes[c]
                val bx0 = (x0 / 8) shr hShift[c]
                val by0 = (y0 / 8) shr vShift[c]
                for (y in 0 until ch.height) {
                    val py = by0 + y
                    if (py >= plane.height) break
                    for (x in 0 until ch.width) {
                        val px = bx0 + x
                        if (px >= plane.width) break
                        plane.data[py * plane.width + px] = ch.data[y * ch.width + x] * scale
                    }
                }
            }
            lfCtx?.let { lc ->
                // ブロックごとに、3 チャンネルの量子化値をしきい値で分けた文脈の番号を求めておく (X・B・Y の順)
                val th = lfGlobal.blockCtx!!.lfThresholds
                val bx0 = x0 / 8
                val by0 = y0 / 8
                val bx1 = minOf(bw, bx0 + gd)
                val by1 = minOf(bh, by0 + gd)
                val ctxOrder = intArrayOf(0, 2, 1)
                // チャンネル c の復号結果は out[order の中の c の位置]
                val chOf = Array(3) { c -> out[order.indexOf(c)] }
                for (by in by0 until by1) for (bx in bx0 until bx1) {
                    var idx = 0
                    for (c in ctxOrder) {
                        val ch = chOf[c]
                        val qx = (bx shr hShift[c]) - (bx0 shr hShift[c])
                        val qy = (by shr vShift[c]) - (by0 shr vShift[c])
                        val t = th[c]
                        idx *= t.size + 1
                        if (qx < ch.width && qy < ch.height) {
                            val q = ch.data[qy * ch.width + qx]
                            for (v in t) if (q > v) idx++
                        }
                    }
                    lc[by * bw + bx] = idx.toByte()
                }
            }
        }
        // Modular の LF グループ
        val parts = groupChannels(true, lgx, lgy, 0, 0)
        decodeModularGroup(br, parts, 1 + numLf + idx)
        if (vardct) decodeHfMetadata(idx, br, x0, y0, lfw, lfh)
    }

    private fun decodeHfMetadata(idx: Int, br: BitReader, x0: Int, y0: Int, lfw: Int, lfh: Int) {
        var lbw = (lfw + 7) / 8
        var lbh = (lfh + 7) / 8
        if (maxHs != 0) lbw = (lbw + 1) / 2 * 2
        if (maxVs != 0) lbh = (lbh + 1) / 2 * 2
        val nbBits = addLog2Ceil(lbw * lbh - 1)
        val nbBlocks = 1 + br.readBits(nbBits)
        val tw = (lfw + 63) / 64
        val th = (lfh + 63) / 64
        val chans = listOf(Channel(tw, th, 0, 0), Channel(tw, th, 0, 0), Channel(nbBlocks, 2, 0, 0), Channel(lbw, lbh, 0, 0))
        val img = ModularImage.read(br, chans, lfGlobal.globalTree, image.bitDepth.bitsPerSample)
        val out = img.decodeAll(br, 1 + 2 * numLf + idx)
        // CfL の係数
        val tx0 = x0 / 64
        val ty0 = y0 / 64
        for (y in 0 until th) for (x in 0 until tw) {
            xFromY[(ty0 + y) * tilesW + tx0 + x] = out[0].data[y * tw + x]
            bFromY[(ty0 + y) * tilesW + tx0 + x] = out[1].data[y * tw + x]
        }
        val info = out[2].data
        val sharp = out[3].data
        val epf = header.epf
        val bx0 = x0 / 8
        val by0 = y0 / 8
        var dataIdx = 0
        val sigmaBase = if (epf != null) epf.quantMul * 65536f / lfGlobal.globalScale else 0f
        for (y in 0 until lbh) {
            var x = 0
            while (x < lbw) {
                val pos = (by0 + y) * bw + bx0 + x
                if (dctSelect[pos].toInt() != -2) {
                    x++
                    continue
                }
                if (dataIdx >= nbBlocks) fail("block info does not fill the LF group")
                val t = info[dataIdx]
                if (t < 0 || t >= Tt.COUNT) fail("invalid transform type")
                val mul = info[nbBlocks + dataIdx] + 1
                if (mul <= 0) fail("non-positive HfMul")
                val dw = Tt.BW[t]
                val dh = Tt.BH[t]
                if ((x % 32) + dw > 32 || (y % 32) + dh > 32) fail("varblock crosses a group border")
                if (x + dw > lbw || y + dh > lbh) fail("varblock does not fit")
                val sigma = sigmaBase / mul
                for (dy in 0 until dh) for (dx in 0 until dw) {
                    val p = (by0 + y + dy) * bw + bx0 + x + dx
                    if (dctSelect[p].toInt() != -2) fail("varblocks overlap")
                    dctSelect[p] = if (dx == 0 && dy == 0) t.toByte() else -1
                    if (dx == 0 && dy == 0) hfMul[p] = mul
                    if (epf != null) {
                        val s = sharp[(y + dy) * lbw + x + dx]
                        if (s < 0 || s > 7) fail("invalid EPF sharpness")
                        epfSigma[p] = sigma * epf.sharpLut[s]
                    }
                }
                dataIdx++
                x += dw
            }
        }
    }

    /** LF 画像の CfL と適応平滑化 */
    private fun finishLf() {
        if (header.useLfFrame) {
            val lf = lfFrames[header.lfLevel] ?: fail("LF frame not available")
            for (c in 0 until 3) {
                val p = lfPlanes[c]
                val src = lf.color[c]
                for (y in 0 until p.height) {
                    val sy = minOf(y, lf.height - 1)
                    for (x in 0 until p.width) {
                        val sx = minOf(x, lf.width - 1)
                        p.data[y * p.width + x] = src[sy * lf.width + sx]
                    }
                }
            }
            return
        }
        val subsampled = maxHs != 0 || maxVs != 0
        if (!subsampled) {
            val kx = lfGlobal.baseCorrX + (lfGlobal.xFactorLf - 128).toFloat() / lfGlobal.colourFactor
            val kb = lfGlobal.baseCorrB + (lfGlobal.bFactorLf - 128).toFloat() / lfGlobal.colourFactor
            val px = lfPlanes[0].data
            val py = lfPlanes[1].data
            val pb = lfPlanes[2].data
            for (i in py.indices) {
                val y = py[i]
                px[i] += kx * y
                pb[i] += kb * y
            }
        }
        if (!header.skipAdaptiveLfSmoothing) {
            if (subsampled) fail("adaptive LF smoothing with chroma subsampling")
            adaptiveLfSmoothing()
        }
    }

    private fun adaptiveLfSmoothing() {
        val scaleInv = lfGlobal.globalScale.toDouble() * lfGlobal.quantLf
        val lfx = (512.0 * lfGlobal.mXLf / scaleInv).toFloat()
        val lfy = (512.0 * lfGlobal.mYLf / scaleInv).toFloat()
        val lfb = (512.0 * lfGlobal.mBLf / scaleInv).toFloat()
        val w = lfPlanes[0].width
        val h = lfPlanes[0].height
        if (w <= 2 || h <= 2) return
        val sSelf = 0.052262735f
        val sSide = 0.2034514f
        val sDiag = 0.03348292f
        // 書き換える前の値は、1 つ上の行と今の行だけ退避しておけばよい (下の行はまだ書き換えていない)
        var up = Array(3) { lfPlanes[it].data.copyOfRange(0, w) }
        var mid = Array(3) { lfPlanes[it].data.copyOfRange(w, 2 * w) }
        val dst = Array(3) { lfPlanes[it].data }
        val mul = floatArrayOf(lfx, lfy, lfb)
        val wa = FloatArray(3)
        val self = FloatArray(3)
        for (y in 1 until h - 1) {
            for (x in 1 until w - 1) {
                var gap = 0.5f
                for (c in 0 until 3) {
                    val u = up[c]
                    val m = mid[c]
                    val d = dst[c]
                    val i = (y + 1) * w + x
                    val v = m[x]
                    val side = m[x - 1] + m[x + 1] + u[x] + d[i]
                    val diag = u[x - 1] + u[x + 1] + d[i - 1] + d[i + 1]
                    val a = v * sSelf + side * sSide + diag * sDiag
                    wa[c] = a
                    self[c] = v
                    val g = kotlin.math.abs(a - v) / mul[c]
                    if (g > gap) gap = g
                }
                val gs = maxOf(0f, 3f - 4f * gap)
                for (c in 0 until 3) dst[c][y * w + x] = (wa[c] - self[c]) * gs + self[c]
            }
            // 次の行のために、書き換えた今の行の元の値を上の行に回し、次の行 (まだ元の値) を退避する
            val t = up
            up = mid
            mid = t
            if (y + 1 < h - 1) for (c in 0 until 3) System.arraycopy(dst[c], (y + 1) * w, mid[c], 0, w)
        }
    }

    // ---------------- HfGlobal ----------------

    private fun readHfGlobal(br: BitReader) {
        quantMatrices = QuantMatrices.parse(br, numLf, lfGlobal.globalTree, image.bitDepth.bitsPerSample)
        numHfPresets = br.readBits(addLog2Ceil(numGroups - 1)) + 1
        val ctx = lfGlobal.blockCtx!!
        hfPasses = Array(header.passes.numPasses) {
            var used = br.readU32(0x5F, 0, 0x13, 0, 0, 0, 0, 13)
            val orders = Array(13) { arrayOfNulls<IntArray>(3) }
            if (used != 0) {
                val code = EntropyCode.parse(br, 8)
                val reader = code.newReader()
                reader.begin(br)
                for (o in 0 until 13) {
                    if (used and 1 != 0) {
                        val size = Tt.ORDER_W[o] * Tt.ORDER_H[o]
                        val nat = Tt.naturalOrder(o)
                        for (c in 0 until 3) {
                            val perm = readPermutation(br, reader, size, size / 64)
                            orders[o][c] = IntArray(size) { nat[perm[it]] }
                        }
                    }
                    used = used ushr 1
                }
                reader.finish()
            }
            val code = EntropyCode.parse(br, 495 * numHfPresets * ctx.numClusters)
            for (o in 0 until 13) for (c in 0 until 3) if (orders[o][c] == null) orders[o][c] = Tt.naturalOrder(o)
            HfPass(orders, code)
        }
    }

    // ---------------- パスグループ ----------------

    /**
     * グループを復号する。target が null なら colorPlanes に、そうでなければ画像の (targetX0, targetY0) を起点とする
     * 行の長さ targetStride の配列に書く
     */
    private fun decodeGroup(g: Int, target: Array<FloatArray>?, targetStride: Int, targetX0: Int, targetY0: Int) {
        val gx = g % header.groupsPerRow
        val gy = g / header.groupsPerRow
        val vardct = !header.isModular
        var gs: GroupScratch? = null
        var gbw = 0
        var gbh = 0
        if (vardct) {
            gbw = minOf(gd / 8, bw - gx * (gd / 8))
            gbh = minOf(gd / 8, bh - gy * (gd / 8))
            gs = GroupScratch.obtain(minOf(gd, bw * 8) * minOf(gd, bh * 8), header.passes.numPasses > 1)
            // 係数は書き込み先 (帯または色平面) のグループの場所に直接置き、その場で逆変換する
            val bx0 = gx * (gd / 8)
            val by0 = gy * (gd / 8)
            for (c in 0 until 3) {
                val cw = (gbw shr hShift[c]) * 8
                val chh = (gbh shr vShift[c]) * 8
                val px0 = (bx0 shr hShift[c]) * 8
                val py0 = (by0 shr vShift[c]) * 8
                if (target != null) {
                    gs.out[c] = target[c]
                    gs.outStride[c] = targetStride
                    gs.outBase[c] = (py0 - targetY0) * targetStride + px0 - targetX0
                } else {
                    val plane = colorPlanes!![c]
                    gs.out[c] = plane.data
                    gs.outStride[c] = plane.width
                    gs.outBase[c] = py0 * plane.width + px0
                }
                val o = gs.out[c]!!
                for (y in 0 until chh) {
                    val r = gs.outBase[c] + y * gs.outStride[c]
                    o.fill(0f, r, r + cw)
                }
                gs.intCoeffs?.get(c)?.fill(0, 0, cw * chh)
                gs.nzRows[c].fill(0, 0, cw * chh / 64)
            }
        }
        try {
            for (p in 0 until header.passes.numPasses) {
                val br = passGroupSection(p, g)
                if (vardct) decodeHfCoefficients(br, p, gx, gy, gbw, gbh, gs!!)
                val mm = header.passes.downsamplingBracket(p)
                val parts = groupChannels(false, gx, gy, mm[0], mm[1])
                decodeModularGroup(br, parts, 1 + 3 * numLf + 17 + p * numGroups + g)
                if (!toc.singleEntry) br.checkOverrun()
            }
            if (vardct) {
                // 複数のパスの係数は整数で足し合わせてから、一度だけ浮動小数点にする
                gs!!.intCoeffs?.let { ints ->
                    for (c in 0 until 3) {
                        val cw = (gbw shr hShift[c]) * 8
                        val chh = (gbh shr vShift[c]) * 8
                        val src = ints[c]
                        val dst = gs.out[c]!!
                        for (y in 0 until chh) {
                            val r = gs.outBase[c] + y * gs.outStride[c]
                            for (x in 0 until cw) dst[r + x] = src[y * cw + x].toFloat()
                        }
                    }
                }
                reconstructGroup(gx, gy, gbw, gbh, gs)
            }
        } finally {
            gs?.recycle()
        }
    }

    /**
     * グループの逆変換の作業領域。グループごと・フレームごとに確保せず、[ScratchPool] に置いて使い回す。
     * 係数は書き込み先 ([out] の [outBase] から、行の長さ [outStride]) のグループの場所に直接置き、その場で逆変換する。
     * パスが 1 つなら、量子化された係数 (整数) を読んだその場で逆量子化して書く
     * (各係数は一度だけ書かれるので、整数から変換するのと同じ値になる)。
     * 複数のパスがあるときは [intCoeffs] (グループの幅を行の長さとする) に整数で足し合わせる
     */
    private class GroupScratch(val size: Int, multiPass: Boolean) {
        val out = arrayOfNulls<FloatArray>(3)
        val outBase = IntArray(3)
        val outStride = IntArray(3)
        val intCoeffs: Array<IntArray>? = if (multiPass) Array(3) { IntArray(size) } else null
        /**
         * ブロック (チャンネルごと、グループ内のブロック座標 y * 幅 + x) の HF 係数のうち 0 でないものがある行の数
         * (係数の配置での上からの行数。0 なら HF 係数はすべて 0)。逆 DCT で 0 だけの行を飛ばすのに使う
         */
        val nzRows = Array(3) { IntArray(maxOf(1, size / 64)) }
        val scratch = FloatArray(512)
        val tmp = FloatArray(64)
        val lfBuf = FloatArray(32 * 32)
        val dct = Dct.Scratch()

        fun recycle() {
            out.fill(null)
            ScratchPool.give(this, size * (if (intCoeffs != null) 12L else 0L) + size * 3L / 16 + 8192L + dct.cols.bytes())
        }

        companion object {
            fun obtain(size: Int, multiPass: Boolean): GroupScratch =
                ScratchPool.take(GroupScratch::class.java) {
                    // 大きすぎるもの (大きいグループ用) を小さい画像に使わない
                    it.size >= size && it.size / 4 <= size && (!multiPass || it.intCoeffs != null)
                } ?: GroupScratch(size, multiPass)
        }
    }

    private fun decodeHfCoefficients(br: BitReader, pass: Int, gx: Int, gy: Int, gbw: Int, gbh: Int, gs: GroupScratch) {
        val outI = gs.intCoeffs
        val nzRows = gs.nzRows
        val ctx = lfGlobal.blockCtx!!
        // パスが 1 つなら、係数を読んだその場で逆量子化する (reconstructGroup で全係数を見直さない)
        val qm = quantMatrices!!
        val qmScale = qmScales()
        val qb = image.quantBias
        val qbn = image.quantBiasNumerator
        val globalScale = lfGlobal.globalScale.toFloat()
        val hfPass = hfPasses[pass]
        val shift = if (pass < header.passes.shift.size) header.passes.shift[pass] else 0
        val hfpBits = addLog2Ceil(numHfPresets - 1)
        val hfp = br.readBits(hfpBits)
        if (hfp >= numHfPresets) fail("HF preset out of range")
        val ctxSize = 495 * ctx.numClusters
        val code = hfPass.code
        val clusterMap = code.clusters
        val ctxBase = ctxSize * hfp
        val reader = code.newReader()
        reader.begin(br)
        // ANS で LZ77 なし (libjxl などの普通の出力) なら、表を 1 つの配列から引く速い経路で読む
        val fastB = code.fastBuckets
        val fastLa = code.fastLogAlphabet
        val fastC = code.fastConfigs
        val lfIdxMul = (ctx.lfThresholds[0].size + 1) * (ctx.lfThresholds[1].size + 1) * (ctx.lfThresholds[2].size + 1)
        val hfIdxMul = ctx.qfThresholds.size + 1
        val nzRow = Array(3) { c -> IntArray(maxOf(1, (gbw + (1 shl hShift[c]) - 1) shr hShift[c])) }
        val bx0 = gx * (gd / 8)
        val by0 = gy * (gd / 8)
        val cOrder = intArrayOf(1, 0, 2)
        val lc = lfCtx
        for (y in 0 until gbh) {
            for (x in 0 until gbw) {
                val pos = (by0 + y) * bw + bx0 + x
                val t = dctSelect[pos].toInt()
                if (t < 0) continue
                val qf = hfMul[pos]
                val w8 = Tt.BW[t]
                val h8 = Tt.BH[t]
                val numBlocks = w8 * h8
                val logNb = Integer.numberOfTrailingZeros(numBlocks)
                val orderId = Tt.ORDER[t]
                val lfIdx = if (lc != null) lc[pos].toInt() else 0
                var hfIdx = 0
                for (v in ctx.qfThresholds) if (qf > v) hfIdx++
                for (ci in 0 until 3) {
                    val c = cOrder[ci]
                    val hs = hShift[c]
                    val vs = vShift[c]
                    val sx = x shr hs
                    val sy = y shr vs
                    if (hs != 0 || vs != 0) {
                        if ((sx shl hs) != x || (sy shl vs) != y) continue
                        // サブサンプリングしたチャンネルでは、その位置にブロックの起点がある場合だけ
                        val p2 = (by0 + sy) * bw + bx0 + sx
                        if (dctSelect[p2] < 0) continue
                    }
                    val chIdx = ci * 13 + orderId
                    val idx = (chIdx * hfIdxMul + hfIdx) * lfIdxMul + lfIdx
                    val blockCtx = ctx.ctxMap[idx]
                    val nzr = nzRow[c]
                    val predicted = if (sy == 0) {
                        if (sx == 0) 32 else nzr[sx - 1]
                    } else if (sx == 0) {
                        nzr[sx]
                    } else {
                        (nzr[sx] + nzr[sx - 1] + 1) shr 1
                    }
                    val nzCtxIdx = if (predicted >= 8) 4 + predicted / 2 else predicted
                    val nzCtx = blockCtx + nzCtxIdx * ctx.numClusters
                    var nonZeros = if (fastB != null) {
                        reader.readAnsFast(br, fastB, fastLa, fastC, clusterMap[ctxBase + nzCtx])
                    } else {
                        reader.readClustered(br, clusterMap[ctxBase + nzCtx], 0)
                    }
                    if (nonZeros > (63 shl logNb)) fail("too many non-zero coefficients")
                    val nzVal = (nonZeros + numBlocks - 1) shr logNb
                    for (dx in 0 until w8) {
                        if (sx + dx < nzr.size) nzr[sx + dx] = nzVal
                    }
                    if (nonZeros == 0) continue
                    val coeffF = gs.out[c]!!
                    val fBase = gs.outBase[c] + sy * 8 * gs.outStride[c] + sx * 8
                    val fStride = gs.outStride[c]
                    val coeffI = outI?.get(c)
                    val cw = (gbw shr hs) * 8
                    var prevNz = if (nonZeros <= numBlocks * 4) 1 else 0
                    val order = hfPass.orders[orderId][c]!!
                    val ow = Tt.ORDER_W[orderId]
                    val transpose = Tt.needTranspose(t)
                    val coeffBase = ctxBase + blockCtx * 458 + 37 * ctx.numClusters
                    val size = order.size
                    val m = qm.get(t, c)
                    val bwv = w8 * 8
                    val scale = 65536f / (globalScale * qf) * qmScale[c]
                    val qbc = qb[c]
                    var rows = 0
                    var k = numBlocks
                    while (k < size) {
                        val i = k - numBlocks
                        val nzc = (nonZeros - 1) shr logNb
                        val cctx = (COEFF_NUM_NONZERO_CONTEXT[nzc] + COEFF_FREQ_CONTEXT[i shr logNb]) * 2 + prevNz
                        if (cctx >= 458) fail("invalid coefficient context")
                        val u = if (fastB != null) {
                            reader.readAnsFast(br, fastB, fastLa, fastC, clusterMap[coeffBase + cctx])
                        } else {
                            reader.readClustered(br, clusterMap[coeffBase + cctx], 0)
                        }
                        if (u == 0) {
                            prevNz = 0
                            k++
                            continue
                        }
                        val v = unpackSigned(u) shl shift
                        val pos2 = order[k]
                        var dx = pos2 % ow
                        var dy = pos2 / ow
                        if (transpose) {
                            val tmp = dx; dx = dy; dy = tmp
                        }
                        if (coeffI != null) {
                            coeffI[(sy * 8 + dy) * cw + sx * 8 + dx] += v
                        } else {
                            // 各係数は一度だけ書かれるので、整数のまま足し合わせてから逆量子化するのと同じ値になる
                            val q = v.toFloat()
                            val dq = if (q == 1f || q == -1f) q * qbc else q - qbn / q
                            coeffF[fBase + dy * fStride + dx] = dq * m[dy * bwv + dx] * scale
                        }
                        if (dy >= rows) rows = dy + 1
                        prevNz = 1
                        nonZeros--
                        if (nonZeros == 0) break
                        k++
                    }
                    if (nonZeros != 0) fail("too many zeros in varblock")
                    val bi = sy * (gbw shr hs) + sx
                    if (rows > nzRows[c][bi]) nzRows[c][bi] = rows
                }
            }
        }
        reader.finish()
    }

    /** 書き込み先に置いた係数を逆量子化 (複数パスのとき)・CfL・逆 DCT する */
    private fun reconstructGroup(gx: Int, gy: Int, gbw: Int, gbh: Int, gs: GroupScratch) {
        val qm = quantMatrices!!
        val bx0 = gx * (gd / 8)
        val by0 = gy * (gd / 8)
        val qmScale = qmScales()
        val qb = image.quantBias
        val qbn = image.quantBiasNumerator
        // パスが 1 つなら係数は読んだときに逆量子化してある
        if (gs.intCoeffs != null) for (c in 0 until 3) {
            val stride = gs.outStride[c]
            val dst = gs.out[c]!!
            forEachVarblock(c, gbw, gbh, bx0, by0) { sx, sy, t, mul ->
                val bwv = Tt.BW[t] * 8
                val bhv = Tt.BH[t] * 8
                val m = qm.get(t, c)
                val scale = 65536f / (lfGlobal.globalScale.toFloat() * mul) * qmScale[c]
                for (yy in 0 until bhv) {
                    val row = gs.outBase[c] + (sy * 8 + yy) * stride + sx * 8
                    val mrow = yy * bwv
                    for (xx in 0 until bwv) {
                        val qf = dst[row + xx]
                        if (qf == 0f) continue
                        val v = if (qf == 1f || qf == -1f) qf * qb[c] else qf - qbn / qf
                        dst[row + xx] = v * m[mrow + xx] * scale
                    }
                }
            }
        }
        // CfL (サブサンプリングなしのとき)
        if (maxHs == 0 && maxVs == 0) {
            val cw = gbw * 8
            val ch = gbh * 8
            val fx = gs.out[0]!!
            val fy = gs.out[1]!!
            val fb = gs.out[2]!!
            // サブサンプリングがなければ 3 チャンネルとも同じ位置・行の長さ
            val o0 = gs.outBase[0]
            val stride = gs.outStride[0]
            if (gs.outBase[1] != o0 || gs.outBase[2] != o0 || gs.outStride[1] != stride || gs.outStride[2] != stride) {
                error("channel layout mismatch")
            }
            val tx0 = gx * gd / 64
            val ty0 = gy * gd / 64
            val cf = lfGlobal.colourFactor.toFloat()
            for (y in 0 until ch) {
                val ty = ty0 + y / 64
                for (x64 in 0 until (cw + 63) / 64) {
                    val ti = ty * tilesW + tx0 + x64
                    val kx = lfGlobal.baseCorrX + xFromY[ti] / cf
                    val kb = lfGlobal.baseCorrB + bFromY[ti] / cf
                    val xEnd = minOf(cw, x64 * 64 + 64)
                    val row = o0 + y * stride
                    for (x in x64 * 64 until xEnd) {
                        val yv = fy[row + x]
                        fx[row + x] += kx * yv
                        fb[row + x] += kb * yv
                    }
                }
            }
        }
        // 逆変換
        val scratch = gs.scratch
        val tmp = gs.tmp
        val lfBuf = gs.lfBuf
        for (c in 0 until 3) {
            val cw = gs.outStride[c]
            val buf = gs.out[c]!!
            val o0 = gs.outBase[c]
            val lf = lfPlanes[c]
            val lbx0 = bx0 shr hShift[c]
            val lby0 = by0 shr vShift[c]
            val nzRows = gs.nzRows[c]
            // CfL で X・B には Y の係数が足されるので、Y の 0 でない行も含める
            val nzY = if (c != 1 && maxHs == 0 && maxVs == 0) gs.nzRows[1] else nzRows
            val bwc = gbw shr hShift[c]
            forEachVarblock(c, gbw, gbh, bx0, by0) { sx, sy, t, _ ->
                val bwb = Tt.BW[t]
                val bhb = Tt.BH[t]
                val base = o0 + (sy * 8) * cw + sx * 8
                val nz = maxOf(nzRows[sy * bwc + sx], nzY[sy * bwc + sx])
                if (bwb == 1 && bhb == 1) {
                    buf[base] = lf.data[(lby0 + sy) * lf.width + lbx0 + sx]
                } else {
                    for (y in 0 until bhb) for (x in 0 until bwb) {
                        lfBuf[y * bwb + x] = lf.data[(lby0 + sy + y) * lf.width + lbx0 + sx + x]
                    }
                    Dct.dct2d(lfBuf, 0, bwb, bwb, bhb, scratch, true)
                    val logbw = Integer.numberOfTrailingZeros(bwb)
                    val logbh = Integer.numberOfTrailingZeros(bhb)
                    for (y in 0 until bhb) for (x in 0 until bwb) {
                        buf[base + y * cw + x] = lfBuf[y * bwb + x] / (Dct.scaleF(y, 5 - logbh) * Dct.scaleF(x, 5 - logbw))
                    }
                }
                when (t) {
                    Tt.DCT2 -> Dct.dct2(buf, base, cw, tmp)
                    Tt.DCT4 -> Dct.dct4(buf, base, cw, tmp, gs.dct)
                    Tt.HORNUSS -> Dct.hornuss(buf, base, cw, tmp, gs.dct)
                    Tt.DCT4X8 -> Dct.dct4x8(buf, base, cw, false, gs.dct)
                    Tt.DCT8X4 -> Dct.dct4x8(buf, base, cw, true, gs.dct)
                    in Tt.AFV0..Tt.AFV3 -> Dct.afv(buf, base, cw, t - Tt.AFV0, gs.dct)
                    // HF 係数がすべて 0 の 8x8 は LF の値で埋まる
                    Tt.DCT8 -> if (nz == 0) Dct.fill8x8(buf, base, cw) else Dct.idct8x8(buf, base, cw, nz, scratch)
                    // 0 でない係数があるのは上から max(LF の行数, nz) 行だけ
                    else -> Dct.idct2d(buf, base, cw, bwb * 8, bhb * 8, maxOf(bhb, nz), scratch, gs.dct.cols)
                }
            }
        }
    }

    /** X・B の量子化行列の倍率 (Y は 1) */
    private fun qmScales(): FloatArray = floatArrayOf(
        Math.pow(0.8, (header.xQmScale - 2).toDouble()).toFloat(), 1f,
        Math.pow(0.8, (header.bQmScale - 2).toDouble()).toFloat(),
    )

    private inline fun forEachVarblock(c: Int, gbw: Int, gbh: Int, bx0: Int, by0: Int, f: (Int, Int, Int, Int) -> Unit) {
        val hs = hShift[c]
        val vs = vShift[c]
        for (y in 0 until gbh) for (x in 0 until gbw) {
            val pos = (by0 + y) * bw + bx0 + x
            val t = dctSelect[pos].toInt()
            if (t < 0) continue
            val sx = x shr hs
            val sy = y shr vs
            if (hs != 0 || vs != 0) {
                if ((sx shl hs) != x || (sy shl vs) != y) continue
                val p2 = (by0 + sy) * bw + bx0 + sx
                if (dctSelect[p2] < 0) continue
            }
            f(sx, sy, t, hfMul[pos])
        }
    }

    // ---------------- 描画 ----------------

    private fun render(modular: List<Channel>?): FrameImage {
        val ecCount = image.extraChannels.size
        val colorCount: Int
        var colors: Array<Plane>
        val ecPlanes = arrayOfNulls<Plane>(ecCount)
        val ecOffset: Int
        if (header.isModular) {
            colorCount = header.encodedColorChannels
            ecOffset = colorCount
            val m = modular!!
            if (image.xybEncoded) {
                val y = m[0]
                val x = m[1]
                val b = m[2]
                val mx = lfGlobal.mXLf / 128f
                val my = lfGlobal.mYLf / 128f
                val mb = lfGlobal.mBLf / 128f
                val px = Plane(colorW, colorH)
                val py = Plane(colorW, colorH)
                val pb = Plane(colorW, colorH)
                for (i in 0 until colorW * colorH) {
                    val yi = y.data[i]
                    px.data[i] = x.data[i] * mx
                    py.data[i] = yi * my
                    pb.data[i] = (b.data[i] + yi) * mb
                }
                x.data = IntArray(0)
                y.data = IntArray(0)
                b.data = IntArray(0)
                colors = arrayOf(px, py, pb)
            } else {
                colors = Array(colorCount) { c ->
                    val ch = m[c]
                    channelToPlane(ch, image.bitDepth).also { ch.data = IntArray(0) }
                }
            }
        } else {
            colorCount = 3
            ecOffset = 0
            val planes = colorPlanes!!
            colorPlanes = null
            colors = Array(3) { c ->
                val src = planes[c]
                val w = (colorW + (1 shl hShift[c]) - 1) shr hShift[c]
                val h = (colorH + (1 shl vShift[c]) - 1) shr vShift[c]
                if (src.width == w && src.height == h) src else {
                    val p = Plane(w, h)
                    for (y in 0 until h) System.arraycopy(src.data, y * src.width, p.data, y * w, w)
                    p
                }
            }
        }
        for (i in 0 until ecCount) {
            val ch = modular!![ecOffset + i]
            ecPlanes[i] = channelToPlane(ch, image.extraChannels[i].bitDepth)
            ch.data = IntArray(0)
        }

        // YCbCr のクロマアップサンプリング
        if (header.doYCbCr && (maxHs != 0 || maxVs != 0)) {
            colors = Array(3) { c -> Filters.upsampleJpeg(colors[c], hShift[c], vShift[c], colorW, colorH) }
        }

        // 復元フィルタ
        val gab = header.gabor
        val epf = header.epf
        if (gab != null || epf != null) {
            if (colors.size < 3) {
                val g = colors[0]
                colors = arrayOf(g, Plane(g.width, g.height).also { System.arraycopy(g.data, 0, it.data, 0, g.data.size) },
                    Plane(g.width, g.height).also { System.arraycopy(g.data, 0, it.data, 0, g.data.size) })
            }
            phase("toPlanes")
            if (gab != null) colors = Filters.gaborish(colors, gab, executor)
            phase("gaborish")
            if (epf != null) {
                val sigma = if (header.isModular) null else epfSigma
                colors = Filters.epf(colors, epf, sigma, bw, executor)
            }
            phase("epf(${epf?.iters})")
            if (colorCount == 1) colors = arrayOf(colors[0])
        }

        val fw = header.sampleWidth(1)
        val fh = header.sampleHeight(1)
        val img = FrameImage(
            header, fw, fh,
            Array(colors.size) { colors[it].data },
            Array(ecCount) { ecPlanes[it]!!.data },
            image.xybEncoded || header.doYCbCr,
        )

        // パッチ・スプライン・ノイズは色のアップサンプリング前 (色サンプル座標) に適用する
        val ups = header.upsampling
        val colorRes = FrameImage(header, colorW, colorH, img.color, img.extra, img.colorTransformPending)
        val ecUpsample = IntArray(ecCount) { Integer.numberOfTrailingZeros(header.ecUpsampling[it]) + image.extraChannels[it].dimShift - Integer.numberOfTrailingZeros(ups) }
        val ecAtColor = Array(ecCount) { i ->
            if (ecUpsample[i] == 0) ecPlanes[i]!! else Filters.upsample(ecPlanes[i]!!, ecUpsample[i], image, colorW, colorH)
        }
        val work = FrameImage(header, colorW, colorH, img.color, Array(ecCount) { ecAtColor[it].data }, img.colorTransformPending)
        lfGlobal.patches?.let { Features.applyPatches(it, work, refs, image) }
        lfGlobal.splines?.let { if (work.color.size == 3) Features.renderSplines(it, work, lfGlobal.baseCorrX, lfGlobal.baseCorrB, ups) }
        lfGlobal.noise?.let { if (work.color.size == 3) Features.renderNoise(it, work, header, lfGlobal.baseCorrX, lfGlobal.baseCorrB, noiseSeed) }

        // アップサンプリング
        var outColor = work.color
        var outExtra = work.extra
        if (ups > 1) {
            val shift = Integer.numberOfTrailingZeros(ups)
            outColor = Array(outColor.size) { c -> Filters.upsample(Plane(colorW, colorH).also { System.arraycopy(work.color[c], 0, it.data, 0, it.data.size) }, shift, image, fw, fh).data }
            outExtra = Array(ecCount) { c -> Filters.upsample(Plane(colorW, colorH).also { System.arraycopy(work.extra[c], 0, it.data, 0, it.data.size) }, shift, image, fw, fh).data }
        }
        return FrameImage(header, fw, fh, outColor, outExtra, img.colorTransformPending)
    }

    // ---------------- Modular の整数サンプルから直接描画 ----------------

    /**
     * 色変換もフィルタもいらない Modular の画像 (sRGB の可逆圧縮など) か。
     * その場合は浮動小数点の平面を作らず、整数のサンプルから ARGB を作る
     */
    private fun canDirectModular(): Boolean {
        if (!header.isModular || image.xybEncoded || header.doYCbCr) return false
        if (header.upsampling != 1 || header.lfLevel != 0) return false
        if (header.gabor != null || header.epf != null) return false
        if (header.hasPatches || header.hasSplines || header.hasNoise) return false
        if (image.bitDepth.isFloat) return false
        for ((i, ec) in image.extraChannels.withIndex()) {
            if (header.ecUpsampling[i] != 1 || ec.dimShift != 0) return false
        }
        // ColorConvert.toSrgb が何もしない色空間だけ
        val enc = image.colourEncoding
        if (!enc.wantIcc) {
            if (enc.transferFunction != TransferFunction.SRGB) return false
            if (enc.colourSpace == ColourSpace.RGB && (!enc.isSrgbPrimaries || !enc.isD65)) return false
        }
        return true
    }

    /** 整数のサンプルから ARGB を作る。出力には最初の色チャンネルの配列を使い回す */
    private fun modularToArgb(m: List<Channel>): IntArray {
        val n = Plane.checkedSize(colorW, colorH)
        // packModular が読むのは色チャンネルとアルファだけ
        val cc = header.encodedColorChannels
        val used = (0 until minOf(cc, 3)) + (if (image.alphaChannel >= 0) listOf(cc + image.alphaChannel) else emptyList())
        for (i in used) if (m[i].width != colorW || m[i].height != colorH) fail("modular channel size mismatch")
        val out = m[0].data
        if (out.size != n) fail("modular channel size mismatch")
        // 出力の画素 i は、どのチャンネルでも同じ位置 i から作るので、最初のチャンネルに上書きしてよい
        parallelRows(executor, colorH) { y0, y1 ->
            packModular(m, out, 0, colorW, 0, colorW, y0, y1)
        }
        for (c in m) if (c.data !== out) c.data = IntArray(0)
        return out
    }

    /**
     * Modular のチャンネル [m] (変換前の並び、幅 w) の行 ya..yb・列 xa..xb を ARGB にして、
     * 画素 (x, y) を out[outBase + y * outStride + x] に書く
     */
    private fun packModular(
        m: List<Channel>, out: IntArray, outBase: Int, outStride: Int, xa: Int, xb: Int, ya: Int, yb: Int,
    ) {
        if (xa >= xb) return
        val cc = header.encodedColorChannels
        val cr = m[0]
        val cg = if (cc >= 3) m[1] else cr
        val cb = if (cc >= 3) m[2] else cr
        val r = cr.data
        val g = cg.data
        val b = cb.data
        val alphaIdx = image.alphaChannel
        val ca = if (alphaIdx >= 0) m[cc + alphaIdx] else null
        val alpha = ca?.data
        val alphaDepth = if (alphaIdx >= 0) image.extraChannels[alphaIdx].bitDepth else null
        val premul = alphaIdx >= 0 && image.extraChannels[alphaIdx].alphaAssociated
        val fast = !premul && image.bitDepth.bitsPerSample == 8 &&
            (alphaDepth == null || (!alphaDepth.isFloat && alphaDepth.bitsPerSample == 8))
        val s = image.bitDepth.scale()
        val alphaScale = alphaDepth?.scale() ?: 0f
        for (y in ya until yb) {
            val ir = cr.offset + y * cr.stride
            val ig = cg.offset + y * cg.stride
            val ib = cb.offset + y * cb.stride
            val ia = if (ca != null) ca.offset + y * ca.stride else 0
            val io = outBase + y * outStride
            if (fast) {
                // 8 ビット整数を正規化してから 8 ビットに丸め直すと、0..255 への切り詰めと同じになる
                for (x in xa until xb) {
                    val rv = r[ir + x]
                    val gv = g[ig + x]
                    val bv = b[ib + x]
                    val av = if (alpha != null) alpha[ia + x] else 255
                    out[io + x] = if (((rv or gv or bv or av) ushr 8) == 0) {
                        (av shl 24) or (rv shl 16) or (gv shl 8) or bv
                    } else {
                        (av.coerceIn(0, 255) shl 24) or (rv.coerceIn(0, 255) shl 16) or
                            (gv.coerceIn(0, 255) shl 8) or bv.coerceIn(0, 255)
                    }
                }
            } else {
                // 高ビット深度や乗算済みアルファは浮動小数点の経路 (channelToPlane → toOutput) と同じ計算にする
                for (x in xa until xb) {
                    var a = 1f
                    if (alpha != null) {
                        val raw = alpha[ia + x]
                        a = if (alphaDepth!!.isFloat) alphaDepth.parseFloatSample(raw) else raw * alphaScale
                        a = a.coerceIn(0f, 1f)
                    }
                    var rv = r[ir + x] * s
                    var gv = g[ig + x] * s
                    var bv = b[ib + x] * s
                    if (premul && a > 0f) {
                        rv /= a; gv /= a; bv /= a
                    }
                    out[io + x] = (to8(a) shl 24) or (to8(rv) shl 16) or (to8(gv) shl 8) or to8(bv)
                }
            }
        }
    }

    /** グローバル Modular 画像のうち、まだ配列のないチャンネルを確保する */
    private fun allocateModular() {
        val g = lfGlobal.gmodular ?: return
        for (i in lfGlobal.gmodularDecoded until g.channels.size) g.channels[i].allocate()
    }

    /**
     * 可逆 (Modular) の画像を、グループごとに小さな配列へ復号して直接 ARGB に書けるか。
     * 全体にかかる変換が RCT だけ (画素ごとの計算) で、全チャンネルがグループに分かれていれば、
     * 画像全体の大きさのチャンネルの配列を作らずに済む
     */
    private fun canStreamModular(): Boolean {
        if (!canDirectModular() || toc.singleEntry || header.passes.numPasses != 1) return false
        val g = lfGlobal.gmodular ?: return false
        if (lfGlobal.gmodularDecoded != 0 || g.nbMetaChannels != 0) return false
        if (!g.header.transforms.all { it is Transform.Rct }) return false
        val cc = header.encodedColorChannels
        if (g.channels.size != cc + image.extraChannels.size) return false
        for (c in g.channels) {
            if (c.hshift != 0 || c.vshift != 0 || c.width != colorW || c.height != colorH) return false
        }
        return true
    }

    private fun decodeModularStreaming(stream: StreamOutput): IntArray {
        val g = lfGlobal.gmodular!!
        val out = stream.target ?: IntArray(Plane.checkedSize(colorW, colorH))
        val outW = if (stream.target != null) stream.targetW else colorW
        val outH = if (stream.target != null) stream.targetH else colorH
        val ox = if (stream.target != null) stream.ox else 0
        val oy = if (stream.target != null) stream.oy else 0
        val nch = g.channels.size
        val rcts = g.header.transforms.reversed().map { it as Transform.Rct }
        parallelFor(executor, numGroups) { gi ->
            val gx = gi % header.groupsPerRow
            val gy = gi / header.groupsPerRow
            val parts = groupChannels(false, gx, gy, 0, 0)
            if (parts.size != nch) fail("modular channel count mismatch")
            val gw = parts[0].ch.width
            val gh = parts[0].ch.height
            val buf = ScratchPool.ints(gw * gh * nch)
            try {
                val targets = parts.mapIndexed { k, p ->
                    if (p.index != k || p.ch.width != gw || p.ch.height != gh) fail("modular channel size mismatch")
                    Channel(gw, gh, 0, 0).also { it.bind(buf, k * gw * gh, gw) }
                }
                val br = passGroupSection(0, gi)
                val img = ModularImage.read(br, parts.map { it.ch }, lfGlobal.globalTree, image.bitDepth.bitsPerSample, targets)
                val result = img.decodeAll(br, 1 + 3 * numLf + 17 + gi).toMutableList()
                br.checkOverrun()
                for (tr in rcts) inverseRct(tr, result)
                // グループの画素 (x, y) はキャンバスの (x0 + ox + x, y0 + oy + y)。はみ出す部分は書かない
                val tx = parts[0].x0 + ox
                val ty = parts[0].y0 + oy
                packModular(result, out, ty * outW + tx, outW,
                    maxOf(0, -tx), minOf(gw, outW - tx), maxOf(0, -ty), minOf(gh, outH - ty))
            } finally {
                ScratchPool.give(buf)
            }
        }
        return out
    }

    // ---------------- 逐次描画 ----------------

    /** グループの行ごとに描画できるか (XYB の静止画で、フレーム全体が必要な処理がない) */
    private fun canStream(): Boolean {
        if (header.isModular || header.doYCbCr || header.upsampling != 1 || header.lfLevel != 0) return false
        if (!image.xybEncoded || header.hasSplines || header.hasNoise) return false
        // パッチは帯ごとに重ねられる (アルファなどの追加チャンネルがない場合)
        if (header.hasPatches && image.extraChannels.isNotEmpty()) return false
        if (image.colourEncoding.wantIcc) return false
        if (toc.singleEntry) return false
        val g = lfGlobal.gmodular
        if (g != null) {
            // 変換 (アルファの Squeeze など) があれば、追加チャンネルは全グループの復号後に逆変換してからアルファを書き込む。
            // 乗算済みアルファは色の計算にアルファが要るので、変換がないときだけ
            if (g.header.transforms.isNotEmpty() && image.alphaChannel >= 0 &&
                image.extraChannels[image.alphaChannel].alphaAssociated
            ) return false
            for ((i, ec) in image.extraChannels.withIndex()) {
                if (header.ecUpsampling[i] != 1 || ec.dimShift != 0) return false
            }
        }
        return true
    }

    private fun decodeStreaming(stream: StreamOutput): IntArray {
        val w = colorW
        val h = colorH
        // 出力先 (キャンバス) と、フレームの置き場所。はみ出す行・列は書かない
        val out = stream.target ?: IntArray(Plane.checkedSize(w, h))
        val outW = if (stream.target != null) stream.targetW else w
        val outH = if (stream.target != null) stream.targetH else h
        val ox = if (stream.target != null) stream.ox else 0
        val oy = if (stream.target != null) stream.oy else 0
        val xa = maxOf(0, -ox)
        val xb = minOf(w, outW - ox)
        val gab = header.gabor
        val epf = header.epf
        val epfSteps = Filters.epfSteps(epf)
        var epfHalo = 0
        for (st in epfSteps) epfHalo += Filters.epfRadius(st)
        val halo = epfHalo + (if (gab != null) 1 else 0)
        // 各グループは、自分の領域を左上に SHIFT 画素ずらした範囲を描く。その計算に要る逆 DCT 後の値は、
        // 自分の領域と、左・上・左上のグループの端の BORDER 画素 (先に復号されている) だけで足りる
        if (halo > SHIFT) error("filter halo too large")
        val alphaIdx = image.alphaChannel
        // 変換のある Modular 画像では、アルファは最後に逆変換してから書き込む (それまでは不透明として描く)
        val lateAlpha = lfGlobal.gmodular?.header?.transforms?.isNotEmpty() == true
        val alpha = if (alphaIdx >= 0 && !lateAlpha) lfGlobal.gmodular?.channels?.get(alphaIdx) else null
        val alphaDepth = if (alphaIdx >= 0) image.extraChannels[alphaIdx].bitDepth else null
        val premul = alphaIdx >= 0 && image.extraChannels[alphaIdx].alphaAssociated

        val ob = image.opsinBias
        val cb0 = Math.cbrt(ob[0].toDouble()).toFloat()
        val cb1 = Math.cbrt(ob[1].toDouble()).toFloat()
        val cb2 = Math.cbrt(ob[2].toDouble()).toFloat()
        val its = 255f / image.toneMapping.intensityTarget
        val m = image.opsinInverseMatrix
        val ob0 = ob[0]
        val ob1 = ob[1]
        val ob2 = ob[2]
        val m0 = m[0]
        val m1 = m[1]
        val m2 = m[2]
        val m3 = m[3]
        val m4 = m[4]
        val m5 = m[5]
        val m6 = m[6]
        val m7 = m[7]
        val m8 = m[8]
        val grey = image.colourEncoding.colourSpace == ColourSpace.GREY
        val alphaScale = alphaDepth?.scale() ?: 0f

        // グループの逆 DCT 後の値は、上に BORDER 行の余白を付けたタイル (行の長さ gd、gd + BORDER 行) に置く。
        // タイルの (0, 0) は画像の (x0, y0 - BORDER)。上の余白は前のグループ行の下端から写す。
        // 各グループは、描く範囲 (自分の領域を左上に SHIFT 画素ずらしたもの) のうち、左のグループとの継ぎ目
        // (x0 - SHIFT .. x0 + SHIFT) 以外を自分のタイルだけで描き、タイルはすぐに返す。継ぎ目は、左右のグループの
        // 端の BORDER 列ずつを写した幅 2 * BORDER の小さな配列から、両方の復号を終えた方のスレッドが描く
        val ts = gd
        val gpr = header.groupsPerRow
        val rowsOfGroups = header.groupsPerColumn
        val sw = bw * 8
        val seamW = 2 * BORDER
        // グループ行の下端 BORDER 行 (次のグループ行が上の余白に使う)。偶数・奇数のグループ行で交互に使う
        val rowStrip = Array(2) { Array(3) { ScratchPool.floats(BORDER * sw) } }
        // 継ぎ目 (グループ gx とその左のグループの間) の配列と、それを描くまでに待つ復号の残りの数
        val seams = java.util.concurrent.atomic.AtomicReferenceArray<Array<FloatArray>>(gpr)
        val pending = java.util.concurrent.atomic.AtomicIntegerArray(gpr)
        val scratchW = minOf(w, gd + BORDER)

        /**
         * 画像の列 [ax0, ax1)・行 [ay0, ay1) を描く。逆 DCT 後の値は src (行の長さ stride、(0, 0) が画像の (sx0, sy0)) の
         * 列 [wx0, wx1) を使う (描く列の左右に SHIFT 画素の余白。画像の端では余白なし)
         */
        fun render(src: Array<FloatArray>, stride: Int, sx0: Int, sy0: Int, wx0: Int, wx1: Int, ax0: Int, ax1: Int, ay0: Int, ay1: Int) {
            if (ax0 >= ax1 || ay0 >= ay1) return
            val ww = wx1 - wx0
            val wo = wx0 - sx0
            // Gaborish 後の行 (窓の列だけ。y は鏡映済みでなくてよい)
            val inputRow = Filters.RowSource { c, yy, dst, off ->
                val y = Filters.mirror(yy, h)
                val sc = src[c]
                val o = (y - sy0) * stride + wo
                if (gab == null) {
                    System.arraycopy(sc, o, dst, off, ww)
                } else {
                    val w0 = gab[2 * c]
                    val w1 = gab[2 * c + 1]
                    val norm = 1f / (1f + 4f * w0 + 4f * w1)
                    val ou = (maxOf(y - 1, 0) - sy0) * stride + wo
                    val od = (minOf(y + 1, h - 1) - sy0) * stride + wo
                    // 左右の端 (隣を鏡映する) だけ別に計算し、内側は分岐なしで回す
                    for (edge in 0 until minOf(ww, 2)) {
                        val x = edge * (ww - 1)
                        val xl = if (x > 0) x - 1 else 0
                        val xr = if (x + 1 < ww) x + 1 else ww - 1
                        val side = sc[ou + x] + sc[od + x] + sc[o + xl] + sc[o + xr]
                        val diag = sc[ou + xl] + sc[ou + xr] + sc[od + xl] + sc[od + xr]
                        dst[off + x] = (sc[o + x] + side * w0 + diag * w1) * norm
                    }
                    // 上下の和 (v) と左右の画素を変数で持ち回り、1 画素あたりの読み込みを 3 回にする
                    if (ww > 2) {
                        var vl = sc[ou] + sc[od]
                        var vm = sc[ou + 1] + sc[od + 1]
                        var sl = sc[o]
                        var sm = sc[o + 1]
                        for (x in 1 until ww - 1) {
                            val vr = sc[ou + x + 1] + sc[od + x + 1]
                            val sr = sc[o + x + 1]
                            dst[off + x] = (sm + (vm + sl + sr) * w0 + (vl + vr) * w1) * norm
                            vl = vm
                            vm = vr
                            sl = sm
                            sm = sr
                        }
                    }
                }
            }
            val ep = Filters.EpfRows.obtain()
            val lin = ScratchPool.floats(3 * scratchW)
            try {
                ep.begin(inputRow, ww, h, epfSteps, epf, epfSigma, bw, ay0, wx0 shr 3)
                val rx = ep.out(0)
                val ry = ep.out(1)
                val rb = ep.out(2)
                val rows = arrayOf(rx, ry, rb)
                // 線形 RGB の 1 行 (r・g・b を lin の 0・scratchW・2 * scratchW から)
                val gOff = scratchW
                val bOff = 2 * scratchW
                // 出力する列 (描く範囲とキャンバスの重なり)
                val cx0 = maxOf(ax0, xa)
                val cx1 = minOf(ax1, xb)
                for (y in ay0 until ay1) {
                    val o = ep.row(y)
                    // パッチは描く範囲の行 y に重ねる (最後の段階の配列は 1 行だけなので、行 y は配列の先頭)
                    lfGlobal.patches?.let { applyPatchesToBand(it, rows, ww + 6, y - o / (ww + 6), y, y + 1, ax0, ax1, wx0, h) }
                    val ty = y + oy
                    if (ty < 0 || ty >= outH || cx0 >= cx1) continue
                    // 窓の中の列 x (画像の列 wx0 + x) の値は rows[c][ro + x]
                    val ro = o + 3
                    val oo = y * w + wx0
                    val to = ty * outW + ox + wx0
                    val la = cx0 - wx0
                    val lz = cx1 - wx0
                    if (premul) {
                        for (x in la until lz) {
                            val i = ro + x
                            val gl = ry[i] + rx[i] - cb0
                            val gm = ry[i] - rx[i] - cb1
                            val gs = rb[i] - cb2
                            val ml = (gl * gl * gl + ob0) * its
                            val mm = (gm * gm * gm + ob1) * its
                            val ms = (gs * gs * gs + ob2) * its
                            var r = m0 * ml + m1 * mm + m2 * ms
                            var g = m3 * ml + m4 * mm + m5 * ms
                            var bb = m6 * ml + m7 * mm + m8 * ms
                            if (grey) {
                                val l = 0.2126f * r + 0.7152f * g + 0.0722f * bb
                                r = l; g = l; bb = l
                            }
                            var av = 1f
                            if (alpha != null) {
                                val raw = alpha.data[oo + x]
                                av = if (alphaDepth!!.isFloat) alphaDepth.parseFloatSample(raw) else raw * alphaScale
                                if (av < 0f) av = 0f else if (av > 1f) av = 1f
                            }
                            var sr = ColorConvert.srgbEncodeFast(r)
                            var sg = ColorConvert.srgbEncodeFast(g)
                            var sb = ColorConvert.srgbEncodeFast(bb)
                            if (av > 0f) {
                                sr /= av; sg /= av; sb /= av
                            }
                            out[to + x] = (to8(av) shl 24) or (to8(sr) shl 16) or (to8(sg) shl 8) or to8(sb)
                        }
                        continue
                    }
                    // XYB から線形 RGB への変換 (分岐のない計算) と、8 ビットへの変換 (表を引く) を別のループにする
                    if (grey) {
                        for (x in la until lz) {
                            val i = ro + x
                            val gl = ry[i] + rx[i] - cb0
                            val gm = ry[i] - rx[i] - cb1
                            val gs = rb[i] - cb2
                            val ml = (gl * gl * gl + ob0) * its
                            val mm = (gm * gm * gm + ob1) * its
                            val ms = (gs * gs * gs + ob2) * its
                            val r = m0 * ml + m1 * mm + m2 * ms
                            val g = m3 * ml + m4 * mm + m5 * ms
                            val bb = m6 * ml + m7 * mm + m8 * ms
                            val l = 0.2126f * r + 0.7152f * g + 0.0722f * bb
                            lin[x] = l
                            lin[gOff + x] = l
                            lin[bOff + x] = l
                        }
                    } else {
                        for (x in la until lz) {
                            val i = ro + x
                            val gl = ry[i] + rx[i] - cb0
                            val gm = ry[i] - rx[i] - cb1
                            val gs = rb[i] - cb2
                            val ml = (gl * gl * gl + ob0) * its
                            val mm = (gm * gm * gm + ob1) * its
                            val ms = (gs * gs * gs + ob2) * its
                            lin[x] = m0 * ml + m1 * mm + m2 * ms
                            lin[gOff + x] = m3 * ml + m4 * mm + m5 * ms
                            lin[bOff + x] = m6 * ml + m7 * mm + m8 * ms
                        }
                    }
                    if (alpha == null) {
                        for (x in la until lz) {
                            out[to + x] = (255 shl 24) or (ColorConvert.linearToSrgb8(lin[x]) shl 16) or
                                (ColorConvert.linearToSrgb8(lin[gOff + x]) shl 8) or ColorConvert.linearToSrgb8(lin[bOff + x])
                        }
                    } else {
                        for (x in la until lz) {
                            val raw = alpha.data[oo + x]
                            var av = if (alphaDepth!!.isFloat) alphaDepth.parseFloatSample(raw) else raw * alphaScale
                            if (av < 0f) av = 0f else if (av > 1f) av = 1f
                            out[to + x] = (to8(av) shl 24) or (ColorConvert.linearToSrgb8(lin[x]) shl 16) or
                                (ColorConvert.linearToSrgb8(lin[gOff + x]) shl 8) or ColorConvert.linearToSrgb8(lin[bOff + x])
                        }
                    }
                }
            } finally {
                ScratchPool.give(lin)
                ep.recycle()
            }
        }

        /** グループ gx とその左のグループの継ぎ目の配列 (無ければ作る) */
        fun seamOf(gx: Int): Array<FloatArray> {
            seams.get(gx)?.let { return it }
            val s = Array(3) { ScratchPool.floats(seamW * (gd + BORDER)) }
            if (seams.compareAndSet(gx, null, s)) return s
            for (c in 0 until 3) ScratchPool.give(s[c])
            return seams.get(gx)
        }

        /** グループ行 gy の、グループ gx とその左のグループの継ぎ目を描く (両方の復号が済んでいること) */
        fun renderSeam(gx: Int, gy: Int) {
            val s = seams.getAndSet(gx, null)
            try {
                val x0 = gx * gd
                val y0 = gy * gd
                val ay0 = if (gy == 0) 0 else y0 - SHIFT
                val ay1 = if (gy == rowsOfGroups - 1) h else minOf(h, y0 + gd - SHIFT)
                render(s, seamW, x0 - BORDER, y0 - BORDER, x0 - BORDER, minOf(w, x0 + BORDER), x0 - SHIFT, minOf(w, x0 + SHIFT), ay0, ay1)
            } finally {
                for (c in 0 until 3) ScratchPool.give(s[c])
            }
        }

        // グループ行 gy のグループ gx を復号し、継ぎ目以外を描く
        fun decodeAndRender(gx: Int, gy: Int) {
            val x0 = gx * gd
            val y0 = gy * gd
            val gw = minOf(gd, sw - x0)
            val gh = minOf(gd, bh * 8 - y0)
            val t = Array(3) { ScratchPool.floats(ts * (gd + BORDER)) }
            try {
                decodeGroup(gy * gpr + gx, t, ts, x0, y0 - BORDER)
                // 左右の継ぎ目に写す列 (上の余白の行を含む)。左の継ぎ目には自分の左端、右の継ぎ目には右端
                val left = if (gx > 0) seamOf(gx) else null
                val right = if (gx < gpr - 1) seamOf(gx + 1) else null
                for (c in 0 until 3) {
                    val d = t[c]
                    // 上の余白は前のグループ行の下端
                    if (gy > 0) {
                        val rs = rowStrip[(gy - 1) and 1][c]
                        for (r in 0 until BORDER) System.arraycopy(rs, r * sw + x0, d, r * ts, gw)
                    }
                    // 継ぎ目に写す列は gw 以上あるとは限らない (右端の狭いグループ) が、写した先の余分な列は使わない
                    if (left != null) {
                        val sc = left[c]
                        for (r in 0 until BORDER + gh) System.arraycopy(d, r * ts, sc, r * seamW + BORDER, BORDER)
                    }
                    if (right != null) {
                        val sc = right[c]
                        for (r in 0 until BORDER + gh) System.arraycopy(d, r * ts + gw - BORDER, sc, r * seamW, BORDER)
                    }
                    // 次のグループ行のために下端の BORDER 行を残す
                    if (gy < rowsOfGroups - 1) {
                        val rs = rowStrip[gy and 1][c]
                        for (r in 0 until BORDER) System.arraycopy(d, (gh + r) * ts, rs, r * sw + x0, gw)
                    }
                }
                // 左右の継ぎ目は、両側の復号を終えた方が描く
                if (gx > 0 && pending.decrementAndGet(gx) == 0) renderSeam(gx, gy)
                if (gx < gpr - 1 && pending.decrementAndGet(gx + 1) == 0) renderSeam(gx + 1, gy)
                val ax0 = if (gx == 0) 0 else x0 + SHIFT
                val ax1 = if (gx == gpr - 1) w else minOf(w, x0 + gd - SHIFT)
                val ay0 = if (gy == 0) 0 else y0 - SHIFT
                val ay1 = if (gy == rowsOfGroups - 1) h else minOf(h, y0 + gd - SHIFT)
                render(t, ts, x0, y0 - BORDER, maxOf(0, ax0 - SHIFT), minOf(w, ax1 + SHIFT), ax0, ax1, ay0, ay1)
            } finally {
                for (c in 0 until 3) ScratchPool.give(t[c])
            }
        }

        try {
            for (gy in 0 until rowsOfGroups) {
                for (gx in 0 until gpr) pending.set(gx, 2)
                parallelFor(executor, gpr) { gx -> decodeAndRender(gx, gy) }
            }
        } finally {
            // 失敗したときに残っている継ぎ目も返す
            for (gx in 0 until gpr) seams.getAndSet(gx, null)?.let { s -> for (c in 0 until 3) ScratchPool.give(s[c]) }
            for (c in 0 until 3) {
                ScratchPool.give(rowStrip[0][c])
                ScratchPool.give(rowStrip[1][c])
            }
        }
        if (lateAlpha) {
            val chans = lfGlobal.gmodular!!.finish()
            if (alphaIdx >= 0) {
                val a = chans[alphaIdx]
                if (a.width != w || a.height != h) fail("modular channel size mismatch")
                val ad = alphaDepth!!
                parallelRows(executor, h) { y0, y1 ->
                    for (y in y0 until y1) {
                        val ty = y + oy
                        if (ty < 0 || ty >= outH) continue
                        val ao = a.offset + y * a.stride
                        val to = ty * outW + ox
                        for (x in xa until xb) {
                            val raw = a.data[ao + x]
                            var av = if (ad.isFloat) ad.parseFloatSample(raw) else raw * alphaScale
                            if (av < 0f) av = 0f else if (av > 1f) av = 1f
                            out[to + x] = (out[to + x] and 0xffffff) or (to8(av) shl 24)
                        }
                    }
                }
            }
        }
        return out
    }

    /**
     * 行 y0..y1 の帯の、列 xMin..xMax にパッチを重ねる (追加チャンネルのない画像)。
     * 帯の配列は、行 lo が先頭で行の長さ rs、画像の列 x は (x - xOff + 3) にある
     */
    private fun applyPatchesToBand(
        patches: Patches, rows: Array<FloatArray>, rs: Int, lo: Int, y0: Int, y1: Int, xMin: Int, xMax: Int, xOff: Int, h: Int,
    ) {
        for (p in patches.refs) {
            val ref = refs[p.refIdx]?.image ?: fail("patch reference not available")
            if (ref.color.size != 3) fail("patch reference channel mismatch")
            for (t in p.targets) {
                val bi = t.blending[0]
                if (bi.mode == 0) continue
                val tx0 = maxOf(t.x, xMin)
                val ty0 = maxOf(t.y, y0)
                var tx1 = minOf(t.x + p.width, xMax)
                var ty1 = minOf(t.y + p.height, minOf(y1, h))
                val rx0 = p.x0 + (tx0 - t.x)
                val ry0 = p.y0 + (ty0 - t.y)
                tx1 = minOf(tx1, tx0 + (ref.width - rx0))
                ty1 = minOf(ty1, ty0 + (ref.height - ry0))
                if (tx1 <= tx0 || ty1 <= ty0 || rx0 < 0 || ry0 < 0) continue
                val mode = when (bi.mode) {
                    1 -> Features.BLEND_REPLACE
                    2 -> Features.BLEND_ADD
                    3 -> Features.BLEND_MUL
                    4, 5 -> Features.BLEND_ALPHA
                    else -> Features.BLEND_MULADD
                }
                for (c in 0 until 3) {
                    Features.blend(
                        mode, bi.clamp, bi.mode == 5 || bi.mode == 7, false,
                        rows[c], rs, tx0 - xOff + 3, ty0 - lo, null,
                        ref.color[c], ref.width, rx0, ry0, null,
                        tx1 - tx0, ty1 - ty0,
                    )
                }
            }
        }
    }

    private fun channelToPlane(ch: Channel, depth: BitDepth): Plane {
        val p = Plane(ch.width, ch.height)
        val d = ch.data
        val o = p.data
        if (depth.isFloat) {
            for (i in d.indices) o[i] = depth.parseFloatSample(d[i])
        } else {
            val s = depth.scale()
            for (i in d.indices) o[i] = d[i] * s
        }
        return p
    }

    companion object {
        val TIMING = System.getenv("JXL_TIMING") != null

        /** 逐次描画で、各グループが描く範囲を左上にずらす画素数 (フィルタの halo 以上で 8 の倍数) */
        private const val SHIFT = 8
        /** 逐次描画で、左・上のグループから受け取る余白の画素数 */
        private const val BORDER = 2 * SHIFT

        fun to8(v: Float): Int {
            val f = v * 255f + 0.5f
            return if (f <= 0f) 0 else if (f >= 255f) 255 else f.toInt()
        }
        private val COEFF_FREQ_CONTEXT = intArrayOf(
            0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 15, 16, 16, 17, 17, 18, 18, 19, 19,
            20, 20, 21, 21, 22, 22, 23, 23, 23, 23, 24, 24, 24, 24, 25, 25, 25, 25, 26, 26, 26, 26, 27,
            27, 27, 27, 28, 28, 28, 28, 29, 29, 29, 29, 30, 30, 30, 30,
        )
        private val COEFF_NUM_NONZERO_CONTEXT = intArrayOf(
            0, 31, 62, 62, 93, 93, 93, 93, 123, 123, 123, 123, 152, 152, 152, 152, 152, 152, 152, 152,
            180, 180, 180, 180, 180, 180, 180, 180, 180, 180, 180, 180, 206, 206, 206, 206, 206, 206,
            206, 206, 206, 206, 206, 206, 206, 206, 206, 206, 206, 206, 206, 206, 206, 206, 206, 206,
            206, 206, 206, 206, 206, 206, 206,
        )
    }
}
