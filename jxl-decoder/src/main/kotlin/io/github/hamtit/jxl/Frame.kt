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
    private var lfQuant: Array<IntArray>? = null
    private var lfQuantW = IntArray(3)
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
        phase("lfGlobal")

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
            decodeGroup(0, null, 0)
        } else {
            parallelFor(executor, numGroups) { g -> decodeGroup(g, null, 0) }
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
            gmod = ModularImage.read(br, chans, tree, image.bitDepth.bitsPerSample)
            val list = gmod.channels
            var end = 0
            while (end < list.size) {
                val c = list[end]
                if (end >= gmod.nbMetaChannels && (c.width > gd || c.height > gd)) break
                end++
            }
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
        if (!header.useLfFrame) {
            lfQuant = Array(3) { c -> IntArray((bw shr hShift[c]) * (bh shr vShift[c])) }
            for (c in 0 until 3) lfQuantW[c] = bw shr hShift[c]
        }
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
                val q = lfQuant!![c]
                val bx0 = (x0 / 8) shr hShift[c]
                val by0 = (y0 / 8) shr vShift[c]
                for (y in 0 until ch.height) {
                    val py = by0 + y
                    if (py >= plane.height) break
                    for (x in 0 until ch.width) {
                        val px = bx0 + x
                        if (px >= plane.width) break
                        val v = ch.data[y * ch.width + x]
                        plane.data[py * plane.width + px] = v * scale
                        q[py * plane.width + px] = v
                    }
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
        val src = Array(3) { lfPlanes[it].data.copyOf() }
        val dst = Array(3) { lfPlanes[it].data }
        val mul = floatArrayOf(lfx, lfy, lfb)
        val wa = FloatArray(3)
        val self = FloatArray(3)
        for (y in 1 until h - 1) {
            for (x in 1 until w - 1) {
                var gap = 0.5f
                for (c in 0 until 3) {
                    val s = src[c]
                    val i = y * w + x
                    val v = s[i]
                    val side = s[i - 1] + s[i + 1] + s[i - w] + s[i + w]
                    val diag = s[i - w - 1] + s[i - w + 1] + s[i + w - 1] + s[i + w + 1]
                    val a = v * sSelf + side * sSide + diag * sDiag
                    wa[c] = a
                    self[c] = v
                    val g = kotlin.math.abs(a - v) / mul[c]
                    if (g > gap) gap = g
                }
                val gs = maxOf(0f, 3f - 4f * gap)
                for (c in 0 until 3) dst[c][y * w + x] = (wa[c] - self[c]) * gs + self[c]
            }
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

    /** グループを復号する。target が null なら colorPlanes に、そうでなければ行 targetY0 から始まる帯に書く */
    private fun decodeGroup(g: Int, target: Array<FloatArray>?, targetY0: Int) {
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
            for (c in 0 until 3) {
                val n = ((gbw shr hShift[c]) * 8) * ((gbh shr vShift[c]) * 8)
                gs.coeffs[c].fill(0f, 0, n)
                gs.intCoeffs?.get(c)?.fill(0, 0, n)
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
                        val n = ((gbw shr hShift[c]) * 8) * ((gbh shr vShift[c]) * 8)
                        val src = ints[c]
                        val dst = gs.coeffs[c]
                        for (i in 0 until n) dst[i] = src[i].toFloat()
                    }
                }
                reconstructGroup(gx, gy, gbw, gbh, gs, target, targetY0)
            }
        } finally {
            gs?.recycle()
        }
    }

    /**
     * グループの係数と逆変換の作業領域。グループごと・フレームごとに確保せず、[ScratchPool] に置いて使い回す。
     * パスが 1 つなら、量子化された係数 (整数) を浮動小数点の配列に直接書き、その場で逆量子化する
     * (各係数は一度だけ書かれるので、整数から変換するのと同じ値になる)。
     * 複数のパスがあるときは [intCoeffs] に整数で足し合わせる
     */
    private class GroupScratch(size: Int, multiPass: Boolean) {
        val coeffs = Array(3) { FloatArray(size) }
        val intCoeffs: Array<IntArray>? = if (multiPass) Array(3) { IntArray(size) } else null
        val scratch = FloatArray(512)
        val tmp = FloatArray(64)
        val lfBuf = FloatArray(32 * 32)
        val dct = Dct.Scratch()

        fun recycle() = ScratchPool.give(this, coeffs[0].size * (if (intCoeffs != null) 24L else 12L) + 8192L)

        companion object {
            fun obtain(size: Int, multiPass: Boolean): GroupScratch =
                ScratchPool.take(GroupScratch::class.java) {
                    // 大きすぎるもの (大きいグループ用) を小さい画像に使わない
                    it.coeffs[0].size >= size && it.coeffs[0].size / 4 <= size && (!multiPass || it.intCoeffs != null)
                } ?: GroupScratch(size, multiPass)
        }
    }

    private fun decodeHfCoefficients(br: BitReader, pass: Int, gx: Int, gy: Int, gbw: Int, gbh: Int, gs: GroupScratch) {
        val outF = gs.coeffs
        val outI = gs.intCoeffs
        val ctx = lfGlobal.blockCtx!!
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
        val lfIdxMul = (ctx.lfThresholds[0].size + 1) * (ctx.lfThresholds[1].size + 1) * (ctx.lfThresholds[2].size + 1)
        val hfIdxMul = ctx.qfThresholds.size + 1
        val nzRow = Array(3) { c -> IntArray(maxOf(1, (gbw + (1 shl hShift[c]) - 1) shr hShift[c])) }
        val bx0 = gx * (gd / 8)
        val by0 = gy * (gd / 8)
        val cOrder = intArrayOf(1, 0, 2)
        val lfOrder = intArrayOf(0, 2, 1)
        val lq = lfQuant
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
                var lfIdx = 0
                if (lq != null) {
                    for (c in lfOrder) {
                        val th = ctx.lfThresholds[c]
                        lfIdx *= th.size + 1
                        val qx = (bx0 + x) shr hShift[c]
                        val qy = (by0 + y) shr vShift[c]
                        val q = lq[c][qy * lfQuantW[c] + qx]
                        for (v in th) if (q > v) lfIdx++
                    }
                }
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
                    var nonZeros = reader.readClustered(br, clusterMap[ctxBase + nzCtx], 0)
                    if (nonZeros > (63 shl logNb)) fail("too many non-zero coefficients")
                    val nzVal = (nonZeros + numBlocks - 1) shr logNb
                    for (dx in 0 until w8) {
                        if (sx + dx < nzr.size) nzr[sx + dx] = nzVal
                    }
                    if (nonZeros == 0) continue
                    val coeffF = outF[c]
                    val coeffI = outI?.get(c)
                    val cw = (gbw shr hs) * 8
                    var prevNz = if (nonZeros <= numBlocks * 4) 1 else 0
                    val order = hfPass.orders[orderId][c]!!
                    val ow = Tt.ORDER_W[orderId]
                    val transpose = Tt.needTranspose(t)
                    val coeffBase = ctxBase + blockCtx * 458 + 37 * ctx.numClusters
                    val size = order.size
                    var k = numBlocks
                    while (k < size) {
                        val i = k - numBlocks
                        val nzc = (nonZeros - 1) shr logNb
                        val cctx = (COEFF_NUM_NONZERO_CONTEXT[nzc] + COEFF_FREQ_CONTEXT[i shr logNb]) * 2 + prevNz
                        if (cctx >= 458) fail("invalid coefficient context")
                        val u = reader.readClustered(br, clusterMap[coeffBase + cctx], 0)
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
                        val ci = (sy * 8 + dy) * cw + sx * 8 + dx
                        if (coeffI != null) coeffI[ci] += v else coeffF[ci] = v.toFloat()
                        prevNz = 1
                        nonZeros--
                        if (nonZeros == 0) break
                        k++
                    }
                    if (nonZeros != 0) fail("too many zeros in varblock")
                }
            }
        }
        reader.finish()
    }

    /** 逆量子化・CfL・逆 DCT をして色平面に書き込む */
    private fun reconstructGroup(
        gx: Int, gy: Int, gbw: Int, gbh: Int, gs: GroupScratch, target: Array<FloatArray>?, targetY0: Int,
    ) {
        val coeffs = gs.coeffs
        val qm = quantMatrices!!
        val bx0 = gx * (gd / 8)
        val by0 = gy * (gd / 8)
        val qmScale = floatArrayOf(
            Math.pow(0.8, (header.xQmScale - 2).toDouble()).toFloat(), 1f,
            Math.pow(0.8, (header.bQmScale - 2).toDouble()).toFloat(),
        )
        val fcoeff = coeffs
        val qb = image.quantBias
        val qbn = image.quantBiasNumerator
        for (c in 0 until 3) {
            val cw = (gbw shr hShift[c]) * 8
            val dst = coeffs[c]
            forEachVarblock(c, gbw, gbh, bx0, by0) { sx, sy, t, mul ->
                val bwv = Tt.BW[t] * 8
                val bhv = Tt.BH[t] * 8
                val m = qm.get(t, c)
                val scale = 65536f / (lfGlobal.globalScale.toFloat() * mul) * qmScale[c]
                for (yy in 0 until bhv) {
                    val row = (sy * 8 + yy) * cw + sx * 8
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
            val fx = fcoeff[0]
            val fy = fcoeff[1]
            val fb = fcoeff[2]
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
                    val row = y * cw
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
            val cw = (gbw shr hShift[c]) * 8
            val buf = fcoeff[c]
            val lf = lfPlanes[c]
            val lbx0 = bx0 shr hShift[c]
            val lby0 = by0 shr vShift[c]
            forEachVarblock(c, gbw, gbh, bx0, by0) { sx, sy, t, _ ->
                val bwb = Tt.BW[t]
                val bhb = Tt.BH[t]
                val base = (sy * 8) * cw + sx * 8
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
                    else -> Dct.dct2d(buf, base, cw, bwb * 8, bhb * 8, scratch, false)
                }
            }
            // 色平面へ
            val px0 = (bx0 shr hShift[c]) * 8
            val py0 = (by0 shr vShift[c]) * 8
            val chh = (gbh shr vShift[c]) * 8
            if (target != null) {
                val stride = bw * 8
                for (y in 0 until chh) System.arraycopy(buf, y * cw, target[c], (py0 - targetY0 + y) * stride + px0, cw)
            } else {
                val plane = colorPlanes!![c]
                for (y in 0 until chh) {
                    System.arraycopy(buf, y * cw, plane.data, (py0 + y) * plane.width + px0, cw)
                }
            }
        }
    }

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
        val cc = header.encodedColorChannels
        val r = m[0].data
        val g = if (cc >= 3) m[1].data else r
        val b = if (cc >= 3) m[2].data else r
        if (r.size != n || g.size != n || b.size != n) fail("modular channel size mismatch")
        val s = image.bitDepth.scale()
        val alphaIdx = image.alphaChannel
        val alpha = if (alphaIdx >= 0) m[cc + alphaIdx].data else null
        if (alpha != null && alpha.size != n) fail("modular channel size mismatch")
        val alphaDepth = if (alphaIdx >= 0) image.extraChannels[alphaIdx].bitDepth else null
        val alphaScale = alphaDepth?.scale() ?: 0f
        val premul = alphaIdx >= 0 && image.extraChannels[alphaIdx].alphaAssociated
        val out = r
        // 色・アルファは浮動小数点の経路 (channelToPlane → toOutput) と同じ計算にする
        parallelRows(executor, colorH) { y0, y1 ->
            for (i in y0 * colorW until y1 * colorW) {
                var a = 1f
                if (alpha != null) {
                    a = if (alphaDepth!!.isFloat) alphaDepth.parseFloatSample(alpha[i]) else alpha[i] * alphaScale
                    a = a.coerceIn(0f, 1f)
                }
                var rv = r[i] * s
                var gv = g[i] * s
                var bv = b[i] * s
                if (premul && a > 0f) {
                    rv /= a; gv /= a; bv /= a
                }
                out[i] = (to8(a) shl 24) or (to8(rv) shl 16) or (to8(gv) shl 8) or to8(bv)
            }
        }
        for (c in m) if (c.data !== out) c.data = IntArray(0)
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
            if (g.header.transforms.isNotEmpty()) return false
            for ((i, ec) in image.extraChannels.withIndex()) {
                if (header.ecUpsampling[i] != 1 || ec.dimShift != 0) return false
            }
        }
        return true
    }

    private fun decodeStreaming(stream: StreamOutput): IntArray {
        val w = colorW
        val h = colorH
        val stride = bw * 8
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
        val epfSteps = when (epf?.iters) {
            null -> IntArray(0)
            1 -> intArrayOf(1)
            2 -> intArrayOf(1, 2)
            else -> intArrayOf(0, 1, 2)
        }
        val epfRadius = IntArray(epfSteps.size) { when (epfSteps[it]) { 0 -> 3; 1 -> 2; else -> 1 } }
        val epfHalo = epfRadius.sum()
        val halo = epfHalo + (if (gab != null) 1 else 0)
        val alphaIdx = image.alphaChannel
        val alpha = if (alphaIdx >= 0) lfGlobal.gmodular?.channels?.get(alphaIdx) else null
        val alphaDepth = if (alphaIdx >= 0) image.extraChannels[alphaIdx].bitDepth else null
        val premul = alphaIdx >= 0 && image.extraChannels[alphaIdx].alphaAssociated
        // 逆 DCT 後の値は、今のグループ行 (行 curY0 から gd 行) と、その直前の carryRows 行だけを持つ。
        // 仕上げ (finalize) が前のグループ行から読むのは、次のグループ行の始まりの 2 * halo + 2 行手前までなので足りる
        val carryRows = minOf(gd, 2 * halo + 2)
        val cur = Array(3) { ScratchPool.floats(stride * gd) }
        val carry = Array(3) { ScratchPool.floats(stride * carryRows) }
        var curY0 = 0

        // 画像の行 y (範囲内) の逆 DCT 後の値
        fun idctRow(c: Int, y: Int): FloatArray {
            if (y >= curY0) return cur[c]
            if (y < curY0 - carryRows) error("row not available")
            return carry[c]
        }
        fun idctOff(y: Int) = if (y >= curY0) (y - curY0) * stride else (y - (curY0 - carryRows)) * stride

        // Gaborish 後の行 (y は鏡映済みでなくてよい)
        val inputRow: (Int, Int, FloatArray, Int) -> Unit = { c, yy, dst, off ->
            val y = Filters.mirror(yy, h)
            val src = idctRow(c, y)
            val o = idctOff(y)
            if (gab == null) {
                System.arraycopy(src, o, dst, off, w)
            } else {
                val w0 = gab[2 * c]
                val w1 = gab[2 * c + 1]
                val norm = 1f / (1f + 4f * w0 + 4f * w1)
                val yu = maxOf(y - 1, 0)
                val yd = minOf(y + 1, h - 1)
                val su = idctRow(c, yu)
                val ou = idctOff(yu)
                val sd = idctRow(c, yd)
                val od = idctOff(yd)
                for (x in 0 until w) {
                    val xl = if (x > 0) x - 1 else 0
                    val xr = if (x + 1 < w) x + 1 else w - 1
                    val side = su[ou + x] + sd[od + x] + src[o + xl] + src[o + xr]
                    val diag = su[ou + xl] + su[ou + xr] + sd[od + xl] + sd[od + xr]
                    dst[off + x] = (src[o + x] + side * w0 + diag * w1) * norm
                }
            }
        }

        val ob = image.opsinBias
        val cb0 = Math.cbrt(ob[0].toDouble()).toFloat()
        val cb1 = Math.cbrt(ob[1].toDouble()).toFloat()
        val cb2 = Math.cbrt(ob[2].toDouble()).toFloat()
        val its = 255f / image.toneMapping.intensityTarget
        val m = image.opsinInverseMatrix
        val grey = image.colourEncoding.colourSpace == ColourSpace.GREY
        val alphaScale = alphaDepth?.scale() ?: 0f

        // 帯の作業領域は ScratchPool から取り、帯ごと・フレームごとに確保し直さない
        val threads = if (executor == null) 1 else Runtime.getRuntime().availableProcessors()
        val rs = w + 6
        val band = (BAND_PIXELS / w).coerceIn(16, 64)

        fun finalizeBand(y0: Int, y1: Int, sc: Filters.EpfScratch) {
            // 帯の行 y0..y1 の色 (XYB)。画素 (x, y) は rows[c][(y - lo) * rs + x + 3]
            val rows: Array<FloatArray>
            val lo: Int
            if (epf != null) {
                Filters.epfBand(inputRow, y0, y1, w, h, epfHalo, epfSteps, epfRadius, epf, epfSigma, bw, sc)
                rows = sc.result
                lo = sc.resultLo
            } else {
                sc.ensure((y1 - y0) * rs)
                rows = sc.a
                lo = y0
                for (c in 0 until 3) for (y in y0 until y1) inputRow(c, y, rows[c], (y - y0) * rs + 3)
            }
            lfGlobal.patches?.let { applyPatchesToBand(it, rows, rs, lo, y0, y1, w, h) }
            val rx = rows[0]
            val ry = rows[1]
            val rb = rows[2]
            for (y in y0 until y1) {
                val ty = y + oy
                if (ty < 0 || ty >= outH) continue
                val ro = (y - lo) * rs + 3
                val oo = y * w
                val to = ty * outW + ox
                for (x in xa until xb) {
                    val i = ro + x
                    val gl = ry[i] + rx[i] - cb0
                    val gm = ry[i] - rx[i] - cb1
                    val gs = rb[i] - cb2
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
                    var av = 1f
                    if (alpha != null) {
                        val raw = alpha.data[oo + x]
                        av = if (alphaDepth!!.isFloat) alphaDepth.parseFloatSample(raw) else raw * alphaScale
                        if (av < 0f) av = 0f else if (av > 1f) av = 1f
                    }
                    if (premul) {
                        var sr = ColorConvert.srgbEncodeFast(r)
                        var sg = ColorConvert.srgbEncodeFast(g)
                        var sb = ColorConvert.srgbEncodeFast(bb)
                        if (av > 0f) {
                            sr /= av; sg /= av; sb /= av
                        }
                        out[to + x] = (to8(av) shl 24) or (to8(sr) shl 16) or (to8(sg) shl 8) or to8(sb)
                    } else {
                        out[to + x] = (to8(av) shl 24) or (ColorConvert.linearToSrgb8(r) shl 16) or
                            (ColorConvert.linearToSrgb8(g) shl 8) or ColorConvert.linearToSrgb8(bb)
                    }
                }
            }
        }

        fun finalize(a: Int, b: Int) {
            if (a >= b) return
            // スレッドごとに連続した行をまとめて受け持ち、その中を最大 band 行の帯に分ける
            val tasks = maxOf(1, minOf(threads, (b - a + 15) / 16))
            val per = (b - a + tasks - 1) / tasks
            parallelFor(executor, tasks) { k ->
                val sc = Filters.EpfScratch.obtain()
                try {
                    val ta = a + k * per
                    val tb = minOf(b, ta + per)
                    var y0 = ta
                    while (y0 < tb) {
                        val y1 = minOf(tb, y0 + band)
                        finalizeBand(y0, y1, sc)
                        y0 = y1
                    }
                } finally {
                    sc.recycle()
                }
            }
        }

        var finalized = 0
        val rowsOfGroups = header.groupsPerColumn
        for (gy in 0 until rowsOfGroups) {
            if (gy > 0) {
                // 前のグループ行の終わりの carryRows 行を残してから、同じ領域に次の行を復号する
                for (c in 0 until 3) System.arraycopy(cur[c], (gd - carryRows) * stride, carry[c], 0, carryRows * stride)
                curY0 = gy * gd
            }
            parallelFor(executor, header.groupsPerRow) { gx ->
                decodeGroup(gy * header.groupsPerRow + gx, cur, gy * gd)
            }
            // 次のグループ行が要らない行まで仕上げる
            val ready = if (gy == rowsOfGroups - 1) h else minOf(h, (gy + 1) * gd - halo - 1)
            finalize(finalized, ready)
            finalized = maxOf(finalized, ready)
        }
        for (c in 0 until 3) {
            ScratchPool.give(cur[c])
            ScratchPool.give(carry[c])
        }
        return out
    }

    /** 行 y0..y1 の帯にパッチを重ねる (追加チャンネルのない画像) */
    private fun applyPatchesToBand(patches: Patches, rows: Array<FloatArray>, rs: Int, lo: Int, y0: Int, y1: Int, w: Int, h: Int) {
        for (p in patches.refs) {
            val ref = refs[p.refIdx]?.image ?: fail("patch reference not available")
            if (ref.color.size != 3) fail("patch reference channel mismatch")
            for (t in p.targets) {
                val bi = t.blending[0]
                if (bi.mode == 0) continue
                val tx0 = maxOf(t.x, 0)
                val ty0 = maxOf(t.y, y0)
                var tx1 = minOf(t.x + p.width, w)
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
                        rows[c], rs, tx0 + 3, ty0 - lo, null,
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

        /**
         * 逐次描画で一度に仕上げる帯の画素数の目安。帯が高いほど EPF の重なりの再計算が減るが、
         * スレッドごとの作業領域が増えるので、幅に応じて 16〜64 行にする
         */
        private const val BAND_PIXELS = 64 * 1024

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
