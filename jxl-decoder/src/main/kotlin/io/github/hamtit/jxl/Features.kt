package io.github.hamtit.jxl

internal class PatchBlending(val mode: Int, val alphaChannel: Int, val clamp: Boolean)

internal class PatchTarget(val x: Int, val y: Int, val blending: Array<PatchBlending>)

internal class PatchRef(val refIdx: Int, val x0: Int, val y0: Int, val width: Int, val height: Int, val targets: List<PatchTarget>)

internal class Patches(val refs: List<PatchRef>) {
    companion object {
        fun parse(br: BitReader, image: ImageHeader, header: FrameHeader): Patches {
            val numExtra = image.extraChannels.size
            val alphaIdx = image.extraChannels.indices.filter { image.extraChannels[it].isAlpha }
            val code = EntropyCode.parse(br, 10)
            val r = code.newReader()
            r.begin(br)
            val maxRefs = minOf(1L shl 24, header.width.toLong() * header.height / 16)
            val numRefs = r.read(br, 0)
            if (numRefs > maxRefs) fail("too many patches")
            var total = 0L
            val refs = ArrayList<PatchRef>()
            for (i in 0 until numRefs) {
                val refIdx = r.read(br, 1)
                if (refIdx >= 4) fail("invalid patch reference")
                val x0 = r.read(br, 3)
                val y0 = r.read(br, 3)
                val w = r.read(br, 2) + 1
                val h = r.read(br, 2) + 1
                val count = r.read(br, 7) + 1
                total += count
                if (total > maxRefs * 4) fail("too many patches")
                val targets = ArrayList<PatchTarget>()
                var px = 0
                var py = 0
                for (k in 0 until count) {
                    val x: Int
                    val y: Int
                    if (k == 0) {
                        x = r.read(br, 4)
                        y = r.read(br, 4)
                    } else {
                        x = px + unpackSigned(r.read(br, 6))
                        y = py + unpackSigned(r.read(br, 6))
                    }
                    px = x
                    py = y
                    val blending = Array(numExtra + 1) {
                        val mode = r.read(br, 5)
                        if (mode > 7) fail("invalid patch blend mode")
                        val alpha = if (mode >= 4 && alphaIdx.size >= 2) r.read(br, 8) else alphaIdx.firstOrNull() ?: 0
                        val clamp = if (mode >= 3) r.read(br, 9) != 0 else false
                        if (mode >= 4 && alpha >= numExtra) fail("invalid patch alpha channel")
                        PatchBlending(mode, alpha, clamp)
                    }
                    targets.add(PatchTarget(x, y, blending))
                }
                refs.add(PatchRef(refIdx, x0, y0, w, h, targets))
            }
            r.finish()
            return Patches(refs)
        }
    }
}

internal class QuantSpline(val points: Array<LongArray>, val xybDct: Array<IntArray>, val sigmaDct: IntArray)

internal class Splines(val splines: List<QuantSpline>, val quantAdjust: Int) {
    companion object {
        fun parse(br: BitReader, header: FrameHeader): Splines {
            val code = EntropyCode.parse(br, 6)
            val r = code.newReader()
            r.begin(br)
            val numPixels = header.width.toLong() * header.height
            val num = r.read(br, 2)
            if (num.toLong() >= minOf(1L shl 24, numPixels / 4)) fail("too many splines")
            val count = num + 1
            val starts = Array(count) { LongArray(2) }
            var sx = r.read(br, 1).toLong()
            var sy = r.read(br, 1).toLong()
            starts[0][0] = sx
            starts[0][1] = sy
            for (i in 1 until count) {
                sx += unpackSigned(r.read(br, 1))
                sy += unpackSigned(r.read(br, 1))
                starts[i][0] = sx
                starts[i][1] = sy
            }
            val quantAdjust = unpackSigned(r.read(br, 0))
            val list = ArrayList<QuantSpline>()
            var acc = 0L
            for (i in 0 until count) {
                val numPoints = r.read(br, 3)
                acc += numPoints
                if (acc > minOf(1L shl 20, numPixels / 2)) fail("too many spline points")
                val pts = ArrayList<LongArray>()
                var cx = starts[i][0]
                var cy = starts[i][1]
                var dx = 0L
                var dy = 0L
                pts.add(longArrayOf(cx, cy))
                for (k in 0 until numPoints) {
                    dx += unpackSigned(r.read(br, 4))
                    dy += unpackSigned(r.read(br, 4))
                    val nx = cx + dx
                    val ny = cy + dy
                    if (nx == cx && ny == cy) fail("duplicate spline control point")
                    cx = nx
                    cy = ny
                    pts.add(longArrayOf(cx, cy))
                }
                val xyb = Array(3) { IntArray(32) { unpackSigned(r.read(br, 5)) } }
                val sigma = IntArray(32) { unpackSigned(r.read(br, 5)) }
                list.add(QuantSpline(pts.toTypedArray(), xyb, sigma))
            }
            r.finish()
            return Splines(list, quantAdjust)
        }
    }
}

internal object Features {

    // ---------------- ブレンド ----------------

    /** パッチのブレンド。base と new は同じ平面配置 */
    fun applyPatches(patches: Patches, img: FrameImage, refs: Array<ReferenceFrame?>, image: ImageHeader) {
        val colorCount = img.color.size
        val numExtra = img.extra.size
        for (p in patches.refs) {
            val ref = refs[p.refIdx]?.image ?: fail("patch reference not available")
            if (ref.color.size != colorCount) fail("patch reference channel mismatch")
            for (t in p.targets) {
                for (ch in 0 until colorCount + numExtra) {
                    val bi = if (ch < colorCount) t.blending[0] else t.blending[1 + ch - colorCount]
                    if (bi.mode == 0) continue
                    // 範囲
                    var tx0 = t.x
                    var ty0 = t.y
                    var tx1 = t.x + p.width
                    var ty1 = t.y + p.height
                    tx0 = maxOf(tx0, 0); ty0 = maxOf(ty0, 0)
                    tx1 = minOf(tx1, img.width); ty1 = minOf(ty1, img.height)
                    var rx0 = p.x0 + (tx0 - t.x)
                    var ry0 = p.y0 + (ty0 - t.y)
                    var w = tx1 - tx0
                    var h = ty1 - ty0
                    w = minOf(w, ref.width - rx0)
                    h = minOf(h, ref.height - ry0)
                    if (w <= 0 || h <= 0 || rx0 < 0 || ry0 < 0) continue
                    val base = if (ch < colorCount) img.color[ch] else img.extra[ch - colorCount]
                    val src = if (ch < colorCount) ref.color[ch] else ref.extra[ch - colorCount]
                    val useAlpha = bi.mode >= 4
                    val alphaCh = colorCount + bi.alphaChannel
                    var mode: Int
                    var swapped = false
                    when (bi.mode) {
                        1 -> mode = BLEND_REPLACE
                        2 -> mode = BLEND_ADD
                        3 -> mode = BLEND_MUL
                        4, 5 -> {
                            swapped = bi.mode == 5
                            mode = if (ch == alphaCh) BLEND_MIX_ALPHA else BLEND_ALPHA
                        }
                        else -> {
                            swapped = bi.mode == 7
                            mode = if (ch == alphaCh) (if (swapped) BLEND_REPLACE else BLEND_SKIP) else BLEND_MULADD
                        }
                    }
                    var baseAlpha: FloatArray? = null
                    var newAlpha: FloatArray? = null
                    var premul = false
                    if (useAlpha && ch != alphaCh) {
                        baseAlpha = img.extra[bi.alphaChannel]
                        newAlpha = ref.extra[bi.alphaChannel]
                        premul = image.extraChannels[bi.alphaChannel].alphaAssociated
                    }
                    blend(
                        mode, bi.clamp, swapped, premul,
                        base, img.width, tx0, ty0, baseAlpha,
                        src, ref.width, rx0, ry0, newAlpha,
                        w, h,
                    )
                }
            }
        }
    }

    const val BLEND_REPLACE = 0
    const val BLEND_ADD = 1
    const val BLEND_MUL = 2
    const val BLEND_ALPHA = 3
    const val BLEND_MULADD = 4
    const val BLEND_MIX_ALPHA = 5
    const val BLEND_SKIP = 6

    fun blend(
        mode: Int, clamp: Boolean, swapped: Boolean, premultiplied: Boolean,
        base: FloatArray, baseW: Int, bx: Int, by: Int, baseAlpha: FloatArray?,
        new: FloatArray, newW: Int, nx: Int, ny: Int, newAlpha: FloatArray?,
        w: Int, h: Int,
    ) {
        if (mode == BLEND_SKIP) return
        for (dy in 0 until h) {
            val bRow = (by + dy) * baseW + bx
            val nRow = (ny + dy) * newW + nx
            when {
                mode == BLEND_REPLACE || (mode == BLEND_ALPHA && newAlpha == null) ->
                    System.arraycopy(new, nRow, base, bRow, w)
                mode == BLEND_ADD || (mode == BLEND_MULADD && newAlpha == null) ->
                    for (x in 0 until w) base[bRow + x] += new[nRow + x]
                mode == BLEND_MUL -> for (x in 0 until w) {
                    var n = new[nRow + x]
                    if (clamp) n = n.coerceIn(0f, 1f)
                    base[bRow + x] *= n
                }
                mode == BLEND_ALPHA -> for (x in 0 until w) {
                    val bs: Float
                    val ns: Float
                    val ba: Float
                    var na: Float
                    if (swapped) {
                        bs = new[nRow + x]
                        ns = base[bRow + x]
                        ba = newAlpha!![nRow + x]
                        na = baseAlpha?.get(bRow + x) ?: 0f
                    } else {
                        bs = base[bRow + x]
                        ns = new[nRow + x]
                        ba = baseAlpha?.get(bRow + x) ?: 0f
                        na = newAlpha!![nRow + x]
                    }
                    if (clamp) na = na.coerceIn(0f, 1f)
                    base[bRow + x] = if (premultiplied) {
                        ns + bs * (1f - na)
                    } else {
                        val mixed = 1f - (1f - na) * (1f - ba)
                        val r = if (mixed > 0f) 1f / mixed else 0f
                        (na * ns + ba * bs * (1f - na)) * r
                    }
                }
                mode == BLEND_MULADD -> for (x in 0 until w) {
                    val bs: Float
                    val ns: Float
                    var na: Float
                    if (swapped) {
                        bs = new[nRow + x]
                        ns = base[bRow + x]
                        na = baseAlpha?.get(bRow + x) ?: 0f
                    } else {
                        bs = base[bRow + x]
                        ns = new[nRow + x]
                        na = newAlpha!![nRow + x]
                    }
                    if (clamp) na = na.coerceIn(0f, 1f)
                    base[bRow + x] = bs + na * ns
                }
                mode == BLEND_MIX_ALPHA -> for (x in 0 until w) {
                    var b = base[bRow + x]
                    var n = new[nRow + x]
                    if (swapped) { val t = b; b = n; n = t }
                    if (clamp) n = n.coerceIn(0f, 1f)
                    base[bRow + x] = b + n * (1f - b)
                }
            }
        }
    }

    // ---------------- スプライン ----------------

    fun renderSplines(splines: Splines, img: FrameImage, corrX: Float, corrB: Float, upsampling: Int) {
        val qa = splines.quantAdjust.toFloat()
        val invQa = if (qa >= 0) 1f / (1f + qa / 8f) else 1f - qa / 8f
        val weights = floatArrayOf(0.0042f, 0.075f, 0.07f, 0.3333f)
        for (qs in splines.splines) {
            val xyb = Array(3) { c -> FloatArray(32) { qs.xybDct[c][it] * weights[c] * invQa } }
            for (i in 0 until 32) {
                xyb[0][i] += corrX * xyb[1][i]
                xyb[2][i] += corrB * xyb[1][i]
            }
            val sigma = FloatArray(32) { qs.sigmaDct[it] * weights[3] * invQa }
            val pts = Array(qs.points.size) { floatArrayOf(qs.points[it][0].toFloat(), qs.points[it][1].toFloat()) }
            val samples = splineSamples(pts)
            val n = samples.size / 3
            val arcLength = n - 2 + samples[(n - 1) * 3 + 2]
            for (i in 0 until n) {
                val px = samples[i * 3]
                val py = samples[i * 3 + 1]
                val len = samples[i * 3 + 2]
                val fromStart = minOf(1f, i / arcLength)
                val t = 31f * fromStart
                val sg = continuousIdct(sigma, t)
                val invSigma = 1f / sg
                val values = FloatArray(3) { continuousIdct(xyb[it], t) * len }
                val maxColor = maxOf(0.01f, maxOf(values[0], maxOf(values[1], values[2])))
                val maxDist = Math.sqrt(2.0 * (Math.log(10.0) * 3.0 + maxColor)).toFloat() * kotlin.math.abs(sg)
                val xb = maxOf(0, Math.floor((px - maxDist + 0.5f).toDouble()).toInt())
                val xe = minOf(img.width, Math.floor((px + maxDist + 1.5f).toDouble()).toInt())
                val yb = maxOf(0, Math.floor((py - maxDist + 0.5f).toDouble()).toInt())
                val ye = minOf(img.height, Math.floor((py + maxDist + 1.5f).toDouble()).toInt())
                val buf0 = img.color[0]
                val buf1 = img.color[1]
                val buf2 = img.color[2]
                val v0 = 0.25f * values[0] * sg
                val v1 = 0.25f * values[1] * sg
                val v2 = 0.25f * values[2] * sg
                val w = img.width
                for (y in yb until ye) {
                    val dy = y - py
                    val dy2 = dy * dy
                    val row = y * w
                    for (x in xb until xe) {
                        val dx = x - px
                        val d = Math.sqrt((dx * dx + dy2).toDouble()).toFloat()
                        val f = erf((0.5f * d + 0.35355338f) * invSigma) - erf((0.5f * d - 0.35355338f) * invSigma)
                        val f2 = f * f
                        val idx = row + x
                        buf0[idx] += v0 * f2
                        buf1[idx] += v1 * f2
                        buf2[idx] += v2 * f2
                    }
                }
            }
        }
        if (upsampling < 0) error("unreachable")
    }

    /** (x, y, 弧長) の並び */
    private fun splineSamples(points: Array<FloatArray>): FloatArray {
        val up = upsampledPoints(points)
        val out = ArrayList<Float>()
        var cx = up[0][0]
        var cy = up[0][1]
        out.add(cx); out.add(cy); out.add(1f)
        var next = 0
        while (next < up.size) {
            var prevX = cx
            var prevY = cy
            var arc = 0f
            while (true) {
                if (next >= up.size) {
                    out.add(prevX); out.add(prevY); out.add(arc)
                    break
                }
                val nx = up[next][0]
                val ny = up[next][1]
                val toNext = Math.sqrt(((nx - prevX) * (nx - prevX) + (ny - prevY) * (ny - prevY)).toDouble()).toFloat()
                if (arc + toNext >= 1f) {
                    val f = (1f - arc) / toNext
                    cx = prevX + (nx - prevX) * f
                    cy = prevY + (ny - prevY) * f
                    out.add(cx); out.add(cy); out.add(1f)
                    break
                }
                arc += toNext
                prevX = nx
                prevY = ny
                next++
            }
        }
        return out.toFloatArray()
    }

    private fun upsampledPoints(s: Array<FloatArray>): Array<FloatArray> {
        if (s.size == 1) return arrayOf(s[0])
        val ext = ArrayList<FloatArray>()
        ext.add(floatArrayOf(2 * s[0][0] - s[1][0], 2 * s[0][1] - s[1][1]))
        ext.addAll(s)
        val l = s.size
        ext.add(floatArrayOf(2 * s[l - 1][0] - s[l - 2][0], 2 * s[l - 1][1] - s[l - 2][1]))
        val out = ArrayList<FloatArray>()
        for (i in 0 until ext.size - 3) {
            val p = Array(4) { ext[i + it] }
            out.add(p[1])
            val t = FloatArray(4)
            for (k in 1 until 4) {
                val dx = p[k][0] - p[k - 1][0]
                val dy = p[k][1] - p[k - 1][1]
                t[k] = t[k - 1] + Math.pow((dx * dx + dy * dy).toDouble(), 0.25).toFloat()
            }
            for (step in 1 until 16) {
                val knot = t[1] + (step / 16f) * (t[2] - t[1])
                val a = Array(3) { k ->
                    val f = (knot - t[k]) / (t[k + 1] - t[k])
                    floatArrayOf(p[k][0] + (p[k + 1][0] - p[k][0]) * f, p[k][1] + (p[k + 1][1] - p[k][1]) * f)
                }
                val b = Array(2) { k ->
                    val f = (knot - t[k]) / (t[k + 2] - t[k])
                    floatArrayOf(a[k][0] + (a[k + 1][0] - a[k][0]) * f, a[k][1] + (a[k + 1][1] - a[k][1]) * f)
                }
                val f = (knot - t[1]) / (t[2] - t[1])
                out.add(floatArrayOf(b[0][0] + (b[1][0] - b[0][0]) * f, b[0][1] + (b[1][1] - b[0][1]) * f))
            }
        }
        out.add(s[l - 1])
        return out.toTypedArray()
    }

    private fun continuousIdct(dct: FloatArray, t: Float): Float {
        var r = dct[0]
        for (i in 1 until 32) {
            val theta = i * (Math.PI.toFloat() / 32f) * (t + 0.5f)
            r += 1.4142135f * dct[i] * Math.cos(theta.toDouble()).toFloat()
        }
        return r
    }

    private fun erf(x: Float): Float {
        val ax = kotlin.math.abs(x)
        val d1 = ax * 7.77394369e-02f + 2.05260015e-04f
        val d2 = d1 * ax + 2.32120216e-01f
        val d3 = d2 * ax + 2.77820801e-01f
        val d4 = d3 * ax + 1.0f
        val d5 = d4 * d4
        val inv = 1f / d5
        val r = -inv * inv + 1f
        return if (x < 0) -r else r
    }

    // ---------------- ノイズ ----------------

    fun renderNoise(lut8: FloatArray, img: FrameImage, header: FrameHeader, corrX: Float, corrB: Float, seed0: Long) {
        val w = img.width
        val h = img.height
        if (w <= 0 || h <= 0) return
        val gd = header.groupDim
        val gpr = (w + gd - 1) / gd
        val lut = FloatArray(9)
        System.arraycopy(lut8, 0, lut, 0, 8)
        lut[8] = lut8[7]
        val px = img.color[0]
        val py = img.color[1]
        val pb = img.color[2]

        // 5 行の環状バッファ (raw と水平 5 点和を 3 チャンネル分だけ持つ)
        val raw = Array(3) { Array(5) { FloatArray(w) } }
        val hsum = Array(3) { Array(5) { FloatArray(w) } }
        val rngs = Array(3) { arrayOfNulls<XorShift128Plus>(gpr) }
        val bits = IntArray(16)
        val maxStride = (minOf(gd, w) + 15) / 16 * 16
        val rowBuf = FloatArray(maxStride)
        var generated = 0

        fun generateRow(ry: Int) {
            if (ry % gd == 0) {
                val gh = minOf(gd, h - ry)
                for (gx in 0 until gpr) {
                    val x0 = gx * gd
                    val gw = minOf(gd, w - x0)
                    val stepsPerRow = (gw + 15) / 16
                    val steps = gh * stepsPerRow
                    val r0 = XorShift128Plus(seed0, (x0.toLong() shl 32) + ry)
                    rngs[0][gx] = r0.copy()
                    r0.skip(steps)
                    rngs[1][gx] = r0.copy()
                    r0.skip(steps)
                    rngs[2][gx] = r0
                }
            }
            val slot = ry % 5
            for (c in 0 until 3) {
                val dst = raw[c][slot]
                val cRngs = rngs[c]
                for (gx in 0 until gpr) {
                    val x0 = gx * gd
                    val gw = minOf(gd, w - x0)
                    val stride = (gw + 15) / 16 * 16
                    val rng = cRngs[gx]!!
                    var k = 0
                    while (k < stride) {
                        rng.fill(bits)
                        for (j in 0 until 16) rowBuf[k + j] = Float.fromBits((bits[j] ushr 9) or 0x3f800000)
                        k += 16
                    }
                    System.arraycopy(rowBuf, 0, dst, x0, gw)
                }
                // 水平 5 点和 (境界は鏡映)
                val hs = hsum[c][slot]
                if (w >= 5) {
                    hs[0] = dst[1] + dst[0] + dst[0] + dst[1] + dst[2]
                    hs[1] = dst[0] + dst[0] + dst[1] + dst[2] + dst[3]
                    for (x in 2 until w - 2) {
                        hs[x] = dst[x - 2] + dst[x - 1] + dst[x] + dst[x + 1] + dst[x + 2]
                    }
                    hs[w - 2] = dst[w - 4] + dst[w - 3] + dst[w - 2] + dst[w - 1] + dst[w - 1]
                    hs[w - 1] = dst[w - 3] + dst[w - 2] + dst[w - 1] + dst[w - 1] + dst[w - 2]
                } else {
                    for (x in 0 until w) {
                        var s = 0f
                        for (dx in -2..2) s += dst[Filters.mirror(x + dx, w)]
                        hs[x] = s
                    }
                }
            }
        }

        for (y in 0 until h) {
            val y0 = Filters.mirror(y - 2, h)
            val y1 = Filters.mirror(y - 1, h)
            val y2 = y
            val y3 = Filters.mirror(y + 1, h)
            val y4 = Filters.mirror(y + 2, h)
            val need = maxOf(y0, maxOf(y1, maxOf(y2, maxOf(y3, y4))))
            while (generated <= need) {
                generateRow(generated)
                generated++
            }
            val s0 = y0 % 5
            val s1 = y1 % 5
            val s2 = y2 % 5
            val s3 = y3 % 5
            val s4 = y4 % 5
            val h00 = hsum[0][s0]; val h01 = hsum[0][s1]; val h02 = hsum[0][s2]; val h03 = hsum[0][s3]; val h04 = hsum[0][s4]
            val h10 = hsum[1][s0]; val h11 = hsum[1][s1]; val h12 = hsum[1][s2]; val h13 = hsum[1][s3]; val h14 = hsum[1][s4]
            val h20 = hsum[2][s0]; val h21 = hsum[2][s1]; val h22 = hsum[2][s2]; val h23 = hsum[2][s3]; val h24 = hsum[2][s4]
            val r0 = raw[0][s2]
            val r1 = raw[1][s2]
            val r2 = raw[2][s2]
            val rowOff = y * w
            for (x in 0 until w) {
                val i = rowOff + x
                val c0 = (h00[x] + h01[x] + h02[x] + h03[x] + h04[x]) * 0.16f - r0[x] * 4f
                val c1 = (h10[x] + h11[x] + h12[x] + h13[x] + h14[x]) * 0.16f - r1[x] * 4f
                val c2 = (h20[x] + h21[x] + h22[x] + h23[x] + h24[x]) * 0.16f - r2[x] * 4f
                val gxv = px[i]
                val gyv = py[i]
                val inX = gxv + gyv
                val inY = gyv - gxv
                val sx3 = maxOf(0f, inX * 3f)
                val sy3 = maxOf(0f, inY * 3f)
                val xi = minOf(sx3.toInt(), 7)
                val xf = sx3 - xi
                val yi = minOf(sy3.toInt(), 7)
                val yf = sy3 - yi
                val sx = (lut[xi + 1] - lut[xi]) * xf + lut[xi]
                val sy = (lut[yi + 1] - lut[yi]) * yf + lut[yi]
                val nx = 0.22f * sx * (0.0078125f * c0 + 0.9921875f * c2)
                val ny = 0.22f * sy * (0.0078125f * c1 + 0.9921875f * c2)
                px[i] = gxv + corrX * (nx + ny) + nx - ny
                py[i] = gyv + nx + ny
                pb[i] += corrB * (nx + ny)
            }
        }
    }

    private class XorShift128Plus {
        private val s0: LongArray
        private val s1: LongArray

        constructor(seed0: Long, seed1: Long) {
            s0 = LongArray(8)
            s1 = LongArray(8)
            s0[0] = splitMix(seed0 + -0x61c8864680b583ebL)
            s1[0] = splitMix(seed1 + -0x61c8864680b583ebL)
            for (i in 1 until 8) {
                s0[i] = splitMix(s0[i - 1])
                s1[i] = splitMix(s1[i - 1])
            }
        }

        private constructor(s0: LongArray, s1: LongArray) {
            this.s0 = s0
            this.s1 = s1
        }

        fun copy(): XorShift128Plus = XorShift128Plus(s0.copyOf(), s1.copyOf())

        private fun splitMix(zIn: Long): Long {
            var z = zIn
            z = (z xor (z ushr 30)) * -0x40a7b892e31b1a47L
            z = (z xor (z ushr 27)) * -0x6b2fb644ecceee15L
            return z xor (z ushr 31)
        }

        fun skip(steps: Int) {
            for (k in 0 until steps) {
                for (i in 0 until 8) {
                    var a = s0[i]
                    val b = s1[i]
                    s0[i] = b
                    a = a xor (a shl 23)
                    s1[i] = a xor (b xor (a ushr 18) xor (b ushr 5))
                }
            }
        }

        fun fill(out: IntArray) {
            for (i in 0 until 8) {
                var a = s0[i]
                val b = s1[i]
                val r = a + b
                s0[i] = b
                a = a xor (a shl 23)
                s1[i] = a xor (b xor (a ushr 18) xor (b ushr 5))
                out[2 * i] = r.toInt()
                out[2 * i + 1] = (r ushr 32).toInt()
            }
        }
    }
}
