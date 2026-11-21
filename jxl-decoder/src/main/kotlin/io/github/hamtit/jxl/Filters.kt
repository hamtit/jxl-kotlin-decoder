package io.github.hamtit.jxl

import java.util.concurrent.ExecutorService

/** 帯の処理に入力の行を渡す (Kotlin の関数型だと Int がボックス化されるので専用の型にする) */
internal fun interface RowSource {
    /** チャンネル c の画像の行 y (範囲外なら鏡映した行) を dst[off until off + w] に書く */
    fun row(c: Int, y: Int, dst: FloatArray, off: Int)
}

internal object Filters {

    /** 境界を鏡映した添字 (-1 -> 0, -2 -> 1, len -> len-1) */
    fun mirror(v: Int, len: Int): Int {
        var o = v
        while (true) {
            o = if (o < 0) -(o + 1) else if (o >= len) 2 * len - 1 - o else return o
        }
    }

    /** 周囲 pad 画素を鏡映した平面を作る */
    private fun padded(src: FloatArray, w: Int, h: Int, pad: Int): FloatArray {
        val pw = w + 2 * pad
        val ph = h + 2 * pad
        val out = FloatArray(pw * ph)
        for (py in 0 until ph) {
            val sy = mirror(py - pad, h)
            val srow = sy * w
            val orow = py * pw
            System.arraycopy(src, srow, out, orow + pad, w)
            for (k in 0 until pad) {
                out[orow + pad - 1 - k] = src[srow + mirror(-1 - k, w)]
                out[orow + pad + w + k] = src[srow + mirror(w + k, w)]
            }
        }
        return out
    }

    /** Gaborish (3x3 の畳み込み)。各平面をその場で書き換える */
    fun gaborish(colors: Array<Plane>, weights: FloatArray, executor: ExecutorService?): Array<Plane> {
        parallelFor(executor, colors.size) { c ->
            val plane = colors[c]
            val w = plane.width
            val h = plane.height
            val w0 = weights[2 * c]
            val w1 = weights[2 * c + 1]
            val norm = 1f / (1f + 4f * w0 + 4f * w1)
            val p = plane.data
            // 上の行と今の行の元の値を残しておき、上から順に書き換える
            var up = FloatArray(w)
            var mid = FloatArray(w)
            val outRow = FloatArray(w)
            System.arraycopy(p, 0, mid, 0, w)
            System.arraycopy(p, 0, up, 0, w)
            for (y in 0 until h) {
                val r2 = minOf(y + 1, h - 1) * w
                for (x in 0 until w) {
                    val xl = if (x > 0) x - 1 else 0
                    val xr = if (x + 1 < w) x + 1 else w - 1
                    val dn = if (y + 1 < h) p[r2 + x] else mid[x]
                    val dnl = if (y + 1 < h) p[r2 + xl] else mid[xl]
                    val dnr = if (y + 1 < h) p[r2 + xr] else mid[xr]
                    val side = up[x] + dn + mid[xl] + mid[xr]
                    val diag = up[xl] + up[xr] + dnl + dnr
                    outRow[x] = (mid[x] + side * w0 + diag * w1) * norm
                }
                val t = up
                up = mid
                mid = t
                if (y + 1 < h) System.arraycopy(p, r2, mid, 0, w)
                System.arraycopy(outRow, 0, p, y * w, w)
            }
        }
        return colors
    }

    /**
     * Gaborish を 1 行分計算する。src の行 up・mid・dn (行の先頭の添字) から dst[off until off + w] に書く。
     * 左右の端は端の画素を繰り返す
     */
    fun gaborishRow(src: FloatArray, up: Int, mid: Int, dn: Int, dst: FloatArray, off: Int, w: Int, w0: Float, w1: Float) {
        val norm = 1f / (1f + 4f * w0 + 4f * w1)
        if (w == 1) {
            val side = src[up] + src[dn] + 2f * src[mid]
            val diag = 2f * (src[up] + src[dn])
            dst[off] = (src[mid] + side * w0 + diag * w1) * norm
            return
        }
        run {
            val side = src[up] + src[dn] + src[mid] + src[mid + 1]
            val diag = src[up] + src[up + 1] + src[dn] + src[dn + 1]
            dst[off] = (src[mid] + side * w0 + diag * w1) * norm
        }
        for (x in 1 until w - 1) {
            val side = src[up + x] + src[dn + x] + src[mid + x - 1] + src[mid + x + 1]
            val diag = src[up + x - 1] + src[up + x + 1] + src[dn + x - 1] + src[dn + x + 1]
            dst[off + x] = (src[mid + x] + side * w0 + diag * w1) * norm
        }
        run {
            val x = w - 1
            val side = src[up + x] + src[dn + x] + src[mid + x - 1] + src[mid + x]
            val diag = src[up + x - 1] + src[up + x] + src[dn + x - 1] + src[dn + x]
            dst[off + x] = (src[mid + x] + side * w0 + diag * w1) * norm
        }
    }

    private val KERNEL0 = intArrayOf(0, -2, -1, -1, 0, -1, 1, -1, -2, 0, -1, 0, 1, 0, 2, 0, -1, 1, 0, 1, 1, 1, 0, 2)
    private val DIST0 = intArrayOf(0, -1, 1, 0, 0, 0, -1, 0, 0, 1)

    /**
     * エッジ保存フィルタ。[sigma] は VarDCT のブロックごとの値 (行の長さ [sigmaStride])、
     * Modular では null (sigma_for_modular を使う)。
     * 全段階を 32 行の帯ごとにまとめて処理するので、作業領域は帯の分だけで済む。
     */
    fun epf(colors: Array<Plane>, params: EpfParams, sigma: FloatArray?, sigmaStride: Int, executor: ExecutorService?): Array<Plane> {
        val w = colors[0].width
        val h = colors[0].height
        val steps = when (params.iters) {
            1 -> intArrayOf(1)
            2 -> intArrayOf(1, 2)
            else -> intArrayOf(0, 1, 2)
        }
        val radius = IntArray(steps.size) { when (steps[it]) { 0 -> 3; 1 -> 2; else -> 1 } }
        val halo = radius.sum()
        val planes = Array(3) { colors[it].data }
        val band = 32
        // スレッドごとの縞 (帯の集まり)。縞の境目の外側の行は、書き換えられる前に退避しておく
        val threads = if (executor == null) 1 else Runtime.getRuntime().availableProcessors()
        val numBands = (h + band - 1) / band
        val stripes = maxOf(1, minOf(threads, numBands))
        val bandsPerStripe = (numBands + stripes - 1) / stripes
        val saved = Array(stripes) { s ->
            val ys = minOf(h, s * bandsPerStripe * band)
            val ye = minOf(h, (s + 1) * bandsPerStripe * band)
            // 上の halo 行と下の halo 行 (画像の外は鏡映)
            Array(3) { c ->
                val buf = FloatArray(2 * halo * w)
                for (k in 0 until halo) {
                    System.arraycopy(planes[c], mirror(ys - halo + k, h) * w, buf, k * w, w)
                    System.arraycopy(planes[c], mirror(ye + k, h) * w, buf, (halo + k) * w, w)
                }
                buf
            }
        }
        parallelFor(executor, stripes) { s ->
            val ys = minOf(h, s * bandsPerStripe * band)
            val ye = minOf(h, (s + 1) * bandsPerStripe * band)
            val sv = saved[s]
            val inputRow = RowSource { c, y, dst, off ->
                if (y < ys) System.arraycopy(sv[c], (y - (ys - halo)) * w, dst, off, w)
                else if (y >= ye) System.arraycopy(sv[c], (halo + y - ye) * w, dst, off, w)
                else System.arraycopy(planes[c], y * w, dst, off, w)
            }
            val scratch = Scratch()
            var pending: Array<FloatArray>? = null
            var pendingY0 = 0
            var pendingY1 = 0
            var y0 = ys
            var outSlot = EPF_OUT_A
            while (y0 < ye) {
                val y1 = minOf(ye, y0 + band)
                val out = epfBand(inputRow, y0, y1, w, h, halo, steps, radius, params, sigma, sigmaStride, scratch, outSlot)
                outSlot = if (outSlot == EPF_OUT_A) EPF_OUT_B else EPF_OUT_A
                // 前の帯の結果は、この帯の入力を読み終えてから書き戻す
                pending?.let { p -> for (c in 0 until 3) System.arraycopy(p[c], 0, planes[c], pendingY0 * w, (pendingY1 - pendingY0) * w) }
                pending = out
                pendingY0 = y0
                pendingY1 = y1
                y0 = y1
            }
            pending?.let { p -> for (c in 0 until 3) System.arraycopy(p[c], 0, planes[c], pendingY0 * w, (pendingY1 - pendingY0) * w) }
        }
        return colors
    }

    /** epfBand の結果を置く [Scratch] の枠 (2 組) */
    const val EPF_OUT_A = 10
    const val EPF_OUT_B = 13

    /**
     * 帯 (行 y0..y1) の全段階を計算する。結果は行の長さ w の 3 平面で、[scratch] の枠 outSlot から 3 つを使う
     * (配列は要求より長いことがある)
     */
    fun epfBand(
        inputRow: RowSource,
        y0: Int, y1: Int, w: Int, h: Int, halo: Int, steps: IntArray, radius: IntArray,
        params: EpfParams, sigma: FloatArray?, sigmaStride: Int, scratch: Scratch, outSlot: Int,
    ): Array<FloatArray> {
        val pw = w + 6
        var lo = y0 - halo
        var hi = y1 + halo
        // 枠 0..2 と 3..5 を交互に入出力に使う
        var curSlot = 0
        var cur = Array(3) { scratch.floats(curSlot + it, (hi - lo) * pw) }
        for (c in 0 until 3) {
            for (y in lo until hi) {
                val row = (y - lo) * pw
                inputRow.row(c, y, cur[c], row + 3)
                padRow(cur[c], row, w)
            }
        }
        for ((k, step) in steps.withIndex()) {
            val r = radius[k]
            val nlo = lo + r
            val nhi = hi - r
            val nextSlot = 3 - curSlot
            val next = Array(3) { scratch.floats(nextSlot + it, (nhi - nlo) * pw) }
            val ya = maxOf(nlo, 0)
            val yb = minOf(nhi, h)
            if (step == 0) epfStep0(cur, lo, next, nlo, ya, yb, w, params, sigma, sigmaStride)
            else epfStep12(cur, lo, hi - lo, next, nlo, ya, yb, w, params, sigma, sigmaStride, step, scratch)
            // 画像の外の行は、計算した行を鏡映して埋める
            for (y in nlo until nhi) {
                if (y in ya until yb) continue
                val sy = mirror(y, h)
                for (c in 0 until 3) System.arraycopy(next[c], (sy - nlo) * pw, next[c], (y - nlo) * pw, pw)
            }
            cur = next
            curSlot = nextSlot
            lo = nlo
            hi = nhi
        }
        return Array(3) { c ->
            val o = scratch.floats(outSlot + c, (y1 - y0) * w)
            for (y in y0 until y1) System.arraycopy(cur[c], (y - lo) * pw + 3, o, (y - y0) * w, w)
            o
        }
    }

    /** 行の左右 3 画素を鏡映で埋める (row は行の先頭、画素は row+3 から) */
    private fun padRow(buf: FloatArray, row: Int, w: Int) {
        for (k in 0 until 3) {
            buf[row + 2 - k] = buf[row + 3 + mirror(-1 - k, w)]
            buf[row + 3 + w + k] = buf[row + 3 + mirror(w + k, w)]
        }
    }

    private const val EPF_C0 = 6.6f * (0.70710677f - 1f)

    /**
     * 十字 (4 近傍) のカーネルの段階 1・2。画像の行 ya..yb を計算する。
     * inp は画像の行 inLo から、out は行 outLo から始まる帯 (行の長さ w + 6)。
     * 隣の画素との重み付き差分 H (右)・V (下) を先に求め、距離はその和で表す。
     */
    private fun epfStep12(
        inp: Array<FloatArray>, inLo: Int, inRows: Int, out: Array<FloatArray>, outLo: Int, ya: Int, yb: Int, w: Int,
        params: EpfParams, sigma: FloatArray?, sigmaStride: Int, step: Int, scratch: Scratch,
    ) {
        val pw = w + 6
        val stepMul = if (step == 2) params.pass2SigmaScale else 1f
        val border = params.borderSadMul
        val s0 = params.channelScale[0]
        val s1 = params.channelScale[1]
        val s2 = params.channelScale[2]
        val i0 = inp[0]
        val i1 = inp[1]
        val i2 = inp[2]
        val o0 = out[0]
        val o1 = out[1]
        val o2 = out[2]
        val plus = step == 1
        // 差分は画像の行 ya-2 .. yb+1 について計算する (局所行 0 = ya-2)
        val rows = (yb - ya) + 4
        val hd = scratch.floats(6, rows * pw)
        val vd = scratch.floats(7, rows * pw)
        for (r in 0 until rows) {
            val iy = ya - 2 + r - inLo
            val lb = r * pw
            if (iy < 0 || iy >= inRows) {
                java.util.Arrays.fill(hd, lb, lb + pw, 0f)
                java.util.Arrays.fill(vd, lb, lb + pw, 0f)
                continue
            }
            val base = iy * pw
            hd[lb + pw - 1] = 0f
            for (x in 0 until pw - 1) {
                val i = base + x
                hd[lb + x] = s0 * kotlin.math.abs(i0[i + 1] - i0[i]) + s1 * kotlin.math.abs(i1[i + 1] - i1[i]) +
                    s2 * kotlin.math.abs(i2[i + 1] - i2[i])
            }
            if (iy + 1 >= inRows) {
                java.util.Arrays.fill(vd, lb, lb + pw, 0f)
            } else {
                for (x in 0 until pw) {
                    val i = base + x
                    vd[lb + x] = s0 * kotlin.math.abs(i0[i + pw] - i0[i]) + s1 * kotlin.math.abs(i1[i + pw] - i1[i]) +
                        s2 * kotlin.math.abs(i2[i + pw] - i2[i])
                }
            }
        }
        var sh = hd
        var sv = vd
        if (plus) {
            // 読むのは行 1..rows-2・列 1..pw-3 だけ (それ以外は書かない)
            sh = scratch.floats(8, rows * pw)
            sv = scratch.floats(9, rows * pw)
            for (r in 1 until rows - 1) {
                val lb = r * pw
                for (x in 1 until pw - 2) {
                    val i = lb + x
                    sh[i] = hd[i] + hd[i - pw] + hd[i + pw] + hd[i - 1] + hd[i + 1]
                    sv[i] = vd[i] + vd[i - pw] + vd[i + pw] + vd[i - 1] + vd[i + 1]
                }
            }
        }
        for (y in ya until yb) {
            val yBorder = ((y + 1) and 6) == 0
            val lrow = (y - ya + 2) * pw
            val prow = (y - inLo) * pw
            val orow = (y - outLo) * pw
            var bx = 0
            while (bx < w) {
                val s = if (sigma != null) sigma[(y shr 3) * sigmaStride + (bx shr 3)] else params.sigmaForModular
                val xEnd = minOf(w, bx + 8)
                if (s < 0.3f) {
                    System.arraycopy(i0, prow + bx + 3, o0, orow + bx + 3, xEnd - bx)
                    System.arraycopy(i1, prow + bx + 3, o1, orow + bx + 3, xEnd - bx)
                    System.arraycopy(i2, prow + bx + 3, o2, orow + bx + 3, xEnd - bx)
                    bx = xEnd
                    continue
                }
                val nisBorder = EPF_C0 / s * stepMul * border
                val nisInner = if (yBorder) nisBorder else EPF_C0 / s * stepMul
                for (x in bx until xEnd) {
                    val pi = prow + x + 3
                    val xb = x and 7
                    val nis = if (xb == 0 || xb == 7) nisBorder else nisInner
                    val li = lrow + x + 3
                    var wr = 1f + sh[li] * nis
                    if (wr < 0f) wr = 0f
                    var wl = 1f + sh[li - 1] * nis
                    if (wl < 0f) wl = 0f
                    var wd = 1f + sv[li] * nis
                    if (wd < 0f) wd = 0f
                    var wu = 1f + sv[li - pw] * nis
                    if (wu < 0f) wu = 0f
                    val inv = 1f / (1f + wr + wl + wd + wu)
                    val oi = orow + x + 3
                    o0[oi] = (i0[pi] + wr * i0[pi + 1] + wl * i0[pi - 1] + wd * i0[pi + pw] + wu * i0[pi - pw]) * inv
                    o1[oi] = (i1[pi] + wr * i1[pi + 1] + wl * i1[pi - 1] + wd * i1[pi + pw] + wu * i1[pi - pw]) * inv
                    o2[oi] = (i2[pi] + wr * i2[pi + 1] + wl * i2[pi - 1] + wd * i2[pi + pw] + wu * i2[pi - pw]) * inv
                }
                bx = xEnd
            }
            for (c in 0 until 3) padRow(out[c], orow, w)
        }
    }

    /** 12 近傍のカーネルの段階 0 (iters = 3 のときだけ使う) */
    private fun epfStep0(
        inp: Array<FloatArray>, inLo: Int, out: Array<FloatArray>, outLo: Int, ya: Int, yb: Int, w: Int,
        params: EpfParams, sigma: FloatArray?, sigmaStride: Int,
    ) {
        val pw = w + 6
        val kernel = KERNEL0
        val dist = DIST0
        val stepMul = params.pass0SigmaScale
        val border = params.borderSadMul
        val scale0 = params.channelScale[0]
        val scale1 = params.channelScale[1]
        val scale2 = params.channelScale[2]
        val nk = kernel.size / 2
        val nd = dist.size / 2
        val kOff = IntArray(nk) { kernel[2 * it + 1] * pw + kernel[2 * it] }
        val dOff = IntArray(nd) { dist[2 * it + 1] * pw + dist[2 * it] }
        val i0 = inp[0]
        val i1 = inp[1]
        val i2 = inp[2]
        val o0 = out[0]
        val o1 = out[1]
        val o2 = out[2]
        for (y in ya until yb) {
            val yBorder = ((y + 1) and 6) == 0
            val prow = (y - inLo) * pw
            val orow = (y - outLo) * pw
            for (x in 0 until w) {
                val s = if (sigma != null) sigma[(y shr 3) * sigmaStride + (x shr 3)] else params.sigmaForModular
                val pi = prow + x + 3
                val oi = orow + x + 3
                if (s < 0.3f) {
                    o0[oi] = i0[pi]
                    o1[oi] = i1[pi]
                    o2[oi] = i2[pi]
                    continue
                }
                val xb = x and 7
                val sm = if (yBorder || xb == 0 || xb == 7) stepMul * border else stepMul
                val negInvSigma = EPF_C0 / s * sm
                var sumW = 1f
                var s0 = i0[pi]
                var s1 = i1[pi]
                var s2 = i2[pi]
                for (k in 0 until nk) {
                    val ki = pi + kOff[k]
                    var a0 = 0f
                    var a1 = 0f
                    var a2 = 0f
                    for (d in 0 until nd) {
                        val od = dOff[d]
                        a0 += kotlin.math.abs(i0[ki + od] - i0[pi + od])
                        a1 += kotlin.math.abs(i1[ki + od] - i1[pi + od])
                        a2 += kotlin.math.abs(i2[ki + od] - i2[pi + od])
                    }
                    val dd = scale0 * a0 + scale1 * a1 + scale2 * a2
                    var wt = 1f + dd * negInvSigma
                    if (wt < 0f) wt = 0f
                    sumW += wt
                    s0 += wt * i0[ki]
                    s1 += wt * i1[ki]
                    s2 += wt * i2[ki]
                }
                val inv = 1f / sumW
                o0[oi] = s0 * inv
                o1[oi] = s1 * inv
                o2[oi] = s2 * inv
            }
            for (c in 0 until 3) padRow(out[c], orow, w)
        }
    }

    /** 非分離型アップサンプリング (2^shift 倍)。出力は outW x outH に切り詰める */
    fun upsample(src: Plane, shift: Int, image: ImageHeader, outW: Int, outH: Int): Plane {
        var cur = src
        var s = shift
        while (s >= 3) {
            cur = upsampleInner(cur, 8, image.up8Weights)
            s -= 3
        }
        if (s == 1) cur = upsampleInner(cur, 2, image.up2Weights)
        if (s == 2) cur = upsampleInner(cur, 4, image.up4Weights)
        if (cur.width == outW && cur.height == outH) return cur
        val out = Plane(outW, outH)
        for (y in 0 until outH) {
            val sy = minOf(y, cur.height - 1)
            System.arraycopy(cur.data, sy * cur.width, out.data, y * outW, minOf(outW, cur.width))
        }
        return out
    }

    private fun upsampleInner(src: Plane, k: Int, weights: FloatArray): Plane {
        val w = src.width
        val h = src.height
        val pad = 2
        val pw = w + 4
        val p = padded(src.data, w, h, pad)
        val matN = k / 2
        val quarter = Array(matN * matN) { FloatArray(25) }
        var wi = 0
        for (y in 0 until 5 * matN) {
            val my = y / 5
            val ky = y % 5
            for (x in y until 5 * matN) {
                val mx = x / 5
                val kx = x % 5
                val v = weights[wi++]
                quarter[my * matN + mx][ky * 5 + kx] = v
                quarter[mx * matN + my][kx * 5 + ky] = v
            }
        }
        val ow = w * k
        val oh = h * k
        val out = Plane(ow, oh)
        val o = out.data
        for (y in 0 until oh) {
            val ry = y / k
            val myy = minOf(y % k, k - y % k - 1)
            val flipV = y % k >= matN
            for (x in 0 until ow) {
                val rx = x / k
                val mxx = minOf(x % k, k - x % k - 1)
                val flipH = x % k >= matN
                val kernel = quarter[myy * matN + mxx]
                var sum = 0f
                var mn = Float.POSITIVE_INFINITY
                var mxv = Float.NEGATIVE_INFINITY
                for (iy in 0 until 5) {
                    val kyy = if (flipV) 4 - iy else iy
                    val row = (ry + iy) * pw + rx
                    for (ix in 0 until 5) {
                        val kxx = if (flipH) 4 - ix else ix
                        val v = p[row + ix]
                        sum += kernel[kyy * 5 + kxx] * v
                        if (v < mn) mn = v
                        if (v > mxv) mxv = v
                    }
                }
                o[y * ow + x] = if (sum < mn) mn else if (sum > mxv) mxv else sum
            }
        }
        return out
    }

    /** JPEG のクロマアップサンプリング (3:1 の線形補間) */
    fun upsampleJpeg(src: Plane, hs: Int, vs: Int, outW: Int, outH: Int): Plane {
        var cur = src
        if (hs != 0) {
            val w = cur.width
            val h = cur.height
            val out = Plane(outW, h)
            for (y in 0 until h) {
                val r = y * w
                for (x in 0 until w) {
                    val c = cur.data[r + x]
                    val l = cur.data[r + maxOf(0, x - 1)]
                    val n = cur.data[r + minOf(w - 1, x + 1)]
                    if (2 * x < outW) out.data[y * outW + 2 * x] = 0.25f * l + 0.75f * c
                    if (2 * x + 1 < outW) out.data[y * outW + 2 * x + 1] = 0.75f * c + 0.25f * n
                }
            }
            cur = out
        } else if (cur.width != outW) {
            val out = Plane(outW, cur.height)
            for (y in 0 until cur.height) System.arraycopy(cur.data, y * cur.width, out.data, y * outW, outW)
            cur = out
        }
        if (vs != 0) {
            val w = cur.width
            val h = cur.height
            val out = Plane(w, outH)
            for (y in 0 until h) {
                val up = maxOf(0, y - 1)
                val dn = minOf(h - 1, y + 1)
                for (x in 0 until w) {
                    val c = cur.data[y * w + x]
                    if (2 * y < outH) out.data[(2 * y) * w + x] = 0.25f * cur.data[up * w + x] + 0.75f * c
                    if (2 * y + 1 < outH) out.data[(2 * y + 1) * w + x] = 0.75f * c + 0.25f * cur.data[dn * w + x]
                }
            }
            cur = out
        } else if (cur.height != outH) {
            val out = Plane(cur.width, outH)
            System.arraycopy(cur.data, 0, out.data, 0, cur.width * outH)
            cur = out
        }
        return cur
    }
}
