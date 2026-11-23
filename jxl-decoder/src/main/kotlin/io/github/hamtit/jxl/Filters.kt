package io.github.hamtit.jxl

import java.util.concurrent.ExecutorService

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
                // 下の行 (最後の行では自分の行の元の値)
                val dn = if (y + 1 < h) p else mid
                val dno = if (y + 1 < h) r2 else 0
                // 左右の端 (隣を鏡映する) だけ別に計算し、内側は分岐なしで回す
                for (edge in 0 until minOf(w, 2)) {
                    val x = edge * (w - 1)
                    val xl = if (x > 0) x - 1 else 0
                    val xr = if (x + 1 < w) x + 1 else w - 1
                    val side = up[x] + dn[dno + x] + mid[xl] + mid[xr]
                    val diag = up[xl] + up[xr] + dn[dno + xl] + dn[dno + xr]
                    outRow[x] = (mid[x] + side * w0 + diag * w1) * norm
                }
                for (x in 1 until w - 1) {
                    val side = up[x] + dn[dno + x] + mid[x - 1] + mid[x + 1]
                    val diag = up[x - 1] + up[x + 1] + dn[dno + x - 1] + dn[dno + x + 1]
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

    private val KERNEL0 = intArrayOf(0, -2, -1, -1, 0, -1, 1, -1, -2, 0, -1, 0, 1, 0, 2, 0, -1, 1, 0, 1, 1, 1, 0, 2)
    private val DIST0 = intArrayOf(0, -1, 1, 0, 0, 0, -1, 0, 0, 1)

    /** EPF の段階 (0〜2) が上下に読む行数 */
    fun epfRadius(step: Int): Int = when (step) { 0 -> 3; 1 -> 2; else -> 1 }

    /** EPF の反復回数から、行う段階の並び */
    fun epfSteps(params: EpfParams?): IntArray = when (params?.iters) {
        null -> IntArray(0)
        1 -> intArrayOf(1)
        2 -> intArrayOf(1, 2)
        else -> intArrayOf(0, 1, 2)
    }

    /**
     * エッジ保存フィルタ (フレーム全体の平面をその場で書き換える)。[sigma] は VarDCT のブロックごとの値
     * (行の長さ [sigmaStride])、Modular では null (sigma_for_modular を使う)。
     * スレッドごとに連続した行 (縞) を受け持ち、[EpfRows] で上から 1 行ずつ計算して書き戻す。
     */
    fun epf(colors: Array<Plane>, params: EpfParams, sigma: FloatArray?, sigmaStride: Int, executor: ExecutorService?): Array<Plane> {
        val w = colors[0].width
        val h = colors[0].height
        val steps = epfSteps(params)
        var halo = 0
        for (st in steps) halo += epfRadius(st)
        val planes = Array(3) { colors[it].data }
        val threads = if (executor == null) 1 else Runtime.getRuntime().availableProcessors()
        val stripes = maxOf(1, minOf(threads, (h + 31) / 32))
        val per = (h + stripes - 1) / stripes
        // 縞の外側の行は、他の縞に書き換えられる前に退避しておく (上の halo 行と下の halo 行、画像の外は鏡映)
        val saved = Array(stripes) { s ->
            val ys = minOf(h, s * per)
            val ye = minOf(h, (s + 1) * per)
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
            val ys = minOf(h, s * per)
            val ye = minOf(h, (s + 1) * per)
            if (ys >= ye) return@parallelFor
            val sv = saved[s]
            // 縞の中の行は、出力の行 y を書き戻す時点で y + halo 行目まで読み終えているので、その場で書き戻してよい
            val src = RowSource { c, y, dst, off ->
                if (y < ys) System.arraycopy(sv[c], (y - (ys - halo)) * w, dst, off, w)
                else if (y >= ye) System.arraycopy(sv[c], (halo + y - ye) * w, dst, off, w)
                else System.arraycopy(planes[c], y * w, dst, off, w)
            }
            val rows = EpfRows.obtain()
            try {
                rows.begin(src, w, h, steps, params, sigma, sigmaStride, ys)
                for (y in ys until ye) {
                    val o = rows.row(y) + 3
                    for (c in 0 until 3) System.arraycopy(rows.out(c), o, planes[c], y * w, w)
                }
            } finally {
                rows.recycle()
            }
        }
        return colors
    }

    /** チャンネル c の行 y (画像の外でもよい) の w 画素を dst[off] から書く */
    fun interface RowSource {
        fun row(c: Int, y: Int, dst: FloatArray, off: Int)
    }

    /**
     * EPF を行ごとに上から順に計算する。入力と各段階の結果は数行の環状バッファにだけ持つので、
     * 帯に分けて重なり (halo) を計算し直すことがなく、作業領域も幅 × 数十行で済む。
     * 1 つのスレッドが連続した行を受け持つ ([begin] の行から順に [row] を呼ぶ)。
     *
     * 各段階の行 y は配列の (y and mask) * (w + 6) から始まり、画素 x は +3+x (左右 3 画素は鏡映)。
     * 画像の外の行は、入力では [RowSource] が、EPF の段階では計算した行を鏡映して返す。
     */
    class EpfRows {
        private class Stage {
            /** -1 は入力、0〜2 は EPF の段階 */
            var step = -1
            var mask = 0
            var rows: Array<FloatArray> = Array(3) { FloatArray(0) }
            /** 次に計算する行 */
            var next = 0

            // 段階 1・2: 右 (hd)・下 (vd) の画素との重み付き差分 (4 行の環状)
            var hd = FloatArray(0)
            var vd = FloatArray(0)
            var rawNext = 0

            // 段階 1: 差分の十字 (5 点) の和。hc は今の行だけ、vc は 4 行の環状
            var hc = FloatArray(0)
            var vc = FloatArray(0)
            var vcNext = 0

            fun ensure(pw: Int) {
                val n = (mask + 1) * pw
                if (rows[0].size < n) rows = Array(3) { FloatArray(n) }
                if (step == 1 || step == 2) {
                    if (hd.size < 4 * pw) {
                        hd = FloatArray(4 * pw)
                        vd = FloatArray(4 * pw)
                    }
                }
                if (step == 1) {
                    if (hc.size < pw) hc = FloatArray(pw)
                    if (vc.size < 4 * pw) vc = FloatArray(4 * pw)
                }
            }

            fun bytes(): Long = 4L * (3L * rows[0].size + hd.size + vd.size + hc.size + vc.size)
        }

        private val stages = Array(4) { Stage() }
        private var n = 1
        private var w = 0
        private var h = 0
        private var pw = 0
        private var src: RowSource? = null
        private var params: EpfParams? = null
        private var sigma: FloatArray? = null
        private var sigmaStride = 0
        private var nis = FloatArray(0)
        private val tapK = IntArray(12)
        private val tapA = IntArray(60)
        private val tapB = IntArray(60)

        /**
         * 行 [y0] から計算を始める。[steps] は EPF の段階の並び (空なら入力をそのまま返す)
         */
        fun begin(
            src: RowSource, w: Int, h: Int, steps: IntArray, params: EpfParams?, sigma: FloatArray?, sigmaStride: Int, y0: Int,
        ) {
            this.src = src
            this.w = w
            this.h = h
            this.params = params
            this.sigma = sigma
            this.sigmaStride = sigmaStride
            pw = w + 6
            if (nis.size < w) nis = FloatArray(w)
            n = steps.size + 1
            var start = y0
            for (k in n - 1 downTo 0) {
                val st = stages[k]
                st.step = if (k == 0) -1 else steps[k - 1]
                // 次の段階が上下に読む行数ぶんを残す (最後の段階は計算した行をすぐ使うので 1 行)
                st.mask = if (k == n - 1) 0 else if (epfRadius(steps[k]) == 1) 3 else 7
                st.ensure(pw)
                st.next = start
                st.rawNext = if (st.step == 1) start - 2 else start - 1
                st.vcNext = start - 1
                if (k > 0) {
                    start -= epfRadius(st.step)
                    // EPF の段階は画像の中の行だけ計算する (外の行は鏡映)
                    if (k > 1) start = maxOf(start, 0)
                }
            }
        }

        /** 最後の段階の配列 (チャンネル c) */
        fun out(c: Int): FloatArray = stages[n - 1].rows[c]

        /** 最後の段階の行 y (0 <= y < h) を計算し、[out] の中の行の先頭の位置を返す */
        fun row(y: Int): Int = ensure(n - 1, y)

        private fun ensure(k: Int, y: Int): Int {
            val st = stages[k]
            val yy = if (k == 0) y else mirror(y, h)
            while (st.next <= yy) {
                produce(k, st.next)
                st.next++
            }
            return (yy and st.mask) * pw
        }

        private fun produce(k: Int, y: Int) {
            val st = stages[k]
            val o = (y and st.mask) * pw
            if (k == 0) {
                val s = src!!
                for (c in 0 until 3) {
                    s.row(c, y, st.rows[c], o + 3)
                    padRow(st.rows[c], o, w)
                }
                return
            }
            if (st.step == 0) step0(k, y, o) else step12(k, y, o)
            for (c in 0 until 3) padRow(st.rows[c], o, w)
        }

        /** 段階 k の差分の行 r (前の段階の行 r と r+1 から) */
        private fun rawDiff(k: Int, r: Int) {
            val st = stages[k]
            val pv = stages[k - 1]
            val a = ensure(k - 1, r)
            val b = ensure(k - 1, r + 1)
            val p = params!!
            val s0 = p.channelScale[0]
            val s1 = p.channelScale[1]
            val s2 = p.channelScale[2]
            val p0 = pv.rows[0]
            val p1 = pv.rows[1]
            val p2 = pv.rows[2]
            val hd = st.hd
            val vd = st.vd
            val lb = (r and 3) * pw
            // 使う列だけ: 段階 1 は十字の和のために左右 1 列多く要る
            val plus = st.step == 1
            val h0 = if (plus) 1 else 2
            val h1 = if (plus) w + 4 else w + 3
            // この行の差分を使うのは、出力の行 r-1 .. r+2 の平滑化するブロック (の左 2・右 1 画素まで) だけ
            forActiveRuns(r - 1, r + 2) { xs, xe ->
                val x0 = maxOf(h0, xs + 1)
                val x1 = minOf(h1, xe + 4)
                for (x in x0 until x1) {
                    hd[lb + x] = s0 * kotlin.math.abs(p0[a + x + 1] - p0[a + x]) + s1 * kotlin.math.abs(p1[a + x + 1] - p1[a + x]) +
                        s2 * kotlin.math.abs(p2[a + x + 1] - p2[a + x])
                }
                for (x in maxOf(h0 + 1, x0) until x1) {
                    vd[lb + x] = s0 * kotlin.math.abs(p0[b + x] - p0[a + x]) + s1 * kotlin.math.abs(p1[b + x] - p1[a + x]) +
                        s2 * kotlin.math.abs(p2[b + x] - p2[a + x])
                }
            }
        }

        /** ブロック行 br・ブロック列 bx を平滑化するか (sigma が小さいブロックはそのまま写す) */
        private fun blockActive(br: Int, bx: Int): Boolean {
            val sig = sigma
            val s = if (sig != null) sig[br * sigmaStride + bx] else params!!.sigmaForModular
            return !(s < 0.3f)
        }

        /**
         * 出力の行 ya..yb のどれかで平滑化するブロックが続く列の範囲 [xs, xe) (画素) ごとに f を呼ぶ。
         * 画像の外の行は含めない
         */
        private inline fun forActiveRuns(ya: Int, yb: Int, f: (Int, Int) -> Unit) {
            val br0 = maxOf(0, ya) shr 3
            val br1 = (minOf(h - 1, yb)) shr 3
            if (br1 < br0) return
            val nb = (w + 7) shr 3
            var bx = 0
            while (bx < nb) {
                if (!blockActive(br0, bx) && (br1 == br0 || !blockActive(br1, bx))) {
                    bx++
                    continue
                }
                val bs = bx
                bx++
                while (bx < nb && (blockActive(br0, bx) || (br1 != br0 && blockActive(br1, bx)))) bx++
                f(bs * 8, minOf(w, bx * 8))
            }
        }

        /**
         * 十字 (4 近傍) のカーネルの段階 1・2。隣の画素との重み付き差分 H (右)・V (下) を先に求め、
         * 段階 1 の距離はその十字の和、段階 2 は差分そのもの
         */
        private fun step12(k: Int, y: Int, o: Int) {
            val st = stages[k]
            val pv = stages[k - 1]
            val p = params!!
            val plus = st.step == 1
            val rawTo = if (plus) y + 1 else y
            while (st.rawNext <= rawTo) {
                rawDiff(k, st.rawNext)
                st.rawNext++
            }
            val sh: FloatArray
            val sho: Int
            val sv: FloatArray
            if (plus) {
                val d = st.hd
                val c = (y and 3) * pw
                val u = ((y - 1) and 3) * pw
                val dn = ((y + 1) and 3) * pw
                val hc = st.hc
                // 十字の和も、平滑化する画素が使う列だけ (hc は左 1 画素まで、vc はその行と次の行で使う)
                forActiveRuns(y, y) { xs, xe ->
                    for (x in maxOf(2, xs + 2) until minOf(w + 3, xe + 3)) {
                        hc[x] = d[c + x] + d[u + x] + d[dn + x] + d[c + x - 1] + d[c + x + 1]
                    }
                }
                while (st.vcNext <= y) {
                    val r = st.vcNext
                    val v = st.vd
                    val vc = st.vc
                    val vr = (r and 3) * pw
                    val vu = ((r - 1) and 3) * pw
                    val vn = ((r + 1) and 3) * pw
                    forActiveRuns(r, r + 1) { xs, xe ->
                        for (x in maxOf(3, xs + 3) until minOf(w + 3, xe + 3)) {
                            vc[vr + x] = v[vr + x] + v[vu + x] + v[vn + x] + v[vr + x - 1] + v[vr + x + 1]
                        }
                    }
                    st.vcNext++
                }
                sh = hc
                sho = 0
                sv = st.vc
            } else {
                sh = st.hd
                sho = (y and 3) * pw
                sv = st.vd
            }
            val svo = (y and 3) * pw
            val svu = ((y - 1) and 3) * pw
            val pm = ensure(k - 1, y)
            val pu = ensure(k - 1, y - 1)
            val pd = ensure(k - 1, y + 1)
            val stepMul = if (st.step == 2) p.pass2SigmaScale else 1f
            val border = p.borderSadMul
            val i0 = pv.rows[0]
            val i1 = pv.rows[1]
            val i2 = pv.rows[2]
            val o0 = st.rows[0]
            val o1 = st.rows[1]
            val o2 = st.rows[2]
            val sig = sigma
            val sRow = (y shr 3) * sigmaStride
            val yBorder = ((y + 1) and 6) == 0
            val nisRow = nis
            var bx = 0
            while (bx < w) {
                // 平滑化するブロックが続く範囲 [bx, runEnd) の距離の倍率を画素ごとに求めてから、まとめて計算する
                var runEnd = bx
                while (runEnd < w) {
                    val s = if (sig != null) sig[sRow + (runEnd shr 3)] else p.sigmaForModular
                    if (s < 0.3f) break
                    val xEnd = minOf(w, runEnd + 8)
                    val nisBorder = EPF_C0 / s * stepMul * border
                    val nisInner = if (yBorder) nisBorder else EPF_C0 / s * stepMul
                    nisRow[runEnd] = nisBorder
                    for (x in runEnd + 1 until xEnd) nisRow[x] = nisInner
                    if (xEnd - runEnd == 8) nisRow[runEnd + 7] = nisBorder
                    runEnd = xEnd
                }
                for (x in bx until runEnd) {
                    val ni = nisRow[x]
                    val li = x + 3
                    val wr = maxOf(0f, 1f + sh[sho + li] * ni)
                    val wl = maxOf(0f, 1f + sh[sho + li - 1] * ni)
                    val wd = maxOf(0f, 1f + sv[svo + li] * ni)
                    val wu = maxOf(0f, 1f + sv[svu + li] * ni)
                    val inv = 1f / (1f + wr + wl + wd + wu)
                    val oi = o + li
                    o0[oi] = (i0[pm + li] + wr * i0[pm + li + 1] + wl * i0[pm + li - 1] + wd * i0[pd + li] + wu * i0[pu + li]) * inv
                    o1[oi] = (i1[pm + li] + wr * i1[pm + li + 1] + wl * i1[pm + li - 1] + wd * i1[pd + li] + wu * i1[pu + li]) * inv
                    o2[oi] = (i2[pm + li] + wr * i2[pm + li + 1] + wl * i2[pm + li - 1] + wd * i2[pd + li] + wu * i2[pu + li]) * inv
                }
                // 平滑化しないブロック (sigma が小さい) はそのまま写す
                bx = runEnd
                while (bx < w) {
                    val s = if (sig != null) sig[sRow + (bx shr 3)] else p.sigmaForModular
                    if (s >= 0.3f) break
                    val xEnd = minOf(w, bx + 8)
                    System.arraycopy(i0, pm + bx + 3, o0, o + bx + 3, xEnd - bx)
                    System.arraycopy(i1, pm + bx + 3, o1, o + bx + 3, xEnd - bx)
                    System.arraycopy(i2, pm + bx + 3, o2, o + bx + 3, xEnd - bx)
                    bx = xEnd
                }
            }
        }

        /** 12 近傍のカーネルの段階 0 (iters = 3 のときだけ使う) */
        private fun step0(k: Int, y: Int, o: Int) {
            val st = stages[k]
            val pv = stages[k - 1]
            val p = params!!
            // 前の段階の行 y-3 .. y+3 の位置から、カーネルの各点と距離の各点の位置を求める
            val rowOff = IntArray(7) { ensure(k - 1, y - 3 + it) }
            val nk = KERNEL0.size / 2
            val nd = DIST0.size / 2
            for (t in 0 until nk) {
                val kx = KERNEL0[2 * t]
                val ky = KERNEL0[2 * t + 1]
                tapK[t] = rowOff[ky + 3] + kx
                for (d in 0 until nd) {
                    val dx = DIST0[2 * d]
                    val dy = DIST0[2 * d + 1]
                    tapA[t * nd + d] = rowOff[ky + dy + 3] + kx + dx
                    tapB[t * nd + d] = rowOff[dy + 3] + dx
                }
            }
            val pm = rowOff[3]
            val stepMul = p.pass0SigmaScale
            val border = p.borderSadMul
            val scale0 = p.channelScale[0]
            val scale1 = p.channelScale[1]
            val scale2 = p.channelScale[2]
            val i0 = pv.rows[0]
            val i1 = pv.rows[1]
            val i2 = pv.rows[2]
            val o0 = st.rows[0]
            val o1 = st.rows[1]
            val o2 = st.rows[2]
            val sig = sigma
            val yBorder = ((y + 1) and 6) == 0
            for (x in 0 until w) {
                val s = if (sig != null) sig[(y shr 3) * sigmaStride + (x shr 3)] else p.sigmaForModular
                val li = x + 3
                val pi = pm + li
                val oi = o + li
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
                for (t in 0 until nk) {
                    val ki = tapK[t] + li
                    var a0 = 0f
                    var a1 = 0f
                    var a2 = 0f
                    for (d in 0 until nd) {
                        val ia = tapA[t * nd + d] + li
                        val ib = tapB[t * nd + d] + li
                        a0 += kotlin.math.abs(i0[ia] - i0[ib])
                        a1 += kotlin.math.abs(i1[ia] - i1[ib])
                        a2 += kotlin.math.abs(i2[ia] - i2[ib])
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
        }

        /** 使い終わったら [ScratchPool] に返す */
        fun recycle() {
            src = null
            sigma = null
            var bytes = 1024L
            for (st in stages) bytes += st.bytes()
            ScratchPool.give(this, bytes)
        }

        companion object {
            /** [ScratchPool] から取り出す。無ければ新しく作る */
            fun obtain(): EpfRows = ScratchPool.take(EpfRows::class.java) { true } ?: EpfRows()
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
