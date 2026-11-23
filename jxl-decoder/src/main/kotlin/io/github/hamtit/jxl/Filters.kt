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
     *
     * 画像の一部の列 (窓) だけを計算することもできる ([begin] の xBlock0)。窓を幅 w の画像とみなして計算するので、
     * 窓の端が画像の端でない側では、端から (入力の Gaborish を含めた) halo 画素以内の結果は正しくない。
     * 使う列の外側に halo 画素以上の余白を付けた窓を渡すこと。
     */
    class EpfRows {
        private class Stage {
            /** -1 は入力、0〜2 は EPF の段階 */
            var step = -1
            var mask = 0
            var rows: Array<FloatArray> = Array(3) { FloatArray(0) }
            /** 次に計算する行 */
            var next = 0

            // 段階 0・1・2: 右 (hd)・下 (vd) の画素との重み付き差分 (4 行の環状)
            var hd = FloatArray(0)
            var vd = FloatArray(0)
            var rawNext = 0

            // 段階 0・1: 差分の十字 (5 点) の和。hc は今の行だけ、vc は 4 行の環状
            var hc = FloatArray(0)
            var vc = FloatArray(0)
            var vcNext = 0

            // 段階 0 (12 近傍、対称な 6 方向の残り 4 方向): 右2 (h2)・下2 (v2)・右下 (dr)・左下 (dl)
            var h2d = FloatArray(0)
            var v2d = FloatArray(0)
            var drd = FloatArray(0)
            var dld = FloatArray(0)
            var h2c = FloatArray(0)
            var v2c = FloatArray(0)
            var drc = FloatArray(0)
            var dlc = FloatArray(0)

            fun ensure(pw: Int) {
                val n = (mask + 1) * pw
                if (rows[0].size < n) rows = Array(3) { FloatArray(n) }
                if (step >= 0) {
                    if (hd.size < 4 * pw) {
                        hd = FloatArray(4 * pw)
                        vd = FloatArray(4 * pw)
                    }
                }
                if (step == 0 || step == 1) {
                    if (hc.size < pw) hc = FloatArray(pw)
                    if (vc.size < 4 * pw) vc = FloatArray(4 * pw)
                }
                if (step == 0) {
                    if (h2d.size < 4 * pw) {
                        h2d = FloatArray(4 * pw)
                        v2d = FloatArray(4 * pw)
                        drd = FloatArray(4 * pw)
                        dld = FloatArray(4 * pw)
                    }
                    if (h2c.size < pw) h2c = FloatArray(pw)
                    if (v2c.size < 4 * pw) {
                        v2c = FloatArray(4 * pw)
                        drc = FloatArray(4 * pw)
                        dlc = FloatArray(4 * pw)
                    }
                }
            }

            fun bytes(): Long = 4L * (
                3L * rows[0].size + hd.size + vd.size + hc.size + vc.size +
                    h2d.size + v2d.size + drd.size + dld.size + h2c.size + v2c.size + drc.size + dlc.size
                )
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
        /** 窓の左端のブロック列 (sigma の列の位置に足す) */
        private var xb0 = 0
        private var nis = FloatArray(0)

        /**
         * 行 [y0] から計算を始める。[steps] は EPF の段階の並び (空なら入力をそのまま返す)。
         * 窓を計算するときは、[w] が窓の幅、[xBlock0] が窓の左端の画像のブロック列 (窓の左端は 8 の倍数の列)
         */
        fun begin(
            src: RowSource, w: Int, h: Int, steps: IntArray, params: EpfParams?, sigma: FloatArray?, sigmaStride: Int, y0: Int,
            xBlock0: Int = 0,
        ) {
            this.src = src
            this.xb0 = xBlock0
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
                st.rawNext = when (st.step) {
                    0 -> start - 3
                    1 -> start - 2
                    else -> start - 1
                }
                st.vcNext = if (st.step == 0) start - 2 else start - 1
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
                    val ax = a + x
                    hd[lb + x] = s0 * kotlin.math.abs(p0[ax + 1] - p0[ax]) +
                        s1 * kotlin.math.abs(p1[ax + 1] - p1[ax]) +
                        s2 * kotlin.math.abs(p2[ax + 1] - p2[ax])
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
            val s = if (sig != null) sig[br * sigmaStride + xb0 + bx] else params!!.sigmaForModular
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
            val sRow = (y shr 3) * sigmaStride + xb0
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
                    // max(0, v)。NaN はそのまま (maxOf と同じ)、-0 は 0 と違うが足し算の結果は変わらない
                    val vr = 1f + sh[sho + li] * ni
                    val vl = 1f + sh[sho + li - 1] * ni
                    val vd = 1f + sv[svo + li] * ni
                    val vu = 1f + sv[svu + li] * ni
                    val wr = if (vr < 0f) 0f else vr
                    val wl = if (vl < 0f) 0f else vl
                    val wd = if (vd < 0f) 0f else vd
                    val wu = if (vu < 0f) 0f else vu
                    val inv = 1f / (1f + wr + wl + wd + wu)
                    val pml = pm + li
                    val oi = o + li
                    o0[oi] = (i0[pml] + wr * i0[pml + 1] + wl * i0[pml - 1] + wd * i0[pd + li] + wu * i0[pu + li]) * inv
                    o1[oi] = (i1[pml] + wr * i1[pml + 1] + wl * i1[pml - 1] + wd * i1[pd + li] + wu * i1[pu + li]) * inv
                    o2[oi] = (i2[pml] + wr * i2[pml + 1] + wl * i2[pml - 1] + wd * i2[pd + li] + wu * i2[pu + li]) * inv
                }
                // 平滑化しないブロック (sigma が小さい) はそのまま写す
                bx = runEnd
                while (bx < w) {
                    val s = if (sig != null) sig[sRow + (bx shr 3)] else p.sigmaForModular
                    if (!(s < 0.3f)) break
                    val xEnd = minOf(w, bx + 8)
                    System.arraycopy(i0, pm + bx + 3, o0, o + bx + 3, xEnd - bx)
                    System.arraycopy(i1, pm + bx + 3, o1, o + bx + 3, xEnd - bx)
                    System.arraycopy(i2, pm + bx + 3, o2, o + bx + 3, xEnd - bx)
                    bx = xEnd
                }
            }
        }

        /** 段階 0 の 6 方向の重み付き差分 (行 r から r+1・r+2 へ) */
        private fun rawDiff0(k: Int, r: Int) {
            val st = stages[k]
            val pv = stages[k - 1]
            val a = ensure(k - 1, r)
            val b = ensure(k - 1, r + 1)
            val c = ensure(k - 1, r + 2)
            val p = params!!
            val s0 = p.channelScale[0]
            val s1 = p.channelScale[1]
            val s2 = p.channelScale[2]
            val p0 = pv.rows[0]
            val p1 = pv.rows[1]
            val p2 = pv.rows[2]
            val hd = st.hd
            val vd = st.vd
            val h2d = st.h2d
            val v2d = st.v2d
            val drd = st.drd
            val dld = st.dld
            val lb = (r and 3) * pw
            // 行 r の差分は出力行 r-1 .. r+3 で使われる
            forActiveRuns(r - 1, r + 3) { xs, xe ->
                val x0 = maxOf(0, xs)
                val x1 = minOf(pw, xe + 5)
                val x1m1 = minOf(pw - 1, x1)
                val x1m2 = minOf(pw - 2, x1)
                for (x in x0 until x1m1) {
                    val ax = a + x
                    hd[lb + x] = s0 * kotlin.math.abs(p0[ax + 1] - p0[ax]) +
                        s1 * kotlin.math.abs(p1[ax + 1] - p1[ax]) +
                        s2 * kotlin.math.abs(p2[ax + 1] - p2[ax])
                }
                for (x in x0 until x1m2) {
                    val ax = a + x
                    h2d[lb + x] = s0 * kotlin.math.abs(p0[ax + 2] - p0[ax]) +
                        s1 * kotlin.math.abs(p1[ax + 2] - p1[ax]) +
                        s2 * kotlin.math.abs(p2[ax + 2] - p2[ax])
                }
                for (x in x0 until x1) {
                    val ax = a + x
                    val bx = b + x
                    val cx = c + x
                    val a0 = p0[ax]; val a1 = p1[ax]; val a2 = p2[ax]
                    vd[lb + x] = s0 * kotlin.math.abs(p0[bx] - a0) +
                        s1 * kotlin.math.abs(p1[bx] - a1) +
                        s2 * kotlin.math.abs(p2[bx] - a2)
                    v2d[lb + x] = s0 * kotlin.math.abs(p0[cx] - a0) +
                        s1 * kotlin.math.abs(p1[cx] - a1) +
                        s2 * kotlin.math.abs(p2[cx] - a2)
                }
                for (x in x0 until x1m1) {
                    val ax = a + x
                    val bx1 = b + x + 1
                    drd[lb + x] = s0 * kotlin.math.abs(p0[bx1] - p0[ax]) +
                        s1 * kotlin.math.abs(p1[bx1] - p1[ax]) +
                        s2 * kotlin.math.abs(p2[bx1] - p2[ax])
                }
                for (x in maxOf(1, x0) until x1) {
                    val ax = a + x
                    val bxm1 = b + x - 1
                    dld[lb + x] = s0 * kotlin.math.abs(p0[bxm1] - p0[ax]) +
                        s1 * kotlin.math.abs(p1[bxm1] - p1[ax]) +
                        s2 * kotlin.math.abs(p2[bxm1] - p2[ax])
                }
            }
        }

        /** 12 近傍のカーネルの段階 0 (iters = 3 のときだけ使う)。対称な 6 方向の差分と 5 点十字和を共有して計算する */
        private fun step0(k: Int, y: Int, o: Int) {
            val st = stages[k]
            val pv = stages[k - 1]
            val p = params!!
            while (st.rawNext <= y + 1) {
                rawDiff0(k, st.rawNext)
                st.rawNext++
            }
            // 縦・斜め 4 方向の 5 点十字和 (4 行の環状バッファ)
            while (st.vcNext <= y) {
                val r = st.vcNext
                val vr = (r and 3) * pw
                val vu = ((r - 1) and 3) * pw
                val vn = ((r + 1) and 3) * pw
                val vd = st.vd; val vc = st.vc
                val v2d = st.v2d; val v2c = st.v2c
                val drd = st.drd; val drc = st.drc
                val dld = st.dld; val dlc = st.dlc
                forActiveRuns(r, r + 2) { xs, xe ->
                    val x0 = maxOf(1, xs + 1)
                    val x1 = minOf(pw - 1, xe + 4)
                    for (x in x0 until x1) {
                        vc[vr + x] = vd[vr + x] + vd[vu + x] + vd[vn + x] + vd[vr + x - 1] + vd[vr + x + 1]
                        v2c[vr + x] = v2d[vr + x] + v2d[vu + x] + v2d[vn + x] + v2d[vr + x - 1] + v2d[vr + x + 1]
                        drc[vr + x] = drd[vr + x] + drd[vu + x] + drd[vn + x] + drd[vr + x - 1] + drd[vr + x + 1]
                        dlc[vr + x] = dld[vr + x] + dld[vu + x] + dld[vn + x] + dld[vr + x - 1] + dld[vr + x + 1]
                    }
                }
                st.vcNext++
            }
            // 水平 2 方向の 5 点十字和 (今の行 y だけ)
            val r0 = (y and 3) * pw
            val r1u = ((y - 1) and 3) * pw
            val r2u = ((y - 2) and 3) * pw
            val r1d = ((y + 1) and 3) * pw
            val hd = st.hd; val hc = st.hc
            val h2d = st.h2d; val h2c = st.h2c
            forActiveRuns(y, y) { xs, xe ->
                val x0 = maxOf(1, xs + 1)
                val x1 = minOf(pw - 1, xe + 3)
                for (x in x0 until x1) {
                    hc[x] = hd[r0 + x] + hd[r1u + x] + hd[r1d + x] + hd[r0 + x - 1] + hd[r0 + x + 1]
                    h2c[x] = h2d[r0 + x] + h2d[r1u + x] + h2d[r1d + x] + h2d[r0 + x - 1] + h2d[r0 + x + 1]
                }
            }
            val pm = ensure(k - 1, y)
            val pu1 = ensure(k - 1, y - 1)
            val pu2 = ensure(k - 1, y - 2)
            val pd1 = ensure(k - 1, y + 1)
            val pd2 = ensure(k - 1, y + 2)
            val vc = st.vc; val v2c = st.v2c; val drc = st.drc; val dlc = st.dlc
            val stepMul = p.pass0SigmaScale
            val border = p.borderSadMul
            val i0 = pv.rows[0]
            val i1 = pv.rows[1]
            val i2 = pv.rows[2]
            val o0 = st.rows[0]
            val o1 = st.rows[1]
            val o2 = st.rows[2]
            val sig = sigma
            val sRow = (y shr 3) * sigmaStride + xb0
            val yBorder = ((y + 1) and 6) == 0
            val nisRow = nis
            var bx = 0
            while (bx < w) {
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
                    // 12 近傍の重み (max(0, 1 + d * ni))
                    val t0 = 1f + v2c[r2u + li] * ni; val w0 = if (t0 < 0f) 0f else t0
                    val t1 = 1f + drc[r1u + li - 1] * ni; val w1 = if (t1 < 0f) 0f else t1
                    val t2 = 1f + vc[r1u + li] * ni; val w2 = if (t2 < 0f) 0f else t2
                    val t3 = 1f + dlc[r1u + li + 1] * ni; val w3 = if (t3 < 0f) 0f else t3
                    val t4 = 1f + h2c[li - 2] * ni; val w4 = if (t4 < 0f) 0f else t4
                    val t5 = 1f + hc[li - 1] * ni; val w5 = if (t5 < 0f) 0f else t5
                    val t6 = 1f + hc[li] * ni; val w6 = if (t6 < 0f) 0f else t6
                    val t7 = 1f + h2c[li] * ni; val w7 = if (t7 < 0f) 0f else t7
                    val t8 = 1f + dlc[r0 + li] * ni; val w8 = if (t8 < 0f) 0f else t8
                    val t9 = 1f + vc[r0 + li] * ni; val w9 = if (t9 < 0f) 0f else t9
                    val t10 = 1f + drc[r0 + li] * ni; val w10 = if (t10 < 0f) 0f else t10
                    val t11 = 1f + v2c[r0 + li] * ni; val w11 = if (t11 < 0f) 0f else t11
                    val inv = 1f / (1f + w0 + w1 + w2 + w3 + w4 + w5 + w6 + w7 + w8 + w9 + w10 + w11)
                    val oi = o + li
                    val pml = pm + li; val pu1l = pu1 + li; val pu2l = pu2 + li; val pd1l = pd1 + li; val pd2l = pd2 + li
                    o0[oi] = (
                        i0[pml] + w0 * i0[pu2l] + w1 * i0[pu1l - 1] + w2 * i0[pu1l] + w3 * i0[pu1l + 1] +
                            w4 * i0[pml - 2] + w5 * i0[pml - 1] + w6 * i0[pml + 1] + w7 * i0[pml + 2] +
                            w8 * i0[pd1l - 1] + w9 * i0[pd1l] + w10 * i0[pd1l + 1] + w11 * i0[pd2l]
                        ) * inv
                    o1[oi] = (
                        i1[pml] + w0 * i1[pu2l] + w1 * i1[pu1l - 1] + w2 * i1[pu1l] + w3 * i1[pu1l + 1] +
                            w4 * i1[pml - 2] + w5 * i1[pml - 1] + w6 * i1[pml + 1] + w7 * i1[pml + 2] +
                            w8 * i1[pd1l - 1] + w9 * i1[pd1l] + w10 * i1[pd1l + 1] + w11 * i1[pd2l]
                        ) * inv
                    o2[oi] = (
                        i2[pml] + w0 * i2[pu2l] + w1 * i2[pu1l - 1] + w2 * i2[pu1l] + w3 * i2[pu1l + 1] +
                            w4 * i2[pml - 2] + w5 * i2[pml - 1] + w6 * i2[pml + 1] + w7 * i2[pml + 2] +
                            w8 * i2[pd1l - 1] + w9 * i2[pd1l] + w10 * i2[pd1l + 1] + w11 * i2[pd2l]
                        ) * inv
                }
                bx = runEnd
                while (bx < w) {
                    val s = if (sig != null) sig[sRow + (bx shr 3)] else p.sigmaForModular
                    if (!(s < 0.3f)) break
                    val xEnd = minOf(w, bx + 8)
                    System.arraycopy(i0, pm + bx + 3, o0, o + bx + 3, xEnd - bx)
                    System.arraycopy(i1, pm + bx + 3, o1, o + bx + 3, xEnd - bx)
                    System.arraycopy(i2, pm + bx + 3, o2, o + bx + 3, xEnd - bx)
                    bx = xEnd
                }
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
