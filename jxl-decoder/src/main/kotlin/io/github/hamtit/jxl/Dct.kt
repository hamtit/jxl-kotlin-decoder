package io.github.hamtit.jxl

/** VarDCT の DCT と小さい変換 (AFV など) */
internal object Dct {
    private const val SQRT2 = 1.4142135f

    private val SEC_HALF = HashMap<Int, FloatArray>()

    private fun secHalf(n: Int): FloatArray {
        synchronized(SEC_HALF) {
            return SEC_HALF.getOrPut(n) {
                FloatArray(n / 2) { k ->
                    val theta = (2 * k + 1).toDouble() / (2 * n).toDouble() * Math.PI
                    (1.0 / Math.cos(theta) / 2.0).toFloat()
                }
            }
        }
    }

    private val SEC8 = floatArrayOf(0.5097955791041592f, 0.6013448869350453f, 0.8999762231364156f, 2.5629154477415055f)
    private val SEC_TABLES = arrayOfNulls<FloatArray>(10).also {
        for (i in 2..9) it[i] = secHalf(1 shl i)
    }

    private val SCALE_F = floatArrayOf(
        1.0000000000000000f, 0.9996047255830407f, 0.9984194528776054f, 0.9964458326264695f,
        0.9936866130906366f, 0.9901456355893141f, 0.9858278282666936f, 0.9807391980963174f,
        0.9748868211368796f, 0.9682788310563117f, 0.9609244059440204f, 0.9528337534340876f,
        0.9440180941651672f, 0.9344896436056892f, 0.9242615922757944f, 0.9133480844001980f,
        0.9017641950288744f, 0.8895259056651056f, 0.8766500784429904f, 0.8631544288990163f,
        0.8490574973847023f, 0.8343786191696513f, 0.8191378932865928f, 0.8033561501721485f,
        0.7870549181591013f, 0.7702563888779096f, 0.7529833816270532f, 0.7352593067735488f,
        0.7171081282466044f, 0.6985543251889097f, 0.6796228528314652f, 0.6603391026591464f,
    )

    fun scaleF(c: Int, logb: Int): Float = SCALE_F[c shl logb]

    /** 1 次元 DCT (順方向は 1/N 正規化、逆方向は係数 0 以外に √2) */
    fun dct1d(io: FloatArray, off: Int, stride: Int, n: Int, scratch: FloatArray, forward: Boolean) {
        when (n) {
            1 -> return
            2 -> {
                val a = io[off]
                val b = io[off + stride]
                if (forward) {
                    io[off] = (a + b) / 2f
                    io[off + stride] = (a - b) / 2f
                } else {
                    io[off] = a + b
                    io[off + stride] = a - b
                }
                return
            }
        }
        for (i in 0 until n) scratch[i] = io[off + i * stride]
        dctInPlace(scratch, 0, n, scratch, n, forward)
        for (i in 0 until n) io[off + i * stride] = scratch[i]
    }

    /** buf[off until off+n] を変換する。tmp[tmpOff until tmpOff+n] を作業領域に使う */
    private fun dctInPlace(buf: FloatArray, off: Int, n: Int, tmp: FloatArray, tmpOff: Int, forward: Boolean) {
        if (n == 1) return
        if (n == 2) {
            val a = buf[off]
            val b = buf[off + 1]
            if (forward) {
                buf[off] = (a + b) / 2f
                buf[off + 1] = (a - b) / 2f
            } else {
                buf[off] = a + b
                buf[off + 1] = a - b
            }
            return
        }
        if (n == 4) {
            dct4(buf, off, forward)
            return
        }
        if (!forward) {
            when (n) {
                8 -> { idct8(buf, off); return }
                16 -> { idct16(buf, off, tmp, tmpOff); return }
                32 -> { idct32(buf, off, tmp, tmpOff); return }
            }
        }
        val half = n / 2
        val sec = SEC_TABLES[Integer.numberOfTrailingZeros(n)] ?: secHalf(n)
        if (forward) {
            for (i in 0 until half) {
                val a = buf[off + i]
                val b = buf[off + n - i - 1]
                tmp[tmpOff + i] = (a + b) / 2f
                tmp[tmpOff + half + i] = (a - b) / 2f * sec[i]
            }
            dctInPlace(tmp, tmpOff, half, buf, off, true)
            dctInPlace(tmp, tmpOff + half, half, buf, off + half, true)
            tmp[tmpOff + half] *= SQRT2
            for (i in 0 until half - 1) tmp[tmpOff + half + i] += tmp[tmpOff + half + i + 1]
            for (i in 0 until half) {
                buf[off + 2 * i] = tmp[tmpOff + i]
                buf[off + 2 * i + 1] = tmp[tmpOff + half + i]
            }
        } else {
            for (i in 0 until half) {
                tmp[tmpOff + i] = buf[off + 2 * i]
                tmp[tmpOff + half + i] = buf[off + 2 * i + 1]
            }
            for (i in 1 until half) tmp[tmpOff + n - i] += tmp[tmpOff + n - i - 1]
            tmp[tmpOff + half] *= SQRT2
            dctInPlace(tmp, tmpOff, half, buf, off, false)
            dctInPlace(tmp, tmpOff + half, half, buf, off + half, false)
            for (i in 0 until half) {
                val a = tmp[tmpOff + i]
                val b = tmp[tmpOff + half + i] * sec[i]
                buf[off + i] = a + b
                buf[off + n - i - 1] = a - b
            }
        }
    }

    private fun dct4(io: FloatArray, off: Int, forward: Boolean) {
        val sec0 = 0.5411961f
        val sec1 = 1.306563f
        val i0 = io[off]
        val i1 = io[off + 1]
        val i2 = io[off + 2]
        val i3 = io[off + 3]
        if (forward) {
            val sum03 = i0 + i3
            val sum12 = i1 + i2
            val tmp0 = (i0 - i3) * sec0
            val tmp1 = (i1 - i2) * sec1
            val out0 = (tmp0 + tmp1) / 4f
            val out1 = (tmp0 - tmp1) / 4f
            io[off] = (sum03 + sum12) / 4f
            io[off + 1] = out0 * SQRT2 + out1
            io[off + 2] = (sum03 - sum12) / 4f
            io[off + 3] = out1
        } else {
            val tmp0 = i1 * SQRT2
            val tmp1 = i1 + i3
            val out0 = (tmp0 + tmp1) * sec0
            val out1 = (tmp0 - tmp1) * sec1
            val sum02 = i0 + i2
            val sub02 = i0 - i2
            io[off] = sum02 + out0
            io[off + 1] = sub02 + out1
            io[off + 2] = sub02 - out1
            io[off + 3] = sum02 - out0
        }
    }

    private fun idct8(io: FloatArray, off: Int) {
        val sec0 = 0.5411961f
        val sec1 = 1.306563f
        // 偶数番目
        val e0 = io[off]
        val e1 = io[off + 2]
        val e2 = io[off + 4]
        val e3 = io[off + 6]
        var t0 = e1 * SQRT2
        var t1 = e1 + e3
        var o0 = (t0 + t1) * sec0
        var o1 = (t0 - t1) * sec1
        var s02 = e0 + e2
        var d02 = e0 - e2
        val a0 = s02 + o0
        val a1 = d02 + o1
        val a2 = d02 - o1
        val a3 = s02 - o0
        // 奇数番目
        val q0 = io[off + 1] * SQRT2
        val q1 = io[off + 3] + io[off + 1]
        val q2 = io[off + 5] + io[off + 3]
        val q3 = io[off + 7] + io[off + 5]
        t0 = q1 * SQRT2
        t1 = q1 + q3
        o0 = (t0 + t1) * sec0
        o1 = (t0 - t1) * sec1
        s02 = q0 + q2
        d02 = q0 - q2
        val b0 = (s02 + o0) * SEC8[0]
        val b1 = (d02 + o1) * SEC8[1]
        val b2 = (d02 - o1) * SEC8[2]
        val b3 = (s02 - o0) * SEC8[3]
        io[off] = a0 + b0
        io[off + 7] = a0 - b0
        io[off + 1] = a1 + b1
        io[off + 6] = a1 - b1
        io[off + 2] = a2 + b2
        io[off + 5] = a2 - b2
        io[off + 3] = a3 + b3
        io[off + 4] = a3 - b3
    }

    /**
     * buf の (x0, y0) から幅 w・高さ h の領域 (行の長さ stride) を 2 次元 DCT する。
     * scratch は 2 * max(w, h) 以上。
     */
    fun dct2d(buf: FloatArray, base: Int, stride: Int, w: Int, h: Int, scratch: FloatArray, forward: Boolean) {
        if (w * h <= 1) return
        if (!forward && w == 8 && h == 8 && scratch.size >= 128) {
            idct8x8(buf, base, stride, scratch)
            return
        }
        // 行方向
        if (w > 1) {
            for (y in 0 until h) {
                val off = base + y * stride
                if (w == 2) {
                    val a = buf[off]
                    val b = buf[off + 1]
                    if (forward) {
                        buf[off] = (a + b) / 2f
                        buf[off + 1] = (a - b) / 2f
                    } else {
                        buf[off] = a + b
                        buf[off + 1] = a - b
                    }
                } else {
                    dctInPlace(buf, off, w, scratch, 0, forward)
                }
            }
        }
        // 列方向
        if (h > 1) {
            val col = scratch
            for (x in 0 until w) {
                for (y in 0 until h) col[y] = buf[base + y * stride + x]
                dctInPlace(col, 0, h, col, h, forward)
                for (y in 0 until h) buf[base + y * stride + x] = col[y]
            }
        }
    }

    private val SEC16 = secHalf(16)
    private val SEC32 = secHalf(32)

    /** 16 点の逆 DCT (tmp[t until t+16] を使う) */
    private fun idct16(buf: FloatArray, off: Int, tmp: FloatArray, t: Int) {
        for (i in 0 until 8) {
            tmp[t + i] = buf[off + 2 * i]
            tmp[t + 8 + i] = buf[off + 2 * i + 1]
        }
        for (i in 15 downTo 9) tmp[t + i] += tmp[t + i - 1]
        tmp[t + 8] *= SQRT2
        idct8(tmp, t)
        idct8(tmp, t + 8)
        for (i in 0 until 8) {
            val a = tmp[t + i]
            val b = tmp[t + 8 + i] * SEC16[i]
            buf[off + i] = a + b
            buf[off + 15 - i] = a - b
        }
    }

    /** 32 点の逆 DCT (tmp[t until t+32] を使う。内側の変換は空いた buf を作業領域にする) */
    private fun idct32(buf: FloatArray, off: Int, tmp: FloatArray, t: Int) {
        for (i in 0 until 16) {
            tmp[t + i] = buf[off + 2 * i]
            tmp[t + 16 + i] = buf[off + 2 * i + 1]
        }
        for (i in 31 downTo 17) tmp[t + i] += tmp[t + i - 1]
        tmp[t + 16] *= SQRT2
        idct16(tmp, t, buf, off)
        idct16(tmp, t + 16, buf, off + 16)
        for (i in 0 until 16) {
            val a = tmp[t + i]
            val b = tmp[t + 16 + i] * SEC32[i]
            buf[off + i] = a + b
            buf[off + 31 - i] = a - b
        }
    }

    /** 8x8 の逆 DCT。行を変換して転置し、もう一度行を変換して転置して戻す */
    private fun idct8x8(buf: FloatArray, base: Int, stride: Int, t: FloatArray) {
        for (y in 0 until 8) {
            val o = base + y * stride
            val r = y * 8
            for (x in 0 until 8) t[r + x] = buf[o + x]
            idct8(t, r)
        }
        for (y in 0 until 8) for (x in 0 until 8) t[64 + x * 8 + y] = t[y * 8 + x]
        for (x in 0 until 8) idct8(t, 64 + x * 8)
        for (y in 0 until 8) {
            val o = base + y * stride
            for (x in 0 until 8) buf[o + x] = t[64 + x * 8 + y]
        }
    }

    // ---- 小さい変換 ----

    private fun auxIdct2(buf: FloatArray, base: Int, stride: Int, size: Int, tmp: FloatArray) {
        val half = size / 2
        for (y in 0 until half) for (x in 0 until half) {
            val c00 = buf[base + y * stride + x]
            val c01 = buf[base + y * stride + x + half]
            val c10 = buf[base + (y + half) * stride + x]
            val c11 = buf[base + (y + half) * stride + x + half]
            tmp[(2 * y) * size + 2 * x] = c00 + c01 + c10 + c11
            tmp[(2 * y) * size + 2 * x + 1] = c00 + c01 - c10 - c11
            tmp[(2 * y + 1) * size + 2 * x] = c00 - c01 + c10 - c11
            tmp[(2 * y + 1) * size + 2 * x + 1] = c00 - c01 - c10 + c11
        }
        for (y in 0 until size) for (x in 0 until size) buf[base + y * stride + x] = tmp[y * size + x]
    }

    fun dct2(buf: FloatArray, base: Int, stride: Int, tmp: FloatArray) {
        auxIdct2(buf, base, stride, 2, tmp)
        auxIdct2(buf, base, stride, 4, tmp)
        auxIdct2(buf, base, stride, 8, tmp)
    }

    /** 小さい変換の作業領域 (ブロックごとに確保しないため) */
    class Scratch {
        val s = FloatArray(64)
        val sc = FloatArray(16)
        val coeff = FloatArray(16)
        val samples = FloatArray(16)
        val s44 = FloatArray(16)
        val s48 = FloatArray(32)
    }

    fun dct4(buf: FloatArray, base: Int, stride: Int, tmp: FloatArray, ds: Scratch) {
        auxIdct2(buf, base, stride, 2, tmp)
        val s = ds.s
        val sc = ds.sc
        for (y in 0 until 2) for (x in 0 until 2) {
            val o = (y * 2 + x) * 16
            for (iy in 0 until 4) for (ix in 0 until 4) {
                // (iy, ix) の位置に (x + ix*2, y + iy*2) を転置して置く
                s[o + ix * 4 + iy] = buf[base + (y + iy * 2) * stride + x + ix * 2]
            }
            dct2d(s, o, 4, 4, 4, sc, false)
        }
        for (y in 0 until 2) for (x in 0 until 2) {
            val o = (y * 2 + x) * 16
            for (iy in 0 until 4) for (ix in 0 until 4) buf[base + (y * 4 + iy) * stride + x * 4 + ix] = s[o + iy * 4 + ix]
        }
    }

    fun hornuss(buf: FloatArray, base: Int, stride: Int, tmp: FloatArray, ds: Scratch) {
        auxIdct2(buf, base, stride, 2, tmp)
        val s = ds.s
        for (y in 0 until 2) for (x in 0 until 2) {
            val o = (y * 2 + x) * 16
            for (iy in 0 until 4) for (ix in 0 until 4) s[o + iy * 4 + ix] = buf[base + (y + iy * 2) * stride + x + ix * 2]
            var sum = 0f
            for (i in 1 until 16) sum += s[o + i]
            val avg = s[o] - sum / 16f
            s[o] = s[o + 5]
            s[o + 5] = 0f
            for (i in 0 until 16) s[o + i] += avg
        }
        for (y in 0 until 2) for (x in 0 until 2) {
            val o = (y * 2 + x) * 16
            for (iy in 0 until 4) for (ix in 0 until 4) buf[base + (y * 4 + iy) * stride + x * 4 + ix] = s[o + iy * 4 + ix]
        }
    }

    fun dct4x8(buf: FloatArray, base: Int, stride: Int, transposed: Boolean, ds: Scratch) {
        val c0 = buf[base]
        val c1 = buf[base + stride]
        buf[base] = c0 + c1
        buf[base + stride] = c0 - c1
        val s = ds.s
        val sc = ds.sc
        for (idx in 0 until 2) {
            val o = idx * 32
            for (iy in 0 until 4) for (ix in 0 until 8) s[o + iy * 8 + ix] = buf[base + (iy * 2 + idx) * stride + ix]
            dct2d(s, o, 8, 8, 4, sc, false)
        }
        if (transposed) {
            for (y in 0 until 8) for (x in 0 until 8) buf[base + x * stride + y] = s[y * 8 + x]
        } else {
            for (y in 0 until 8) for (x in 0 until 8) buf[base + y * stride + x] = s[y * 8 + x]
        }
    }

    fun afv(buf: FloatArray, base: Int, stride: Int, n: Int, ds: Scratch) {
        val flipX = n % 2
        val flipY = n / 2
        fun g(x: Int, y: Int) = buf[base + y * stride + x]
        val coeff = ds.coeff
        coeff[0] = (g(0, 0) + g(1, 0) + g(0, 1)) * 4f
        for (idx in 1 until 16) coeff[idx] = g(2 * (idx % 4), 2 * (idx / 4))
        val samples = ds.samples
        samples.fill(0f)
        for (k in 0 until 16) {
            val c = coeff[k]
            val basis = AFV_BASIS[k]
            for (i in 0 until 16) samples[i] += c * basis[i]
        }
        val s44 = ds.s44
        val s48 = ds.s48
        val sc = ds.sc
        s44[0] = g(0, 0) - g(1, 0) + g(0, 1)
        for (iy in 0 until 4) for (ix in 0 until 4) {
            if (ix or iy == 0) continue
            s44[ix * 4 + iy] = g(2 * ix + 1, 2 * iy)
        }
        dct2d(s44, 0, 4, 4, 4, sc, false)
        s48[0] = g(0, 0) - g(0, 1)
        for (iy in 0 until 4) for (ix in 0 until 8) {
            if (ix or iy == 0) continue
            s48[iy * 8 + ix] = g(ix, 2 * iy + 1)
        }
        dct2d(s48, 0, 8, 8, 4, sc, false)
        for (iy in 0 until 4) {
            val ay = if (flipY == 0) iy else 3 - iy
            for (ix in 0 until 4) {
                val ax = if (flipX == 0) ix else 3 - ix
                buf[base + (flipY * 4 + iy) * stride + flipX * 4 + ix] = samples[ay * 4 + ax]
            }
        }
        for (iy in 0 until 4) {
            val y = flipY * 4 + iy
            for (ix in 0 until 4) buf[base + y * stride + (1 - flipX) * 4 + ix] = s44[iy * 4 + ix]
        }
        for (iy in 0 until 4) {
            val y = (1 - flipY) * 4 + iy
            for (ix in 0 until 8) buf[base + y * stride + ix] = s48[iy * 8 + ix]
        }
    }

    private val AFV_BASIS: Array<FloatArray> = arrayOf(
        FloatArray(16) { 0.25f },
        floatArrayOf(
            0.876902929799142f, 0.2206518106944235f, -0.10140050393753763f, -0.1014005039375375f,
            0.2206518106944236f, -0.10140050393753777f, -0.10140050393753772f, -0.10140050393753763f,
            -0.10140050393753758f, -0.10140050393753769f, -0.1014005039375375f, -0.10140050393753768f,
            -0.10140050393753768f, -0.10140050393753759f, -0.10140050393753763f, -0.10140050393753741f,
        ),
        floatArrayOf(
            0.0f, 0.0f, 0.40670075830260755f, 0.44444816619734445f,
            0.0f, 0.0f, 0.19574399372042936f, 0.2929100136981264f,
            -0.40670075830260716f, -0.19574399372042872f, 0.0f, 0.11379074460448091f,
            -0.44444816619734384f, -0.29291001369812636f, -0.1137907446044814f, 0.0f,
        ),
        floatArrayOf(
            0.0f, 0.0f, -0.21255748058288748f, 0.3085497062849767f,
            0.0f, 0.4706702258572536f, -0.1621205195722993f, 0.0f,
            -0.21255748058287047f, -0.16212051957228327f, -0.47067022585725277f, -0.1464291867126764f,
            0.3085497062849487f, 0.0f, -0.14642918671266536f, 0.4251149611657548f,
        ),
        floatArrayOf(
            0.0f, -0.7071067811865474f, 0.0f, 0.0f,
            0.7071067811865476f, 0.0f, 0.0f, 0.0f,
            0.0f, 0.0f, 0.0f, 0.0f,
            0.0f, 0.0f, 0.0f, 0.0f,
        ),
        floatArrayOf(
            -0.4105377591765233f, 0.6235485373547691f, -0.06435071657946274f, -0.06435071657946266f,
            0.6235485373547694f, -0.06435071657946284f, -0.0643507165794628f, -0.06435071657946274f,
            -0.06435071657946272f, -0.06435071657946279f, -0.06435071657946266f, -0.06435071657946277f,
            -0.06435071657946277f, -0.06435071657946273f, -0.06435071657946274f, -0.0643507165794626f,
        ),
        floatArrayOf(
            0.0f, 0.0f, -0.4517556589999482f, 0.15854503551840063f,
            0.0f, -0.04038515160822202f, 0.0074182263792423875f, 0.39351034269210167f,
            -0.45175565899994635f, 0.007418226379244351f, 0.1107416575309343f, 0.08298163094882051f,
            0.15854503551839705f, 0.3935103426921022f, 0.0829816309488214f, -0.45175565899994796f,
        ),
        floatArrayOf(
            0.0f, 0.0f, -0.304684750724869f, 0.5112616136591823f,
            0.0f, 0.0f, -0.290480129728998f, -0.06578701549142804f,
            0.304684750724884f, 0.2904801297290076f, 0.0f, -0.23889773523344604f,
            -0.5112616136592012f, 0.06578701549142545f, 0.23889773523345467f, 0.0f,
        ),
        floatArrayOf(
            0.0f, 0.0f, 0.3017929516615495f, 0.25792362796341184f,
            0.0f, 0.16272340142866204f, 0.09520022653475037f, 0.0f,
            0.3017929516615503f, 0.09520022653475055f, -0.16272340142866173f, -0.35312385449816297f,
            0.25792362796341295f, 0.0f, -0.3531238544981624f, -0.6035859033230976f,
        ),
        floatArrayOf(
            0.0f, 0.0f, 0.40824829046386274f, 0.0f,
            0.0f, 0.0f, 0.0f, -0.4082482904638628f,
            -0.4082482904638635f, 0.0f, 0.0f, -0.40824829046386296f,
            0.0f, 0.4082482904638634f, 0.408248290463863f, 0.0f,
        ),
        floatArrayOf(
            0.0f, 0.0f, 0.1747866975480809f, 0.0812611176717539f,
            0.0f, 0.0f, -0.3675398009862027f, -0.307882213957909f,
            -0.17478669754808135f, 0.3675398009862011f, 0.0f, 0.4826689115059883f,
            -0.08126111767175039f, 0.30788221395790305f, -0.48266891150598584f, 0.0f,
        ),
        floatArrayOf(
            0.0f, 0.0f, -0.21105601049335784f, 0.18567180916109802f,
            0.0f, 0.0f, 0.49215859013738733f, -0.38525013709251915f,
            0.21105601049335806f, -0.49215859013738905f, 0.0f, 0.17419412659916217f,
            -0.18567180916109904f, 0.3852501370925211f, -0.1741941265991621f, 0.0f,
        ),
        floatArrayOf(
            0.0f, 0.0f, -0.14266084808807264f, -0.3416446842253372f,
            0.0f, 0.7367497537172237f, 0.24627107722075148f, -0.08574019035519306f,
            -0.14266084808807344f, 0.24627107722075137f, 0.14883399227113567f, -0.04768680350229251f,
            -0.3416446842253373f, -0.08574019035519267f, -0.047686803502292804f, -0.14266084808807242f,
        ),
        floatArrayOf(
            0.0f, 0.0f, -0.13813540350758585f, 0.3302282550303788f,
            0.0f, 0.08755115000587084f, -0.07946706605909573f, -0.4613374887461511f,
            -0.13813540350758294f, -0.07946706605910261f, 0.49724647109535086f, 0.12538059448563663f,
            0.3302282550303805f, -0.4613374887461554f, 0.12538059448564315f, -0.13813540350758452f,
        ),
        floatArrayOf(
            0.0f, 0.0f, -0.17437602599651067f, 0.0702790691196284f,
            0.0f, -0.2921026642334881f, 0.3623817333531167f, 0.0f,
            -0.1743760259965108f, 0.36238173335311646f, 0.29210266423348785f, -0.4326608024727445f,
            0.07027906911962818f, 0.0f, -0.4326608024727457f, 0.34875205199302267f,
        ),
        floatArrayOf(
            0.0f, 0.0f, 0.11354987314994337f, -0.07417504595810355f,
            0.0f, 0.19402893032594343f, -0.435190496523228f, 0.21918684838857466f,
            0.11354987314994257f, -0.4351904965232251f, 0.5550443808910661f, -0.25468277124066463f,
            -0.07417504595810233f, 0.2191868483885728f, -0.25468277124066413f, 0.1135498731499429f,
        ),
    )
}
