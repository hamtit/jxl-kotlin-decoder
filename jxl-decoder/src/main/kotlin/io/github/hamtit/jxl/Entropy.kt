package io.github.hamtit.jxl

/** hybrid-uint の設定 */
internal class IntegerConfig(val splitExponent: Int, val msbInToken: Int, val lsbInToken: Int) {
    val split: Int = 1 shl splitExponent

    companion object {
        fun parse(br: BitReader, logAlphabetSize: Int): IntegerConfig {
            val splitExponent = br.readBits(addLog2Ceil(logAlphabetSize))
            var msb = 0
            var lsb = 0
            if (splitExponent != logAlphabetSize) {
                msb = br.readBits(addLog2Ceil(splitExponent))
                if (msb > splitExponent) fail("invalid integer config")
                lsb = br.readBits(addLog2Ceil(splitExponent - msb))
            }
            if (lsb + msb > splitExponent) fail("invalid integer config")
            return IntegerConfig(splitExponent, msb, lsb)
        }
    }
}

/** ceil(log2(x + 1)) */
internal fun addLog2Ceil(x: Int): Int {
    val v = x + 1
    return if (v <= 1) 0 else 32 - Integer.numberOfLeadingZeros(v - 1)
}

/** ANS のヒストグラム構築の作業領域 (クラスタ間で使い回す) */
internal class AnsParseScratch(tableSize: Int) {
    val dist = IntArray(tableSize)
    val repeatStart = IntArray(tableSize)
    val repeatEnd = IntArray(tableSize)
    val aliasSymbol = IntArray(tableSize)
    val aliasOffset = IntArray(tableSize)
    val aliasCutoff = IntArray(tableSize)
    val underfull = IntArray(tableSize)
    val overfull = IntArray(tableSize)
}

/** ANS のヒストグラム (エイリアス表) のパーサー */
internal object AnsHistogram {
    private fun readU8(br: BitReader): Int =
        if (br.readBool()) {
            val n = br.readBits(3)
            (1 shl n) + br.readBits(n)
        } else 0

    private fun readPrefix(br: BitReader): Int = when (br.readBits(3)) {
        0 -> 10
        1 -> {
            if (br.readBool()) 4
            else if (br.readBool()) 0
            else if (br.readBool()) 11
            else if (br.readBool()) 13
            else 12
        }
        2 -> 7
        3 -> if (br.readBool()) 1 else 3
        4 -> 6
        5 -> 8
        6 -> 9
        else -> if (br.readBool()) 2 else 5
    }

    /**
     * 1 クラスタ分のエイリアス表を [buckets] の [base] から書き込む。
     * 単一シンボルのヒストグラムならそのシンボル、そうでなければ -1 を返す。
     */
    fun parseInto(br: BitReader, logAlphabetSize: Int, buckets: LongArray, base: Int, sc: AnsParseScratch): Int {
        val tableSize = 1 shl logAlphabetSize
        val logBucketSize = 12 - logAlphabetSize
        val bucketSize = 1 shl logBucketSize
        val dist = sc.dist
        dist.fill(0, 0, tableSize)
        val alphabetSize: Int
        if (br.readBool()) {
            if (br.readBool()) {
                val v0 = readU8(br)
                val v1 = readU8(br)
                if (v0 == v1) fail("invalid ANS histogram")
                alphabetSize = maxOf(v0, v1) + 1
                if (alphabetSize > tableSize) fail("invalid ANS histogram")
                val prob = br.readBits(12)
                dist[v0] = prob
                dist[v1] = 4096 - prob
            } else {
                val v = readU8(br)
                alphabetSize = v + 1
                if (alphabetSize > tableSize) fail("invalid ANS histogram")
                dist[v] = 4096
            }
        } else if (br.readBool()) {
            alphabetSize = readU8(br) + 1
            if (alphabetSize > tableSize) fail("invalid ANS histogram")
            val baseVal = 4096 / alphabetSize
            val leftover = 4096 % alphabetSize
            for (i in 0 until alphabetSize) dist[i] = if (i < leftover) baseVal + 1 else baseVal
        } else {
            var len = 0
            while (len < 3) {
                if (br.readBool()) len++ else break
            }
            val shift = br.readBits(len) + (1 shl len) - 1
            if (shift > 13) fail("invalid ANS histogram")
            alphabetSize = readU8(br) + 3
            if (alphabetSize > tableSize) fail("invalid ANS histogram")
            // 13 は直前の値の繰り返し
            val repeatStart = sc.repeatStart
            val repeatEnd = sc.repeatEnd
            var numRepeat = 0
            var omitLog = -1
            var omitPos = -1
            var idx = 0
            while (idx < alphabetSize) {
                dist[idx] = readPrefix(br)
                if (dist[idx] == 13) {
                    val repeatCount = readU8(br) + 4
                    if (idx + repeatCount > alphabetSize) fail("invalid ANS histogram")
                    repeatStart[numRepeat] = idx
                    repeatEnd[numRepeat] = idx + repeatCount
                    numRepeat++
                    idx += repeatCount
                    continue
                }
                if (omitPos < 0 || dist[idx] > omitLog) {
                    omitLog = dist[idx]
                    omitPos = idx
                }
                idx++
            }
            if (omitPos < 0) fail("invalid ANS histogram")
            if (omitPos + 1 < tableSize && dist[omitPos + 1] == 13) fail("invalid ANS histogram")
            var rIdx = 0
            var acc = 0
            var prev = 0
            for (i in 0 until tableSize) {
                if (rIdx < numRepeat && repeatStart[rIdx] <= i) {
                    if (repeatEnd[rIdx] == i) {
                        rIdx++
                    } else {
                        dist[i] = prev
                        acc += prev
                        if (acc > 4096) fail("invalid ANS histogram")
                        continue
                    }
                }
                val code = dist[i]
                if (code == 0 || i == omitPos) {
                    prev = 0
                    continue
                }
                var value = code
                if (code > 1) {
                    val zeros = code - 1
                    val bitcount = (shift - ((12 - zeros) shr 1)).coerceIn(0, zeros)
                    value = (1 shl zeros) + (br.readBits(bitcount) shl (zeros - bitcount))
                }
                dist[i] = value
                prev = value
                acc += value
                if (acc > 4096) fail("invalid ANS histogram")
            }
            dist[omitPos] = 4096 - acc
            if (dist[omitPos] < 0) fail("invalid ANS histogram")
        }

        var single = -1
        for (i in 0 until tableSize) if (dist[i] == 4096) { single = i; break }
        if (single >= 0) {
            for (i in 0 until tableSize) {
                val d = dist[i]
                buckets[base + i] = pack(single, 0, d, bucketSize * i, d xor 4096)
            }
            return single
        }

        val aliasSymbol = sc.aliasSymbol
        val aliasOffset = sc.aliasOffset
        val aliasCutoff = sc.aliasCutoff
        for (i in 0 until tableSize) {
            aliasSymbol[i] = if (i < alphabetSize) i else 0
            aliasOffset[i] = 0
            aliasCutoff[i] = dist[i]
        }
        val underfull = sc.underfull
        var nu = 0
        val overfull = sc.overfull
        var no = 0
        for (i in 0 until tableSize) {
            if (dist[i] < bucketSize) underfull[nu++] = i
            else if (dist[i] > bucketSize) overfull[no++] = i
        }
        while (no > 0 && nu > 0) {
            val o = overfull[--no]
            val u = underfull[--nu]
            val by = bucketSize - aliasCutoff[u]
            aliasCutoff[o] -= by
            aliasSymbol[u] = o
            aliasOffset[u] = aliasCutoff[o]
            if (aliasCutoff[o] < bucketSize) underfull[nu++] = o
            else if (aliasCutoff[o] > bucketSize) overfull[no++] = o
        }
        for (i in 0 until tableSize) {
            buckets[base + i] = if (aliasCutoff[i] == bucketSize) {
                pack(i, 0, dist[i], 0, 0)
            } else {
                pack(
                    aliasSymbol[i], aliasCutoff[i], dist[i],
                    aliasOffset[i] - aliasCutoff[i], dist[i] xor dist[aliasSymbol[i]],
                )
            }
        }
        return -1
    }

    private fun pack(sym: Int, cutoff: Int, dist: Int, offset: Int, xor: Int): Long =
        (sym.toLong() and 0xff) or ((cutoff.toLong() and 0xff) shl 8) or
            ((dist.toLong() and 0xffff) shl 16) or ((offset.toLong() and 0xffff) shl 32) or
            ((xor.toLong() and 0xffff) shl 48)
}

/** プレフィックス符号 (Brotli 形式) のヒストグラム */
internal class PrefixHistogram private constructor(
    private val toplevelBits: Int,
    private val toplevel: IntArray,
    private val second: IntArray,
) {
    // エントリ: bit 31 = nested, bits 16..23 = ビット数 または マスク, bits 0..15 = シンボル または オフセット
    val singleSymbol: Int
        get() = if (toplevel.size == 1 && (toplevel[0] ushr 16) and 0xff == 0 && toplevel[0] >= 0) {
            toplevel[0] and 0xffff
        } else -1

    fun readSymbol(br: BitReader): Int {
        val peeked = br.peek(15)
        val e = toplevel[peeked and ((1 shl toplevelBits) - 1)]
        if (e < 0) {
            val mask = (e ushr 16) and 0xff
            val off = (e and 0xffff) + ((peeked ushr toplevelBits) and mask)
            val s = second[off]
            br.consume((s ushr 16) and 0xff)
            return s and 0xffff
        }
        br.consume((e ushr 16) and 0xff)
        return e and 0xffff
    }

    companion object {
        private const val MAX_BITS = 15
        private const val MAX_TOPLEVEL_BITS = 10

        private fun entry(nested: Boolean, bits: Int, sym: Int): Int =
            (if (nested) Int.MIN_VALUE else 0) or ((bits and 0xff) shl 16) or (sym and 0xffff)

        private fun single(sym: Int) = PrefixHistogram(0, intArrayOf(entry(false, 0, sym)), IntArray(0))

        private fun reverseInto(v: IntArray, len: Int, out: IntArray, outOff: Int) {
            val bits = Integer.numberOfTrailingZeros(len)
            for (i in 0 until len) {
                val rev = if (bits == 0) 0 else Integer.reverse(i) ushr (32 - bits)
                out[outOff + i] = v[rev]
            }
        }

        private fun withCodeLengths(lengths: IntArray): PrefixHistogram {
            var maxLen = 0
            for (l in lengths) if (l > maxLen) maxLen = l
            val symsForLength = Array(maxLen) { ArrayList<Int>() }
            for (sym in lengths.indices) {
                val l = lengths[sym]
                if (l > 0) symsForLength[l - 1].add(sym)
            }
            val toplevelBits = minOf(maxLen, MAX_TOPLEVEL_BITS)
            val entries = IntArray(1 shl toplevelBits)
            var current = 0
            for (idx in 0 until toplevelBits) {
                val shifts = toplevelBits - 1 - idx
                for (sym in symsForLength[idx]) {
                    val e = entry(false, idx + 1, sym)
                    val n = 1 shl shifts
                    if (current + n > entries.size) fail("invalid prefix code")
                    java.util.Arrays.fill(entries, current, current + n, e)
                    current += n
                }
            }
            val second = ArrayList<Int>()
            if (toplevelBits < maxLen) {
                var remaining = IntArray(0)
                var remainingBits = 0
                for (idx in toplevelBits until maxLen) {
                    val syms = symsForLength[idx]
                    if (syms.isEmpty()) continue
                    val chunkBits = idx + 1 - toplevelBits
                    val chunkSize = 1 shl chunkBits
                    var chunk = IntArray(chunkSize)
                    var chunkLen = 0
                    if (remaining.isNotEmpty()) {
                        val mult = 1 shl (chunkBits - remainingBits)
                        for (e in remaining) {
                            for (k in 0 until mult) chunk[chunkLen++] = e
                        }
                    }
                    for (sym in syms) {
                        chunk[chunkLen++] = entry(false, idx + 1, sym)
                        if (chunkLen == chunkSize) {
                            if (current >= entries.size) fail("invalid prefix code")
                            entries[current] = entry(true, chunkSize - 1, second.size)
                            val tmp = IntArray(chunkSize)
                            reverseInto(chunk, chunkSize, tmp, 0)
                            for (t in tmp) second.add(t)
                            current++
                            chunk = IntArray(chunkSize)
                            chunkLen = 0
                        }
                    }
                    remaining = chunk.copyOf(chunkLen)
                    remainingBits = chunkBits
                }
                if (remaining.isNotEmpty()) fail("invalid prefix code")
            }
            if (current != (1 shl toplevelBits)) fail("invalid prefix code")
            val top = IntArray(entries.size)
            reverseInto(entries, entries.size, top, 0)
            return PrefixHistogram(toplevelBits, top, second.toIntArray())
        }

        fun parse(br: BitReader, alphabetSize: Int): PrefixHistogram {
            if (alphabetSize == 1) return single(0)
            if (alphabetSize > (1 shl MAX_BITS)) fail("prefix alphabet too large")
            val hskip = br.readBits(2)
            return if (hskip == 1) parseSimple(br, alphabetSize) else parseComplex(br, alphabetSize, hskip)
        }

        private fun parseSimple(br: BitReader, alphabetSize: Int): PrefixHistogram {
            val alphabetBits = 32 - Integer.numberOfLeadingZeros(alphabetSize - 1)
            val nsym = br.readBits(2) + 1
            val syms: IntArray
            val lens: IntArray
            when (nsym) {
                1 -> {
                    val s = br.readBits(alphabetBits)
                    if (s >= alphabetSize) fail("invalid prefix code")
                    return single(s)
                }
                2 -> {
                    syms = intArrayOf(br.readBits(alphabetBits), br.readBits(alphabetBits))
                    lens = intArrayOf(1, 1)
                }
                3 -> {
                    syms = intArrayOf(br.readBits(alphabetBits), br.readBits(alphabetBits), br.readBits(alphabetBits))
                    lens = intArrayOf(1, 2, 2)
                }
                else -> {
                    syms = IntArray(4) { br.readBits(alphabetBits) }
                    lens = if (br.readBool()) intArrayOf(1, 2, 3, 3) else intArrayOf(2, 2, 2, 2)
                }
            }
            val codeLengths = IntArray(alphabetSize)
            for (i in syms.indices) {
                if (syms[i] >= alphabetSize) fail("invalid prefix code")
                codeLengths[syms[i]] = lens[i]
            }
            return withCodeLengths(codeLengths)
        }

        private val CODE_LENGTH_ORDER = intArrayOf(1, 2, 3, 4, 0, 5, 17, 6, 16, 7, 8, 9, 10, 11, 12, 13, 14, 15)

        private fun parseComplex(br: BitReader, alphabetSize: Int, hskip: Int): PrefixHistogram {
            val clcl = IntArray(18)
            var bitacc = 0
            var nonzeroCount = 0
            var nonzeroSym = 0
            for (i in hskip until 18) {
                val idx = CODE_LENGTH_ORDER[i]
                // U32(0, 4, 3, 8)、8 のときは追加ビットで 2 / 1 / 5
                val len = when (br.peek(4)) {
                    0, 4, 8, 12 -> { br.consume(2); 0 }
                    1, 5, 9, 13 -> { br.consume(2); 4 }
                    2, 6, 10, 14 -> { br.consume(2); 3 }
                    3, 11 -> { br.consume(3); 2 }
                    7 -> { br.consume(4); 1 }
                    else -> { br.consume(4); 5 }
                }
                clcl[idx] = len
                if (len != 0) {
                    nonzeroCount++
                    nonzeroSym = idx
                    bitacc += 32 shr len
                    if (bitacc == 32) break
                    if (bitacc > 32) fail("invalid prefix code")
                }
            }
            val clHist = if (nonzeroCount == 1) single(nonzeroSym)
            else if (bitacc != 32) fail("invalid prefix code")
            else withCodeLengths(clcl)

            val codeLengths = IntArray(alphabetSize)
            var acc = 0
            var prevSym = 8
            var lastNonzero = 8
            var lastRepeatCount = 0
            var repeatCount = 0
            var repeatSym = 0
            var i = 0
            while (i < alphabetSize) {
                var len: Int
                if (repeatCount > 0) {
                    len = repeatSym
                    repeatCount--
                } else {
                    val sym = clHist.readSymbol(br)
                    when {
                        sym == 0 -> len = 0
                        sym <= 15 -> {
                            len = sym
                            lastNonzero = sym
                        }
                        sym == 16 -> {
                            repeatCount = br.readBits(2) + 3
                            if (prevSym == 16) {
                                repeatCount += lastRepeatCount * 3 - 8
                                lastRepeatCount += repeatCount
                            } else {
                                lastRepeatCount = repeatCount
                            }
                            repeatSym = lastNonzero
                            len = repeatSym
                            repeatCount--
                        }
                        else -> {
                            repeatCount = br.readBits(3) + 3
                            if (prevSym == 17) {
                                repeatCount += lastRepeatCount * 7 - 16
                                lastRepeatCount += repeatCount
                            } else {
                                lastRepeatCount = repeatCount
                            }
                            repeatSym = 0
                            len = 0
                            repeatCount--
                        }
                    }
                    prevSym = sym
                }
                codeLengths[i] = len
                i++
                if (len != 0) {
                    acc += 1 shl (MAX_BITS - len)
                    if (acc > (1 shl MAX_BITS)) fail("invalid prefix code")
                    if (acc == (1 shl MAX_BITS) && repeatCount == 0) break
                }
            }
            if (acc != (1 shl MAX_BITS) || repeatCount > 0) fail("invalid prefix code")
            return withCodeLengths(codeLengths)
        }
    }
}

/**
 * エントロピー符号のパラメータ (クラスタ・ヒストグラム・LZ77)。
 * 読み出し状態は [EntropyReader] が持つので、同じ符号を複数スレッドから使える。
 */
internal class EntropyCode private constructor(
    val clusters: ByteArray,
    val configs: Array<IntegerConfig>,
    private val ansBuckets: LongArray?,
    private val ansLogAlphabet: Int,
    private val ansSingles: ShortArray?,
    private val prefix: Array<PrefixHistogram>?,
    val lz77Enabled: Boolean,
    val lz77MinSymbol: Int,
    val lz77MinLength: Int,
    val lz77LenConfig: IntegerConfig?,
) {
    val usesAns: Boolean get() = ansBuckets != null

    /**
     * ANS で LZ77 なしのときの、全クラスタのエイリアス表を並べた表 (クラスタ c の表は c shl [fastLogAlphabet] から)。
     * [EntropyReader.readAnsFast] で読むのに使う。それ以外は null
     */
    val fastBuckets: LongArray? = if (ansBuckets != null && !lz77Enabled) ansBuckets else null
    val fastLogAlphabet: Int = if (fastBuckets != null) ansLogAlphabet else 0
    /** クラスタ c の hybrid-uint の設定 (split・splitExponent・msbInToken・lsbInToken) を c * 4 から並べた表 */
    val fastConfigs: IntArray = IntArray(configs.size * 4).also { fc ->
        for ((i, c) in configs.withIndex()) {
            fc[4 * i] = c.split
            fc[4 * i + 1] = c.splitExponent
            fc[4 * i + 2] = c.msbInToken
            fc[4 * i + 3] = c.lsbInToken
        }
    }

    fun newReader(): EntropyReader = EntropyReader(this)

    fun readSymbol(br: BitReader, cluster: Int, reader: EntropyReader): Int {
        val ab = ansBuckets
        if (ab != null) {
            if (reader.ansInitial) {
                reader.ansState = br.readBits(32).toLong() and 0xffffffffL
                reader.ansInitial = false
            }
            val state = reader.ansState
            val idx = (state and 0xfff).toInt()
            val logBucket = 12 - ansLogAlphabet
            val i = idx ushr logBucket
            val pos = idx and ((1 shl logBucket) - 1)
            val b = ab[(cluster shl ansLogAlphabet) + i]
            val cutoff = ((b ushr 8) and 0xff).toInt()
            var dist = ((b ushr 16) and 0xffff).toInt()
            val symbol: Int
            val offset: Int
            if (pos >= cutoff) {
                offset = ((b ushr 32) and 0xffff).toInt() + pos
                dist = dist xor ((b ushr 48) and 0xffff).toInt()
                symbol = (b and 0xff).toInt()
            } else {
                offset = pos
                symbol = i
            }
            var next = (state ushr 12) * dist + offset
            if (next < 0x10000) {
                next = (next shl 16) or br.readBits16().toLong()
            }
            reader.ansState = next
            return symbol
        }
        return prefix!![cluster].readSymbol(br)
    }

    fun singleSymbol(cluster: Int): Int =
        ansSingles?.get(cluster)?.toInt() ?: prefix!![cluster].singleSymbol

    /** LZ77 なしで、常に同じ値を返すクラスタならその値 */
    fun singleToken(cluster: Int): Int {
        if (lz77Enabled) return -1
        val s = singleSymbol(cluster)
        if (s < 0) return -1
        return if (s < configs[cluster].split) s else -1
    }

    val lzDistCluster: Int get() = clusters[clusters.size - 1].toInt() and 0xff

    companion object {
        fun parse(br: BitReader, numDist: Int): EntropyCode = parseInner(br, numDist, true)

        private fun parseInner(br: BitReader, numDistIn: Int, allowLz77: Boolean): EntropyCode {
            var numDist = numDistIn
            val lz77 = br.readBool()
            var minSymbol = 0
            var minLength = 0
            var lenConf: IntegerConfig? = null
            if (lz77) {
                if (!allowLz77) fail("LZ77 not allowed here")
                minSymbol = br.readU32(224, 0, 512, 0, 4096, 0, 8, 15)
                minLength = br.readU32(3, 0, 4, 0, 5, 2, 9, 8)
                lenConf = IntegerConfig.parse(br, 8)
                numDist++
            }
            val clusters = readClusters(br, numDist)
            var numClusters = 0
            for (c in clusters) {
                val ci = c.toInt() and 0xff
                if (ci + 1 > numClusters) numClusters = ci + 1
            }
            val usePrefix = br.readBool()
            val logAlphabetSize = if (usePrefix) 15 else br.readBits(2) + 5
            val configs = Array(numClusters) { IntegerConfig.parse(br, logAlphabetSize) }
            var ansBuckets: LongArray? = null
            var ansSingles: ShortArray? = null
            var prefix: Array<PrefixHistogram>? = null
            if (usePrefix) {
                val counts = IntArray(numClusters) {
                    if (br.readBool()) {
                        val n = br.readBits(4)
                        1 + (1 shl n) + br.readBits(n)
                    } else 1
                }
                for (c in counts) if (c > (1 shl 15)) fail("invalid prefix histogram")
                prefix = Array(numClusters) { PrefixHistogram.parse(br, counts[it]) }
            } else {
                val buckets = LongArray(numClusters shl logAlphabetSize)
                val singles = ShortArray(numClusters)
                val scratch = AnsParseScratch(1 shl logAlphabetSize)
                for (c in 0 until numClusters) {
                    singles[c] = AnsHistogram.parseInto(br, logAlphabetSize, buckets, c shl logAlphabetSize, scratch).toShort()
                }
                ansBuckets = buckets
                ansSingles = singles
            }
            return EntropyCode(clusters, configs, ansBuckets, logAlphabetSize, ansSingles, prefix, lz77, minSymbol, minLength, lenConf)
        }

        fun readClusters(br: BitReader, numDist: Int): ByteArray {
            if (numDist == 1) return ByteArray(1)
            val out: ByteArray
            if (br.readBool()) {
                val nbits = br.readBits(2)
                out = ByteArray(numDist) { br.readBits(nbits).toByte() }
            } else {
                val useMtf = br.readBool()
                val code = parseInner(br, 1, numDist > 2)
                val reader = code.newReader()
                reader.begin(br)
                out = ByteArray(numDist) {
                    val v = reader.read(br, 0)
                    if (v !in 0..255) fail("invalid cluster index")
                    v.toByte()
                }
                reader.finish()
                if (useMtf) {
                    val mtf = ByteArray(256) { it.toByte() }
                    for (i in out.indices) {
                        val idx = out[i].toInt() and 0xff
                        val v = mtf[idx]
                        out[i] = v
                        System.arraycopy(mtf, 0, mtf, 1, idx)
                        mtf[0] = v
                    }
                }
            }
            var numClusters = 0
            for (c in out) {
                val ci = c.toInt() and 0xff
                if (ci + 1 > numClusters) numClusters = ci + 1
            }
            val seen = BooleanArray(numClusters)
            for (c in out) seen[c.toInt() and 0xff] = true
            for (s in seen) if (!s) fail("cluster map has a hole")
            return out
        }
    }
}

private val SPECIAL_DISTANCES = intArrayOf(
    0, 1, 1, 0, 1, 1, -1, 1, 0, 2, 2, 0, 1, 2, -1, 2, 2, 1, -2, 1,
    2, 2, -2, 2, 0, 3, 3, 0, 1, 3, -1, 3, 3, 1, -3, 1, 2, 3, -2, 3,
    3, 2, -3, 2, 0, 4, 4, 0, 1, 4, -1, 4, 4, 1, -4, 1, 3, 3, -3, 3,
    2, 4, -2, 4, 4, 2, -4, 2, 0, 5, 3, 4, -3, 4, 4, 3, -4, 3, 5, 0,
    1, 5, -1, 5, 5, 1, -5, 1, 2, 5, -2, 5, 5, 2, -5, 2, 4, 4, -4, 4,
    3, 5, -3, 5, 5, 3, -5, 3, 0, 6, 6, 0, 1, 6, -1, 6, 6, 1, -6, 1,
    2, 6, -2, 6, 6, 2, -6, 2, 4, 5, -4, 5, 5, 4, -5, 4, 3, 6, -3, 6,
    6, 3, -6, 3, 0, 7, 7, 0, 1, 7, -1, 7, 5, 5, -5, 5, 7, 1, -7, 1,
    4, 6, -4, 6, 6, 4, -6, 4, 2, 7, -2, 7, 7, 2, -7, 2, 3, 7, -3, 7,
    7, 3, -7, 3, 5, 6, -5, 6, 6, 5, -6, 5, 8, 0, 4, 7, -4, 7, 7, 4,
    -7, 4, 8, 1, 8, 2, 6, 6, -6, 6, 8, 3, 5, 7, -5, 7, 7, 5, -7, 5,
    8, 4, 6, 7, -6, 7, 7, 6, -7, 6, 8, 5, 7, 7, -7, 7, 8, 6, 8, 7,
)

/** [EntropyCode] を使って値を読むときの状態 (ANS の状態・LZ77 の窓) */
internal class EntropyReader(val code: EntropyCode) {
    var ansState: Long = 0
    var ansInitial: Boolean = true
    private var window: IntArray = IntArray(0)
    private var numToCopy = 0
    private var copyPos = 0
    private var numDecoded = 0

    fun begin(br: BitReader) {
        if (code.usesAns) {
            ansState = br.readBits(32).toLong() and 0xffffffffL
            ansInitial = false
        }
    }

    fun finish() {
        if (code.usesAns && ansState != 0x130000L) fail("ANS stream checksum mismatch")
    }

    fun readUint(br: BitReader, config: IntegerConfig, token: Int): Int {
        if (token < config.split) return token
        val n = (config.splitExponent - (config.msbInToken + config.lsbInToken) +
            ((token - config.split) ushr (config.msbInToken + config.lsbInToken))) and 31
        val rest = br.readBits(n).toLong() and 0xffffffffL
        val lowBits = token and ((1 shl config.lsbInToken) - 1)
        var t = token ushr config.lsbInToken
        t = t and ((1 shl config.msbInToken) - 1)
        t = t or (1 shl config.msbInToken)
        val result = ((((t.toLong() shl n) or rest) shl config.lsbInToken) or lowBits.toLong())
        return result.toInt()
    }

    fun read(br: BitReader, ctx: Int): Int = readClustered(br, code.clusters[ctx].toInt() and 0xff, 0)

    /**
     * [EntropyCode.fastBuckets] がある (ANS で LZ77 なし) ときに、クラスタ [cluster] の値を読む。
     * [readClustered] と同じ結果を、表を 1 つの配列から引いて求める ([begin] で状態を読んでおくこと)
     */
    fun readAnsFast(br: BitReader, buckets: LongArray, logAlpha: Int, cfg: IntArray, cluster: Int): Int {
        // ANS のシンボル
        val state = ansState
        val idx = (state and 0xfff).toInt()
        val logBucket = 12 - logAlpha
        val i = idx ushr logBucket
        val pos = idx and ((1 shl logBucket) - 1)
        val b = buckets[(cluster shl logAlpha) + i]
        val cutoff = ((b ushr 8) and 0xff).toInt()
        var dist = ((b ushr 16) and 0xffff).toInt()
        val token: Int
        val offset: Int
        if (pos >= cutoff) {
            offset = ((b ushr 32) and 0xffff).toInt() + pos
            dist = dist xor ((b ushr 48) and 0xffff).toInt()
            token = (b and 0xff).toInt()
        } else {
            offset = pos
            token = i
        }
        var next = (state ushr 12) * dist + offset
        if (next < 0x10000) next = (next shl 16) or br.readBits16().toLong()
        ansState = next
        // hybrid-uint (readUint と同じ)
        val c4 = cluster * 4
        if (token < cfg[c4]) return token
        val msb = cfg[c4 + 2]
        val lsb = cfg[c4 + 3]
        val n = (cfg[c4 + 1] - (msb + lsb) + ((token - cfg[c4]) ushr (msb + lsb))) and 31
        val rest = br.readBits(n).toLong() and 0xffffffffL
        val lowBits = token and ((1 shl lsb) - 1)
        val t = ((token ushr lsb) and ((1 shl msb) - 1)) or (1 shl msb)
        return (((t.toLong() shl n) or rest) shl lsb or lowBits.toLong()).toInt()
    }

    fun readClustered(br: BitReader, cluster: Int, distMultiplier: Int): Int {
        if (!code.lz77Enabled) {
            val token = code.readSymbol(br, cluster, this)
            return readUint(br, code.configs[cluster], token)
        }
        return readLz77(br, cluster, distMultiplier)
    }

    private fun pushWindow(v: Int) {
        val off = numDecoded and 0xfffff
        if (off >= window.size) {
            val newSize = if (window.isEmpty()) 4096 else minOf(window.size * 2, 1 shl 20)
            window = window.copyOf(newSize)
        }
        window[off] = v
        numDecoded++
    }

    private fun readLz77(br: BitReader, cluster: Int, distMultiplier: Int): Int {
        val r: Int
        if (numToCopy > 0) {
            r = window[copyPos and 0xfffff]
            copyPos++
            numToCopy--
        } else {
            val token = code.readSymbol(br, cluster, this)
            if (token >= code.lz77MinSymbol) {
                if (numDecoded == 0) fail("LZ77 repeat without preceding symbols")
                val lzCluster = code.lzDistCluster
                val n = readUint(br, code.lz77LenConfig!!, token - code.lz77MinSymbol)
                val copyLen = n.toLong() + code.lz77MinLength
                if (copyLen > Int.MAX_VALUE) fail("invalid LZ77 length")
                numToCopy = copyLen.toInt()
                val dtoken = code.readSymbol(br, lzCluster, this)
                var distance = readUint(br, code.configs[lzCluster], dtoken)
                distance = if (distMultiplier == 0) {
                    distance
                } else if (distance in 0 until 120) {
                    val offset = SPECIAL_DISTANCES[distance * 2]
                    val dist = SPECIAL_DISTANCES[distance * 2 + 1]
                    maxOf(offset + distMultiplier * dist - 1, 0)
                } else {
                    distance - 120
                }
                val d = minOf(minOf(distance.toLong() and 0xffffffffL, (1L shl 20) - 1) + 1, numDecoded.toLong()).toInt()
                copyPos = numDecoded - d
                r = window[copyPos and 0xfffff]
                copyPos++
                numToCopy--
            } else {
                r = readUint(br, code.configs[cluster], token)
            }
        }
        pushWindow(r)
        return r
    }
}

/** レーマー符号で符号化された順列を読む */
internal fun readPermutation(br: BitReader, reader: EntropyReader, size: Int, skip: Int): IntArray {
    fun ctx(x: Int) = minOf(addLog2Ceil(x), 7)
    val end = reader.read(br, ctx(size))
    if (end > size - skip) fail("invalid permutation")
    val lehmer = IntArray(end)
    var prev = 0
    for (i in 0 until end) {
        val v = reader.read(br, ctx(prev))
        if (v >= size - skip - i) fail("invalid permutation")
        lehmer[i] = v
        prev = v
    }
    val temp = ArrayList<Int>(size - skip)
    for (i in skip until size) temp.add(i)
    val out = IntArray(size)
    for (i in 0 until skip) out[i] = i
    var idx = skip
    for (l in lehmer) out[idx++] = temp.removeAt(l)
    for (t in temp) out[idx++] = t
    return out
}
