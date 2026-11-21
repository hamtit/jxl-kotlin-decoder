package io.github.hamtit.jxl

import kotlin.test.Test
import kotlin.test.assertEquals
import kotlin.test.assertNotNull

/** 1 種類のプロパティだけで分岐する木の表 (MaTree.singlePropertyLut) が、木をたどった結果と一致するか */
class MaTreeLutTest {
    /** 木の形だけを使うので、符号は何でもよい */
    private val code = EntropyCode.parse(BitReader(ByteArray(64)), 1)

    /** 分割値の並びから、右 (v <= split) に伸びる 1 本の深い木を作る。葉 i の predictor に i を入れて区別する */
    private fun chain(property: Int, splits: IntArray): MaTree {
        val n = splits.size * 2 + 1
        val prop = IntArray(n) { -1 }
        val split = IntArray(n)
        val left = IntArray(n)
        val right = IntArray(n)
        val pred = IntArray(n)
        // ノード 2k が分岐、2k + 1 がその左の葉、最後が右端の葉
        for (k in splits.indices) {
            prop[2 * k] = property
            split[2 * k] = splits[k]
            left[2 * k] = 2 * k + 1
            right[2 * k] = 2 * k + 2
            pred[2 * k + 1] = k
        }
        pred[n - 1] = splits.size
        return MaTree(prop, split, left, right, pred, IntArray(n), IntArray(n) { 1 }, IntArray(n), code)
    }

    private fun walk(t: MaTree, v: Int): Int {
        var node = 0
        while (t.property[node] >= 0) node = if (v > t.splitValue[node]) t.left[node] else t.right[node]
        return node
    }

    private fun lookup(t: MaTree, v: Int): Int {
        val lut = assertNotNull(t.singlePropertyLut)
        val c = if (v < lut.min) lut.min else if (v > lut.max) lut.max else v
        return lut.leaf[c - lut.min]
    }

    private fun check(t: MaTree, values: List<Int>) {
        for (v in values) assertEquals(walk(t, v), lookup(t, v), "v = $v")
    }

    @Test
    fun splitAtIntMax() {
        val m = Int.MAX_VALUE
        val t = chain(9, intArrayOf(m, m - 1, m - 5, m - 3))
        check(t, listOf(Int.MIN_VALUE, 0, m - 6, m - 5, m - 4, m - 3, m - 2, m - 1, m))
    }

    @Test
    fun splitAtIntMin() {
        val m = Int.MIN_VALUE
        val t = chain(15, intArrayOf(m + 2, m, m + 1))
        check(t, listOf(m, m + 1, m + 2, m + 3, 0, Int.MAX_VALUE))
    }

    @Test
    fun deepTree() {
        // 分割値が狭い範囲に並ぶ深い木 (表を 1 回の走査で埋める)
        val splits = IntArray(8000) { (it * 7919) % 8000 - 4000 }
        val t = chain(9, splits)
        check(t, (-4100..4100).toList())
    }
}
