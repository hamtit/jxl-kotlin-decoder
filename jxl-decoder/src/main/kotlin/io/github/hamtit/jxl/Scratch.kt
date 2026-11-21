package io.github.hamtit.jxl

import java.util.concurrent.ConcurrentLinkedQueue

/**
 * 作業用の配列の置き場。1 つのタスクが使っている間は他から使わない。
 * 配列は使い回すので中身は前の値が残っている (必要なら呼び出し側で初期化する)。
 * 要求より大きい配列を返すことがあるので、大きさは配列の size ではなく要求した値を使う。
 */
internal class Scratch {
    private val floats = arrayOfNulls<FloatArray>(SLOTS)
    private val ints = arrayOfNulls<IntArray>(SLOTS)

    fun floats(slot: Int, size: Int): FloatArray {
        val a = floats[slot]
        if (a != null && a.size >= size) return a
        return FloatArray(size).also { floats[slot] = it }
    }

    fun ints(slot: Int, size: Int): IntArray {
        val a = ints[slot]
        if (a != null && a.size >= size) return a
        return IntArray(size).also { ints[slot] = it }
    }

    companion object {
        const val SLOTS = 24
    }
}

/** フレームの復号中に使い回す [Scratch] の集まり (同時に動くタスクの数だけ作られる) */
internal class ScratchPool {
    val free = ConcurrentLinkedQueue<Scratch>()

    inline fun <R> use(block: (Scratch) -> R): R {
        val s = free.poll() ?: Scratch()
        try {
            return block(s)
        } finally {
            free.offer(s)
        }
    }
}
