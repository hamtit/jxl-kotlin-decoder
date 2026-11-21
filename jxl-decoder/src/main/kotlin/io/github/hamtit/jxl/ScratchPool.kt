package io.github.hamtit.jxl

import java.lang.ref.SoftReference

/**
 * 復号の作業領域 (EPF の帯・グループの係数・逐次描画の行など) を、フレームや画像をまたいで使い回す置き場。
 * アニメーションではフレームごとに数 MB の作業領域を確保し直していたので、それをなくす。
 *
 * 置いておく量は全体で [LIMIT_BYTES] までで、超える分は古いものから捨てる。
 * また SoftReference で持つので、メモリが足りなくなれば GC が回収する。
 * 複数のデコーダー・スレッドから同時に使ってよい。
 */
internal object ScratchPool {
    private const val LIMIT_BYTES = 24L shl 20

    private class Entry(val ref: SoftReference<Any>, val bytes: Long)

    // 古いものが先頭。取り出しは新しい方から探す
    private val entries = ArrayList<Entry>()
    private var pooledBytes = 0L

    /** 条件に合う作業領域を取り出す。無ければ null */
    @Suppress("UNCHECKED_CAST")
    fun <T : Any> take(cls: Class<T>, accept: (T) -> Boolean): T? = synchronized(this) {
        var i = entries.size - 1
        while (i >= 0) {
            val e = entries[i]
            val v = e.ref.get()
            if (v == null) {
                removeAt(i)
            } else if (cls.isInstance(v) && accept(v as T)) {
                removeAt(i)
                return v
            }
            i--
        }
        null
    }

    /** 作業領域を返す (bytes はおおよその大きさ)。上限を超える大きさなら置かない */
    fun give(v: Any, bytes: Long) {
        if (bytes > LIMIT_BYTES) return
        synchronized(this) {
            while (entries.isNotEmpty() && pooledBytes + bytes > LIMIT_BYTES) removeAt(0)
            entries.add(Entry(SoftReference(v), bytes))
            pooledBytes += bytes
        }
    }

    private fun removeAt(i: Int) {
        pooledBytes -= entries.removeAt(i).bytes
    }

    /** 長さ n 以上の FloatArray (中身は不定) */
    fun floats(n: Int): FloatArray {
        // 大きすぎるものを小さい用途に使うと、後の大きい用途で確保し直すことになるので 2 倍までにする
        return take(FloatArray::class.java) { it.size >= n && it.size / 2 <= n } ?: FloatArray(n)
    }

    fun give(a: FloatArray) = give(a, a.size * 4L)
}
