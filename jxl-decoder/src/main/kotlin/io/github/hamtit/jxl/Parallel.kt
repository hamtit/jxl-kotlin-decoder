package io.github.hamtit.jxl

import java.util.concurrent.ExecutorService
import java.util.concurrent.Future

/** 0 until n を並列に処理する。executor が null なら直列 */
internal fun parallelFor(executor: ExecutorService?, n: Int, body: (Int) -> Unit) {
    if (executor == null || n <= 1) {
        for (i in 0 until n) body(i)
        return
    }
    val futures = ArrayList<Future<*>>(n)
    for (i in 0 until n) futures.add(executor.submit { body(i) })
    var error: Throwable? = null
    for (f in futures) {
        try {
            f.get()
        } catch (e: java.util.concurrent.ExecutionException) {
            if (error == null) error = e.cause ?: e
        }
    }
    error?.let { throw it }
}

/** 行の範囲を 32 行ずつの帯に分けて (並列に) 処理する */
internal fun parallelRows(executor: ExecutorService?, height: Int, body: (Int, Int) -> Unit) {
    val chunk = 32
    val n = (height + chunk - 1) / chunk
    if (executor == null || n <= 1) {
        var y0 = 0
        while (y0 < height) {
            body(y0, minOf(height, y0 + chunk))
            y0 += chunk
        }
        return
    }
    parallelFor(executor, n) { i ->
        val y0 = i * chunk
        body(y0, minOf(height, y0 + chunk))
    }
}
