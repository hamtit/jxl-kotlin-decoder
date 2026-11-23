# jxl-kotlin-decoder

純 Kotlin で書いた JPEG XL デコーダーです。ネイティブコード (JNI) を使わないので、Android アプリに入れても ABI ごとの `.so` が増えません。JVM 17 のバイトコードなので、Android 12 (API 31) 以降と JVM で動きます。外部の依存ライブラリはありません (Kotlin 標準ライブラリのみ)。

JPEG XL の仕様と、参照実装の [libjxl](https://github.com/libjxl/libjxl) および [jxl-oxide](https://github.com/tirr-c/jxl-oxide) を元に実装しています。

## 使い方

```kotlin
import io.github.hamtit.jxl.JxlDecoder

// 静止画 (最初のフレーム)
val frame = JxlDecoder.decode(bytes)
val bitmap = Bitmap.createBitmap(frame.pixels, frame.width, frame.height, Bitmap.Config.ARGB_8888)

// アニメーション: フレームを先頭から順に取り出す (全フレームを一度に展開しない)
val decoder = JxlDecoder(bytes)
println("${decoder.info.width}x${decoder.info.height} animated=${decoder.info.isAnimated}")
while (true) {
    val f = decoder.nextFrame() ?: break   // f.pixels は ARGB (アルファは非乗算)、f.durationMillis は表示時間
}
decoder.rewind()                           // 最初から読み直す (ループ再生)
```

- `JxlDecoder.isJxl(head)` で先頭のバイト列が JPEG XL (コードストリーム `FF 0A` またはコンテナ) か判定できます。
- `JxlDecodeOptions(threads = 1)` で呼び出したスレッドだけで復号します。既定では CPU 数のスレッドを持つ共有のスレッドプールでグループを並列に処理します。独自の `ExecutorService` も渡せますが、そのプールのスレッドからデコーダーを呼ばないでください (待ち合わせで止まります)。
- 出力は 8 ビット sRGB の ARGB です。向き (orientation) は既定で反映します。
- 壊れたデータは `JxlException` になります。`maxPixels` (既定 2^28) より大きい画像は復号しません。
- `JxlDecoder` はスレッドセーフではありません。

## 依存関係への追加

[リリース](https://github.com/hamtit/jxl-kotlin-decoder/releases) に jar を添付しています。Gradle なら ivy リポジトリとしてリリースから直接取れます (Kotlin 標準ライブラリ以外の依存はありません)。

```kotlin
// settings.gradle.kts
dependencyResolutionManagement {
    repositories {
        ivy("https://github.com/hamtit/jxl-kotlin-decoder/releases/download") {
            patternLayout { artifact("v[revision]/[artifact]-[revision](-[classifier]).[ext]") }
            metadataSources { artifact() }
            content { includeModule("io.github.hamtit", "jxl-decoder") }
        }
    }
}

// build.gradle.kts
dependencies { implementation("io.github.hamtit:jxl-decoder:<バージョン>") }
```

ソースから使う場合は、このリポジトリを clone して Gradle の composite build (`includeBuild`) で `:jxl-decoder` を参照してください。

## 対応している機能

- コードストリームとコンテナ (`jxlc`・`jxlp` ボックス)
- VarDCT: 全ての変換 (DCT2〜DCT256、非正方形、Hornuss、AFV)、適応量子化、独自の量子化行列、係数の順序、CfL、適応 LF 平滑化、JPEG からの再圧縮 (YCbCr・クロマサブサンプリング)
- Modular: MA 木、全ての予測器 (自己補正予測器を含む)、RCT・パレット (デルタパレットを含む)・Squeeze、LZ77、16 ビット超・浮動小数点のサンプル
- Gaborish・エッジ保存フィルタ (EPF)、アップサンプリング、パッチ、スプライン、ノイズ
- プログレッシブ (複数パス、LF フレーム、TOC の並べ替え)
- アニメーション、フレームのブレンド (全モード)、参照フレーム、切り出したフレーム
- アルファ (乗算済みを含む)、グレースケール、向き
- 色: XYB は sRGB に変換、列挙型で指定された色空間 (Display P3・BT.2020・線形・ガンマ・PQ・HLG など) は sRGB に変換

### 対応していないもの

- ICC プロファイルによる色変換 (ICC プロファイル付きの画像は、プロファイルを読み飛ばして変換せずに出力します)。CMYK もこのため正しい色になりません
- スポットカラーの合成、JPEG のバイト列の復元 (`jbrd`)、途中までしか受信していないデータの表示、HDR のトーンマッピング (PQ・HLG は SDR の白を 1 として切り詰めます)

## 正確さ

libjxl の `djxl` で 8 ビットに復号した結果と比べ、全画素の差が非可逆で ±1 以内 (丸めの違い)、可逆で完全一致です。`jxl-decoder/src/test` に cjxl などで作った画像と djxl の出力を入れてあり、`./gradlew test` で確かめます。テスト画像は [Kodak の画像セット](https://r0k.us/graphics/kodak/) の一部を切り抜いたものと、スクリプトで描いた画像 (文字・グラデーション・アニメーション) です。

[libjxl の conformance テスト](https://github.com/libjxl/conformance) 40 件では、次の 7 件以外が djxl と ±1 以内で一致しました。

| テスト | 理由 |
|---|---|
| bench_oriented_brg (×2)・cafe (×2)・cmyk_layers | ICC プロファイルによる色変換が必要 |
| spot | ICC プロファイル + スポットカラーの合成 |
| alpha_premultiplied | 差は出力形式だけ (djxl は乗算済みのまま書き出し、このライブラリは非乗算に戻す) |

## 性能

Linux x86-64 (4 コア) で、JVM (OpenJDK 21) 上のこのライブラリと、同じマシンの libjxl 0.14 (`djxl`、main の開発版 12e4bf4 をビルド、AVX2) を比べた結果です。画像は libjxl 0.14 の cjxl と、別の軽量な JPEG XL エンコーダー (品質 65、Gaborish なし・EPF 2 回) で作りました。時間はどちらもファイルの読み込みを除いた値で、このライブラリは 30〜60 回の中央値、libjxl は 60 回の平均です。測ったマシンは負荷の揺れが大きく、このライブラリの値は 10% 程度ぶれます。

| 画像 | libjxl 1 スレッド | Kotlin 1 スレッド | libjxl 4 スレッド | Kotlin 4 スレッド |
|---|---:|---:|---:|---:|
| 768x512 非可逆 (d1) | 19 ms | 30 ms | 9 ms | 18 ms |
| 1920x1280 非可逆 (d1) | 81 ms | 175 ms | 33 ms | 93 ms |
| 1920x1280 非可逆 (軽量エンコーダー、品質 65) | 80 ms | 190 ms | 40 ms | 98 ms |
| 768x512 可逆 (e3) | 56 ms | 80 ms | 20 ms | 27 ms |
| 1920x1280 可逆 (e2) | 179 ms | 230 ms | 49 ms | 72 ms |
| 256x170 6 フレームのアニメーション (d1) | 4 ms | 8 ms | 7 ms | 10 ms |

SIMD を使う libjxl に比べ、非可逆で約 2〜2.5 倍、可逆で約 1.3〜1.5 倍の時間がかかります。0.1.3 と比べると (同じマシンで交互に測定)、非可逆の画像は 1 スレッドで 15〜23%、4 スレッドで 12〜16% 速くなっています。Android の ART では JVM (HotSpot) より遅くなる傾向があります (HotSpot の自動ベクトル化を止めて測っても、0.1.3 からの短縮はほぼ同じでした)。

メモリは、VarDCT (非可逆) で置き換え (REPLACE) で重ねるフレーム (普通の静止画と、多くのアニメーション) ならグループ (256x256 画素) ごとに復号してすぐに描画し、直接 ARGB に書き出します。係数はスレッドごとのグループの大きさの配列 (0.8 MB 程度) に直接置いてその場で逆変換し、復元フィルタ (Gaborish・EPF) は段階ごとに数行だけを持って上から 1 行ずつ計算します。各グループは自分の領域を左上に 8 画素ずらした範囲を描くので、要るのは上のグループの下端 16 行 (画像の幅ぶん) と、左右のグループとの継ぎ目の 32 列だけで、画像の幅に比例するグループ 1 行分の配列を持ちません。1920x1280 の非可逆画像で Java ヒープ 17 MB 程度 (出力の ARGB 配列 10 MB を含む。1 スレッドで復号できる最小のヒープ。グループの行ごとに描画していた前の版は 23 MB 程度)、4032x3024 で 59 MB 程度 (出力の 49 MB を含む。前の版は 69 MB 程度) です。スレッドを増やすと、同時に復号するグループの数だけ作業領域が増えます (1920x1280 の 4 スレッドで 23 MB 程度)。アルファ付きの非可逆画像 (アルファが Squeeze などで変換されている場合も) も同じように描画し、アルファは最後に書き込みます。アニメーションで後のフレームが置き換えで重ねるだけの参照フレームも、浮動小数点ではなく ARGB (1 画素 4 バイト) で持ちます。色変換やフィルタが不要な可逆画像は、整数サンプルから直接 ARGB を作ります。画像全体にかかる変換が RCT だけ (cjxl の可逆の多くの設定) なら、グループごとに小さな配列へ復号してそのまま ARGB に書くので、画像の大きさのチャンネルの配列を作りません (1920x1280 の可逆画像で Java ヒープ 16 MB 程度、出力の 10 MB を含む)。それ以外の可逆画像は、最初の色チャンネルの配列を出力に使い、別の ARGB 配列を確保しません。前のフレームとアルファブレンドなどで重ねるアニメーションや、直接描画できない画像は、フレーム全体を浮動小数点で持ちます。復号の作業領域 (EPF の行・逆変換の作業領域・帯など) は合わせて 24 MB までをフレームや画像をまたいで使い回します (メモリが足りなくなれば GC が回収します)。

## 開発

```sh
./gradlew test                           # テスト (libjxl の出力との比較を含む)
./gradlew :cli:installDist
cli/build/install/cli/bin/cli decode in.jxl out.png
cli/build/install/cli/bin/cli bench in.jxl --iterations 20 --threads 1
cli/build/install/cli/bin/cli compare in.jxl ref.png   # PSNR と最大誤差
```

環境変数 `JXL_TIMING=1` を付けると、フレームごとの処理時間の内訳を標準エラーに出します。

### バージョンとリリース

バージョンは `gradle.properties` の `version` だけで決まります。上げて `main` にマージすると、CI (`.github/workflows/ci.yml`) が `v<version>` のタグと GitHub のリリース (jar・ソースの jar・SHA256SUMS) を作ります。同じバージョンのタグを手で push しても同じようにリリースされます (タグと `gradle.properties` が食い違うと失敗します)。

CI は PR と `main` への push のたびに、警告をエラーにしたコンパイル、Linux arm64 でのビルドとテスト (jar は CPU によらないので、Android 端末と同じ arm64 で 1 回だけ)、jar の中身 (ライセンス・Java 17 のクラス) と jar での復号の確認を行います。

## 謝辞

- [libjxl](https://github.com/libjxl/libjxl) (BSD-3-Clause, Copyright (c) the JPEG XL Project Authors): 仕様の参照、定数・表 (量子化行列、AFV の基底、アップサンプリングの重み、sRGB の近似式など)
- [jxl-oxide](https://github.com/tirr-c/jxl-oxide) (MIT OR Apache-2.0, Copyright (c) Wonwoo Choi): デコード処理の構成と各処理の実装の参考

## ライセンス

BSD-3-Clause ([LICENSE](LICENSE))。jxl-oxide を元にした部分の MIT ライセンスの表示は [LICENSE-jxl-oxide](LICENSE-jxl-oxide) にあります。どちらも jar の `META-INF` に入っています。
