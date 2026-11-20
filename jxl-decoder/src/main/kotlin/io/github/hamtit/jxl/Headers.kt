package io.github.hamtit.jxl

internal class BitDepth(val bitsPerSample: Int, val expBits: Int) {
    val isFloat: Boolean get() = expBits != 0

    /** 整数サンプルを 0..1 の浮動小数点に直す (浮動小数点サンプルはビット列として解釈する) */
    fun scale(): Float = 1.0f / ((1L shl bitsPerSample) - 1).toFloat()

    fun parseFloatSample(sample: Int): Float {
        val mantissaBits = bitsPerSample - expBits - 1
        val mantissaMask = (1 shl mantissaBits) - 1
        val expMask = ((1 shl (bitsPerSample - 1)) - 1) xor mantissaMask
        val signed = sample and (1 shl (bitsPerSample - 1)) != 0
        var mantissa = sample and mantissaMask
        var exp = (sample and expMask) ushr mantissaBits
        if (exp == 0 && mantissa == 0) return if (signed) -0.0f else 0.0f
        exp -= (1 shl (expBits - 1)) - 1
        mantissa = if (mantissaBits < 23) mantissa shl (23 - mantissaBits) else mantissa ushr (mantissaBits - 23)
        val bits = ((if (signed) 1 else 0) shl 31) or ((exp + 127) shl 23) or mantissa
        return Float.fromBits(bits)
    }

    companion object {
        val DEFAULT = BitDepth(8, 0)

        fun parse(br: BitReader): BitDepth {
            return if (br.readBool()) {
                val bps = br.readU32(32, 0, 16, 0, 24, 0, 1, 6)
                val exp = br.readBits(4) + 1
                if (exp !in 2..8) fail("invalid exp_bits")
                val mant = bps - exp - 1
                if (mant !in 2..23) fail("invalid mantissa bits")
                BitDepth(bps, exp)
            } else {
                val bps = br.readU32(8, 0, 10, 0, 12, 0, 1, 6)
                if (bps > 31) fail("invalid bits_per_sample")
                BitDepth(bps, 0)
            }
        }
    }
}

internal object ExtraChannelType {
    const val ALPHA = 0
    const val DEPTH = 1
    const val SPOT_COLOUR = 2
    const val SELECTION_MASK = 3
    const val BLACK = 4
    const val CFA = 5
    const val THERMAL = 6
    const val NON_OPTIONAL = 15
    const val OPTIONAL = 16
}

internal class ExtraChannelInfo(
    val type: Int,
    val bitDepth: BitDepth,
    val dimShift: Int,
    val alphaAssociated: Boolean,
    val spotColour: FloatArray?,
) {
    val isAlpha: Boolean get() = type == ExtraChannelType.ALPHA

    companion object {
        fun parse(br: BitReader): ExtraChannelInfo {
            if (br.readBool()) return ExtraChannelInfo(ExtraChannelType.ALPHA, BitDepth.DEFAULT, 0, false, null)
            val type = br.readEnum()
            if (type !in 0..6 && type != 15 && type != 16) fail("invalid extra channel type")
            val bitDepth = BitDepth.parse(br)
            val dimShift = br.readU32(0, 0, 3, 0, 4, 0, 1, 3)
            val nameLen = br.readU32(0, 0, 0, 4, 16, 5, 48, 10)
            for (i in 0 until nameLen) br.readBits(8)
            var assoc = false
            var spot: FloatArray? = null
            when (type) {
                ExtraChannelType.ALPHA -> assoc = br.readBool()
                ExtraChannelType.SPOT_COLOUR -> spot = FloatArray(4) { br.readF16() }
                ExtraChannelType.CFA -> br.readU32(1, 0, 0, 2, 3, 4, 19, 8)
            }
            return ExtraChannelInfo(type, bitDepth, dimShift, assoc, spot)
        }
    }
}

internal object ColourSpace {
    const val RGB = 0
    const val GREY = 1
    const val XYB = 2
    const val UNKNOWN = 3
}

internal object TransferFunction {
    const val BT709 = 1
    const val UNKNOWN = 2
    const val LINEAR = 8
    const val SRGB = 13
    const val PQ = 16
    const val DCI = 17
    const val HLG = 18
    const val GAMMA = -1
}

internal class ColourEncoding(
    val wantIcc: Boolean,
    val colourSpace: Int,
    /** CIE xy: white(2), red(2), green(2), blue(2) */
    val whitePoint: FloatArray,
    val primaries: FloatArray,
    val transferFunction: Int,
    /** transferFunction == GAMMA のときの gamma (例 1/2.2 = 0.4545) */
    val gamma: Float,
    val renderingIntent: Int,
) {
    val isSrgb: Boolean
        get() = !wantIcc && colourSpace == ColourSpace.RGB && transferFunction == TransferFunction.SRGB &&
            isD65 && isSrgbPrimaries

    val isD65: Boolean get() = whitePoint[0] == D65[0] && whitePoint[1] == D65[1]
    val isSrgbPrimaries: Boolean get() = primaries.contentEquals(SRGB_PRIMARIES)

    companion object {
        val D65 = floatArrayOf(0.3127f, 0.3290f)
        val SRGB_PRIMARIES = floatArrayOf(0.639998686f, 0.330010138f, 0.300003784f, 0.600003357f, 0.150002046f, 0.059997204f)
        private val BT2100 = floatArrayOf(0.708f, 0.292f, 0.170f, 0.797f, 0.131f, 0.046f)
        private val P3 = floatArrayOf(0.680f, 0.320f, 0.265f, 0.690f, 0.150f, 0.060f)
        private val DCI_WHITE = floatArrayOf(0.314f, 0.351f)
        private val E_WHITE = floatArrayOf(1f / 3f, 1f / 3f)

        val SRGB = ColourEncoding(false, ColourSpace.RGB, D65, SRGB_PRIMARIES, TransferFunction.SRGB, 0f, 1)

        private fun customXy(br: BitReader): FloatArray {
            val x = unpackSigned(br.readU32(0, 19, 524288, 19, 1048576, 20, 2097152, 21))
            val y = unpackSigned(br.readU32(0, 19, 524288, 19, 1048576, 20, 2097152, 21))
            return floatArrayOf(x / 1e6f, y / 1e6f)
        }

        fun parse(br: BitReader): ColourEncoding {
            if (br.readBool()) return SRGB
            val wantIcc = br.readBool()
            val cs = br.readEnum()
            if (cs > 3) fail("invalid colour space")
            if (wantIcc) return ColourEncoding(true, cs, D65, SRGB_PRIMARIES, TransferFunction.SRGB, 0f, 0)
            var white = D65
            if (cs != ColourSpace.XYB) {
                white = when (br.readEnum()) {
                    1 -> D65
                    2 -> customXy(br)
                    10 -> E_WHITE
                    11 -> DCI_WHITE
                    else -> fail("invalid white point")
                }
            }
            var prim = SRGB_PRIMARIES
            if (cs != ColourSpace.XYB && cs != ColourSpace.GREY) {
                prim = when (br.readEnum()) {
                    1 -> SRGB_PRIMARIES
                    2 -> customXy(br) + customXy(br) + customXy(br)
                    9 -> BT2100
                    11 -> P3
                    else -> fail("invalid primaries")
                }
            }
            var tf: Int
            var gamma = 0f
            if (br.readBool()) {
                tf = TransferFunction.GAMMA
                gamma = br.readBits(24) / 1e7f
            } else {
                tf = br.readEnum()
                if (tf !in intArrayOf(1, 2, 8, 13, 16, 17, 18)) fail("invalid transfer function")
            }
            val intent = br.readEnum()
            if (intent > 3) fail("invalid rendering intent")
            return ColourEncoding(false, cs, white, prim, tf, gamma, intent)
        }
    }
}

internal class ToneMapping(val intensityTarget: Float, val minNits: Float, val relativeToMaxDisplay: Boolean, val linearBelow: Float) {
    companion object {
        val DEFAULT = ToneMapping(255f, 0f, false, 0f)
        fun parse(br: BitReader): ToneMapping {
            if (br.readBool()) return DEFAULT
            return ToneMapping(br.readF16(), br.readF16(), br.readBool(), br.readF16())
        }
    }
}

internal class AnimationHeader(val tpsNumerator: Int, val tpsDenominator: Int, val numLoops: Int, val haveTimecodes: Boolean)

internal class ImageHeader(
    val width: Int,
    val height: Int,
    val orientation: Int,
    val animation: AnimationHeader?,
    val bitDepth: BitDepth,
    val modular16BitBuffers: Boolean,
    val extraChannels: List<ExtraChannelInfo>,
    val xybEncoded: Boolean,
    val colourEncoding: ColourEncoding,
    val toneMapping: ToneMapping,
    val opsinInverseMatrix: FloatArray,
    val opsinBias: FloatArray,
    val quantBias: FloatArray,
    val quantBiasNumerator: Float,
    val up2Weights: FloatArray,
    val up4Weights: FloatArray,
    val up8Weights: FloatArray,
    val previewWidth: Int,
    val previewHeight: Int,
) {
    val alphaChannel: Int = extraChannels.indexOfFirst { it.isAlpha }
    val isGrey: Boolean get() = colourEncoding.colourSpace == ColourSpace.GREY

    companion object {
        private fun sizeU32(br: BitReader) = br.readU32(1, 9, 1, 13, 1, 18, 1, 30)

        private fun widthFromRatio(ratio: Int, height: Int): Int {
            val h = height.toLong()
            return when (ratio) {
                1 -> h
                2 -> h * 12 / 10
                3 -> h * 4 / 3
                4 -> h * 3 / 2
                5 -> h * 16 / 9
                6 -> h * 5 / 4
                7 -> h * 2
                else -> 0L
            }.toInt()
        }

        /** 署名 (0xFF0A) の直後から読む */
        fun parse(br: BitReader): ImageHeader {
            // SizeHeader
            val div8 = br.readBool()
            val height = if (div8) (br.readBits(5) + 1) * 8 else sizeU32(br)
            val ratio = br.readBits(3)
            val width = if (ratio == 0) {
                if (div8) (br.readBits(5) + 1) * 8 else sizeU32(br)
            } else widthFromRatio(ratio, height)

            // ImageMetadata
            val allDefault = br.readBool()
            var orientation = 1
            var animation: AnimationHeader? = null
            var bitDepth = BitDepth.DEFAULT
            var modular16 = true
            var extra = emptyList<ExtraChannelInfo>()
            var xyb = true
            var colour = ColourEncoding.SRGB
            var tone = ToneMapping.DEFAULT
            var previewW = 0
            var previewH = 0
            if (!allDefault) {
                val extraFields = br.readBool()
                if (extraFields) {
                    orientation = br.readBits(3) + 1
                    if (br.readBool()) {
                        // intrinsic size
                        val d8 = br.readBool()
                        val h = if (d8) (br.readBits(5) + 1) * 8 else sizeU32(br)
                        val r = br.readBits(3)
                        if (r == 0) { if (d8) br.readBits(5) else sizeU32(br) }
                        h.hashCode()
                    }
                    if (br.readBool()) {
                        // preview
                        val d8 = br.readBool()
                        previewH = if (d8) br.readU32(16, 0, 32, 0, 1, 5, 33, 9) * 8 else br.readU32(1, 6, 65, 8, 321, 10, 1345, 12)
                        val r = br.readBits(3)
                        previewW = if (r == 0) {
                            if (d8) br.readU32(16, 0, 32, 0, 1, 5, 33, 9) * 8 else br.readU32(1, 6, 65, 8, 321, 10, 1345, 12)
                        } else widthFromRatio(r, previewH)
                    }
                    if (br.readBool()) {
                        val num = br.readU32(100, 0, 1000, 0, 1, 10, 1, 30)
                        val den = br.readU32(1, 0, 1001, 0, 1, 8, 1, 10)
                        val loops = br.readU32(0, 0, 0, 3, 0, 16, 0, 32)
                        val tc = br.readBool()
                        animation = AnimationHeader(num, den, loops, tc)
                    }
                }
                bitDepth = BitDepth.parse(br)
                modular16 = br.readBool()
                val numExtra = br.readU32(0, 0, 1, 0, 2, 4, 1, 12)
                if (numExtra > 256) fail("too many extra channels")
                extra = List(numExtra) { ExtraChannelInfo.parse(br) }
                xyb = br.readBool()
                colour = ColourEncoding.parse(br)
                if (extraFields) tone = ToneMapping.parse(br)
                br.extensions()
            }
            var invMat = DEFAULT_INV_MAT
            var opsinBias = DEFAULT_OPSIN_BIAS
            var quantBias = DEFAULT_QUANT_BIAS
            var quantBiasNum = 0.145f
            var up2 = UpsamplingWeights.UP2
            var up4 = UpsamplingWeights.UP4
            var up8 = UpsamplingWeights.UP8
            val defaultM = br.readBool()
            if (!defaultM) {
                if (xyb) {
                    if (!br.readBool()) {
                        invMat = FloatArray(9) { br.readF16() }
                        opsinBias = FloatArray(3) { br.readF16() }
                        quantBias = FloatArray(3) { br.readF16() }
                        quantBiasNum = br.readF16()
                    }
                }
                val cwMask = br.readBits(3)
                if (cwMask and 1 != 0) up2 = FloatArray(15) { br.readF16() }
                if (cwMask and 2 != 0) up4 = FloatArray(55) { br.readF16() }
                if (cwMask and 4 != 0) up8 = FloatArray(210) { br.readF16() }
            }
            if (width <= 0 || height <= 0) fail("invalid image size")
            if (tone.intensityTarget <= 0f) fail("invalid intensity target")
            return ImageHeader(
                width, height, orientation, animation, bitDepth, modular16, extra, xyb, colour, tone,
                invMat, opsinBias, quantBias, quantBiasNum, up2, up4, up8, previewW, previewH,
            )
        }

        private val DEFAULT_INV_MAT = floatArrayOf(
            11.031566901960783f, -9.866943921568629f, -0.16462299647058826f,
            -3.254147380392157f, 4.418770392156863f, -0.16462299647058826f,
            -3.6588512862745097f, 2.7129230470588235f, 1.9459282392156863f,
        )
        private val DEFAULT_OPSIN_BIAS = floatArrayOf(-0.0037930732552754493f, -0.0037930732552754493f, -0.0037930732552754493f)
        private val DEFAULT_QUANT_BIAS = floatArrayOf(1.0f - 0.05465007330715401f, 1.0f - 0.07005449891748593f, 1.0f - 0.049935103337343655f)
    }
}


internal object FrameType {
    const val REGULAR = 0
    const val LF = 1
    const val REFERENCE_ONLY = 2
    const val SKIP_PROGRESSIVE = 3
}

internal object BlendMode {
    const val REPLACE = 0
    const val ADD = 1
    const val BLEND = 2
    const val MUL_ADD = 3
    const val MUL = 4
}

internal class BlendingInfo(val mode: Int, val alphaChannel: Int, val clamp: Boolean, val source: Int)

internal class EpfParams(
    val iters: Int,
    val sharpLut: FloatArray,
    val channelScale: FloatArray,
    val quantMul: Float,
    val pass0SigmaScale: Float,
    val pass2SigmaScale: Float,
    val borderSadMul: Float,
    val sigmaForModular: Float,
)

internal class Passes(val numPasses: Int, val shift: IntArray, val downsample: IntArray, val lastPass: IntArray) {
    /** パス [pass] で復号する Modular チャンネルのシフト範囲 [minShift, maxShift] */
    fun downsamplingBracket(pass: Int): IntArray {
        var maxShift = 2
        var minShift = 3
        var i = 0
        while (true) {
            for (j in downsample.indices) {
                if (i == lastPass[j]) {
                    when (downsample[j]) {
                        8 -> minShift = 3
                        4 -> minShift = 2
                        2 -> minShift = 1
                        1 -> minShift = 0
                    }
                }
            }
            if (i == numPasses - 1) minShift = 0
            if (i == pass) return intArrayOf(minShift, maxShift)
            maxShift = minShift - 1
            i++
        }
    }
}

internal class FrameHeader(
    val frameType: Int,
    val isModular: Boolean,
    val flags: Long,
    val doYCbCr: Boolean,
    val encodedColorChannels: Int,
    val jpegUpsampling: IntArray,
    val upsampling: Int,
    val ecUpsampling: IntArray,
    val groupSizeShift: Int,
    val xQmScale: Int,
    val bQmScale: Int,
    val passes: Passes,
    val lfLevel: Int,
    val x0: Int,
    val y0: Int,
    val width: Int,
    val height: Int,
    val blending: BlendingInfo,
    val ecBlending: List<BlendingInfo>,
    val duration: Int,
    val isLast: Boolean,
    val saveAsReference: Int,
    val saveBeforeCt: Boolean,
    val resetsCanvas: Boolean,
    val gabor: FloatArray?,
    val epf: EpfParams?,
) {
    val hasNoise: Boolean get() = flags and 1L != 0L
    val hasPatches: Boolean get() = flags and 2L != 0L
    val hasSplines: Boolean get() = flags and 0x10L != 0L
    val useLfFrame: Boolean get() = flags and 0x20L != 0L
    val skipAdaptiveLfSmoothing: Boolean get() = flags and 0x80L != 0L

    val isNormalFrame: Boolean get() = frameType == FrameType.REGULAR || frameType == FrameType.SKIP_PROGRESSIVE
    val isKeyframe: Boolean get() = isNormalFrame && (isLast || duration != 0)
    val canReference: Boolean get() = !isLast && (duration == 0 || saveAsReference != 0) && frameType != FrameType.LF

    fun sampleWidth(ups: Int): Int {
        var w = width
        if (ups > 1) w = (w + ups - 1) / ups
        if (lfLevel > 0) {
            val div = 1 shl (3 * lfLevel)
            w = (w + div - 1) shr (3 * lfLevel)
        }
        return w
    }

    fun sampleHeight(ups: Int): Int {
        var h = height
        if (ups > 1) h = (h + ups - 1) / ups
        if (lfLevel > 0) {
            val div = 1 shl (3 * lfLevel)
            h = (h + div - 1) shr (3 * lfLevel)
        }
        return h
    }

    val colorWidth: Int get() = sampleWidth(upsampling)
    val colorHeight: Int get() = sampleHeight(upsampling)
    val groupDim: Int get() = 128 shl groupSizeShift
    val lfGroupDim: Int get() = groupDim * 8
    val groupsPerRow: Int get() = (colorWidth + groupDim - 1) / groupDim
    val groupsPerColumn: Int get() = (colorHeight + groupDim - 1) / groupDim
    val numGroups: Int get() = groupsPerRow * groupsPerColumn
    val lfGroupsPerRow: Int get() = (colorWidth + lfGroupDim - 1) / lfGroupDim
    val lfGroupsPerColumn: Int get() = (colorHeight + lfGroupDim - 1) / lfGroupDim
    val numLfGroups: Int get() = lfGroupsPerRow * lfGroupsPerColumn

    fun lfGroupIdxFromGroupIdx(groupIdx: Int): Int {
        val gx = (groupIdx % groupsPerRow) / 8
        val gy = (groupIdx / groupsPerRow) / 8
        return gx + gy * lfGroupsPerRow
    }

    companion object {
        fun parse(br: BitReader, image: ImageHeader): FrameHeader {
            val allDefault = br.readBool()
            var frameType = FrameType.REGULAR
            var modular = false
            var flags = 0L
            var doYCbCr = false
            if (!allDefault) {
                frameType = br.readBits(2)
                modular = br.readBits(1) == 1
                flags = br.readU64()
                if (!image.xybEncoded) doYCbCr = br.readBool()
            }
            val useLfFrame = flags and 0x20L != 0L
            val jpegUpsampling = IntArray(3)
            if (doYCbCr && !useLfFrame) {
                for (i in 0 until 3) jpegUpsampling[i] = br.readBits(2)
            }
            var upsampling = 1
            val numExtra = image.extraChannels.size
            val ecUpsampling = IntArray(numExtra) { 1 }
            if (!allDefault && !useLfFrame) {
                upsampling = 1 shl br.readBits(2)
                for (i in 0 until numExtra) ecUpsampling[i] = 1 shl br.readBits(2)
            }
            var groupSizeShift = 1
            if (modular) groupSizeShift = br.readBits(2)
            var xQm = if (image.xybEncoded && !modular) 3 else 2
            var bQm = 2
            if (!allDefault && image.xybEncoded && !modular) {
                xQm = br.readBits(3)
                bQm = br.readBits(3)
            }
            var passes = Passes(1, IntArray(0), IntArray(0), IntArray(0))
            if (!allDefault && frameType != FrameType.REFERENCE_ONLY) {
                val numPasses = br.readU32(1, 0, 2, 0, 3, 0, 4, 3)
                if (numPasses != 1) {
                    val numDs = br.readU32(0, 0, 1, 0, 2, 0, 3, 1)
                    if (numDs >= numPasses) fail("invalid passes")
                    val shift = IntArray(numPasses - 1) { br.readBits(2) }
                    val down = IntArray(numDs) { 1 shl br.readBits(2) }
                    val last = IntArray(numDs) { br.readU32(0, 0, 1, 0, 2, 0, 0, 3) }
                    passes = Passes(numPasses, shift, down, last)
                } else {
                    passes = Passes(1, IntArray(0), IntArray(0), IntArray(0))
                }
            }
            var lfLevel = 0
            if (frameType == FrameType.LF) lfLevel = br.readBits(2) + 1
            var haveCrop = false
            if (!allDefault && frameType != FrameType.LF) haveCrop = br.readBool()
            var x0 = 0
            var y0 = 0
            var width = image.width
            var height = image.height
            if (haveCrop) {
                if (frameType != FrameType.REFERENCE_ONLY) {
                    x0 = unpackSigned(br.readU32(0, 8, 256, 11, 2304, 14, 18688, 30))
                    y0 = unpackSigned(br.readU32(0, 8, 256, 11, 2304, 14, 18688, 30))
                }
                width = br.readU32(0, 8, 256, 11, 2304, 14, 18688, 30)
                height = br.readU32(0, 8, 256, 11, 2304, 14, 18688, 30)
            }
            val isNormal = frameType == FrameType.REGULAR || frameType == FrameType.SKIP_PROGRESSIVE
            val fullFrame = x0 <= 0 && y0 <= 0 &&
                x0.toLong() + width >= image.width && y0.toLong() + height >= image.height

            fun readBlending(): BlendingInfo {
                val mode = br.readU32(0, 0, 1, 0, 2, 0, 3, 2)
                if (mode > 4) fail("invalid blend mode")
                var alpha = 0
                var clamp = false
                if (numExtra > 0 && (mode == BlendMode.BLEND || mode == BlendMode.MUL_ADD)) {
                    alpha = br.readU32(0, 0, 1, 0, 2, 0, 3, 3)
                }
                if ((numExtra > 0 && (mode == BlendMode.BLEND || mode == BlendMode.MUL_ADD)) || mode == BlendMode.MUL) {
                    clamp = br.readBool()
                }
                // libjxl と同じく、各チャンネル自身のモードで判断する
                var source = 0
                if (mode != BlendMode.REPLACE || (haveCrop && !fullFrame)) source = br.readBits(2)
                return BlendingInfo(mode, alpha, clamp, source)
            }

            var blending = BlendingInfo(BlendMode.REPLACE, 0, false, 0)
            var ecBlending = List(numExtra) { BlendingInfo(BlendMode.REPLACE, 0, false, 0) }
            if (!allDefault && isNormal) {
                blending = readBlending()
                ecBlending = List(numExtra) { readBlending() }
            }
            var duration = 0
            if (!allDefault && isNormal && image.animation != null) {
                duration = br.readU32(0, 0, 1, 0, 0, 8, 0, 32)
            }
            if (!allDefault && isNormal && image.animation?.haveTimecodes == true) br.readBits(32)
            var isLast = frameType == FrameType.REGULAR
            if (!allDefault && isNormal) isLast = br.readBool()
            var saveAsRef = 0
            if (!allDefault && frameType != FrameType.LF && !isLast) saveAsRef = br.readBits(2)
            val resetsCanvas = blending.mode == BlendMode.REPLACE && (!haveCrop || fullFrame)
            val canRef = !isLast && (duration == 0 || saveAsRef != 0) && frameType != FrameType.LF
            var saveBeforeCt = !isNormal
            if (!allDefault && (frameType == FrameType.REFERENCE_ONLY || (resetsCanvas && canRef))) {
                saveBeforeCt = br.readBool()
            }
            var gabor: FloatArray? = DEFAULT_GABOR
            var epf: EpfParams? = defaultEpf()
            if (!allDefault) {
                // name
                val nameLen = br.readU32(0, 0, 0, 4, 16, 5, 48, 10)
                for (i in 0 until nameLen) br.readBits(8)
                // RestorationFilter
                if (!br.readBool()) {
                    if (br.readBool()) {
                        if (br.readBool()) {
                            val w = FloatArray(6) { br.readF16() }
                            for (c in 0 until 3) {
                                if (kotlin.math.abs(1f + (w[2 * c] + w[2 * c + 1]) * 4f) < 1.1920929E-7f) fail("invalid gaborish weights")
                            }
                            gabor = w
                        }
                    } else {
                        gabor = null
                    }
                    val iters = br.readBits(2)
                    if (iters == 0) {
                        epf = null
                    } else {
                        var sharp = DEFAULT_SHARP_LUT
                        if (!modular && br.readBool()) sharp = FloatArray(8) { br.readF16() }
                        var scale = DEFAULT_CHANNEL_SCALE
                        if (br.readBool()) {
                            scale = FloatArray(3) { br.readF16() }
                            br.readBits(32)
                        }
                        var quantMul = 0.46f
                        var s0 = 0.9f
                        var s2 = 6.5f
                        var border = 2.0f / 3.0f
                        if (br.readBool()) {
                            if (!modular) quantMul = br.readF16()
                            s0 = br.readF16()
                            s2 = br.readF16()
                            border = br.readF16()
                        }
                        var sigmaModular = 1.0f
                        if (modular) sigmaModular = br.readF16()
                        epf = EpfParams(iters, sharp, scale, quantMul, s0, s2, border, sigmaModular)
                    }
                    br.extensions()
                }
                br.extensions()
            }
            val grey = image.isGrey
            val encodedColor = if (modular && !doYCbCr && !image.xybEncoded && grey) 1 else 3
            if (upsampling > 1) {
                for (u in ecUpsampling) if (u < upsampling) fail("ec_upsampling < upsampling")
            }
            return FrameHeader(
                frameType, modular, flags, doYCbCr, encodedColor, jpegUpsampling, upsampling, ecUpsampling,
                groupSizeShift, xQm, bQm, passes, lfLevel, x0, y0, width, height, blending, ecBlending, duration,
                isLast, saveAsRef, saveBeforeCt, resetsCanvas, gabor, epf,
            )
        }

        private val DEFAULT_GABOR = floatArrayOf(0.115169525f, 0.061248592f, 0.115169525f, 0.061248592f, 0.115169525f, 0.061248592f)
        private val DEFAULT_SHARP_LUT = floatArrayOf(0f, 1f / 7f, 2f / 7f, 3f / 7f, 4f / 7f, 5f / 7f, 6f / 7f, 1f)
        private val DEFAULT_CHANNEL_SCALE = floatArrayOf(40.0f, 5.0f, 3.5f)

        private fun defaultEpf() = EpfParams(2, DEFAULT_SHARP_LUT, DEFAULT_CHANNEL_SCALE, 0.46f, 0.9f, 6.5f, 2.0f / 3.0f, 1.0f)
    }
}

/** TOC: 各セクションのバイト範囲 (data 内の絶対オフセット) */
internal class Toc(val offsets: IntArray, val sizes: IntArray, val singleEntry: Boolean) {
    companion object {
        fun parse(br: BitReader, header: FrameHeader): Toc {
            val numGroups = header.numGroups
            val numPasses = header.passes.numPasses
            val single = numGroups == 1 && numPasses == 1
            val count = if (single) 1 else 1 + header.numLfGroups + 1 + numGroups * numPasses
            if (count > 65536 * 16) fail("too many TOC entries")
            var permutation: IntArray? = null
            if (br.readBool()) {
                val code = EntropyCode.parse(br, 8)
                val reader = code.newReader()
                reader.begin(br)
                permutation = readPermutation(br, reader, count, 0)
                reader.finish()
            }
            br.zeroPadToByte()
            val sizes = IntArray(count) { br.readU32(0, 10, 1024, 14, 17408, 22, 4211712, 30) }
            br.zeroPadToByte()
            val base = br.bytePosition()
            val offsetsBitstream = IntArray(count)
            var acc = base.toLong()
            for (i in 0 until count) {
                offsetsBitstream[i] = acc.toInt()
                acc += sizes[i]
                if (acc > Int.MAX_VALUE) fail("TOC too large")
            }
            // セクション i の中身は、ビットストリーム上 permutation[i] 番目にある
            val offsets = IntArray(count)
            val outSizes = IntArray(count)
            if (permutation != null) {
                for (i in 0 until count) {
                    offsets[i] = offsetsBitstream[permutation[i]]
                    outSizes[i] = sizes[permutation[i]]
                }
            } else {
                System.arraycopy(offsetsBitstream, 0, offsets, 0, count)
                System.arraycopy(sizes, 0, outSizes, 0, count)
            }
            return Toc(offsets, outSizes, single)
        }
    }

    val end: Long get() {
        var max = 0L
        for (i in offsets.indices) max = maxOf(max, offsets[i].toLong() + sizes[i])
        return max
    }
}
