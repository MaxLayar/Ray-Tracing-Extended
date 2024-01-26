struct DDGIVolume
{
    int3 probeCounts;
    int3 logProbeCounts;
    float3 probeGridOrigin;
    float3 probeSpacing;
    float3 invProbeSpacing; // 1 / probeSpacing
    int3 phaseOffsets;
    uint2 invIrradianceTextureSize;
    uint2 invVisibilityTextureSize;
    // probeOffsetLimit on [0,0.5] where max probe 
    // offset = probeOffsetLimit * probeSpacing
    // Usually 0.4, controllable from GUI.
    float probeOffsetLimit;
    int irradianceProbeSideLength;
    int visibilityProbeSideLength;
    float selfShadowBias;
    float irradianceGamma;
    float invIrradianceGamma;

    //DEBUG
    float debugMeanBias; // 0 in production code
    float debugVarianceBias; // 0 in production code
    float debugChebyshevBias; // 0 in production code
    float debugChebyshevNormalize; // 1/(1-debugChebyshevBias) = 1 in production code
    int cameraLocked;

    float energyConservation;
    float depthSharpness;
    float hysteresis;
    float maxDistance;
};

StructuredBuffer<DDGIVolume> DDGIVolumes;

RWStructuredBuffer<float4> irradianceTexture;
int2 irradianceTextureSize;
RWStructuredBuffer<float4> visibilityTexture;
int2 visibilityTextureSize;
RWStructuredBuffer<float4> probeOffsetsTexture;
int2 probeOffsetsTextureSize;
RWStructuredBuffer<float4> probeOffsetsImage;
int2 probeOffsetsImageSize;