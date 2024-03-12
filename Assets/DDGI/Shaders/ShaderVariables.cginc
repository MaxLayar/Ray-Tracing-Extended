struct DDGIVolume
{
    uint3 probeCounts;
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

    float depthSharpness;
    float hysteresis;
    float energyConservation;
    float maxDistance;
};

StructuredBuffer<DDGIVolume> DDGIVolumes;

Texture2D<float4> irradianceTexture;
SamplerState sampler_irradianceTexture;

Texture2D<float4> visibilityTexture;
SamplerState sampler_visibilityTexture;

RWTexture2D<float4> Result;

RWStructuredBuffer<float4> probeOffsetsTexture;
uint2 probeOffsetsTextureSize;

RWStructuredBuffer<float4> probeOffsetsImage;
uint2 probeOffsetsImageSize;