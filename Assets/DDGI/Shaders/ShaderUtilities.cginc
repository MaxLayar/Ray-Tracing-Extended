float4 Load(RWStructuredBuffer<float4> buffer, uint2 coord, uint2 size)
{
    return buffer[coord.x + coord.y * size.x];
}

float4 LoadBilinear(RWStructuredBuffer<float4> buffer, float2 coord, float2 size)
{
    // Clamp coordinates to valid range
    coord = clamp(coord, 0.5, size - 0.5);

    // Calculate integer and fractional parts
    float2 flooredCoord = floor(coord - 0.5) + 0.5;
    float2 fraction = coord - flooredCoord;

    // Sample corner values
    uint2 coord00 = flooredCoord;
    uint2 coord10 = flooredCoord + uint2(1, 0);
    uint2 coord01 = flooredCoord + uint2(0, 1);
    uint2 coord11 = flooredCoord + uint2(1, 1);

    float4 p00 = Load(buffer, coord00, size);
    float4 p10 = Load(buffer, coord10, size);
    float4 p01 = Load(buffer, coord01, size);
    float4 p11 = Load(buffer, coord11, size);

    // Perform bilinear interpolation
    float4 lerpX0 = lerp(p00, p10, fraction.x);
    float4 lerpX1 = lerp(p01, p11, fraction.x);
    return lerp(lerpX0, lerpX1, fraction.y);
}

void Write(RWStructuredBuffer<float4> buffer, uint2 coord, uint2 size, float4 value)
{
    buffer[coord.x + coord.y * size.x] = value;
}