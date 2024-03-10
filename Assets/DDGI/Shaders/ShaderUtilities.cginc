float4 Load(RWStructuredBuffer<float4> buffer, uint2 coord, uint2 size)
{
    return buffer[coord.x + coord.y * size.x];
}

float4 LoadBilinear(RWStructuredBuffer<float4> buffer, float2 coord, float2 size)
{
    //bilinear interpolation
    float2 flooredCoord = floor(coord);
    float2 fractionalCoord = frac(coord);

    // Calculate adjacent integer coordinates for bilinear interpolation
    uint2 coord00 = uint2(flooredCoord);
    uint2 coord10 = uint2(flooredCoord.x + 1, flooredCoord.y);
    uint2 coord01 = uint2(flooredCoord.x, flooredCoord.y + 1);
    uint2 coord11 = uint2(flooredCoord.x + 1, flooredCoord.y + 1);

    // Clamp coordinates to texture bounds
    coord00 = clamp(coord00, uint2(0, 0), size - uint2(1, 1));
    coord10 = clamp(coord10, uint2(0, 0), size - uint2(1, 1));
    coord01 = clamp(coord01, uint2(0, 0), size - uint2(1, 1));
    coord11 = clamp(coord11, uint2(0, 0), size - uint2(1, 1));

    // Sample corner values
    float4 p00 = Load(buffer, coord00, size);
    float4 p10 = Load(buffer, coord10, size);
    float4 p01 = Load(buffer, coord01, size);
    float4 p11 = Load(buffer, coord11, size);

    // Perform bilinear interpolation
    float4 lerpX0 = lerp(p00, p10, fractionalCoord.x);
    float4 lerpX1 = lerp(p01, p11, fractionalCoord.x);
    return lerp(lerpX0, lerpX1, fractionalCoord.y);
}

void Write(RWStructuredBuffer<float4> buffer, uint2 coord, uint2 size, float4 value)
{
    buffer[coord.x + coord.y * size.x] = value;
}