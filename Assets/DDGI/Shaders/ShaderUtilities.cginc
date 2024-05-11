float4 Load(RWStructuredBuffer<float4> buffer, uint2 coord, uint2 size)
{
    return buffer[coord.x + coord.y * size.x];
}

void Write(RWStructuredBuffer<float4> buffer, uint2 coord, uint2 size, float4 value)
{
    buffer[coord.x + coord.y * size.x] = value;
}