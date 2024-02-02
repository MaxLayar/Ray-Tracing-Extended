// Returns ±1
float2 SignNotZero(float2 v) {
    return float2((v.x >= 0.0) ? +1.0 : -1.0, (v.y >= 0.0) ? +1.0 : -1.0);
}


uint CoordToIndex(uint3 coord, uint3 size)
{
    return coord.z * size.x * size.y + coord.y * size.x + coord.x;
}

uint3 IndexToCoord(uint index, uint3 size)
{
    int3 result;
    int a = (size.x * size.y);
    result.z = index / a;
    int b = index - a * result.z;
    result.y = b / size.x;
    result.x = b % size.x;
    return result;
}

float4 Load(RWStructuredBuffer<float4> buffer, uint2 coord, uint2 size)
{
    return buffer[coord.x + coord.y * size.x];
}

void Write(RWStructuredBuffer<float4> buffer, uint2 coord, uint2 size, float4 value)
{
    buffer[coord.x + coord.y * size.x] = value;
}