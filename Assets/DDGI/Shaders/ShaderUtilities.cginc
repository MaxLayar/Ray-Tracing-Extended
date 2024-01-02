// Returns ±1
float2 SignNotZero(float2 v) {
    return float2((v.x >= 0.0) ? +1.0 : -1.0, (v.y >= 0.0) ? +1.0 : -1.0);
}

// Assume normalized input. Output is on [-1, 1] for each component.
float2 SphereToOctProjected(float3 v) {
    // Project the sphere onto the octahedron, and then onto the xy plane
    float2 p = v.xy * (1.0 / (abs(v.x) + abs(v.y) + abs(v.z)));

    // Reflect the folds of the lower hemisphere over the diagonals
    return (v.z <= 0.0) ? ((1.0 - abs(p.yx)) * SignNotZero(p)) : p;
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