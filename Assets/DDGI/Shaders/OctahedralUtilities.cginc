/** Efficient GPU implementation of the octahedral unit vector encoding from
    Cigolle, Donow, Evangelakos, Mara, McGuire, Meyer,
    A Survey of Efficient Representations for Independent Unit Vectors, Journal of Computer Graphics Techniques (JCGT), vol. 3, no. 2, 1-30, 2014
    Available online http://jcgt.org/published/0003/02/01/
*/

float2 signNotZero(float2 v)
{
    return float2(
        v.x >= 0 ? 1 : -1,
        v.y >= 0 ? 1 : -1
    );
}

/** Assumes that v is a unit vector. The result is an octahedral vector on the [-1, +1] square. */
float2 octEncode(in float3 v) {
    float l1norm = abs(v.x) + abs(v.y) + abs(v.z);
    float2 result = v.xy * (1.0 / l1norm);
    if (v.z < 0.0) {
        result = (1.0 - abs(result.yx)) * signNotZero(result.xy);
    }
    return result;
}


/** Returns a unit vector. Argument o is an octahedral vector packed via octEncode,
    on the [-1, +1] square*/
float3 octDecode(float2 o) {
    float3 v = float3(o.x, o.y, 1.0 - abs(o.x) - abs(o.y));
    if (v.z < 0.0) {
        v.xy = (1.0 - abs(v.yx)) * signNotZero(v.xy);
    }
    return normalize(v);
}


// Compute normalized oct coord, mapping top left of top left pixel to (-1,-1)
float2 normalizedOctCoord(uint2 fragCoord, int PROBE_SIDE_LENGTH) {
    int probeWithBorderSide = PROBE_SIDE_LENGTH + 2;

    float2 octFragCoord = int2(fragCoord.x % probeWithBorderSide, fragCoord.y % probeWithBorderSide) - int2(1, 1);
    // Add back the half pixel to get pixel center normalized coordinates
    return (octFragCoord / (float)PROBE_SIDE_LENGTH) * 2.0 - 1.0;
}