int PROBE_SIDE_LENGTH;

int texelToIndex(float2 texelXY, uint width) 
{
    int probeWithBorderSide = PROBE_SIDE_LENGTH + 2;
    int probesPerSide = width / probeWithBorderSide;
    return int(texelXY.x / probeWithBorderSide) + probesPerSide * int(texelXY.y / probeWithBorderSide);
}

/** 
 \param probeCoords Integer (stored in float) coordinates of the probe on the probe grid 
 */
 int gridCoordToProbeIndex(in DDGIVolume ddgiVolume, in float3 probeCoords) {
    //probeCoords = int3(mod(probeCoords + ddgiVolume.phaseOffsets, ddgiVolume.probeCounts));
    return int(probeCoords.x + probeCoords.y * ddgiVolume.probeCounts.x + probeCoords.z * ddgiVolume.probeCounts.x * ddgiVolume.probeCounts.y);
}


int3 baseGridCoord(in DDGIVolume ddgiVolume, float3 X) {
    // Implicit floor in the convert to int
   int3 unOffsetGridCoord = clamp(int3((X - ddgiVolume.probeGridOrigin) * ddgiVolume.invProbeSpacing),
                int3(0, 0, 0), 
                ddgiVolume.probeCounts - 1);
   return int3(fmod((unOffsetGridCoord - ddgiVolume.phaseOffsets) , ddgiVolume.probeCounts));
}


/** Returns the index of the probe at the floor along each dimension. */
int baseProbeIndex(in DDGIVolume ddgiVolume, float3 X) {
    return gridCoordToProbeIndex(ddgiVolume, baseGridCoord(ddgiVolume, X));
}


/** Matches code in LightFieldModel::debugDraw() */
float3 gridCoordToColor(int3 gridCoord) {
    gridCoord.x &= 1;
    gridCoord.y &= 1;
    gridCoord.z &= 1;

    if (gridCoord.x + gridCoord.y + gridCoord.z == 0) {
        return float3(0.1, 0.1, 0.1);
    } else {
        return float3(gridCoord) * 0.9;
    }
}

int findMSB(int x)
{
  int i;
  int mask;
  int res = -1;
  if (x < 0)
  {
    int not = ~x;
    x = not;
  }
  for(i = 0; i < 32; i++)
  {
    mask = 0x80000000 >> i;
    if (x & mask) {
      res = 31 - i;
      break;
    }
  }
  return res;
}

int3 probeIndexToGridCoord(in DDGIVolume ddgiVolume, int index) {    
    /* Works for any # of probes */
    int3 iPos;
    iPos.x = index % ddgiVolume.probeCounts.x;
    iPos.y = (index % (ddgiVolume.probeCounts.x * ddgiVolume.probeCounts.y)) / ddgiVolume.probeCounts.x;
    iPos.z = index / (ddgiVolume.probeCounts.x * ddgiVolume.probeCounts.y);

    // // Assumes probeCounts are powers of two.
    // // Saves ~10ms compared to the divisions above
    // // Precomputing the MSB actually slows this code down substantially
    // int3 iPos;
    // iPos.x = index & (ddgiVolume.probeCounts.x - 1);
    // iPos.y = (index & ((ddgiVolume.probeCounts.x * ddgiVolume.probeCounts.y) - 1)) >> findMSB(ddgiVolume.probeCounts.x);
    // iPos.z = index >> findMSB(ddgiVolume.probeCounts.x * ddgiVolume.probeCounts.y);

    return iPos;
}


float3 probeIndexToColor(in DDGIVolume ddgiVolume, int index) {
    return gridCoordToColor(probeIndexToGridCoord(ddgiVolume, index));
}


/** probeCoords Coordinates of the probe, computed as part of the process. */
int nearestProbeIndex(in DDGIVolume ddgiVolume, float3 X, out float3 probeCoords) {
    probeCoords = clamp(round((X - ddgiVolume.probeGridOrigin) / ddgiVolume.probeSpacing),
                    float3(0, 0, 0), 
                    float3(ddgiVolume.probeCounts - 1));

    return gridCoordToProbeIndex(ddgiVolume, probeCoords);
}


float highestSignedValue;

float4 readProbeOffset(in DDGIVolume ddgiVolume, int2 texelCoord) {
    float4 v = Load(probeOffsetsTexture, texelCoord, probeOffsetsTextureSize);
    return float4(v.xyz * ddgiVolume.probeOffsetLimit * ddgiVolume.probeSpacing / highestSignedValue, v.w);
}

void writeProbeOffset(in DDGIVolume ddgiVolume, in int2 texelCoord, in float4 offsetAndFlags) {
    Write(probeOffsetsImage, texelCoord, probeOffsetsImageSize, int4(int3(ceil(highestSignedValue * offsetAndFlags.xyz * ddgiVolume.invProbeSpacing / ddgiVolume.probeOffsetLimit)), offsetAndFlags.w));
}

// Apply the per-axis phase offset to derive the correct location for each probe.
float3 gridCoordToPositionNoOffset(in DDGIVolume ddgiVolume, int3 c) {
    // Phase offset may be negative, which is fine for modular arithmetic.
    int3 phaseOffsetGridCoord = int3(fmod(c + ddgiVolume.phaseOffsets, ddgiVolume.probeCounts));
    return ddgiVolume.probeSpacing * float3(phaseOffsetGridCoord) + ddgiVolume.probeGridOrigin;
}

float3 gridCoordToPosition(in DDGIVolume ddgiVolume, int3 c) {

    //Add per-probe offset
    uint idx = gridCoordToProbeIndex(ddgiVolume, c);
    float probeXY = ddgiVolume.probeCounts.x * ddgiVolume.probeCounts.y;
    int2 C = probeXY != 0 ? int2(fmod(idx, probeXY), idx / probeXY) : int2(0, 0);

    float3 offset =
#if FIRST_FRAME == 1
        int3(0, 0, 0);
#else
        readProbeOffset(ddgiVolume, C).xyz; // readProbeOffset multiplies by probe step.
#endif
    return gridCoordToPositionNoOffset(ddgiVolume, c) + offset;
}

float3 probeLocation(in DDGIVolume ddgiVolume, int index) {
    return gridCoordToPosition(ddgiVolume, probeIndexToGridCoord(ddgiVolume, index));
}

/**
   \param baseProbeIndex Index into ddgiVolume.radianceProbeGrid's TEXTURE_2D_ARRAY. This is the probe
   at the floor of the current ray sampling position.

   \param relativeIndex on [0, 7]. This is used as a set of three 1-bit offsets

   Returns a probe index into ddgiVolume.radianceProbeGrid. It may be the *same* index as 
   baseProbeIndex.

   This will wrap in crazy ways when the camera is outside of the bounding box
   of the probes...but that's ok. If that case arises, then the trace is likely to 
   be poor quality anyway. Regardless, this function will still return the index 
   of some valid probe, and that probe can either be used or fail because it does not 
   have visibility to the location desired.

   \see nextCycleIndex, baseProbeIndex
 */
int relativeProbeIndex(in DDGIVolume ddgiVolume, int baseProbeIndex, int relativeIndex) {
    int numProbes = ddgiVolume.probeCounts.x * ddgiVolume.probeCounts.y * ddgiVolume.probeCounts.z;

    // Use the bits of 0 <= relativeIndex < 8 to enumerate the +1 or +0 offsets along each axis.
    //
    // relativeIndex bit 0 = x offset
    // relativeIndex bit 1 = y offset
    // relativeIndex bit 2 = z offset
    int3 offset = int3(relativeIndex & 1, (relativeIndex >> 1) & 1, (relativeIndex >> 2) & 1);
    int3 stride = int3(1, ddgiVolume.probeCounts.x, ddgiVolume.probeCounts.x * ddgiVolume.probeCounts.y);
    
    // If the probe is outside of the grid, return *some* probe so that the code 
    // doesn't crash. With cascades implemented (even one), this case is never needed because
    // the cascade will fade out before the sample point leaves the grid.
    //
    // (numProbes is guaranteed to be a power of 2 in the current implementation, 
    // which allows us to use a bitand instead of a modulo operation.)
    return (baseProbeIndex + dot(offset, stride)) & (numProbes - 1);
}

float2 probeTextureCoordFromDirection
   (float3             dir, 
    int3           probeGridCoord,
    const in bool       useIrradiance,
    DDGIVolume          ddgiVolume) {
    int probeSideLength = useIrradiance ? ddgiVolume.irradianceProbeSideLength : ddgiVolume.visibilityProbeSideLength;

    float2 signedOct = octEncode(dir);
    float2 unsignedOct = (signedOct + 1.0f) * 0.5f;
    int probeWithBorderSide = probeSideLength + 2;

    float2 coordInProbe = unsignedOct * probeSideLength + float2(1,1);
    uint2 probeTexCoordStart = uint2((probeGridCoord.x + probeGridCoord.y * ddgiVolume.probeCounts.x) * probeWithBorderSide,
        probeGridCoord.z * probeWithBorderSide);
    
    return probeTexCoordStart + (uint2) coordInProbe;
}