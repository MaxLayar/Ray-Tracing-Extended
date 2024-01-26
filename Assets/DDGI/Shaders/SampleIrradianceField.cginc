#include "UnityShaderVariables.cginc"
#include "UnityShaderUtilities.cginc"
#include "UnityInstancing.cginc"

#include "OctahedralUtilities.cginc"
#include "ShaderUtilities.cginc"
#include "ShaderVariables.cginc"

#define M_PI 3.1415926535897932384626433832795
#define NUM_DDGIVOLUMES 1

int OFFSET_BITS_PER_CHANNEL;
float highestSignedValue;

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
   return int3(modf((unOffsetGridCoord - ddgiVolume.phaseOffsets) , ddgiVolume.probeCounts));
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
    /*
    iPos.x = index % ddgiVolume.probeCounts.x;
    iPos.y = (index % (ddgiVolume.probeCounts.x * ddgiVolume.probeCounts.y)) / ddgiVolume.probeCounts.x;
    iPos.z = index / (ddgiVolume.probeCounts.x * ddgiVolume.probeCounts.y);
    */

    // Assumes probeCounts are powers of two.
    // Saves ~10ms compared to the divisions above
    // Precomputing the MSB actually slows this code down substantially
    int3 iPos;
    iPos.x = index & (ddgiVolume.probeCounts.x - 1);
    iPos.y = (index & ((ddgiVolume.probeCounts.x * ddgiVolume.probeCounts.y) - 1)) >> findMSB(ddgiVolume.probeCounts.x);
    iPos.z = index >> findMSB(ddgiVolume.probeCounts.x * ddgiVolume.probeCounts.y);

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
    int3 phaseOffsetGridCoord = int3(modf(c + ddgiVolume.phaseOffsets, ddgiVolume.probeCounts));
    return ddgiVolume.probeSpacing * float3(phaseOffsetGridCoord) + ddgiVolume.probeGridOrigin;
}

float3 gridCoordToPosition(in DDGIVolume ddgiVolume, int3 c) {

    //Add per-probe offset
    uint idx = gridCoordToProbeIndex(ddgiVolume, c);
    float probeXY = ddgiVolume.probeCounts.x * ddgiVolume.probeCounts.y;
    int2 C = probeXY != 0 ? int2(modf(idx, probeXY), idx / probeXY) : int2(0, 0);

    float3 offset =
#if FIRST_FRAME
        int3(0);
#else
        readProbeOffset(ddgiVolume, C).xyz; // readProbeOffset multiplies by probe step.
#endif
    return gridCoordToPositionNoOffset(ddgiVolume, c) + offset;
}

float3 probeLocation(in DDGIVolume ddgiVolume, int index) {
    return gridCoordToPosition(ddgiVolume, probeIndexToGridCoord(ddgiVolume, index));
}


/** GLSL's dot on ivec3 returns a float. This is an all-integer version */
int idot(int3 a, int3 b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

float square(float f)
{
    return f * f;
}

float3 square(float3 f)
{
    return float3(f.x * f.x, f.y * f.y, f.z * f.z);
}

float pow3(float f)
{
    return f * f * f;
}

float pow(float a, float b)
{
    return exp(a * log(b));
}

float3 pow(float3 f, float3 x)
{
    return float3(pow(f.x, x.x), pow(f.y, x.y), pow(f.z, x.z));
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
    return (baseProbeIndex + idot(offset, stride)) & (numProbes - 1);
}

float2 probeTextureCoordFromDirection
   (float3             dir, 
    int3           probeGridCoord,
    const in bool       useIrradiance,
    DDGIVolume          ddgiVolume) {

    float2 invTextureSize = useIrradiance ? ddgiVolume.invIrradianceTextureSize : ddgiVolume.invVisibilityTextureSize;
    int probeSideLength = useIrradiance ? ddgiVolume.irradianceProbeSideLength : ddgiVolume.visibilityProbeSideLength;

    float2 signedOct = octEncode(dir);
    float2 unsignedOct = (signedOct + 1.0f) * 0.5f;
    float2 octCoordNormalizedToTextureDimensions = (unsignedOct * (float)probeSideLength) * invTextureSize;
    
    int probeWithBorderSide = probeSideLength + 2;

    float2 probeTopLeftPosition = float2((probeGridCoord.x + probeGridCoord.y * ddgiVolume.probeCounts.x) * probeWithBorderSide,
        probeGridCoord.z * probeWithBorderSide) + float2(1,1);

    float2 normalizedProbeTopLeftPosition = float2(probeTopLeftPosition) * invTextureSize;

    return float2(normalizedProbeTopLeftPosition + octCoordNormalizedToTextureDimensions);
}


/**
  Result.rgb = Irradiance3
  Result.a   = weight based on wsPosition vs. probe bounds
*/
float4 sampleOneDDGIVolume
   (DDGIVolume             ddgiVolume,
    float3                 wsPosition,
    float3                offsetPos,
    float3                sampleDirection,
	float3				   cameraPos,
    // Set to false in production 
    const bool             debugDisableBackface,
    // Set to false in production 
    const bool             debugDisableChebyshev,
    // Set to -1 in production
    const int              debugProbeIndex,

    // Can we skip this volume if it has zero weight?
    bool                   skippable) {

    highestSignedValue = float(1 << (OFFSET_BITS_PER_CHANNEL - 1));
    // Compute the weight for this volume relative to other volumes
    float volumeWeight = 1.0;
    if (skippable) {
        // Compute the non-integer baseGridCoord. Use the unshifted position so that weights are consistent between
        // volumes. Use the geometric mean across all axes.
		float3 shiftedOrigin = ddgiVolume.probeGridOrigin;

		if (ddgiVolume.cameraLocked) {
			shiftedOrigin = cameraPos - (ddgiVolume.probeSpacing * (ddgiVolume.probeCounts - float3(1, 1, 1)) * 0.5);
		}
        float3 realGridCoord = (wsPosition - shiftedOrigin) * ddgiVolume.invProbeSpacing;
        for (int axis = 0; axis < 3; ++axis) {
            float a = realGridCoord[axis];
            if (a < 1.0) {
                volumeWeight *= clamp(a, 0.0, 1.0);
            } else if (a > float(ddgiVolume.probeCounts[axis]) - 2.0 - (ddgiVolume.cameraLocked ? 1.0 : 0.0)) {
				volumeWeight *= clamp(float(ddgiVolume.probeCounts[axis]) - 1.0 - (ddgiVolume.cameraLocked ? 1.0 : 0.0) - a, 0.0, 1.0);
            }
        }
		// Blending is improved without logarithmic fallof
        //volumeWeight = pow(volumeWeight, 1.0 / 3.0);
    }

    if (volumeWeight == 0.0) {
        // Don't bother sampling, this volume won't be used
        return float4(0,0,0,0);
    }

    offsetPos *= ddgiVolume.selfShadowBias;
    
    const bool debugDisableNonlinear = debugDisableChebyshev;

    // We're sampling at (wsPosition + offsetPos). This is inside of some grid cell, which
    // is bounded by 8 probes. Find the coordinate of the corner for the LOWEST probe (i.e.,
    // floor along x,y,z) and call that baseGridCoord. The other seven probes are offset by
    // +0 or +1 in grid space along each axis from this. We'll process them all in the main
    // loop below.
    //
    // This is all analogous to bilinear interpolation for texture maps, but we're doing it
    // in 3D, with visibility, and nonlinearly.
    int3 anchorGridCoord = baseGridCoord(ddgiVolume, wsPosition + offsetPos);

    // Don't use the offset to compute trilinear.
    float3 baseProbePos = //gridCoordToPosition(ddgiVolume, baseGridCoord);
        gridCoordToPositionNoOffset(ddgiVolume, anchorGridCoord);

    // Weighted irradiance accumulated in RGB across probes. The Alpha channel contains the
    // sum of the weights, which is used for normalization at the end.
    float4 irradiance = float4(0,0,0,0);
    
    // `alpha` is how far from the floor(currentVertex) position. On [0, 1] for each axis.
    // Used for trilinear weighting. 
    float3 alpha = clamp((wsPosition + offsetPos - baseProbePos) * ddgiVolume.invProbeSpacing, float3(0,0,0), float3(1,1,1));
    
    // This term is experimental and not in use in the current implementation
    float chebWeightSum = 0;

    // Iterate over adjacent probe cage
    for (int i = 0; i < 8; ++i) {
        // Compute the offset grid coord and clamp to the probe grid boundary
        // Offset = 0 or 1 along each axis. Pull the offsets from the bits of the 
        // loop index: x = bit 0, y = bit 1, z = bit 2
        int3  offset = int3(i, i >> 1, i >> 2) & int3(1,1,1);

        // Compute the trilinear weights based on the grid cell vertex to smoothly
        // transition between probes. Offset is binary, so we're
        // using 1-a when offset = 0 and a when offset = 1.
        float3 trilinear3 = max(float3(0.001,0.001,0.001), lerp(1.0 - alpha, alpha, offset));
        float trilinear = trilinear3.x * trilinear3.y * trilinear3.z;

		// Because of the phase offset applied for camera locked volumes,
		// we need to add the computed offset modulo the probecounts.
        int3 probeGridCoord = int3(modf((anchorGridCoord + offset), ddgiVolume.probeCounts));

        // Make cosine falloff in tangent plane with respect to the angle from the surface to the probe so that we never
        // test a probe that is *behind* the surface.
        // It doesn't have to be cosine, but that is efficient to compute and we must clip to the tangent plane.
        float3 probePos = gridCoordToPosition(ddgiVolume, probeGridCoord);

        float weight = 1.0;
        // Clamp all of the multiplies. We can't let the weight go to zero because then it would be 
        // possible for *all* weights to be equally low and get normalized
        // up to 1/n. We want to distinguish between weights that are 
        // low because of different factors.

        if (! debugDisableBackface) {
            // Computed without the biasing applied to the "dir" variable. 
            // This test can cause reflection-map looking errors in the image
            // (stuff looks shiny) if the transition is poor.
            float3 trueDirectionToProbe = normalize(probePos - wsPosition);

            // The naive soft backface weight would ignore a probe when
            // it is behind the surface. That's good for walls. But for small details inside of a
            // room, the normals on the details might rule out all of the probes that have mutual
            // visibility to the point. So, we instead use a "wrap shading" test below inspired by
            // NPR work.

            // The small offset at the end reduces the "going to zero" impact
            // where this is really close to exactly opposite
#if SHOW_CHEBYSHEV_WEIGHTS == 0
            weight *= square((dot(trueDirectionToProbe, sampleDirection) + 1.0) * 0.5) + 0.2;
#endif
        }
        
        // Bias the position at which visibility is computed; this avoids performing a shadow 
        // test *at* a surface, which is a dangerous location because that is exactly the line
        // between shadowed and unshadowed. If the normal bias is too small, there will be
        // light and dark leaks. If it is too large, then samples can pass through thin occluders to
        // the other side (this can only happen if there are MULTIPLE occluders near each other, a wall surface
        // won't pass through itself.)
        float3 probeToBiasedPointDirection = (wsPosition + offsetPos) - probePos;
        float distanceToBiasedPoint = length(probeToBiasedPointDirection);
        probeToBiasedPointDirection *= 1.0 / distanceToBiasedPoint;
                 
        if (! debugDisableChebyshev) {
            float2 texCoord = probeTextureCoordFromDirection(probeToBiasedPointDirection, probeGridCoord, false, ddgiVolume);

            float2 temp = Load(visibilityTexture, texCoord, visibilityTextureSize).xy;
            float meanDistanceToOccluder = temp.x;
            float variance = abs(square(temp.x) - temp.y);

            variance += ddgiVolume.debugVarianceBias;
            meanDistanceToOccluder += ddgiVolume.debugMeanBias;

            float chebyshevWeight = 1.0;
            if (distanceToBiasedPoint > meanDistanceToOccluder) {
                // In "shadow"

                // http://www.punkuser.net/vsm/vsm_paper.pdf; equation 5
                // Need the max in the denominator because biasing can cause a negative displacement
                chebyshevWeight = variance / (variance + square(distanceToBiasedPoint - meanDistanceToOccluder));
                
                // Increase contrast in the weight
                chebyshevWeight = max(pow3(chebyshevWeight) - ddgiVolume.debugChebyshevBias, 0.0) * ddgiVolume.debugChebyshevNormalize;

            }


            // Avoid visibility weights ever going all of the way to zero because when *no* probe has
            // visibility we need some fallback value.
            chebyshevWeight = max(0.05, chebyshevWeight);
            weight *= chebyshevWeight;
        }

        // Avoid zero weight
        weight = max(0.000001, weight);

        float2 texCoord = probeTextureCoordFromDirection(sampleDirection, probeGridCoord, true, ddgiVolume);

        if (! debugDisableNonlinear) {
            // A tiny bit of light is really visible due to log perception, so
            // crush tiny weights but keep the curve continuous.
            const float crushThreshold = 0.2;
            if (weight < crushThreshold) {
                weight *= square(weight) * (1.0 / square(crushThreshold)); 
            }
        }
      
		weight = weight * trilinear;
               
        float3 probeIrradiance = Load(irradianceTexture, texCoord, irradianceTextureSize).rgb;

        // Decode the tone curve, but leave a gamma = 2 curve (=sqrt here) to approximate sRGB blending for the trilinear
        float comp = ddgiVolume.irradianceGamma * 0.5;
        probeIrradiance = pow(probeIrradiance, float3(comp,comp,comp));

        if (debugDisableNonlinear) {
            // Undo the gamma encoding for debugging
            probeIrradiance = square(probeIrradiance);
        }

#       if DEBUG_VISUALIZATION_MODE == DebugVisualizationMode_IRRADIANCE_PROBE_CONTRIBUTIONS
        {
            int p = gridCoordToProbeIndex(ddgiVolume, probeGridCoord);
			probeIrradiance = (p == debugProbeIndex) ? float3(1,1,1) : float3(0,0,0);
        }
#       endif
#       if SHOW_CHEBYSHEV_WEIGHTS
		{
			int p = gridCoordToProbeIndex(ddgiVolume, probeGridCoord);
			probeIrradiance = (p == debugProbeIndex) ? float3(1,1,1) : float3(0,0,0);
		}
#       endif

        irradiance += float4(weight * probeIrradiance, weight);
    }

    // Normalize by the sum of the weights
    irradiance.xyz *= 1.0 / irradiance.a;

    // Go back to linear irradiance
    if (! debugDisableNonlinear) {   
        irradiance.xyz = square(irradiance.xyz);
    }

    // Was factored out of probes
    irradiance.xyz *= 0.5 * M_PI;

    return float4(irradiance.xyz, volumeWeight);
}

#if NUM_DDGIVOLUMES > 0
float3 sampleIrradiance
   (DDGIVolume             ddgiVolumeArray[NUM_DDGIVOLUMES],
    float3                 wsPosition,
    float3                offsetPos,
    float3                sampleDirection,
	float3			       cameraPos,
    // Set to false in production 
    const bool             debugDisableBackface,
    // Set to false in production 
    const bool             debugDisableChebyshev,
    // Set to -1 in production
    const int              debugProbeIndex) {

    float4 sum = float4(0, 0, 0, 0);
    
    // Sample until we have "100%" weight covered, and then *stop* looking at lower-resolution
    // volumes because they aren't needed.
    for (int v = 0; (v < NUM_DDGIVOLUMES) && (sum.a < 1.0); ++v) {
        // Can skip if not the last volume or some other volume has already contributed
        bool skippable = (v < NUM_DDGIVOLUMES - 1) || (sum.a > 0.9f);
        float4 irradiance = sampleOneDDGIVolume(ddgiVolumeArray[v], wsPosition, offsetPos, sampleDirection, cameraPos, debugDisableBackface, debugDisableChebyshev, debugProbeIndex, skippable);

        // Visualize weights per volume:
		//irradiance.rgb = float3(0); irradiance[v] = irradiance.a;// 1.0;

        // Max contribution of other probes should be limited by how much more weight is required
        irradiance.a *= saturate(1.0 - sum.a);

        // Premultiply
        irradiance.rgb *= irradiance.a;
        sum += irradiance;
    }

    // Normalize
    return sum.rgb / max(0.001, sum.a);
}
#endif