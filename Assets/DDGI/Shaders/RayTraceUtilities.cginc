struct Ray
{
    float3 origin;
    float3 direction;
    float tMin;
    float tMax;
};

struct RayTracingMaterial
{
    float4 colour;
    float4 emissionColour;
    float4 specularColour;
    float emissionStrength;
    float smoothness;
    float specularProbability;
    int flag;
};

struct HitInfo
{
    bool didHit;
    float dst;
    float3 hitPoint;
    float3 normal;
    RayTracingMaterial material;
    bool isProbe;
};

struct Triangle
{
    float3 posA, posB, posC;
    float3 normalA, normalB, normalC;
};

struct MeshInfo
{
    uint firstTriangleIndex;
    uint numTriangles;
    RayTracingMaterial material;
    float3 boundsMin;
    float3 boundsMax;
};

struct Sun
{
    float3 direction;
    float padding;
    float3 color;
    float padding2;
};


//Raytracing settings
uint NumRaysPerProbe;

// Mesh Scene
StructuredBuffer<Triangle> Triangles;
StructuredBuffer<MeshInfo> AllMeshInfo;
int NumMeshes;


// Thanks to https://gist.github.com/DomNomNom/46bb1ce47f68d255fd5d
bool RayBoundingBox(Ray ray, float3 boxMin, float3 boxMax)
{
    float3 invDir = 1 / ray.direction;
    float3 tMin = (boxMin - ray.origin) * invDir;
    float3 tMax = (boxMax - ray.origin) * invDir;
    float3 t1 = min(tMin, tMax);
    float3 t2 = max(tMin, tMax);
    float tNear = max(max(t1.x, t1.y), t1.z);
    float tFar = min(min(t2.x, t2.y), t2.z);
    return tNear <= tFar;
};

// Calculate the intersection of a ray with a triangle using Möller–Trumbore algorithm
// Thanks to https://stackoverflow.com/a/42752998
HitInfo RayTriangle(Ray ray, Triangle tri)
{
    float3 edgeAB = tri.posB - tri.posA;
    float3 edgeAC = tri.posC - tri.posA;
    float3 normalVector = cross(edgeAB, edgeAC);
    float3 ao = ray.origin - tri.posA;
    float3 dao = cross(ao, ray.direction);

    float determinant = -dot(ray.direction, normalVector);
    float invDet = 1 / determinant;
    
    // Calculate dst to triangle & barycentric coordinates of intersection point
    float dst = dot(ao, normalVector) * invDet;
    float u = dot(edgeAC, dao) * invDet;
    float v = -dot(edgeAB, dao) * invDet;
    float w = 1 - u - v;
    
    // Initialize hit info
    HitInfo hitInfo;
    hitInfo.didHit = determinant >= 1E-6 && dst >= 0 && u >= 0 && v >= 0 && w >= 0;
    hitInfo.hitPoint = ray.origin + ray.direction * dst;
    hitInfo.normal = normalize(tri.normalA * w + tri.normalB * u + tri.normalC * v);
    hitInfo.isProbe = false;
    hitInfo.dst = dst;
    return hitInfo;
}

// Calculate the intersection of a ray with a sphere
HitInfo RaySphere(Ray ray, float3 sphereCentre, float sphereRadius)
{
    HitInfo hitInfo = (HitInfo)0;
    float3 offsetRayOrigin = ray.origin - sphereCentre;
    // From the equation: sqrLength(rayOrigin + rayDir * dst) = radius^2
    // Solving for dst results in a quadratic equation with coefficients:
    float a = dot(ray.direction, ray.direction); // a = 1 (assuming unit vector)
    float b = 2 * dot(offsetRayOrigin, ray.direction);
    float c = dot(offsetRayOrigin, offsetRayOrigin) - sphereRadius * sphereRadius;
    // Quadratic discriminant
    float discriminant = b * b - 4 * a * c; 

    // No solution when d < 0 (ray misses sphere)
    if (discriminant >= 0) {
        // Distance to nearest intersection point (from quadratic formula)
        float dst = (-b - sqrt(discriminant)) / (2 * a);

        // Ignore intersections that occur behind the ray
        if (dst >= 0) {
            hitInfo.didHit = true;
            hitInfo.dst = dst;
            hitInfo.hitPoint = ray.origin + ray.direction * dst;
            hitInfo.normal = normalize(hitInfo.hitPoint - sphereCentre);
            hitInfo.isProbe = true;
        }
    }
    return hitInfo;
}

#if defined(SHOW_PROBES)
StructuredBuffer<float3> ProbesPositions;
float DebugProbesRadius;
bool debugIrradiance;
#endif

// Find the first point that the given ray collides with, and return hit info
bool CalculateRayCollision(in Ray ray, out HitInfo info, bool detectProbes = false)
{
    HitInfo closestHit = (HitInfo)0;
    // We haven't hit anything yet, so 'closest' hit is infinitely far away
    closestHit.dst = 1.#INF;

    #if defined(SHOW_PROBES)
    if (detectProbes)
    {
        int numProbes = DDGIVolumes[0].probeCounts.x * DDGIVolumes[0].probeCounts.y * DDGIVolumes[0].probeCounts.z;
        // Raycast against all spheres and keep info about the closest hit
        for (int i = 0; i < numProbes; i ++)
        {
            HitInfo hitInfo = RaySphere(ray, ProbesPositions[i], DebugProbesRadius);

            if (hitInfo.didHit && hitInfo.dst < closestHit.dst)
            {
                closestHit = hitInfo;
            
                RayTracingMaterial m = (RayTracingMaterial) 0;

                int3 probeGridCoord = probeIndexToGridCoord(DDGIVolumes[0], i);
                float2 texCoord = probeTextureCoordFromDirection(hitInfo.normal, probeGridCoord, debugIrradiance, DDGIVolumes[0]);

                if (debugIrradiance)
                    m.colour = irradianceTexture.SampleLevel(sampler_irradianceTexture, texCoord * DDGIVolumes[0].invIrradianceTextureSize, 0);
                else
                    m.colour = visibilityTexture.SampleLevel(sampler_visibilityTexture, texCoord * DDGIVolumes[0].invVisibilityTextureSize, 0);

                closestHit.material = m;
            }
        }
    }
    #endif
    
    // Raycast against all meshes and keep info about the closest hit
    for (int meshIndex = 0; meshIndex < NumMeshes; meshIndex ++)
    {
        MeshInfo meshInfo = AllMeshInfo[meshIndex];
        if (!RayBoundingBox(ray, meshInfo.boundsMin, meshInfo.boundsMax)) {
            continue;
        }

        for (uint i = 0; i < meshInfo.numTriangles; i ++) {
            int triIndex = meshInfo.firstTriangleIndex + i;
            Triangle tri = Triangles[triIndex];
            info = RayTriangle(ray, tri);

            if (info.didHit && info.dst < closestHit.dst)
            {
                closestHit = info;
                closestHit.material = meshInfo.material;
            }
        }
    }
    info = closestHit;

    return info.dst < ray.tMax;
}

float4 ComputeShadingAt(HitInfo info, float3 viewVec, float3 sunDirection, float4 sunColor, float maxDistance, bool useIndirect, bool debugDisableBackface, bool debugDisableChebyshev)
{
    float4 indirectL = float4(sampleIrradiance(DDGIVolumes, info.hitPoint, info.normal, info.normal, _WorldSpaceCameraPos, debugDisableBackface, debugDisableChebyshev, -1), 1);

    HitInfo shadowHit;
    Ray shadowRay;
    shadowRay.direction = -sunDirection;
    shadowRay.origin = info.hitPoint + (info.normal * DDGIVolumes[0].selfShadowBias);
    shadowRay.tMin = 0.01;
    shadowRay.tMax = maxDistance;
    int lit = !CalculateRayCollision(shadowRay, shadowHit); //sHit ? 0 : 1;

    SurfaceOutputStandard s;
    s.Albedo = info.material.colour.rgb; // base (diffuse or specular) color
    s.Normal = normalize(info.normal);
    s.Emission = info.material.emissionColour * info.material.emissionStrength;
    s.Metallic = info.material.specularProbability; // 0=non-metal, 1=metal
    s.Smoothness = info.material.smoothness; // 0=rough, 1=smooth
    s.Alpha = info.material.colour.a; // alpha for transparencies

    UnityLight light;
    light.color = sunColor * lit;
    light.dir = -sunDirection;

    UnityIndirect indirect;
    indirect.diffuse = useIndirect ? indirectL.rgb : 0;
    indirect.specular = 0;

    half oneMinusReflectivity;
    half3 specColor;
    s.Albedo = DiffuseAndSpecularFromMetallic (s.Albedo, s.Metallic, /*out*/ specColor, /*out*/ oneMinusReflectivity);

    // shader relies on pre-multiply alpha-blend (_SrcBlend = One, _DstBlend = OneMinusSrcAlpha)
    // this is necessary to handle transparency in physically correct way - only diffuse component gets affected by alpha
    half outputAlpha;
    s.Albedo = PreMultiplyAlpha (s.Albedo, s.Alpha, oneMinusReflectivity, /*out*/ outputAlpha);

    half4 c = BRDF1_Unity_PBS (s.Albedo, specColor, oneMinusReflectivity, s.Smoothness, s.Normal, viewVec, light, indirect);
    c.a = outputAlpha;
    return c;
}