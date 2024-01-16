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
    hitInfo.dst = dst;
    return hitInfo;
}

// Find the first point that the given ray collides with, and return hit info
bool CalculateRayCollision(in Ray ray, out HitInfo info)
{
    HitInfo closestHit = (HitInfo)0;
    // We haven't hit anything yet, so 'closest' hit is infinitely far away
    closestHit.dst = 1.#INF;

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

bool TraceRay(in Ray ray, out HitInfo info)
{
    return CalculateRayCollision(ray, info);
}

bool TraceRaySimple(in Ray ray)
{
    HitInfo info;
    return CalculateRayCollision(ray, info);
}