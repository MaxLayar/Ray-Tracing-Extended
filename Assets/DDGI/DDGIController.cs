using UnityEditor;
using UnityEngine;
using System.Collections.Generic;
using System.Linq;

[ExecuteInEditMode]
public class DDGIController : MonoBehaviour
{
    public ComputeShader computeRays;

    public Bounds volume;
    public float minProbesSpacing = 1f;

    public bool isRealtimeRaytracing = false;

    public bool debugShowProbes = false;

    public int bufferDimension = 6;

    int debugFrame = 1;

    [Header("Raytracing")]
    [Range(1,256)]
    public int numRaysPerProbe = 10;

    [Range(0,24)]
    public int maxBounceCount = 1;

    private float realProbesSpacing;

    Vector3[] probesPositions;

    Vector3Int numberOfProbes;

    public Light sun;


    //Buffers
    RenderTexture rayHitLocationsBuffer;
    RenderTexture rayHitRadianceBuffer;
    RenderTexture rayHitNormalsBuffer;
    RenderTexture rayDirectionsBuffer;
    RenderTexture rayOriginsBuffer;



    // Prepare Scene
    ComputeBuffer probesPositionsBuffer;
    ComputeBuffer triangleBuffer;
	ComputeBuffer meshInfoBuffer;
    List<Triangle> allTriangles;
    List<MeshInfo> allMeshInfo;
    int NumMeshes;

    void Update()
    {
        if (isRealtimeRaytracing)
        {
            PrepareScene();
            FirstPass();
            //DispatchDDGIRayTracing();
            debugFrame++;
        }
    }

    public void FirstPass()
    {
	    int surfelWidth = numRaysPerProbe;
	    int surfelHeight = numberOfProbes.x * numberOfProbes.y * numberOfProbes.z;

		RefreshBufferIfNeeded(ref rayHitLocationsBuffer, "rayHitLocations", surfelWidth, surfelHeight);
		RefreshBufferIfNeeded(ref rayHitRadianceBuffer, "rayHitRadiance", surfelWidth, surfelHeight);
		RefreshBufferIfNeeded(ref rayHitNormalsBuffer, "rayHitNormals", surfelWidth, surfelHeight);
		RefreshBufferIfNeeded(ref rayDirectionsBuffer, "rayDirections", surfelWidth, surfelHeight);
		RefreshBufferIfNeeded(ref rayOriginsBuffer, "rayOrigins", surfelWidth, surfelHeight);

		// Set the buffers to your compute shader material
		computeRays.SetTexture(0, "rayHitLocations", rayHitLocationsBuffer);
        computeRays.SetTexture(0, "rayHitRadiance", rayHitRadianceBuffer);
        computeRays.SetTexture(0, "rayHitNormals", rayHitNormalsBuffer);
        computeRays.SetTexture(0, "rayDirections", rayDirectionsBuffer);
        computeRays.SetTexture(0, "rayOrigins", rayOriginsBuffer);

        CreateStructuredBuffer(ref probesPositionsBuffer, probesPositions.ToList());
        computeRays.SetBuffer(0, "ProbesPositions", probesPositionsBuffer);

        computeRays.SetVector("NumProbes", new Vector4(numberOfProbes.x, numberOfProbes.z, numberOfProbes.y, 0));
        computeRays.SetInt("BufferDimension", bufferDimension);
        computeRays.SetInt("MaxBounceCount", maxBounceCount);
        computeRays.SetInt("Frame", debugFrame);
        computeRays.SetInt("NumRaysPerProbe", numRaysPerProbe);

        // Get a random direction in spherical coordinates
        float theta = Mathf.Acos(2f * Random.value - 1f);
        float phi = 2f * Mathf.PI * Random.value;

        Vector3 axis = new Vector3(Mathf.Sin(theta) * Mathf.Cos(phi), Mathf.Sin(theta) * Mathf.Sin(phi), Mathf.Cos(theta));

        Matrix4x4 randomOrientationMatrix = Matrix4x4.TRS(Vector3.zero, Quaternion.AngleAxis(Random.value * 360f, axis), Vector3.one);

        computeRays.SetMatrix("randomOrientation", randomOrientationMatrix);

        // Dispatch your compute shader
        int threadGroupsX = Mathf.CeilToInt(surfelWidth / 8.0f);
        int threadGroupsY = Mathf.CeilToInt(surfelHeight / 8.0f);
        computeRays.Dispatch(0, threadGroupsX, threadGroupsY, 1);
    }

	private void RefreshBufferIfNeeded(ref RenderTexture buffer, string bufferName, int width, int heigh)
	{
		if (!buffer || buffer.width != width || buffer.height != heigh)
		{
			buffer?.Release();

			buffer = new RenderTexture(width, heigh, 16);
			buffer.enableRandomWrite = true;
			buffer.name = bufferName;
		}
	}

	public void PrepareScene()
    {
        CreateMeshes();
        CreateLights();
    }

	void CreateMeshes()
	{
		RayTracedMesh[] meshObjects = FindObjectsByType<RayTracedMesh>(FindObjectsInactive.Exclude, FindObjectsSortMode.None);

		allTriangles ??= new List<Triangle>();
		allMeshInfo ??= new List<MeshInfo>();
		allTriangles.Clear();
		allMeshInfo.Clear();

		for (int i = 0; i < meshObjects.Length; i++)
		{
			MeshChunk[] chunks = meshObjects[i].GetSubMeshes();
			foreach (MeshChunk chunk in chunks)
			{
				RayTracingMaterial material = meshObjects[i].GetMaterial(chunk.subMeshIndex);
				allMeshInfo.Add(new MeshInfo(allTriangles.Count, chunk.triangles.Length, material, chunk.bounds));
				allTriangles.AddRange(chunk.triangles);
			}
		}

		CreateStructuredBuffer(ref triangleBuffer, allTriangles);
		CreateStructuredBuffer(ref meshInfoBuffer, allMeshInfo);
		computeRays.SetBuffer(0, "Triangles", triangleBuffer);
		computeRays.SetBuffer(0, "AllMeshInfo", meshInfoBuffer);
		computeRays.SetInt("NumMeshes", allMeshInfo.Count);
	}

    void CreateLights()
    {
        computeRays.SetVector("sunColor", sun.isActiveAndEnabled ? sun.color : Color.black);
        computeRays.SetVector("sunDirection", sun.transform.forward);
    }

	// Create a compute buffer containing the given data (Note: data must be blittable)
	public static void CreateStructuredBuffer<T>(ref ComputeBuffer buffer, List<T> data) where T : struct
	{
		// Cannot create 0 length buffer (not sure why?)
		int length = Mathf.Max(1, data.Count);
		// The size (in bytes) of the given data type
		int stride = System.Runtime.InteropServices.Marshal.SizeOf(typeof(T));

		// If buffer is null, wrong size, etc., then we'll need to create a new one
		if (buffer == null || !buffer.IsValid() || buffer.count != length || buffer.stride != stride)
		{
			if (buffer != null) { buffer.Release(); }
			buffer = new ComputeBuffer(length, stride, ComputeBufferType.Structured);
		}

		buffer.SetData(data);
	}
    public static int GetStride<T>() => System.Runtime.InteropServices.Marshal.SizeOf(typeof(T));


    private void OnValidate()
    {
        if (minProbesSpacing != realProbesSpacing)
        {
            RefreshProbesPlacement();
        }
    }

    public void RefreshProbesPlacement()
    {
        Debug.Log($"Refresh probes placement with a distance of {minProbesSpacing}.");

        numberOfProbes = new(
            Mathf.CeilToInt(volume.size.x / minProbesSpacing),
            Mathf.CeilToInt(volume.size.y / minProbesSpacing),
            Mathf.CeilToInt(volume.size.z / minProbesSpacing)
        );

        Vector3 spacing = new(
            volume.size.x / (numberOfProbes.x - 1),
            volume.size.y / (numberOfProbes.y - 1),
            volume.size.z / (numberOfProbes.z - 1)
        );

        spacing.x = (float.IsInfinity(spacing.x) || float.IsNaN(spacing.x)) ? 0 : spacing.x;
        spacing.y = (float.IsInfinity(spacing.y) || float.IsNaN(spacing.y)) ? 0 : spacing.y;
        spacing.z = (float.IsInfinity(spacing.z) || float.IsNaN(spacing.z)) ? 0 : spacing.z;

        probesPositions = new Vector3[numberOfProbes.x * numberOfProbes.y * numberOfProbes.z];

        int idx = 0;
		for (int x = 0; x < numberOfProbes.x; x++)
		{
			for (int y = 0; y < numberOfProbes.y; y++)
			{
				for (int z = 0; z < numberOfProbes.z; z++)
				{
					probesPositions[idx++] = volume.center - volume.extents + new Vector3(x * spacing.x, y * spacing.y, z * spacing.z);
				}
			}
		}

        realProbesSpacing = minProbesSpacing;
    }

    private void OnDrawGizmos()
    {
        if (debugShowProbes)
        {
            Gizmos.color = Color.yellow;
            foreach (var probe in probesPositions)
            {
                Gizmos.DrawSphere(probe, GizmoUtility.iconSize);
            }
        }
    }
}