using UnityEditor;
using UnityEngine;
using System.Collections.Generic;
using System;
using System.Linq;

[ExecuteInEditMode]
public class DDGIController : MonoBehaviour
{
    public ComputeShader ddgiComputeShader;
    public RenderTexture renderTexture;

    public Bounds volume;
    public float minProbesSpacing = 1f;

    public bool isRealtimeRaytracing = false;

    public bool debugShowProbes = false;

    public int bufferDimension = 6;

    int debugFrame = 1;

    [Header("Raytracing")]
    [Range(1,256)]
    public int numRaysPerPixel = 10;

    [Range(0,24)]
    public int maxBounceCount = 1;

    private float realProbesSpacing;

    Vector3[] probesPositions;

    Vector3Int numberOfProbes;

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
            DispatchDDGIRayTracing();
            debugFrame++;
        }
    }

    public void DispatchDDGIRayTracing()
    {
        if (!renderTexture)
            RecomputeBuffers();

        int threadGroupsX = renderTexture.width / 8;
        int threadGroupsY = renderTexture.height / 8;

        SetShaderParams();
        ddgiComputeShader.Dispatch(0, threadGroupsX, threadGroupsY, 2);
    }

    void SetShaderParams()
	{
		CreateStructuredBuffer(ref probesPositionsBuffer, probesPositions.ToList());
        ddgiComputeShader.SetBuffer(0, "ProbesPositions", probesPositionsBuffer);

        ddgiComputeShader.SetVector("NumProbes", new Vector4(numberOfProbes.x, numberOfProbes.z, numberOfProbes.y, 0));
        ddgiComputeShader.SetInt("BufferDimension", bufferDimension);
        ddgiComputeShader.SetInt("MaxBounceCount", maxBounceCount);
        ddgiComputeShader.SetInt("Frame", debugFrame);
        ddgiComputeShader.SetTexture(0, "Result", renderTexture);
        ddgiComputeShader.SetInt("NumRaysPerPixel", numRaysPerPixel);
	}

    public void PrepareScene()
    {
        CreateMeshes();
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
		ddgiComputeShader.SetBuffer(0, "Triangles", triangleBuffer);
		ddgiComputeShader.SetBuffer(0, "AllMeshInfo", meshInfoBuffer);
		ddgiComputeShader.SetInt("NumMeshes", allMeshInfo.Count);
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
            for (int y = 0; y < numberOfProbes.y; y++)
                for (int z = 0; z < numberOfProbes.z; z++)
                    probesPositions[idx++] = volume.center - volume.extents + new Vector3(x * spacing.x, y * spacing.y, z * spacing.z);

        realProbesSpacing = minProbesSpacing;

        RecomputeBuffers();
    }

    void RecomputeBuffers()
    {
        if (renderTexture)
            renderTexture.Release();

        renderTexture = new RenderTexture(bufferDimension * numberOfProbes.x * numberOfProbes.y, bufferDimension * numberOfProbes.z, 16);
        renderTexture.enableRandomWrite = true;
        renderTexture.Create();
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