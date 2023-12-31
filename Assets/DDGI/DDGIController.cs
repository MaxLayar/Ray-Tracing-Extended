using UnityEditor;
using UnityEngine;
using System.Collections.Generic;

public class DDGIController : MonoBehaviour
{
    public ComputeShader ddgiComputeShader;
    public RenderTexture renderTexture;

    public Bounds volume;
    public float minProbesSpacing = 1f;

    public bool isRealtimeRaytracing = false;

    public bool debugShowProbes = false;
    public Vector3 debugRayDir = Vector3.right;

    public int debugFrame = 1;

    [Header("Raytracing")]
    [Range(1,256)]
    public int numRaysPerPixel = 10;

    [Range(0,24)]
    public int maxBounceCount = 1;

    private float realProbesSpacing;

    Vector3[] probePositions;

    // Prepare Scene
    ComputeBuffer triangleBuffer;
	ComputeBuffer meshInfoBuffer;
    List<Triangle> allTriangles;
    List<MeshInfo> allMeshInfo;
    int NumMeshes;

    public void DispatchDDGIRayTracing()
    {
        if (!renderTexture)
        {
            renderTexture = new RenderTexture(256, 256, 24);
            renderTexture.enableRandomWrite = true;
            renderTexture.Create();
        }

        int threadGroupsX = renderTexture.width / 8;
        int threadGroupsY = renderTexture.height / 8;

        SetShaderParams();
        ddgiComputeShader.Dispatch(0, threadGroupsX, threadGroupsY, 1);
    }

    void SetShaderParams()
	{
        ddgiComputeShader.SetInt("MaxBounceCount", maxBounceCount);
        ddgiComputeShader.SetInt("Frame", debugFrame);
        ddgiComputeShader.SetVector("DebugRayDir", debugRayDir);
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

        int probesInX = Mathf.CeilToInt(volume.size.x / minProbesSpacing);
        int probesInY = Mathf.CeilToInt(volume.size.y / minProbesSpacing);
        int probesInZ = Mathf.CeilToInt(volume.size.z / minProbesSpacing);

        float spacingX = volume.size.x / (probesInX-1);
        float spacingY = volume.size.y / (probesInY-1);
        float spacingZ = volume.size.z / (probesInZ-1);

        probePositions = new Vector3[probesInX * probesInY * probesInZ];

        int idx = 0;
        for (int x = 0; x < probesInX; x++)
            for (int y = 0; y < probesInY; y++)
                for (int z = 0; z < probesInZ; z++)
                    probePositions[idx++] = volume.center - volume.extents + new Vector3(x * spacingX, y * spacingY, z * spacingZ);

        realProbesSpacing = minProbesSpacing;
    }

    private void OnDrawGizmosSelected()
    {
        if (debugShowProbes)
        {
            Gizmos.color = Color.yellow;
            foreach (var probe in probePositions)
            {
                Gizmos.DrawSphere(probe, GizmoUtility.iconSize);
            }
        }
        Gizmos.color = Color.white;
        Gizmos.DrawRay(Vector3.zero, debugRayDir*10f);
    }
}