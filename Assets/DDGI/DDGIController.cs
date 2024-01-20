using System.Collections.Generic;
using System.Linq;

using UnityEditor;
using UnityEngine;

[ExecuteAlways, ImageEffectAllowedInSceneView]
public class DDGIController : MonoBehaviour
{
	[Header("     Probes"), Space(5)]
	public Bounds volume;
    public float minProbesSpacing = 1f;

	Vector3Int numberOfProbes;
    Vector3[] probesPositions;
    float realProbesSpacing;
	int bufferDimension = 6;

	ComputeBuffer probesPositionsBuffer;

	LightField L;
	struct LightField
	{
		public Vector3Int probeCounts;
		public int raysPerProbe;
		public Vector3 probeStartPosition;
		public float normalBias;
		public Vector3 probeStep;
		public int irradianceTextureWidth;
		public int irradianceTextureHeight;
		public int irradianceProbeSideLength;
		public int depthTextureWidth;
		public int depthTextureHeight;
		public int depthProbeSideLength;
		public float chebBias, minRayDst, energyConservation;
	};

	[Space(20), Header("     RayTracing"), Space(5)]
	public ComputeShader computeRays;
    [Range(1,300)]
    public int numRaysPerProbe = 10;

	//Buffers
	ComputeBuffer lightFieldBuffer;
	RenderTexture rayHitLocationsBuffer;
    RenderTexture rayHitRadianceBuffer;
    RenderTexture rayHitNormalsBuffer;
    RenderTexture rayDirectionsBuffer;
    RenderTexture rayOriginsBuffer;

	int randomSeed = 1;
	Matrix4x4 randomOrientationMatrix;


	[Space(20), Header("     Irradiance Settings"), Space(5)]
	public ComputeShader computeIrradiance;
	[Range(0f, 1f)]
	public float energyConservation = 0.85f;
	[Range(0f, 1f)]
	public float hysteresis = 0.965f;

	public float depthSharpness = 50.0f;

	[SerializeField] RenderTexture irradianceTex;
	[SerializeField] RenderTexture weightTex;

	ComputeBuffer uniformsBuffer;
	Values V;
	struct Values
	{
		public float depthSharpness;
		public float hysteresis;
		public float maxDistance;
	}


	[Space(20), Header("     Scene Objects"), Space(5)]
	ComputeBuffer triangleBuffer;
	ComputeBuffer meshInfoBuffer;
    List<Triangle> allTriangles;
    List<MeshInfo> allMeshInfo;
	public Light sun;


	[Space(20), Header("     DEBUG"), Space(5)]

	public bool debugShowProbes = false;
	public bool isRealtimeRaytracing = false;
	public bool isRandomDirection = true;
	[Header("DANGEROUS")]
	public bool showHitPoints = false;
	Texture2D hitPointsDebugTexture;
	Texture2D hitNormalsDebugTexture;
	Texture2D hitRadianceDebugTexture;


	void Update()
	{
		if (isRealtimeRaytracing)
		{
			PrepareScene();
			ComputeProbesRays();
			UpdateProbes(true);
			UpdateProbes(false);
			randomSeed++;
		}
	}

	public void ComputeProbesRays()
	{
		(Vector3 minScene, Vector3 maxScene) = (volume.min, volume.max);
		L = new();
		L.probeCounts = numberOfProbes;
		L.depthProbeSideLength = 14;
		L.irradianceProbeSideLength = 4;
		L.normalBias = 0.10f;
		L.minRayDst = 0.00f;
		L.irradianceTextureWidth = (L.irradianceProbeSideLength + 2) /* 1px Border around probe left and right */ * L.probeCounts.x * L.probeCounts.y + 2 /* 1px Border around whole texture left and right*/;
		L.irradianceTextureHeight = (L.irradianceProbeSideLength + 2) * L.probeCounts.z + 2;
		L.depthTextureWidth = (L.depthProbeSideLength + 2) * L.probeCounts.x * L.probeCounts.y + 2;
		L.depthTextureHeight = (L.depthProbeSideLength + 2) * L.probeCounts.z + 2;
		L.probeStartPosition = minScene - Vector3.one;
		L.probeStep = DivideVectors((maxScene - minScene + (Vector3.one * 1.3f)),(L.probeCounts - Vector3.one));
		L.raysPerProbe = numRaysPerProbe;
		L.energyConservation = energyConservation;

		CreateStructuredBuffer(ref lightFieldBuffer, new List<LightField>() { L });
		computeRays.SetBuffer(0, "LBuffer", lightFieldBuffer);

		int surfelWidth = numRaysPerProbe;
	    int surfelHeight = numberOfProbes.x * numberOfProbes.y * numberOfProbes.z;

		RefreshBufferIfNeeded(ref rayHitLocationsBuffer, "rayHitLocations", surfelWidth, surfelHeight);
		RefreshBufferIfNeeded(ref rayHitRadianceBuffer, "rayHitRadiance", surfelWidth, surfelHeight);
		RefreshBufferIfNeeded(ref rayHitNormalsBuffer, "rayHitNormals", surfelWidth, surfelHeight);
		RefreshBufferIfNeeded(ref rayDirectionsBuffer, "rayDirections", surfelWidth, surfelHeight);
		RefreshBufferIfNeeded(ref rayOriginsBuffer, "rayOrigins", surfelWidth, surfelHeight);

		RefreshBufferIfNeeded(ref irradianceTex, "irradianceTex", L.irradianceTextureWidth, L.irradianceTextureHeight);
		RefreshBufferIfNeeded(ref weightTex, "weightTex", L.depthTextureWidth, L.depthTextureHeight);

		// Set the buffers to your compute shader material
		computeRays.SetTexture(0, "rayHitLocations", rayHitLocationsBuffer);
        computeRays.SetTexture(0, "rayHitRadiance", rayHitRadianceBuffer);
        computeRays.SetTexture(0, "rayHitNormals", rayHitNormalsBuffer);
        computeRays.SetTexture(0, "rayDirections", rayDirectionsBuffer);
        computeRays.SetTexture(0, "rayOrigins", rayOriginsBuffer);

		computeRays.SetTexture(0, "irradianceTex", irradianceTex);
		computeRays.SetTexture(0, "weightTex", weightTex);

		CreateStructuredBuffer(ref probesPositionsBuffer, probesPositions.ToList());
        computeRays.SetBuffer(0, "ProbesPositions", probesPositionsBuffer);

        computeRays.SetVector("NumProbes", new Vector4(numberOfProbes.x, numberOfProbes.z, numberOfProbes.y, 0));
        computeRays.SetInt("BufferDimension", bufferDimension);
        computeRays.SetInt("Frame", randomSeed);
        computeRays.SetInt("NumRaysPerProbe", numRaysPerProbe);

		if (isRandomDirection)
		{
			// Get a random direction in spherical coordinates
			float theta = Mathf.Acos(2f * Random.value - 1f);
			float phi = 2f * Mathf.PI * Random.value;

			Vector3 axis = new Vector3(Mathf.Sin(theta) * Mathf.Cos(phi), Mathf.Sin(theta) * Mathf.Sin(phi), Mathf.Cos(theta));

			randomOrientationMatrix = Matrix4x4.TRS(Vector3.zero, Quaternion.AngleAxis(Random.value * 360f, axis), Vector3.one);
		}
		
        computeRays.SetMatrix("randomOrientation", randomOrientationMatrix);

        // Dispatch your compute shader
        int threadGroupsX = Mathf.CeilToInt(surfelWidth / 8.0f);
        int threadGroupsY = Mathf.CeilToInt(surfelHeight / 8.0f);
        computeRays.Dispatch(0, threadGroupsX, threadGroupsY, 1);
	}


	public void UpdateProbes(bool isOutputIrradiance)
	{
		computeIrradiance.SetFloat("PROBE_SIDE_LENGTH", isOutputIrradiance ? L.irradianceProbeSideLength : L.depthProbeSideLength);
		computeIrradiance.SetInt("RAYS_PER_PROBE", numRaysPerProbe);

		if (isOutputIrradiance)
		{
			computeIrradiance.EnableKeyword("OUTPUT_IRRADIANCE");
			computeIrradiance.SetTexture(0, "rayHitRadiance", rayHitRadianceBuffer);
		}
		else
		{
			computeIrradiance.DisableKeyword("OUTPUT_IRRADIANCE");
			computeIrradiance.SetTexture(0, "rayHitLocations", rayHitLocationsBuffer);
			computeIrradiance.SetTexture(0, "rayOrigins", rayOriginsBuffer);
			computeIrradiance.SetTexture(0, "rayHitNormals", rayHitNormalsBuffer);
		}

		computeIrradiance.SetTexture(0, "rayDirections", rayDirectionsBuffer);

		V = new();
		V.depthSharpness = depthSharpness;
		V.hysteresis = hysteresis;
		Vector3 probeEnd = L.probeStartPosition + MultiplyVectors(new Vector3(L.probeCounts.x - 1, L.probeCounts.y - 1, L.probeCounts.z -1), L.probeStep);
		Vector3 probeSpan = probeEnd - L.probeStartPosition;
		V.maxDistance = Vector3.Magnitude(DivideVectors(probeSpan, L.probeCounts)) * 1.5f;

		CreateStructuredBuffer(ref uniformsBuffer, new List<Values>() { V });
		computeIrradiance.SetBuffer(0, "uniforms", uniformsBuffer);
		computeIrradiance.SetTexture(0, "tex", isOutputIrradiance ? irradianceTex : weightTex);

		// Dispatch your compute shader
		int threadGroupsX = Mathf.CeilToInt((isOutputIrradiance ? L.irradianceTextureWidth : L.depthTextureWidth) / 8f);
		int threadGroupsY = Mathf.CeilToInt((isOutputIrradiance ? L.irradianceTextureHeight : L.depthTextureHeight) / 8f);
		computeIrradiance.Dispatch(0, threadGroupsX, threadGroupsY, 1);
	}

	private void RefreshBufferIfNeeded(ref RenderTexture buffer, string bufferName, int width, int heigh)
	{
		if (!buffer)
		{
			buffer = new RenderTexture(width, heigh, 16, UnityEngine.Experimental.Rendering.DefaultFormat.HDR);
			buffer.enableRandomWrite = true;
			buffer.name = bufferName;
		}

		if (buffer.width != width || buffer.height != heigh)
		{
			buffer.Release();
			buffer.width = width;
			buffer.height = heigh;
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
		if (sun && sun.isActiveAndEnabled)
		{
			computeRays.SetVector("sunColor", sun.color);
			computeRays.SetVector("sunDirection", sun.transform.forward);
		}
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

		numRaysPerPixel = Mathf.Max(1, numRaysPerPixel);
    }

    public void RefreshProbesPlacement()
    {
        UnityEngine.Debug.Log($"Refresh probes placement with a distance of {minProbesSpacing}.");

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
			int probeID = 0;
            foreach (var probe in probesPositions)
			{
				Gizmos.color = GetProbeColor(probeID++);
				Gizmos.DrawSphere(probe, GizmoUtility.iconSize);

				Gizmos.color = Color.red;
				Gizmos.DrawRay(probe, randomOrientationMatrix.GetColumn(2));
            }
        }

		if (showHitPoints)
		{
			hitPointsDebugTexture = GetRTPixels(rayHitLocationsBuffer);
			hitNormalsDebugTexture = GetRTPixels(rayHitNormalsBuffer);
			hitRadianceDebugTexture = GetRTPixels(rayHitRadianceBuffer);

			for (int y = 0; y < rayHitLocationsBuffer.height; y++)
			{
				for (int x = 0; x < rayHitLocationsBuffer.width; x++)
				{
					//Gizmos.color = GetProbeColor(y);
					Gizmos.color = hitRadianceDebugTexture.GetPixel(x, y);

					Color positionColor = hitPointsDebugTexture.GetPixel(x, y);
					Vector3 position = new Vector3(positionColor.r, positionColor.g, positionColor.b);
					Gizmos.DrawSphere(position, GizmoUtility.iconSize / 4);

					Color normalColor = hitNormalsDebugTexture.GetPixel(x, y);
					Vector3 normal = new Vector3(normalColor.r, normalColor.g, normalColor.b);
					Gizmos.DrawRay(position, DivideVectors(normal, new Vector3(10,10,10)));

					Gizmos.color *= new Color(1,1,1,.1f);
					Gizmos.DrawLine(position, probesPositions[y]);
				}
			}
		}
	}

	static public Texture2D GetRTPixels(RenderTexture rt)
	{
		// Remember currently active render texture
		RenderTexture currentActiveRT = RenderTexture.active;

		// Set the supplied RenderTexture as the active one
		RenderTexture.active = rt;

		// Create a new Texture2D and read the RenderTexture image into it
		Texture2D tex = new Texture2D(rt.width, rt.height, UnityEngine.Experimental.Rendering.DefaultFormat.HDR, UnityEngine.Experimental.Rendering.TextureCreationFlags.DontInitializePixels);
		tex.ReadPixels(new Rect(0, 0, tex.width, tex.height), 0, 0);
		tex.Apply();

		// Restorie previously active render texture
		RenderTexture.active = currentActiveRT;
		return tex;
	}

	private Color GetProbeColor(int probeID)
	{
		Random.State state = Random.state;

		Random.InitState(probeID);
		float r = Random.value;
		Random.InitState(probeID + 1);
		float g = Random.value;
		Random.InitState(probeID + 2);
		float b = Random.value;

		Random.state = state;

		return new Color(r, g, b);
	}

	private static Vector3 MultiplyVectors(Vector3 a, Vector3 b)
	{
		return new Vector3(a.x * b.x, a.y * b.y, a.z * b.z);
	}

	private static Vector3 DivideVectors(Vector3 a, Vector3 b)
	{
		return new Vector3(a.x / b.x, a.y / b.y, a.z / b.z);
	}

	public void ResetBuffers()
	{
		ClearOutRenderTexture(weightTex);
		ClearOutRenderTexture(irradianceTex);

		ClearOutRenderTexture(rayHitLocationsBuffer);
		ClearOutRenderTexture(rayHitRadianceBuffer);
		ClearOutRenderTexture(rayHitNormalsBuffer);
		ClearOutRenderTexture(rayDirectionsBuffer);
		ClearOutRenderTexture(rayOriginsBuffer);
	}

	public void ClearOutRenderTexture(RenderTexture renderTexture)
	{
		RenderTexture rt = RenderTexture.active;
		RenderTexture.active = renderTexture;
		GL.Clear(true, true, Color.clear);
		RenderTexture.active = rt;
	}



	   [Header("Ray Tracing Settings")]
	[SerializeField, Range(0, 64)] int numRaysPerPixel = 2;


	[Header("View Settings")]
	[SerializeField] bool useShaderInSceneView;

	[Header("References")]
	[SerializeField] Shader rayTracingShader;

	[Header("Info")]
	[SerializeField] int numRenderedFrames;



	// Materials and render textures
	Material rayTracingMaterial;
	RenderTexture resultTexture;


	void Start()
	{
		numRenderedFrames = 0;
	}

    void OnRenderImage(RenderTexture src, RenderTexture target)
    {
        bool isSceneCam = Camera.current.name == "SceneCamera";

        if (isSceneCam)
        {
            if (useShaderInSceneView)
            {
                InitFrame();
                Graphics.Blit(src, target, rayTracingMaterial);
            }
            else
            {
                Graphics.Blit(src, target); // Draw the unaltered camera render to the screen
            }
        }
        else
        {
            InitFrame();

            // Create copy of prev frame
            //RenderTexture prevFrameCopy = RenderTexture.GetTemporary(src.width, src.height, 0, ShaderHelper.RGBA_SFloat);
            //Graphics.Blit(resultTexture, prevFrameCopy);


            Graphics.Blit(src, target, rayTracingMaterial);

            // Draw result to screen
            //Graphics.Blit(resultTexture, target);

            // Release temps
            //RenderTexture.ReleaseTemporary(prevFrameCopy);

            numRenderedFrames += Application.isPlaying ? 1 : 0;
        }
    }


	void InitFrame()
	{
		// Create materials used in blits
		ShaderHelper.InitMaterial(rayTracingShader, ref rayTracingMaterial);
		
		Camera.current.depthTextureMode = DepthTextureMode.Depth;

		// Run the ray tracing shader and draw the result to a temp texture
		rayTracingMaterial.SetInt("Frame", numRenderedFrames);

		rayTracingMaterial.SetBuffer("LBuffer", lightFieldBuffer);

		rayTracingMaterial.SetTexture("irradianceTex", irradianceTex);
		rayTracingMaterial.SetTexture("weightTex", weightTex);

		rayTracingMaterial.SetMatrix("_InverseView", Camera.current.cameraToWorldMatrix);

		// Create result render texture
		ShaderHelper.CreateRenderTexture(ref resultTexture, Screen.width, Screen.height, FilterMode.Bilinear, ShaderHelper.RGBA_SFloat, "Result");
	}


	void OnDisable()
	{
		ShaderHelper.Release(resultTexture);
	}
}



/* EXECUTION ORDER : 
 * 
 * 
 * 1 Write borders with 1s (in irr & weight buffers)
 * 2 Ray pass (STEP 1 from Morgan's explaination)
 * 3 Update weight probe pass (STEP 2)
 * 3.5 copy border weight pass
 * 4 Update irradiance probe pass
 * 4.5 copy border irr pass
 * 5 samples the probes per pixel, per frame (pixel shader?) (STEP 3)
 * 
 * 
 */