using System.Collections.Generic;
using System.Linq;

using UnityEditor;
using UnityEngine;

[ExecuteAlways, ImageEffectAllowedInSceneView]
public class DDGIController : MonoBehaviour
{
	[Header("     Probes"), Space(5)]
	public Bounds boundsVolume;
    public Vector3 minProbesSpacing = Vector3.one;
	Vector3 realProbeSpacing = Vector3.one;

    Vector3[] probesPositions;
	ComputeBuffer probesPositionsBuffer;
	bool isWriteOnesDone = false;

	[SerializeField] List<DDGIVolume> ddgiVolumes;

	[System.Serializable]
	struct DDGIVolume {
		public Vector3Int              probeCounts;
		public Vector3Int              logProbeCounts;
		public Vector3          	   probeGridOrigin;
		public Vector3          	   probeSpacing;
		public Vector3          	   invProbeSpacing; // 1 / probeSpacing
		public Vector3Int              phaseOffsets;
		public Vector2Int              invIrradianceTextureSize;
		public Vector2Int              invVisibilityTextureSize;
		// probeOffsetLimit on [0,0.5] where max probe 
		// offset = probeOffsetLimit * probeSpacing
		// Usually 0.4, controllable from GUI.
		[Range(0f,0.5f)]
		public float                   probeOffsetLimit;
		[Range(2,64)]
		public int                     irradianceProbeSideLength;
		[Range(2,64)]
		public int                     visibilityProbeSideLength;
		public float                   selfShadowBias;
		[Range(0.01f, 10f)]
		public float                   irradianceGamma;
		public float                   invIrradianceGamma;
		public float                   debugMeanBias; // 0 in production code
		public float                   debugVarianceBias; // 0 in production code
		public float                   debugChebyshevBias; // 0 in production code
		public float                   debugChebyshevNormalize; // 1/(1-debugChebyshevBias) = 1 in production code
		public int					   isCameraLocked;

		[Range(0, 10000)]
		public float 					depthSharpness;
		[Range(0f,1f)]
		public float 					hysteresis;
		[Range(0f,1f)]
		public float					energyConservation;
		public float 					maxDistance;
	};
	[SerializeField] RenderTexture irradianceTexture;
	[SerializeField] RenderTexture visibilityTexture;
	[SerializeField] RenderTexture probeOffsetsTexture;
	[SerializeField] RenderTexture probeOffsetsImage;

	[Space(20), Header("     RayTracing"), Space(5)]
	public ComputeShader computeRays;
    [Range(1,300)]
    public int numRaysPerProbe = 10;

	//Buffers
	ComputeBuffer ddgiVolumesBuffer;
	RenderTexture rayHitLocationsBuffer;
    RenderTexture rayHitRadianceBuffer;
    RenderTexture rayHitNormalsBuffer;
    RenderTexture rayDirectionsBuffer;
    RenderTexture rayOriginsBuffer;

	int randomSeed = 1;
	Matrix4x4 randomOrientationMatrix;


	[Space(20), Header("     Irradiance Settings"), Space(5)]
	public ComputeShader computeIrradiance;
	public ComputeShader computeBorders;

	public float chebBias = 0;

	[Space(20), Header("     Scene Objects"), Space(5)]
	ComputeBuffer triangleBuffer;
	ComputeBuffer meshInfoBuffer;
    List<Triangle> allTriangles;
    List<MeshInfo> allMeshInfo;
	public Light sun;

	[Space(20), Header("     DEBUG"), Space(5)]

	public int irradianceProbeSideLength = 6;
	public int depthProbeSideLength = 14;

	public bool debugShowProbes = false;
	public bool isRealtimeRaytracing = false;
	public bool isRandomDirection = true;
	[Header("DANGEROUS")]
	public bool showHitPoints = false;
	Texture2D hitPointsDebugTexture;
	Texture2D hitNormalsDebugTexture;
	Texture2D hitRadianceDebugTexture;

	void OnEnable()
	{
		if (ddgiVolumes == null)
			ddgiVolumes = new List<DDGIVolume>();

		if (ddgiVolumes.Count == 0)
			ddgiVolumes.Add(new DDGIVolume());
		
		for(int i = 0; i < ddgiVolumes.Count; i++)
		{
			if (!irradianceTexture)
				RefreshBufferIfNeeded(ref irradianceTexture, "Irradiance", 1, 1);
			if (!visibilityTexture)
				RefreshBufferIfNeeded(ref visibilityTexture, "Visibility", 1, 1);
			if (!probeOffsetsTexture)
				RefreshBufferIfNeeded(ref probeOffsetsTexture, "Probe Offsets Texture", 1, 1);
			if (!probeOffsetsImage)
				RefreshBufferIfNeeded(ref probeOffsetsImage, "Probe Offsets Image", 1, 1);
		}

		NotifyOfCameraPosition(Camera.current ? Camera.current.transform.position : Camera.main.transform.position);
	}

	void Update()
	{
		if (isRealtimeRaytracing)
		{
			if (!isWriteOnesDone)
			{
				foreach(DDGIVolume ddgiVolume in ddgiVolumes)
				{
					if (visibilityTexture)
						UpdateProbesBorders(visibilityTexture, ddgiVolume.visibilityProbeSideLength, 0); //set ones
					if (irradianceTexture)
						UpdateProbesBorders(irradianceTexture, ddgiVolume.irradianceProbeSideLength, 0); //set ones
				}
				isWriteOnesDone = true;
			}

			PrepareScene(computeRays);
			ComputeProbesRays();
			UpdateProbes(false); //weight
			UpdateProbes(true); //irradiance

			foreach(DDGIVolume ddgiVolume in ddgiVolumes)
			{
				if (visibilityTexture)
					UpdateProbesBorders(visibilityTexture, ddgiVolume.visibilityProbeSideLength, 1); //copy borders
				if (irradianceTexture)
					UpdateProbesBorders(irradianceTexture, ddgiVolume.irradianceProbeSideLength, 1); //copy borders
			}

			randomSeed++;
		}
	}

	public void ComputeProbesRays()
	{
		DDGIVolume copy = ddgiVolumes[0];

		CreateStructuredBuffer(ref ddgiVolumesBuffer, ddgiVolumes);
		computeRays.SetBuffer(0, "DDGIVolumes", ddgiVolumesBuffer);

		int irrW = (copy.irradianceProbeSideLength + 2) /* 1px Border around probe left and right */ * copy.probeCounts.x * copy.probeCounts.y + 2; /* 1px Border around whole texture left and right*/
		int irrH = (copy.irradianceProbeSideLength + 2) * copy.probeCounts.z + 2;
		int visW = (copy.visibilityProbeSideLength + 2) * copy.probeCounts.x * copy.probeCounts.y + 2;
		int visH = (copy.visibilityProbeSideLength + 2) * copy.probeCounts.z + 2;

		RefreshBufferIfNeeded(ref irradianceTexture, "irradianceTexture", irrW, irrH);
		RefreshBufferIfNeeded(ref visibilityTexture, "visibilityTexture", visW, visH);
		RefreshBufferIfNeeded(ref probeOffsetsTexture, "probeOffsetsTexture", visW, visH);
		RefreshBufferIfNeeded(ref probeOffsetsImage, "probeOffsetsImage", visW, visH);

		computeRays.SetTexture(0, "irradianceTexture", irradianceTexture);
		computeRays.SetTexture(0, "visibilityTexture", visibilityTexture);
		computeRays.SetTexture(0, "probeOffsetsTexture", probeOffsetsTexture);
		computeRays.SetTexture(0, "probeOffsetsImage", probeOffsetsImage);

		computeRays.SetInt("OFFSET_BITS_PER_CHANNEL", 8); //probeOffsetsTexture->format()->redBits


		int surfelWidth = numRaysPerProbe;
	    int surfelHeight = Mathf.Max(copy.probeCounts.x * copy.probeCounts.y * copy.probeCounts.z, 1);

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
		DDGIVolume copy = ddgiVolumes[0];
		
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

		computeIrradiance.SetTexture(0, "tex", isOutputIrradiance ? irradianceTexture : visibilityTexture);

		// Dispatch your compute shader
		int threadGroupsX = Mathf.CeilToInt((isOutputIrradiance ? irradianceTexture.width : visibilityTexture.width) / 8f);
		int threadGroupsY = Mathf.CeilToInt((isOutputIrradiance ? irradianceTexture.height : visibilityTexture.height) / 8f);
		computeIrradiance.Dispatch(0, threadGroupsX, threadGroupsY, 1);
	}

	void UpdateProbesBorders(RenderTexture probesTex, float probeSideLength, int kernelIndex)
	{
		computeBorders.SetTexture(kernelIndex, "tex", probesTex);
		computeBorders.SetFloat("PROBE_SIDE_LENGTH", probeSideLength);

		int threadGroupsX = Mathf.CeilToInt(probesTex.width / 8f);
		int threadGroupsY = Mathf.CeilToInt(probesTex.height / 8f);
		computeBorders.Dispatch(kernelIndex, threadGroupsX, threadGroupsY, 1);
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

	public void PrepareScene(ComputeShader cs)
    {
		CreateMeshBuffers();
        SetMeshesBuffer(cs);
        SetLightsValues(cs);
    }

	void CreateMeshBuffers()
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
	}

	void SetMeshesBuffer(ComputeShader cs)
	{
		cs.SetBuffer(0, "Triangles", triangleBuffer);
		cs.SetBuffer(0, "AllMeshInfo", meshInfoBuffer);
		cs.SetInt("NumMeshes", allMeshInfo.Count);
	}

    void SetLightsValues(ComputeShader cs)
    {
		if (sun && sun.isActiveAndEnabled)
		{
			cs.SetVector("sunColor", sun.color);
			cs.SetVector("sunDirection", sun.transform.forward);
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
		minProbesSpacing = Vector3.Max(minProbesSpacing, new Vector3(.05f,.05f,.05f));
        //if (minProbesSpacing != realProbeSpacing)
		RefreshProbesPlacement();

		DDGIVolume copy = ddgiVolumes[0];
		copy.invIrradianceTextureSize = DivideVectors(Vector2Int.one, new Vector2Int(irradianceTexture.width, irradianceTexture.height));
		copy.invVisibilityTextureSize = DivideVectors(Vector2Int.one, new Vector2Int(visibilityTexture.width, visibilityTexture.height));;
		copy.invIrradianceGamma = 1.0f / copy.irradianceGamma; 

		Vector3 probeEnd = copy.probeGridOrigin + MultiplyVectors(copy.probeCounts - Vector3Int.one, copy.probeSpacing);
		Vector3 probeSpan = probeEnd - copy.probeGridOrigin;
		copy.maxDistance = Vector3.Magnitude(DivideVectors(probeSpan, copy.probeCounts)) * 1.5f;
		ddgiVolumes[0] = copy;

		numRaysPerPixel = Mathf.Max(1, numRaysPerPixel);
		isWriteOnesDone = false;
    }

    public void RefreshProbesPlacement()
    {
        Debug.Log($"Refresh probes placement with a distance of {minProbesSpacing}.");

		DDGIVolume copy = ddgiVolumes[0];

        copy.probeCounts = new(
            Mathf.CeilToInt(boundsVolume.size.x / minProbesSpacing.x),
            Mathf.CeilToInt(boundsVolume.size.y / minProbesSpacing.y),
            Mathf.CeilToInt(boundsVolume.size.z / minProbesSpacing.z)
        );

        copy.probeSpacing = new(
            boundsVolume.size.x / (copy.probeCounts.x - 1),
            boundsVolume.size.y / (copy.probeCounts.y - 1),
            boundsVolume.size.z / (copy.probeCounts.z - 1)
        );

        copy.probeSpacing.x *= IsValid(copy.probeSpacing.x);
        copy.probeSpacing.y *= IsValid(copy.probeSpacing.y);
        copy.probeSpacing.z *= IsValid(copy.probeSpacing.z);

        probesPositions = new Vector3[copy.probeCounts.x * copy.probeCounts.y * copy.probeCounts.z];

        int idx = 0;
		for (int z = 0; z < copy.probeCounts.z; z++)
		{
			for (int y = 0; y < copy.probeCounts.y; y++)
			{
				for (int x = 0; x < copy.probeCounts.x; x++)
				{
					probesPositions[idx++] = boundsVolume.center - boundsVolume.extents + new Vector3(x * copy.probeSpacing.x, y * copy.probeSpacing.y, z * copy.probeSpacing.z);
				}
			}
		}

		copy.logProbeCounts = new Vector3Int(
			Mathf.FloorToInt(Mathf.Log(copy.probeCounts.x)),
			Mathf.FloorToInt(Mathf.Log(copy.probeCounts.y)),
			Mathf.FloorToInt(Mathf.Log(copy.probeCounts.z))
		);

		copy.probeGridOrigin = boundsVolume.min;
		copy.invProbeSpacing = Vector3.Max(DivideVectors(Vector3.one, copy.probeSpacing), Vector3.zero);

        realProbeSpacing = minProbesSpacing;
		ddgiVolumes[0] = copy;

		int IsValid(float spacing)
		{
			return !(float.IsInfinity(spacing) || float.IsNaN(spacing)) ? 1 : 0;
		}
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

	private static Vector2Int DivideVectors(Vector2Int a, Vector2Int b)
	{
		return new Vector2Int(a.x / b.x, a.y / b.y);
	}

	public void ResetBuffers()
	{
		ClearOutRenderTexture(visibilityTexture);
		ClearOutRenderTexture(irradianceTexture);

		ClearOutRenderTexture(rayHitLocationsBuffer);
		ClearOutRenderTexture(rayHitRadianceBuffer);
		ClearOutRenderTexture(rayHitNormalsBuffer);
		ClearOutRenderTexture(rayDirectionsBuffer);
		ClearOutRenderTexture(rayOriginsBuffer);

		ClearOutRenderTexture(resultTexture);
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
	[SerializeField] bool isRaytracingRendering;

	[Header("References")]
	[SerializeField] Shader rasterizedDDGIShader;
	[SerializeField] ComputeShader raytracedDDGIShader;

	[Header("Info")]
	[SerializeField] int numRenderedFrames;

	// Materials and render textures
	Material blitMaterial;
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
				if (isRaytracingRendering)
					Graphics.Blit(resultTexture, target);
				else
                	Graphics.Blit(src, target, blitMaterial);
			}
            else
                Graphics.Blit(src, target); // Draw the unaltered camera render to the screen
        }
        else
        {
			InitFrame();
			if (isRaytracingRendering)
				Graphics.Blit(resultTexture, target);
			else
				Graphics.Blit(src, target, blitMaterial);

            numRenderedFrames += Application.isPlaying ? 1 : 0;
        }
    }


	void InitFrame()
	{
		if (!isRaytracingRendering)
		{
			Camera cam = Camera.current;
			cam.depthTextureMode = DepthTextureMode.Depth;

			// Create materials used in blits
			ShaderHelper.InitMaterial(rasterizedDDGIShader, ref blitMaterial);

			blitMaterial.SetBuffer("DDGIVolumes", ddgiVolumesBuffer);

			blitMaterial.SetTexture("irradianceTex", irradianceTexture);
			blitMaterial.SetTexture("weightTex", visibilityTexture);

			blitMaterial.SetMatrix("_InverseView", cam.cameraToWorldMatrix);
			//blitMaterial.SetMatrix("", (cam.projectionMatrix * cam.worldToCameraMatrix).inverse.transpose);
		}
		else
		{
			Camera cam = Camera.current;
			RefreshBufferIfNeeded(ref resultTexture, "Result", cam.scaledPixelWidth, cam.scaledPixelHeight);
			raytracedDDGIShader.SetTexture(0, "Result", resultTexture);

			raytracedDDGIShader.SetBuffer(0, "DDGIVolumes", ddgiVolumesBuffer);

			raytracedDDGIShader.SetTexture(0, "irradianceTex", irradianceTexture);
			raytracedDDGIShader.SetTexture(0, "weightTex", visibilityTexture);

			// Set camera parameters
			float planeHeight = cam.farClipPlane * Mathf.Tan(cam.fieldOfView * 0.5f * Mathf.Deg2Rad) * 2;
			float planeWidth = planeHeight * cam.aspect;
			raytracedDDGIShader.SetVector("ViewParams", new Vector3(planeWidth, planeHeight, cam.farClipPlane));
			raytracedDDGIShader.SetMatrix("_CameraLocalToWorld", cam.transform.localToWorldMatrix);

			SetMeshesBuffer(raytracedDDGIShader);
			SetLightsValues(raytracedDDGIShader);

			// Dispatch your compute shader
			int threadGroupsX = Mathf.CeilToInt(cam.scaledPixelWidth/8.0f);
			int threadGroupsY = Mathf.CeilToInt(cam.scaledPixelHeight/8.0f);
			raytracedDDGIShader.Dispatch(0, threadGroupsX, threadGroupsY, 1);
		}
	}


	void OnDisable()
	{
		ShaderHelper.Release(resultTexture);
	}

	

	public bool NotifyOfCameraPosition(Vector3 cameraWSPosition)
	{
		DDGIVolume copy = ddgiVolumes[0];
		if (copy.isCameraLocked == 0)
			return false;

		// Centerpoint of the volume.
		Vector3 volumeCenter = copy.probeGridOrigin + MultiplyVectors(copy.probeSpacing, (copy.probeCounts - new Vector3(1, 1, 1)) / 2.0f);
		Bounds volumeCenterRegion = new Bounds(volumeCenter, copy.probeSpacing * 2.0f);

		// If the camera is within our center volume, don't need to move.
		if (volumeCenterRegion.Contains(cameraWSPosition))
			return false;

		// Move 2 probeSpacing lengths along the axis that would get the camera closest to the new centerpoint

		// Distance to center, pointing towards the camera to make computation easier.
		Vector3 vectorToCenter = cameraWSPosition - volumeCenter;

		Vector3Int[] movementOptions = {
			new Vector3Int(1, 0, 0),
			new Vector3Int(-1, 0, 0),
			new Vector3Int(0, 1, 0),
			new Vector3Int(0, -1, 0),
			new Vector3Int(0, 0, 1),
			new Vector3Int(0, 0, -1)
		};

		int index = -1;
		float length = vectorToCenter.magnitude;
		for (int i = 0; i < movementOptions.Length; ++i)
		{
			float newLength = (vectorToCenter - MultiplyVectors(movementOptions[i], copy.probeSpacing)).magnitude;
			if (newLength < length)
			{
				length = newLength;
				index = i;
			}
		}

		copy.probeGridOrigin += MultiplyVectors(movementOptions[index], copy.probeSpacing);

		// Phase offset works opposite the motion of the camera.
		copy.phaseOffsets -= movementOptions[index];

		ddgiVolumes[0] = copy;
		return true;
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