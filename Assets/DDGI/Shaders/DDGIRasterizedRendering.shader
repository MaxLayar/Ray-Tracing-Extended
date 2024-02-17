Shader "DDGI/ScreenSpaceApplyDDGI"
{
    HLSLINCLUDE
    #pragma multi_compile_local __ SHOW_PROBES

	#include <UnityPBSLighting.cginc>

	#include "UnityCG.cginc"
	#include "SampleIrradianceField.cginc"

	#pragma target 5.0

    sampler2D _MainTex;
    sampler2D _CameraDepthTexture;
    sampler2D _CameraGBufferTexture2;

    float4x4 _InverseView;
    float4x4 _ViewProjInv;

#if defined(SHOW_PROBES)
    StructuredBuffer<float3> ProbesPositions;
    float DebugProbesRadius;
#endif

    struct Varyings
    {
        float4 position : SV_Position;
        float2 texcoord : TEXCOORD0;
        float3 ray : TEXCOORD1;
    };

    // Vertex shader that procedurally outputs a full screen triangle
    Varyings Vertex(uint vertexID : SV_VertexID)
    {
        // Render settings
        float far = _ProjectionParams.z;
        float2 orthoSize = unity_OrthoParams.xy;
        float isOrtho = unity_OrthoParams.w; // 0: perspective, 1: orthographic

        // Vertex ID -> clip space vertex position
        float x = (vertexID != 1) ? -1 : 3;
        float y = (vertexID == 2) ? -3 : 1;
        float3 vpos = float3(x, y, 1.0);

        // Perspective: view space vertex position of the far plane
        float3 rayPers = mul(unity_CameraInvProjection, vpos.xyzz * far).xyz;

        // Orthographic: view space vertex position
        float3 rayOrtho = float3(orthoSize * vpos.xy, 0);

        Varyings o;
        o.position = float4(vpos.x, -vpos.y, 1, 1);
        o.texcoord = (vpos.xy + 1) / 2;
        o.ray = lerp(rayPers, rayOrtho, isOrtho);
        return o;
    }

    float3 ComputeViewSpacePosition(Varyings input)
    {
        // Render settings
        float near = _ProjectionParams.y;
        float far = _ProjectionParams.z;
        float isOrtho = unity_OrthoParams.w; // 0: perspective, 1: orthographic

        // Z buffer sample
        float z = tex2D(_CameraDepthTexture, input.texcoord);

        // Far plane exclusion
        #if !defined(EXCLUDE_FAR_PLANE)
        float mask = 1;
        #elif defined(UNITY_REVERSED_Z)
        float mask = z > 0;
        #else
        float mask = z < 1;
        #endif

        // Perspective: view space position = ray * depth
        float3 vposPers = input.ray * Linear01Depth(z);

        // Orthographic: linear depth (with reverse-Z support)
        #if defined(UNITY_REVERSED_Z)
        float depthOrtho = -lerp(far, near, z);
        #else
        float depthOrtho = -lerp(near, far, z);
        #endif

        // Orthographic: view space position
        float3 vposOrtho = float4(input.ray.xy, depthOrtho, 1);

        // Result: view space position
        return lerp(vposPers, vposOrtho, isOrtho) * mask;
    }

    // https://github.com/DiligentGraphics/DiligentSamples/blob/master/Tutorials/Tutorial22_HybridRendering/assets/Utils.fxh
    float3 ScreenPosToWorldPos(float2 ScreenSpaceUV)
    {
        float4 PosClipSpace;
        PosClipSpace.xy = ScreenSpaceUV * float2(2.0, -2.0) + float2(-1.0, 1.0);
        PosClipSpace.z = tex2D(_CameraDepthTexture, ScreenSpaceUV);
        PosClipSpace.w = 1.0;
        float4 WorldPos = mul(PosClipSpace, _ViewProjInv);
        return WorldPos.xyz / WorldPos.w;
    }

    // https://www.iquilezles.org/www/articles/spherefunctions/spherefunctions.htm
    float sphIntersect(float3 ro, float3 rd, float4 sph)
    {
        float3 oc = ro - sph.xyz;
        float b = dot( oc, rd );
        float c = dot( oc, oc ) - sph.w*sph.w;
        float h = b*b - c;
        if( h<0.0 ) return -1.0;
        h = sqrt( h );
        return -b - h;
    }

    ENDHLSL
    Properties
    {
        _MainTex ("Texture", 2D) = "white" {}
    }
    SubShader
    {
        Cull Off ZWrite Off ZTest Always
        Pass
        {
            HLSLPROGRAM

            #pragma vertex Vertex
            #pragma fragment Fragment

            half4 Fragment(Varyings input) : SV_Target
            {
                float4 vpos = float4(ComputeViewSpacePosition(input), 1);
                float4 wpos = mul(_InverseView, vpos);

				float3 wNormal = tex2D(_CameraGBufferTexture2, input.texcoord).rgb * 2. - 1.;
				wNormal = normalize(wNormal);

                int notSphereHit = 1;

#if defined(SHOW_PROBES)
                // ray origin
                float3 rayOrigin = _WorldSpaceCameraPos.xyz;

                // normalize ray vector
                float3 rayDir = normalize(wpos - _WorldSpaceCameraPos.xyz);

                float dist = distance(wpos, _WorldSpaceCameraPos);

                int numProbes = DDGIVolumes[0].probeCounts.x * DDGIVolumes[0].probeCounts.y * DDGIVolumes[0].probeCounts.z;

                float minDist = dist;

                for(int i = 0; i < numProbes; i++)
                {
                    float3 spherePos = ProbesPositions[i];

                    // ray box intersection
                    float rayHit = sphIntersect(rayOrigin, rayDir, float4(spherePos, DebugProbesRadius));
                    
                    if (rayHit < 0)
                        continue;
                    
                    if (rayHit >= minDist)
                        continue;
                    
                    minDist = rayHit;

                    notSphereHit = 0;

                    // override world space position from ray, front hit ray length, and ray origin
                    wpos = float4(rayDir * rayHit + rayOrigin, 0);

                    // override world space surface normal
                    wNormal = normalize(wpos - spherePos);
                }
#endif
                if (distance(wpos, _WorldSpaceCameraPos) > 50.0)
                    return float4(0,0,0.1,1);
                    
                float3 viewVec = normalize(UnityWorldSpaceViewDir(wpos));

                float3 irradiance = sampleIrradiance(DDGIVolumes, wpos, wNormal*.2+viewVec*.8, wNormal, _WorldSpaceCameraPos, false, false, -1);

                float3 color = tex2D(_MainTex, input.texcoord);

                return float4(irradiance + color * notSphereHit, 1);
            }

            ENDHLSL
        }
    }
}