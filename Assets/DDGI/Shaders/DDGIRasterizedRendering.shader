Shader "DDGI/ScreenSpaceApplyDDGI"
{
    HLSLINCLUDE
    #include "UnityCG.cginc"
    #include "SampleIrradianceField.cginc"

    sampler2D _MainTex;
    sampler2D _CameraDepthTexture;
    sampler2D _CameraGBufferTexture2;

    float4x4 _InverseView;

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
                float3 normal = tex2D(_CameraGBufferTexture2, input.texcoord);
                //normal = UnityObjectToWorldNormal(normal);

                float4 vpos = float4(ComputeViewSpacePosition(input), 1);
                float4 wpos = mul(_InverseView, vpos);

                float3 viewDir = WorldSpaceViewDir(vpos);
                
                float3 color = tex2D(_MainTex, input.texcoord);
                float3 irradiance = SampleIrradianceField(wpos.xyz, normal, LBuffer[0].energyConservation, viewDir);
                return float4(color + irradiance, 0);
            }

            ENDHLSL
        }
    }
}