using UnityEditor;
using UnityEditor.IMGUI.Controls;
using UnityEngine;
using UnityEngine.Rendering;

[CustomEditor(typeof(DDGIController))]
public class DDGIControllerEditor : Editor
{
    private BoxBoundsHandle m_BoundsHandle = new BoxBoundsHandle();
    private DDGIController ddgiController;

    void Awake()
    {
        ddgiController = (DDGIController)target;
    }

    // the OnSceneGUI callback uses the Scene view camera for drawing handles by default
    protected virtual void OnSceneGUI()
    {
        // copy the target object's data to the handle
        m_BoundsHandle.center = ddgiController.volume.center;
        m_BoundsHandle.size = ddgiController.volume.size;

        // draw the handle
        EditorGUI.BeginChangeCheck();
        m_BoundsHandle.DrawHandle();
        if (EditorGUI.EndChangeCheck())
        {
            // record the target object before setting new values so changes can be undone/redone
            Undo.RecordObject(ddgiController, "Change Bounds");

            // copy the handle's updated data back to the target object
            Bounds newBounds = new()
            {
                center = m_BoundsHandle.center,
                size = m_BoundsHandle.size
            };
            ddgiController.volume = newBounds;

            ddgiController.RefreshProbesPlacement();
        }
    }


    public override void OnInspectorGUI()
    {
        base.OnInspectorGUI();

        if (!ddgiController.isRealtimeRaytracing){
            if (GUILayout.Button("Prepare Scene"))
            {
                ddgiController.PrepareScene();
            }

            if (GUILayout.Button("Trace Rays"))
            {
                ddgiController.DispatchDDGIRayTracing();
            }
        }
    }
}
