using System;
using System.Diagnostics;
using System.Linq;
using System.Runtime.InteropServices;
using TMPro;
using UnityEngine;
using UnityEngine.InputSystem.iOS;

public class GameController : MonoBehaviour
{
    public static GameController main; // static self reference object

    private int atomCount;
    private int electronCount;

    public GameObject atomPrefab;
    public GameObject electronPrefab;
    public GameObject infoBoxPrefab_e;
    public GameObject infoBoxPrefab_a;
    public GameObject infoBoxAnchor;
    public GameObject[] particles;
    public TextMeshProUGUI runningStatusText;

    public InputFieldManager inputFields;

    public const string PATH_TO_EFF_APP = "/home/perry/FireCore/cpp/Build/apps/EFF/libEFFapp_console.so";
    
[DllImport(PATH_TO_EFF_APP, CallingConvention = CallingConvention.Cdecl, CharSet = CharSet.Ansi)]
public static extern IntPtr unityInit(string fileName);

[DllImport(PATH_TO_EFF_APP, CallingConvention = CallingConvention.Cdecl, CharSet = CharSet.Ansi)]
public static extern IntPtr unityNextFrame(out int size);

[DllImport(PATH_TO_EFF_APP, CallingConvention = CallingConvention.Cdecl, CharSet = CharSet.Ansi)]
public static extern void cleanupPositions(IntPtr positions);

    [DllImport(PATH_TO_EFF_APP, CallingConvention = CallingConvention.Cdecl, CharSet = CharSet.Ansi)]
    public static extern void unitySetElectronPosition(int electronIndex, float x, float y, float z, float size);
    
    [DllImport(PATH_TO_EFF_APP, CallingConvention = CallingConvention.Cdecl, CharSet = CharSet.Ansi)]
    public static extern void unitySetAtomPosition(int atomIndex, float x, float y, float z);

    public bool isRunning = false;

    public Vector3[] positions { get; private set; }
    public float[] sizes { get; private set; }

    private (Vector3[] positions, float[] sizes) GetPositionsFromNative()
    {
        int size;
        IntPtr positionsPtr = unityNextFrame(out size);
        
        // Convert to managed array
        float[] rawPositions = new float[size];
    
        Marshal.Copy(positionsPtr, rawPositions, 0, size);
        
        
        // Convert to Vector3 array
        Vector3[] positions = new Vector3[(size - electronCount) / 3];
        float[] sizes = new float[electronCount];
        for (int i = 0; i < positions.Length; i++)
        {
            positions[i] = new Vector3(
                rawPositions[i * 3],
                rawPositions[i * 3 + 1],
                rawPositions[i * 3 + 2]
            );
            //UnityEngine.Debug.Log(positions[i]);
        }
        for (int i = 0; i < electronCount; i++) {
            sizes[i] = rawPositions[i + (positions.Length * 3)];
            //UnityEngine.Debug.Log(sizes[i]);
        }
        
        
        // Clean up native memory
        cleanupPositions(positionsPtr);
        
        return (positions, sizes);
    }


    void Awake() {
        main = this;
        enabled = false;

        // unityInit(@"../cpp/sketches_SDL/Molecular/data/H2.fgo");
        // int[] counts = new int[2];
        // Marshal.Copy(unityInit(@"../cpp/sketches_SDL/Molecular/data/H2O_shoot.fgo"), counts, 0, 2);
        // electronCount = counts[0];
        // atomCount = counts[1];


        // UnityEngine.Debug.Log($"{electronCount}, {atomCount}");
    }
    // Start is called once before the first execution of Update after the MonoBehaviour is created
    void Start()
    {
        particles = new GameObject[atomCount + electronCount];

        for (int i = 0; i < electronCount; i++) {
            particles[i] = Instantiate(electronPrefab, new Vector3(0, 0, 0), Quaternion.identity);

            var info = Instantiate(infoBoxPrefab_e, infoBoxAnchor.transform);
            info.GetComponent<InfoBox>().SetConnector(i, i, ObjectType.ELECTRON);
            info.GetComponent<RectTransform>().anchoredPosition = new Vector2(0, -63 * i);
        }
        for (int i = 0; i < atomCount; i++) {
            particles[electronCount + i] = Instantiate(atomPrefab, new Vector3(0, 0, 0), Quaternion.identity);

            var info = Instantiate(infoBoxPrefab_a, infoBoxAnchor.transform);
            info.GetComponent<InfoBox>().SetConnector(i, electronCount + i, ObjectType.ATOM);
            info.GetComponent<RectTransform>().anchoredPosition = new Vector2(0, (-63 * electronCount) + (-50 * i));
        }

        
        DoSimUpdateCycle();
        Array.ForEach(infoBoxAnchor.GetComponentsInChildren<InfoBox>(), x => x.ForceUpdate());
    }

    // Update is called once per frame

    public void StartSimulation(string fileName) {
        int[] counts = new int[2];
        Marshal.Copy(unityInit(@"../cpp/sketches_SDL/Molecular/data/" + fileName), counts, 0, 2);
        electronCount = counts[0];
        atomCount = counts[1];
        enabled = true;
    }

    void Update()
    {
        if(isRunning){
            DoSimUpdateCycle();
        }
    }

    private void DoSimUpdateCycle() {
        var data = GetPositionsFromNative();
        positions = data.positions;
        sizes = data.sizes;

        for (int i = 0; i < particles.Length; i++) {
            particles[i].transform.position = positions[i];
            if(i < electronCount){
                particles[i].transform.localScale = new Vector3(sizes[i], sizes[i], sizes[i]);
                //UnityEngine.Debug.Log(sizes[i]);
            }
            //UnityEngine.Debug.Log(positions[i].x + " " + positions[i].y + " " + positions[i].z);

        }
    }

    public void Continue() {
        isRunning = true;
        runningStatusText.SetText("RUNNING");
    }

    public void Stop() {
        isRunning = false;
        runningStatusText.SetText("PAUSED");
    }
    public void ToggleRunning() {
        if(isRunning) {
            Stop();
        }
        else {
            Continue();
        }
    }

    public Vector3 GetPosition(ObjectType type, int displayId) {
        return positions[type == ObjectType.ELECTRON ? displayId : displayId + electronCount];
    }

    public void SetParticlePosition(ObjectType type, int id, Vector3 pos, float size) {
        try {
            if(type == ObjectType.ELECTRON) {
                unitySetElectronPosition(id, pos.x, pos.y, pos.z, size);
                particles[id].transform.position = pos;
            } else {
                unitySetAtomPosition(id, pos.x, pos.y, pos.z);
                particles[id + electronCount].transform.localScale = new Vector3(size, size, size);
            }
        } catch (Exception e) {
            UnityEngine.Debug.LogError($"Exception in SetParticlePosition: {e.GetType()} : {e.Message}\n{e.StackTrace}");
        }
    }
}
