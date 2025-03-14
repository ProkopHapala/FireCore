using System;
using System.Diagnostics;
using System.Runtime.InteropServices;
using UnityEngine;
using UnityEngine.InputSystem.iOS;

public class CppConnector : MonoBehaviour
{
    private int atomCount;
    private int electronCount;

    public GameObject atomPrefab;
    public GameObject electronPrefab;

    public GameObject[] particles;

    public const string PATH_TO_EFF_APP = "/home/perry/FireCore/cpp/Build/apps/EFF/libEFFapp_console.so";
    
    [DllImport(PATH_TO_EFF_APP)]
    public static extern IntPtr unityInit(string fileName);
    
    [DllImport(PATH_TO_EFF_APP)]
    public static extern IntPtr unityNextFrame(out int size);
    
    [DllImport(PATH_TO_EFF_APP)]
    public static extern void cleanupPositions(IntPtr positions);

    public bool isRunning = false;


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
        }
        for (int i = 0; i < atomCount; i++) {
            particles[electronCount + i] = Instantiate(atomPrefab, new Vector3(0, 0, 0), Quaternion.identity);
        }
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
            var data = GetPositionsFromNative();
            Vector3[] positions = data.positions;
            float[] sizes = data.sizes;

            for (int i = 0; i < particles.Length; i++) {
                particles[i].transform.position = positions[i];
                if(i < electronCount){
                    particles[i].transform.localScale = new Vector3(sizes[i], sizes[i], sizes[i]);
                    //UnityEngine.Debug.Log(sizes[i]);
                }
                UnityEngine.Debug.Log(positions[i].x + " " + positions[i].y + " " + positions[i].z);

            }
        }
    }
}
