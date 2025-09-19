using System;
using System.Collections;
using System.Linq;
using System.Runtime.InteropServices;
using TMPro;
using Unity.VisualScripting;
using UnityEngine;
using UnityEngine.InputSystem.iOS;

public class GameController : MonoBehaviour
{
    public static GameController main; // static self reference object

    private int atomCount;
    private int electronCount;

    public GameObject atomPrefab;
    public GameObject electronPrefabPlusSpin;
    public GameObject electronPrefabMinusSpin;
    public GameObject infoBoxPrefab_e;
    public GameObject infoBoxPrefab_a;
    public GameObject infoBoxAnchor;
    public Electron[] electrons;
    public Atom[] atoms;
    public TextMeshProUGUI runningStatusText;
    public CameraControl cameraControl;
    public GameObject keybinds;
    public GameObject canvas;

    public InputFieldManager inputFields;

    public const string PATH_TO_EFF_APP = "/home/perry/FireCore/cpp/Build/apps/EFF/libEFFapp_console.so";
    
[DllImport(PATH_TO_EFF_APP, CallingConvention = CallingConvention.Cdecl, CharSet = CharSet.Ansi)]
public static extern IntPtr unityInit(string fileName);

[DllImport(PATH_TO_EFF_APP, CallingConvention = CallingConvention.Cdecl, CharSet = CharSet.Ansi)]
public static extern IntPtr unityNextFrame(out int size);

[DllImport(PATH_TO_EFF_APP, CallingConvention = CallingConvention.Cdecl, CharSet = CharSet.Ansi)]
public static extern void cleanupPositions(IntPtr positions);

[DllImport(PATH_TO_EFF_APP, CallingConvention = CallingConvention.Cdecl, CharSet = CharSet.Ansi)]
public static extern void cleanupInitData(IntPtr data);

    [DllImport(PATH_TO_EFF_APP, CallingConvention = CallingConvention.Cdecl, CharSet = CharSet.Ansi)]
    public static extern void unitySetElectronPosition(int electronIndex, float x, float y, float z, float size);
    
    [DllImport(PATH_TO_EFF_APP, CallingConvention = CallingConvention.Cdecl, CharSet = CharSet.Ansi)]
    public static extern void unitySetAtomPosition(int atomIndex, float x, float y, float z);

    [DllImport(PATH_TO_EFF_APP, CallingConvention = CallingConvention.Cdecl, CharSet = CharSet.Ansi)]
    public static extern void fixAtom(int atomIndex);

    [DllImport(PATH_TO_EFF_APP, CallingConvention = CallingConvention.Cdecl, CharSet = CharSet.Ansi)]
    public static extern void fixElectron(int atomIndex);

    [DllImport(PATH_TO_EFF_APP, CallingConvention = CallingConvention.Cdecl, CharSet = CharSet.Ansi)]
    public static extern void unfixAtom(int atomIndex);

    [DllImport(PATH_TO_EFF_APP, CallingConvention = CallingConvention.Cdecl, CharSet = CharSet.Ansi)]
    public static extern void unfixElectron(int atomIndex);



    public bool isRunning = false;

    private Vector3[] positions;
    private float[] sizes;
    private int[] espins;


    public IParticle GetParticle(int id) {
        if(id >= electronCount) {
            return atoms[id - electronCount];
        }
        return electrons[id];
    }
    public IParticle GetParticle(int id, ObjectType type) {
        if(type == ObjectType.ATOM) {
            return atoms[id];
        }
        return electrons[id];
    }

    private (Vector3[] positions, float[] sizes) GetPositionsFromNative()
    {
        int size;
        IntPtr positionsPtr = unityNextFrame(out size);
        
        // Convert to managed array
        float[] rawPositions = new float[size];
    
        Marshal.Copy(positionsPtr, rawPositions, 0, size);
        
        
        // BUG: This should be (atomCount + electronCount), not (size - electronCount) / 3
        Vector3[] positions = new Vector3[atomCount + electronCount];
        float[] sizes = new float[electronCount];
        
        int idx = 0;
        // First, extract electron positions
        for (int i = 0; i < electronCount; i++)
        {
            positions[i] = new Vector3(
                rawPositions[idx++],
                rawPositions[idx++],
                rawPositions[idx++]
            );
        }
        
        // Then extract atom positions
        for (int i = 0; i < atomCount; i++)
        {
            positions[electronCount + i] = new Vector3(
                rawPositions[idx++],
                rawPositions[idx++],
                rawPositions[idx++]
            );
        }
        
        // Finally extract electron sizes
        for (int i = 0; i < electronCount; i++) {
            sizes[i] = rawPositions[idx++];
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
        electrons = new Electron[electronCount];
        atoms = new Atom[atomCount];
        //particles = new GameObject[atomCount + electronCount];

        for (int i = 0; i < electronCount; i++) {
            // electrons[i] = 
            // (
            //     espins[i] == 1 ? 
            //     Instantiate(electronPrefabPlusSpin, new Vector3(0, 0, 0), Quaternion.identity) : 
            //     Instantiate(electronPrefabMinusSpin, new Vector3(0, 0, 0), Quaternion.identity)
            // ).GetComponent<Electron>();
            // Debug.Log(string.Join(", ", espins));
            electrons[i] = Electron.CreateNew(espins[i]);

            // var info = Instantiate(infoBoxPrefab_e, infoBoxAnchor.transform);
            // info.GetComponent<InfoBox>().SetConnector(i, i, ObjectType.ELECTRON);
            // info.GetComponent<RectTransform>().anchoredPosition = new Vector2(0, -63 * i);
        }
        for (int i = 0; i < atomCount; i++) {
            // particles[electronCount + i] = Instantiate(atomPrefab, new Vector3(0, 0, 0), Quaternion.identity);

            Atom.CreateNew();
            Debug.Log(atoms[i]);
            // var info = Instantiate(infoBoxPrefab_a, infoBoxAnchor.transform);
            // info.GetComponent<InfoBox>().SetConnector(i, electronCount + i, ObjectType.ATOM);
            // info.GetComponent<RectTransform>().anchoredPosition = new Vector2(0, (-63 * electronCount) + (-50 * i));
        }
        SpawnGUI();
        
        DoSimUpdateCycle();
        Array.ForEach(infoBoxAnchor.GetComponentsInChildren<InfoBox>(), x => x.ForceUpdate());
    }

    public void SpawnGUI() {
        Array.ForEach(GameObject.FindGameObjectsWithTag("InfoBox"), Destroy);

        int column = 0;
        int row = 0;
        float worldHeight = Camera.main.orthographicSize * 2;
        for (int i = 0; i < electronCount; i++) {
            var info = Instantiate(infoBoxPrefab_e, infoBoxAnchor.transform);
            info.GetComponent<InfoBox>().SetConnector(i, ObjectType.ELECTRON);
            //info.GetComponent<RectTransform>().anchoredPosition = new Vector2(340, -63 * row);
    
            row++;
            if(63*(row + 2) > Screen.height) {
                row = 0;
                column++;
            }
        }

        int offsetRows = row;
        for (int i = 0; i < atomCount; i++) {
            var info = Instantiate(infoBoxPrefab_a, infoBoxAnchor.transform);
            info.GetComponent<InfoBox>().SetConnector(i, ObjectType.ATOM);
            //info.GetComponent<RectTransform>().anchoredPosition = new Vector2(340 * column, (-50 * (row - offsetRows)) + (-63 * offsetRows));
            
            row++;
            if(63*(row + 2) > Screen.height) {
                offsetRows = 0;
                row = 0;
                column++;
            }
        }
    }

    // Update is called once per frame

    public void StartSimulation(string fileName) {
        // Get the pointer to the data returned from unityInit
        IntPtr dataPtr = unityInit(@"../cpp/sketches_SDL/Molecular/data/" + fileName);
        
        // First get the electron and atom counts
        electronCount = Marshal.ReadInt32(dataPtr, 0);
        atomCount = Marshal.ReadInt32(dataPtr, sizeof(int));
        
        // Create and populate the espins array
        espins = new int[electronCount];
        for (int i = 0; i < electronCount; i++) {
            espins[i] = Marshal.ReadInt32(dataPtr, (2 + i) * sizeof(int));
        }
        
        // Free the unmanaged memory
        cleanupInitData(dataPtr);
        
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

        // Now the positions array contains electrons first (indices 0 to electronCount-1)
        // then atoms (indices electronCount to electronCount+atomCount-1)
        
        for (int i = 0; i < electronCount; i++) {
            electrons[i].Position = positions[i]; // Electrons are at indices 0 to electronCount-1
            electrons[i].Size = sizes[i];
        }
        
        for (int i = 0; i < atomCount; i++) {
            atoms[i].Position = positions[electronCount + i]; // Atoms are at indices electronCount onwards
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

    public void SetParticlePosition(ObjectType type, IParticle particle, Vector3 pos, float size) {
        try {
            particle.Position = pos;
            if(type == ObjectType.ELECTRON) {
                unitySetElectronPosition(particle.Id, pos.x, pos.y, pos.z, size);
                ((Electron)particle).Size = size;
            } else {
                unitySetAtomPosition(particle.Id, pos.x, pos.y, pos.z);
            }
        } catch (Exception e) {
            UnityEngine.Debug.LogError($"Exception in SetParticlePosition: {e.GetType()} : {e.Message}\n{e.StackTrace}");
        }
    }

    public void ToggleKeybindList() {
        keybinds.SetActive(!keybinds.activeInHierarchy);
    }
}
