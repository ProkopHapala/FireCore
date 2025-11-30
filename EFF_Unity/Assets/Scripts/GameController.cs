using System;
using System.Collections;
using System.IO;
using System.Linq;
using System.Runtime.InteropServices;
using particles;
using TMPro;
using Unity.VisualScripting;
using UnityEngine;
using UnityEngine.InputSystem.iOS;
using Directory = UnityEngine.Windows.Directory;
using File = UnityEngine.Windows.File;


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
    public GameObject fileNameText;

    public InputFieldManager inputFields;

    private IDataFeeder feeder;
    
    public bool isRunning = false;
    public bool isPlayback = false;
    
    private Vector3[] positions;
    private float[] sizes;
    private int[] espins;
    private int[] protonNumbers;

    public static string rootPath;


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


    void Awake() {
        main = this;
        enabled = false;
    }
    
    
    // Start is called once before the first execution of Update after the MonoBehaviour is created
    void Start()
    {
        electrons = new Electron[electronCount];
        atoms = new Atom[atomCount];

        for (int i = 0; i < electronCount; i++) {
            Electron.CreateNew(espins[i]);
        }
        for (int i = 0; i < atomCount; i++) {

            Atom atom = Atom.CreateNew();
            // atom.ProtonNumber = new[] { 6, 1, 1, 1, 1 }[i];
            atom.ProtonNumber = protonNumbers[i];
            if (atom.Style == RenderingStyle.SPHERICAL) {
                var visualData = AtomVisualParams.Get(atom.ProtonNumber);
                atom.Sphere.transform.localScale += visualData.size * 0.25f * Vector3.one;
                atom.Sphere.GetComponent<Renderer>().material.color = visualData.color;
            }
            
            Debug.Log(atoms[i]);
        }
        SpawnGUI();
        
        
        
        DoVisualUpdateCycle();
        Array.ForEach(infoBoxAnchor.GetComponentsInChildren<InfoBox>(), x => x.ForceUpdate());
    }

    public void SpawnGUI() {
        Array.ForEach(GameObject.FindGameObjectsWithTag("InfoBox"), Destroy);

        int column = 0;
        int row = 0;
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
            
            row++;
            if(63*(row + 2) > Screen.height) {
                row = 0;
                column++;
            }
        }
    }

    // Update is called once per frame

    public void StartSimulation(string fileName)
    {
        isPlayback = System.IO.File.ReadLines(fileName).First().Contains("#playback");
        // isPlayback = false;
        feeder = isPlayback ? new PlaybackDataFeeder() : new SimulationDataFeeder();
        fileNameText.GetComponent<TextMeshProUGUI>().text = Path.GetFileName(fileName) + " " + (isPlayback ? "[playback]" : "[simulated]");
        fileNameText.SetActive(true);
        if (isPlayback) {
            SetParticlePosition = (_, _, _, _) => { }; // positions cannot be set in playback mode
        }
        else {
            SetParticlePosition = _simSetParticlePosition;
        }
        
        var data = feeder.Init(fileName, () => (electronCount, atomCount));
        electronCount = data.electronCount;
        atomCount = data.atomCount;
        espins = data.spins;
        
        // Get proton numbers from C++
        if (!isPlayback) {
            IntPtr protonNumbersPtr = SimulationManager.unityGetProtonNumbers();
            protonNumbers = new int[atomCount];
            Marshal.Copy(protonNumbersPtr, protonNumbers, 0, atomCount);
            SimulationManager.cleanupProtonNumbers(protonNumbersPtr);
        } else {
            protonNumbers = new int[atomCount]; // Empty array for playback mode
        }
        
        enabled = true;
    }

    public void FindAnchor()
    {
        // Start from current directory (represented by ".")
        string currentPath = Path.GetFullPath(".");
        
        while (currentPath != null)
        {
            // Check if .anchor file exists in current directory
            string anchorFilePath = Path.Combine(currentPath, ".anchor");
            
            if (System.IO.File.Exists(anchorFilePath))
            {
                // Found the anchor file, return the absolute path of its containing folder
                rootPath = currentPath;
                return;
            }
            
            // Move to parent directory
            DirectoryInfo parentDir = System.IO.Directory.GetParent(currentPath);
            
            // If parent is null, we've reached the root directory
            if (parentDir == null)
            {
                break;
            }

            // Debug.Log(Path.Combine(currentPath, ".anchor"));
            currentPath = parentDir.FullName;
        }
        
        throw new Exception("Anchor file not found");
    }

    void Update()
    {
        if(isRunning){
            DoVisualUpdateCycle();
        }
    }

    private void DoVisualUpdateCycle() {
        var data = feeder.GetNextFrame();
        positions = data.positions;
        sizes = data.sizes;
        
        for (int i = 0; i < electronCount; i++) {
            electrons[i].Position = positions[i]; // Electrons are at indices 0 to electronCount-1
            electrons[i].Size = sizes[i];
        }
        
        for (int i = 0; 
             i < atomCount; i++) {
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

    public void ToggleKeybindList() {
        keybinds.SetActive(!keybinds.activeInHierarchy);
    }

    // Switching functions
    
    // Init
    // private delegate (int electronCount, int atomCount, int[] spins) InitDel(string fileName);
    // private readonly InitDel _simInit = fileName =>
    // {
    //     var dataPtr = SimulationManager.unityInit(fileName);
    //     // First get the electron and atom counts
    //     int ne = Marshal.ReadInt32(dataPtr, 0);
    //     int na = Marshal.ReadInt32(dataPtr, sizeof(int));
    //     
    //     // Create and populate the espins array
    //     var spins = new int[ne];
    //     for (int i = 0; i < ne; i++) {
    //         spins[i] = Marshal.ReadInt32(dataPtr, (2 + i) * sizeof(int));
    //     }
    //     
    //     // Free the unmanaged memory
    //     SimulationManager.cleanupInitData(dataPtr);
    //
    //     return (ne, na, spins);
    // };
    //
    // private readonly InitDel _playbackInit;
    //
    // private static InitDel Init;
    
    
    // // GetPositions
    // private delegate (Vector3[] positions, float[] sizes) GetNextFrameDel(int electronCount, int atomCount);
    // private readonly GetNextFrameDel _simGetNextFrame = (electronCount, atomCount) =>
    // {
    //     int size;
    //     IntPtr positionsPtr = SimulationManager.unityNextFrame(out size);
    //     
    //     // Convert to managed array
    //     float[] rawPositions = new float[size];
    //
    //     Marshal.Copy(positionsPtr, rawPositions, 0, size);
    //     
    //     // Clean up native memory
    //     SimulationManager.cleanupPositions(positionsPtr);
    //     
    //     Vector3[] positions = new Vector3[atomCount + electronCount];
    //     float[] sizes = new float[electronCount];
    //     
    //     int idx = 0;
    //     // First, extract electron positions
    //     for (int i = 0; i < electronCount; i++)
    //     {
    //         positions[i] = new Vector3(
    //         rawPositions[idx++],
    //         rawPositions[idx++],
    //         rawPositions[idx++]
    //         );
    //     }
    //     
    //     // Then extract atom positions
    //     for (int i = 0; i < atomCount; i++)
    //     {
    //         positions[electronCount + i] = new Vector3(
    //         rawPositions[idx++],
    //         rawPositions[idx++],
    //         rawPositions[idx++]
    //         );
    //     }
    //     
    //     // Finally extract electron sizes
    //     for (int i = 0; i < electronCount; i++) {
    //         sizes[i] = rawPositions[idx++];
    //     }
    //     
    //     
    //     return (positions, sizes);
    // };
    //
    // private readonly GetNextFrameDel _playbackGetNextFrame;
    //
    // private static GetNextFrameDel GetNextFrame;
    
    // SetParticlePosition
    public delegate void SetParticlePositionDel(ObjectType type, IParticle particle, Vector3 pos, float size);
    private static readonly SetParticlePositionDel _simSetParticlePosition = (type, particle, pos, size) =>
    {
        try {
            particle.Position = pos;
            if(type == ObjectType.ELECTRON) {
                SimulationManager.unitySetElectronPosition(particle.Id, pos.x, pos.y, pos.z, size);
                ((Electron)particle).Size = size;
            } else {
                SimulationManager.unitySetAtomPosition(particle.Id, pos.x, pos.y, pos.z);
            }
        } catch (Exception e) {
            Debug.LogError($"Exception in SetParticlePosition: {e.GetType()} : {e.Message}\n{e.StackTrace}");
        }
    };

    public SetParticlePositionDel SetParticlePosition = _simSetParticlePosition;
}