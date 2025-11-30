using System.IO;
using System.Linq;
using TMPro;
using UnityEngine;
using UnityEngine.UI;
using SFB; // StandaloneFileBrowser namespace

public class FgoPicker : MonoBehaviour
{
    [SerializeField]
    private GameObject scrollBoxContent;
    
    [SerializeField]
    private GameObject entryPrefab;
    
    [SerializeField]
    private GameObject importButton;

    // Start is called once before the first execution of Update after the MonoBehaviour is created
    void Start()
    {
        GameController.main.FindAnchor();

        string[] fileNames = 
        Directory.GetFiles(Path.Combine(GameController.rootPath,  "cpp/sketches_SDL/Molecular/data/"))
        .Select(x => x.Split("/").Last()).ToArray();
        
        for (int i = 0; i < fileNames.Length; i++)
        {
            if(!fileNames[i].EndsWith(".fgo")) {
                continue;
            }
            GameObject entry = Instantiate(entryPrefab, scrollBoxContent.transform);
            string iName = fileNames[i];
            entry.GetComponentInChildren<TextMeshProUGUI>().SetText(iName);
            entry.GetComponent<Button>().onClick.AddListener(() => ProcessButtonPress(iName));
            //entry.GetComponent<RectTransform>().pivot = new Vector2(1.18f, -3.62f + (0.5f * i));
        }
        
        // Add listener to the import button
        if (importButton != null)
        {
            importButton.GetComponent<Button>().onClick.AddListener(OnImportButtonPressed);
        }
    }

    public void ProcessButtonPress(string name) {
        GetComponent<GameController>().StartSimulation(Path.Combine(GameController.rootPath, "cpp/sketches_SDL/Molecular/data/", name));
        scrollBoxContent.transform.parent.parent.gameObject.SetActive(false);
    }

    private void OnImportButtonPressed()
    {
        // Create extension filter for .fgo files
        var extensions = new ExtensionFilter[] {
            new("FGO Files", "fgo"),
            new("All Files", "*")
        };

        // Open native file dialog with .fgo filter
        var paths = StandaloneFileBrowser.OpenFilePanel("Select .fgo file", "", extensions, false);
        
        if (paths != null && paths.Length > 0 && !string.IsNullOrEmpty(paths[0]))
        {
            string selectedFilePath = paths[0];
            
            // Call StartSimulation with the selected file
            GameController.main.StartSimulation(selectedFilePath);
            
            // Hide the picker UI
            scrollBoxContent.transform.parent.parent.gameObject.SetActive(false);
        }
    }
}
