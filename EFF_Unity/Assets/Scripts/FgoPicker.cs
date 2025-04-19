using System.IO;
using System.Linq;
using TMPro;
using UnityEngine;
using UnityEngine.UI;

public class FgoPicker : MonoBehaviour
{
    public GameObject scrollBoxContent;
    public GameObject entryPrefab;

    // Start is called once before the first execution of Update after the MonoBehaviour is created
    void Start()
    {
        string[] fileNames = 
        Directory.GetFiles("/home/perry/FireCore/cpp/sketches_SDL/Molecular/data/")
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
        
    }

    public void ProcessButtonPress(string name) {
        GetComponent<GameController>().StartSimulation(name);
        scrollBoxContent.transform.parent.parent.gameObject.SetActive(false);
    }

    // Update is called once per frame
    void Update()
    {
        
    }

    
}
