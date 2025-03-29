using System;
using System.Security;
using TMPro;
using UnityEngine;
using UnityEngine.UI;

public class InfoBox : MonoBehaviour
{
    public GameObject button;
    public TextMeshProUGUI objText;
    public TextMeshProUGUI posText;
    public TextMeshProUGUI sizeText;

    private CppConnector connector;
    private delegate void del();
    private del SetStats;
    
    // Start is called once before the first execution of Update after the MonoBehaviour is created
    void Awake() {
        button.GetComponent<Button>().onClick.AddListener(OnButtonClick);
    }

    // Update is called once per frame
    void Update()
    {
        if(!connector.isRunning) {
            return;
        }
        // var pos = connector.positions[id];
        // if(type == ObjectType.ATOM) {
        //     posText.SetText($"pos: x {pos.x} | y {pos.y} | z {pos.z}");
        // }
        // else {
        //     posText.SetText($"pos: x {pos.x} | y {pos.y} | z {pos.z}");
        // }

        // sizeText.SetText($"siz: {connector.sizes[id]}");

        SetStats();
    }

    // id: to be displayed
    // posId: index in positions to read from
    public void SetConnector(CppConnector connector, int id, int posId, ObjectType type) {
        objText.SetText($"{(type == ObjectType.ATOM ? "ATOM" : "ELECTRON")} {id}");
        this.connector = connector;

        if(type == ObjectType.ELECTRON){
            SetStats = () => {
                var pos = connector.positions[posId];
                posText.SetText($"pos: x {pos.x:F3} | y {pos.y:F3} | z {pos.z:F3}");
                sizeText.SetText($"siz: {connector.sizes[id]:F5}");
            };
        }
        else {
            SetStats = () => {
                var pos = connector.positions[posId];
                posText.SetText($"pos: x {pos.x:F3} | y {pos.y:F3} | z {pos.z:F3}");
            };
        }
        gameObject.SetActive(true);
    }

    public void ForceUpdate() {
        SetStats?.Invoke();
    }

    public void OnButtonClick() {
        
    }
}
