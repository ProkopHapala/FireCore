using TMPro;
using UnityEngine;
using UnityEngine.InputSystem;
using UnityEngine.UIElements;

public class InputTracker : MonoBehaviour
{

    public CppConnector connector;
    public TextMeshProUGUI statusText;

    private const string RUNNING_MSG = "RUNNING";
    private const string PAUSED_MSG = "PAUSED";

        // Start is called once before the first execution of Update after the MonoBehaviour is created
    void Start()
    {
        statusText.SetText(PAUSED_MSG);
        // this.RegisterCallback<KeyDownEvent>(OnKeyDown, TrickleDown.TrickleDown);
    }

    // Update is called once per frame
    void Update()
    {
        if(Input.GetKeyDown(KeyCode.Space)) {
            connector.isRunning = !connector.isRunning;
            statusText.SetText(connector.isRunning ? RUNNING_MSG : PAUSED_MSG);
        }
    }

    // void OnKeyDown(KeyDownEvent ev) {
    //     if(ev.keyCode == KeyCode.Space) {
    //         connector.isRunning = !connector.isRunning;
    //         statusText.SetText(connector.isRunning ? RUNNING_MSG : PAUSED_MSG);
    //     }
    // }
}
