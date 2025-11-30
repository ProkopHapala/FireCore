using System;
using System.Linq;
using System.Security;
using NUnit.Framework;
using TMPro;
using Unity.PlasticSCM.Editor.WebApi;
using UnityEngine;
using UnityEngine.EventSystems;
using UnityEngine.UI;

public class InfoBox : MonoBehaviour {
    public GameObject button;
    public TextMeshProUGUI objText;
    public TextMeshProUGUI posText;
    public TextMeshProUGUI sizeText;

    //public GameObject inputFieldAnchor;

    // inputFields can't be assigned directly from the editor, idk why
    //public GameObject inputFieldAnchor;
    public InputFieldManager inputFields;

    //private GameController connector;
    private int id;
    private IParticle particle;
    // private Outline particleOutline;
    private CameraControl cameraControl;
    private delegate void del();
    private del SetStats;
    private bool cursorOn;
    public bool IsFixed = false;

    public ObjectType Type {get; private set;}
    
    // Start is called once before the first execution of Update after the MonoBehaviour is created
    void Awake() {
        button.GetComponent<InfoButton>().AddListener(MouseButtons.LEFT, () => OnButtonClick());
        button.GetComponent<InfoButton>().AddListener(MouseButtons.RIGHT, () => OnButtonRightClick());
        cameraControl = Camera.main.GetComponent<CameraControl>();

        EventTrigger trigger = button.GetComponent<EventTrigger>();
        EventTrigger.Entry entry = new EventTrigger.Entry
        {
            eventID = EventTriggerType.PointerEnter
        };
        entry.callback.AddListener( eventData => OnPointerEnter((PointerEventData)eventData) );
        trigger.triggers.Add(entry);

        entry = new EventTrigger.Entry
        {
            eventID = EventTriggerType.PointerExit
        };
        entry.callback.AddListener( eventData => OnPointerExit((PointerEventData)eventData) );
        trigger.triggers.Add(entry);

    }


    // Update is called once per frame
    void Update()
    {
        if(cursorOn){
            if(Input.GetKeyDown(KeyCode.S)) {
                if(IsFixed) {
                    if(Type == ObjectType.ATOM){
                        SimulationManager.unfixAtom(id);
                        Debug.Log("Unfixed atom " + id);
                    }
                    else {
                        SimulationManager.unfixElectron(id);
                        Debug.Log("Unfixed electron " + id);
                    }
                    objText.SetText(objText.text[..^4]);
                    IsFixed = false;
                }
                else {
                    if(Type == ObjectType.ATOM){
                        SimulationManager.fixAtom(id);
                        Debug.Log("Fixed atom " + id);
                    }
                    else {
                        SimulationManager.fixElectron(id);
                        Debug.Log("fixed electron " + id);
                    }
                    objText.SetText(objText.text + " [S]");
                    IsFixed = true;
                }
            }

            if(Input.GetKeyDown(KeyCode.F) && !Input.GetKey(KeyCode.LeftShift)) {
                GameController.main.cameraControl.FollowParticle(((MonoBehaviour)particle).gameObject);
            }

            if(!GameController.main.isRunning) {
                return;
            }
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
    public void SetConnector(int id, ObjectType type) {
        inputFields = GameController.main.inputFields;

        objText.SetText($"{(type == ObjectType.ATOM ? "ATOM" : "ELECTRON")} {id} {(type == ObjectType.ELECTRON ? (GameController.main.electrons[id].Spin == 1 ? "(+1/2)" : "(-1/2)") : "")}");
        this.id = id;
        Type = type;
        particle = GameController.main.GetParticle(id, type);
        if(type == ObjectType.ELECTRON){
            
            SetStats = () => {
                var pos = particle.Position;
                posText.SetText($"pos: x {pos.x:F3} | y {pos.y:F3} | z {pos.z:F3}");
                sizeText.SetText($"siz: {(particle as Electron).Size:F5}");
            };
        }
        else {
            SetStats = () => {
                var pos = particle.Position;
                posText.SetText($"pos: x {pos.x:F3} | y {pos.y:F3} | z {pos.z:F3}");
            };
        }
        gameObject.SetActive(true);
    }



    public void ForceUpdate() {
        SetStats?.Invoke();
    }

    public void OnButtonClick() {
        cameraControl.SetPivot(particle.Position);
    }
    public void OnButtonRightClick() {
        // inputFieldAnchor.GetComponent<RectTransform>().anchoredPosition = GetComponent<RectTransform>().anchoredPosition + new Vector2(-50, 0);
        // inputFieldAnchor.SetActive(true);
        if(Type == ObjectType.ATOM){
            inputFields.SetDefault(particle.Position);
        }
        else {
            inputFields.SetDefault(particle.Position, (particle as Electron).Size);
        }

        // Vector3 screenPoint = Camera.main.WorldToScreenPoint(transform.position);
        inputFields.Spawn(Type, particle, new Vector2(transform.position.x, transform.position.y) + new Vector2(350, -550));

    }

    public void OnPointerEnter(PointerEventData eventData)
    {
        particle.SetOutline(true);
        cursorOn = true;
    }

    public void OnPointerExit(PointerEventData eventData)
    {
        particle.SetOutline(false);
        cursorOn = false;
    }


}
