using System;
using System.Security;
using TMPro;
using UnityEngine;
using UnityEngine.EventSystems;
using UnityEngine.UI;

public class InfoBox : MonoBehaviour
{
    public GameObject button;
    public TextMeshProUGUI objText;
    public TextMeshProUGUI posText;
    public TextMeshProUGUI sizeText;

    //public GameObject inputFieldAnchor;

    // inputFields can't be assigned directly from the editor, idk why
    //public GameObject inputFieldAnchor;
    public InputFieldManager inputFields;

    //private GameController connector;
    private int posId;
    private int id;
    private GameObject particle;
    private Outline particleOutline;
    private CameraControl cameraControl;
    private delegate void del();
    private del SetStats;
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
        if(!GameController.main.isRunning) {
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
    public void SetConnector(int id, int posId, ObjectType type) {
        inputFields = GameController.main.inputFields;

        objText.SetText($"{(type == ObjectType.ATOM ? "ATOM" : "ELECTRON")} {id}");
        this.posId = posId;
        this.id = id;
        Type = type;
        particle = GameController.main.particles[posId];
        particleOutline = particle.GetComponent<Outline>();

        if(type == ObjectType.ELECTRON){
            SetStats = () => {
                var pos = GameController.main.positions[posId];
                posText.SetText($"pos: x {pos.x:F3} | y {pos.y:F3} | z {pos.z:F3}");
                sizeText.SetText($"siz: {GameController.main.sizes[id]:F5}");
            };
        }
        else {
            SetStats = () => {
                var pos = GameController.main.positions[posId];
                posText.SetText($"pos: x {pos.x:F3} | y {pos.y:F3} | z {pos.z:F3}");
            };
        }
        gameObject.SetActive(true);
    }



    public void ForceUpdate() {
        SetStats?.Invoke();
    }

    public void OnButtonClick() {
        cameraControl.SetPivot(particle.transform.position);
    }
    public void OnButtonRightClick() {
        // inputFieldAnchor.GetComponent<RectTransform>().anchoredPosition = GetComponent<RectTransform>().anchoredPosition + new Vector2(-50, 0);
        // inputFieldAnchor.SetActive(true);
        if(Type == ObjectType.ATOM){
            inputFields.SetDefault(GameController.main.positions[posId]);
        }
        else {
            inputFields.SetDefault(GameController.main.positions[posId], GameController.main.sizes[posId]);
        }
        inputFields.Spawn(Type, id, GetComponent<RectTransform>().anchoredPosition + new Vector2(350, 0));

    }

    public void OnPointerEnter(PointerEventData eventData)
    {
        particleOutline.OutlineMode = Outline.Mode.OutlineAndSilhouette;
    }

    public void OnPointerExit(PointerEventData eventData)
    {
        particleOutline.OutlineMode = Outline.Mode.OutlineHidden;
    }


}
