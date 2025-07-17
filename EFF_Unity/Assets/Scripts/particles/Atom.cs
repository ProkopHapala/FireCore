using NUnit.Framework;
using TMPro;
using Unity.VisualScripting;
using UnityEngine;
using UnityEngine.InputSystem.iOS;
using UnityEngine.Networking;

public class Atom : MonoBehaviour, IParticle
{
    private static int atomCount = 0;
    private Vector3 _position;
    public Vector3 Position {get => _position; set {transform.position = value; _position = value;}}

    [SerializeField]
    private GameObject sphere;

    [SerializeField]
    private GameObject spritePrefab;

    private GameObject sprite;

    private Transform canvas;

    private Outline outline;

    private RenderingStyle _style;
    public RenderingStyle Style {get => _style; 
    set {
        _style = value;
        // sphere.SetActive(value == RenderingStyle.SPHERICAL);
        sphere.SetActive(value == RenderingStyle.SPHERICAL);
        sprite.SetActive(value == RenderingStyle.FOG_AND_POINT);
        Debug.Log("Setting rendering style to: " + value);
        switch(value) {
            case RenderingStyle.SPHERICAL:
                _SetOutline = _SetOutlineSphere;
                break;
            case RenderingStyle.FOG_AND_POINT:
                _SetOutline = _SetOutlineSprite;
                break;
            default: break;
        }
    }}


    private int _id = -1;
    public int Id { get => _id; set {_id = _id == -1 ? value : _id;} }

    public static Atom CreateNew() {
        Atom atom = Instantiate(GameController.main.atomPrefab, Vector3.zero, Quaternion.identity).GetComponent<Atom>();
        atom.outline = atom.sphere.GetComponent<Outline>();
        atom.canvas = Camera.main.transform.GetChild(0);
        atom.sprite = Instantiate(atom.spritePrefab, atom.canvas.transform);
        atom.Style = RenderingStyle.FOG_AND_POINT;
        atom.Id = atomCount++;
        atom.sprite.GetComponentInChildren<TextMeshProUGUI>().SetText(atom.Id.ToString());
        GameController.main.atoms[atom.Id] = atom;
        //atom.sprite.SetActive(atom._style == RenderingStyle.FOG_AND_POINT);
        return atom;
    }

    // void Start() {
    //     outline = sphere.GetComponent<Outline>();
    //     canvas = Camera.main.transform.GetChild(0);
    //     sprite = Instantiate(spritePrefab, canvas.transform);
    //     Style = RenderingStyle.FOG_AND_POINT;
    //     sprite.SetActive(_style == RenderingStyle.FOG_AND_POINT);
    // }

    public void UpdateSpritePositions() {
        if(_style == RenderingStyle.SPHERICAL) {
            return;
        }

        // RectTransformUtility.ScreenPointToLocalPointInRectangle((RectTransform)GameController.main.canvas.transform, Camera.main.WorldToScreenPoint(transform.position), null, out var CanvasPoint);
        // sprite.GetComponent<RectTransform>().position = CanvasPoint;

        // RectTransform rectTransform = sprite.GetComponent<RectTransform>();
        // Vector2 viewportPosition = Camera.main.WorldToViewportPoint(transform.position);
        // sprite.GetComponent<RectTransform>().position = new Vector2(viewportPosition.x * 1920, viewportPosition.y * 1080);

        sprite.GetComponent<RectTransform>().position = Camera.main.WorldToScreenPoint(transform.position);
    }

    private delegate void outlineSetter(bool value, Outline outline);
    private outlineSetter _SetOutline = _SetOutlineSphere;
    private static void _SetOutlineSphere(bool value, Outline outline) {
        outline.OutlineMode = value ? Outline.Mode.OutlineAndSilhouette : Outline.Mode.OutlineHidden;
    }
    private static void _SetOutlineSprite(bool value, Outline outline) {
        
    }
    public void SetOutline(bool value)
    {
        _SetOutline(value, outline);
    }

    public override string ToString()
    {
        return $"Atom: (id: {_id}, pos: {_position}, style: {Style})";
    }
}
