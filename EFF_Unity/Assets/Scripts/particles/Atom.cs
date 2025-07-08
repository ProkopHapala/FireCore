using Unity.VisualScripting;
using UnityEngine;
using UnityEngine.InputSystem.iOS;
using UnityEngine.Networking;

public class Atom : MonoBehaviour, IParticle
{
    private Vector3 _position;
    public Vector3 Position {get => _position; set {transform.position = value; _position = value;}}

    [SerializeField]
    private GameObject sphere;

    [SerializeField]
    private GameObject sprite;

    private Outline outline;

    private RenderingStyle _style;
    public RenderingStyle Style {get => _style; 
    set {
        _style = value;
        sphere.SetActive(value == RenderingStyle.SPHERICAL);
        sprite.SetActive(value == RenderingStyle.FOG_AND_POINT);
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
        return Instantiate(GameController.main.atomPrefab, Vector3.zero, Quaternion.identity).GetComponent<Atom>();
    }

    void Start() {
        outline = sphere.GetComponent<Outline>();

    }

    private delegate void outlineSetter(bool value);
    private outlineSetter _SetOutline;
    private void _SetOutlineSphere(bool value) {
        outline.OutlineMode = value ? Outline.Mode.OutlineAndSilhouette : Outline.Mode.OutlineHidden;
    }
    private void _SetOutlineSprite(bool value) {
        
    }
    public void SetOutline(bool value)
    {

    }

    public override string ToString()
    {
        return $"Atom: (id: {_id}, pos: {_position}, style: {Style})";
    }
}
