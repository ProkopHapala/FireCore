// using System.Numerics;

using System;
using System.Drawing;
using System.Threading;
using UnityEngine;
using UnityEngine.UIElements;

public class Electron : MonoBehaviour, IParticle
{
    [SerializeField] private GameObject sphere;
    [SerializeField] private GameObject fog;
    [SerializeField] private ParticleSystem fogPS;

    private Outline outline;
    
    private static int electronCount = 0;
  
    //public bool UseFancyRendering {get; private set;} = true;
    public int Spin {get; set;}
    private Vector3 _position;
    public Vector3 Position {get => _position; set {transform.position = value; _position = value;}}
    private RenderingStyle _style;
    public RenderingStyle Style {get => _style; 
    set {
        _style = value;
        sphere.SetActive(value == RenderingStyle.SPHERICAL);
        fog.SetActive(value == RenderingStyle.FOG_AND_POINT);
        switch(value) {
            case RenderingStyle.SPHERICAL:
                _SetSize = SetSizeSphere;
                _SetOutline = SetOutlineSphere;
                break;
            case RenderingStyle.FOG_AND_POINT:
                _SetSize = SetSizeFog;
                _SetOutline = SetOutlineFog;
                break;
            default: break;
        }
    }}
    private float _size;
    public float Size {get => _size; set {_SetSize(value, this); _size = value;}}

    private int _id = -1;
    public int Id { get => _id; set {_id = _id == -1 ? value : _id;} }

    private delegate void sizeSetter(float size, Electron particle);
    

    private static void SetSizeSphere(float size, Electron particle) { particle.sphere.transform.localScale = new Vector3(size, size, size); }
    private static void SetSizeFog(float size, Electron particle) {
        ParticleSystem.MainModule main = particle.fogPS.main;
        ParticleSystem.ShapeModule shape = particle.fogPS.shape;
        main.startSize = size;
        shape.radius = size / 2;
        
        //ps.startSpeed = new ParticleSystem.MinMaxCurve(0, 2 * (size / 10));
        }
        
    private sizeSetter _SetSize = SetSizeFog;

    void Start() {
        outline = sphere.GetComponent<Outline>();
    }
    public static Electron CreateNew(int spin) {
        Electron inst;
        if(spin == 1) {
            inst = Instantiate(GameController.main.electronPrefabPlusSpin, Vector3.zero, Quaternion.identity).GetComponent<Electron>();
        }
        else {
            inst = Instantiate(GameController.main.electronPrefabMinusSpin, Vector3.zero, Quaternion.identity).GetComponent<Electron>();
        }
        
        inst.Id = electronCount++;
        inst.Spin = spin;
        inst.Style = RenderingStyle.FOG_AND_POINT;
        return inst;
    }



    private delegate void outlineSetter(bool value, Outline outline);

    private outlineSetter _SetOutline = SetOutlineFog;

    private static void SetOutlineSphere(bool value, Outline outline) {
        outline.OutlineMode = value ? Outline.Mode.OutlineAndSilhouette : Outline.Mode.OutlineHidden;
    }
    private static void SetOutlineFog(bool value, Outline outline) {
        return;
    }

    public void SetOutline(bool value) { _SetOutline(value, outline); }


    public override string ToString()
    {
        return $"Electron: (id: {_id}, pos: {_position}, size: {_size}, style: {Style})";
    }
    
    public static void Clear()
    {
        electronCount = 0;
        GameController.main.atoms = Array.Empty<Atom>();
    }
}
