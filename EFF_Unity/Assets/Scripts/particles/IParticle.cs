

using UnityEngine;

public interface IParticle {
    int Id {get; set;}
    Vector3 Position {get; set;}
    RenderingStyle Style {get; set;}
    void SetOutline(bool value);
    
}