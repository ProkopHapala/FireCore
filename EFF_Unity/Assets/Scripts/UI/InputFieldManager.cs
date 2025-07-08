using System;
using System.Linq;
using TMPro;
using Unity.VisualScripting;
using UnityEngine;
using UnityEngine.Networking;
using UnityEngine.Rendering.Universal;

public class InputFieldManager : MonoBehaviour
{
    public GameObject[] fieldObjects;
    private TMP_InputField[] fields;
    private IParticle particle;
    private ObjectType particleType;

    void Awake() {
        fields = Array.ConvertAll(fieldObjects, o => o.GetComponent<TMP_InputField>());
        for (int i = 0; i < fields.Length; i++)
        {
            int fieldIndex = i; // Capture the current value of i
            fields[i].onSubmit.AddListener(str => SetPositions(str, fieldIndex));
        }
    }

    public void SetDefault(Vector3 position) {
        fields[0].text = position.x.ToString();
        fields[1].text = position.y.ToString();
        fields[2].text = position.z.ToString();
    }

    public void SetDefault(Vector3 position, float size) {
        fields[0].text = position.x.ToString();
        fields[1].text = position.y.ToString();
        fields[2].text = position.z.ToString();
        fields[3].text = size.ToString();

    }

    public void SetPositions(string input, int fieldId) {
        try {
            float value = float.Parse(input);
            Vector3 pos = particle.Position;
            float size = particleType == ObjectType.ELECTRON ? ((Electron)particle).Size : 0;
            
            // Update the appropriate field
            switch(fieldId) {
                case 0:
                    pos.x = value;
                    break;
                case 1:
                    pos.y = value;
                    break;
                case 2:
                    pos.z = value;
                    break;
                case 3:
                    size = value;
                    break;
                default:
                    UnityEngine.Debug.LogWarning($"Unknown field ID: {fieldId}");
                    return;
            }

            // Apply the changes
            UnityEngine.Debug.Log($"Setting {particleType} {particle.Id} position to {pos}, size: {size}");
            GameController.main.SetParticlePosition(particleType, particle, pos, size);
        }
        catch (Exception e) {
            UnityEngine.Debug.LogError($"Error in SetPositions: {e.Message}\n{e.StackTrace}");
        }
    }

    public void Spawn(ObjectType type, IParticle particle, Vector2 posUI) {
        this.particleType = type;
        this.particle = particle;
        GetComponent<RectTransform>().anchoredPosition = posUI;

        Array.ForEach(fieldObjects, o => o.SetActive(false));

        if(type == ObjectType.ATOM) {
            foreach(var field in fieldObjects.Take(3)) {
                field.SetActive(true);
            }
        }
        else {
            foreach(var field in fieldObjects) {
                field.SetActive(true);
            }
        }
    }
    public void Despawn() {
        foreach(var field in fieldObjects) {
            field.SetActive(false);
        }
    }
}
