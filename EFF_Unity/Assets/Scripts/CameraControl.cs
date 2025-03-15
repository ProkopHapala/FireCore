using UnityEngine;

public class CameraControl : MonoBehaviour
{
    // Update is called once per frame

    public float rotationSpeed = 100;
    void Update()
    {
        if(Input.GetKey(KeyCode.RightArrow)) {
            transform.RotateAround(Vector3.zero, Vector3.up, -rotationSpeed * Time.deltaTime);
            return;
        }
        if(Input.GetKey(KeyCode.LeftArrow)) {
            transform.RotateAround(Vector3.zero, Vector3.up, rotationSpeed * Time.deltaTime);
            return;
        }
        if(Input.GetKey(KeyCode.UpArrow)) {
            transform.RotateAround(Vector3.zero, new Vector3(transform.position.z, 0, -transform.position.x), -rotationSpeed * Time.deltaTime);
            return;
        }
        if(Input.GetKey(KeyCode.DownArrow)) {
            transform.RotateAround(Vector3.zero, new Vector3(transform.position.z, 0, -transform.position.x), rotationSpeed * Time.deltaTime);
            return;
        }
    }
}
