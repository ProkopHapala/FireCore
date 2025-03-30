using UnityEditor;
using UnityEngine;

public class CameraControl : MonoBehaviour
{
    // Update is called once per frame

    private new Camera camera;
    private Vector3 pivotPoint;
    public float rotationSpeed = 100;
    public float zoomSpeed = 3;
    public float moveSpeed = 1;

    void Awake() {
        camera = GetComponent<Camera>();
        pivotPoint = Vector3.zero;
        //transform.Rotate(new Vector3(1.95f, 0, 0), Space.Self);
    }
    void Update()
    {
        if(Input.GetKey(KeyCode.LeftShift)) {
            Vector3 moveVector = Vector3.zero;
            if(Input.GetKey(KeyCode.RightArrow)) {
                moveVector += new Vector3(-transform.position.z, 0, transform.position.x) * (moveSpeed * Time.deltaTime);
            }
            else if(Input.GetKey(KeyCode.LeftArrow)) {
                moveVector += new Vector3(transform.position.z, 0, -transform.position.x) * (moveSpeed * Time.deltaTime);
            }
            if(Input.GetKey(KeyCode.UpArrow)) {
                moveVector += new Vector3(-transform.position.x, 0, -transform.position.z) * (moveSpeed * Time.deltaTime);              
            }
            else if(Input.GetKey(KeyCode.DownArrow)) {
                moveVector += new Vector3(transform.position.x, 0, transform.position.z) * (moveSpeed * Time.deltaTime);
            }
            // * 30 because that is the distance of the camera from origin -> speed will be synced with horizontal movement
            if(Input.GetKey(KeyCode.Comma)) {
                moveVector += new Vector3(0, -30, 0) * (moveSpeed * Time.deltaTime);
            }
            else if(Input.GetKey(KeyCode.Period)) {
                moveVector += new Vector3(0, 30, 0) * (moveSpeed * Time.deltaTime);
            }

            pivotPoint += moveVector;
            transform.position += moveVector;
            return;
        }

        if(Input.GetKey(KeyCode.RightArrow)) {
            transform.RotateAround(pivotPoint, Vector3.up, -rotationSpeed * Time.deltaTime);
        }
        else if(Input.GetKey(KeyCode.LeftArrow)) {
            transform.RotateAround(pivotPoint, Vector3.up, rotationSpeed * Time.deltaTime);
        }
        if(Input.GetKey(KeyCode.UpArrow)) {
            transform.RotateAround(pivotPoint, new Vector3(transform.position.z, 0, -transform.position.x), -rotationSpeed * Time.deltaTime);
        }
        else if(Input.GetKey(KeyCode.DownArrow)) {
            transform.RotateAround(pivotPoint, new Vector3(transform.position.z, 0, -transform.position.x), rotationSpeed * Time.deltaTime);
        }
        if(Input.GetKey(KeyCode.Comma)) {
            camera.fieldOfView += zoomSpeed * Time.deltaTime;
        }
        else if(Input.GetKey(KeyCode.Period)) {
            camera.fieldOfView -= zoomSpeed * Time.deltaTime;
        }

        
    }

    public void SetPivot(Vector3 newPivot) {
        transform.position += newPivot - pivotPoint;
        pivotPoint = newPivot;
    }
}
