using System;
using System.Linq.Expressions;
using UnityEditor;
using UnityEngine;

public class CameraControl : MonoBehaviour
{
    // Update is called once per frame

    private new Camera camera;
    private Vector3 pivotPoint;
    [SerializeField] private float rotationSpeed;
    [SerializeField] private float zoomSpeed;
    [SerializeField] private float moveSpeed;

    private Vector2 rotation = new Vector2(90, 0);
    private float distanceToPivot = 30;

    private GameObject trackingObject;
    private bool isFollowing = false;

    void Awake() {
        camera = GetComponent<Camera>();
        pivotPoint = Vector3.zero;
        //transform.Rotate(new Vector3(1.95f, 0, 0), Space.Self);
    }
    void Update()
    {
        if(isFollowing && GameController.main.isRunning) {
            SetPivot(trackingObject.transform.position);
        }

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
            if(Input.GetKey(KeyCode.K)) {
                moveVector += new Vector3(0, -30, 0) * (moveSpeed * Time.deltaTime);
            }
            else if(Input.GetKey(KeyCode.L)) {
                moveVector += new Vector3(0, 30, 0) * (moveSpeed * Time.deltaTime);
            }

            pivotPoint += moveVector;
            transform.position += moveVector;
            return;
        }

        // if(Input.GetKey(KeyCode.RightArrow)) {
        //     transform.RotateAround(pivotPoint, Vector3.up, -rotationSpeed * Time.deltaTime);
        // }
        // else if(Input.GetKey(KeyCode.LeftArrow)) {
        //     transform.RotateAround(pivotPoint, Vector3.up, rotationSpeed * Time.deltaTime);
        // }
        // if(Input.GetKey(KeyCode.UpArrow)) {
            
        //     transform.RotateAround(pivotPoint, new Vector3(transform.position.z, 0, -transform.position.x), -rotationSpeed * Time.deltaTime);
        // }
        // else if(Input.GetKey(KeyCode.DownArrow)) {
        //     transform.RotateAround(pivotPoint, new Vector3(transform.position.z, 0, -transform.position.x), rotationSpeed * Time.deltaTime);
        // }
        // if(Input.GetKey(KeyCode.K)) {
        //     camera.fieldOfView += zoomSpeed * Time.deltaTime;
        // }
        // else if(Input.GetKey(KeyCode.L)) {
        //     camera.fieldOfView -= zoomSpeed * Time.deltaTime;
        // }

        if(Input.GetKey(KeyCode.RightArrow)) {
            rotation.x += rotationSpeed * Time.deltaTime;
        }
        else if(Input.GetKey(KeyCode.LeftArrow)) {
            rotation.x -= rotationSpeed * Time.deltaTime;
        }
        if(Input.GetKey(KeyCode.UpArrow)) {
            rotation.y += rotationSpeed * Time.deltaTime;
        }
        else if(Input.GetKey(KeyCode.DownArrow)) {
            rotation.y -= rotationSpeed * Time.deltaTime;
        }
        if(Input.GetKey(KeyCode.K)) {
            camera.fieldOfView += zoomSpeed * Time.deltaTime;
        }
        else if(Input.GetKey(KeyCode.L)) {
            camera.fieldOfView -= zoomSpeed * Time.deltaTime;
        }

        UpdateRotation();
    }

    private void UpdateRotation() {
        // transform.rotation = Quaternion.Euler(rotation.x, rotation.y, 0);
        //float distanceToPivot = Vector3.Distance(transform.position, pivotPoint);
        transform.position = pivotPoint + new Vector3(
            (distanceToPivot * (float)Math.Cos(rotation.x)) * (float)Math.Cos(rotation.y), 
            distanceToPivot * (float)Math.Sin(rotation.y), 
            (distanceToPivot * (float)Math.Sin(rotation.x)) * (float)Math.Cos(rotation.y));

        transform.LookAt(pivotPoint);
    }

    public void FollowParticle(GameObject particleObject) {
        isFollowing = true;
        trackingObject = particleObject;
    }

    public void StopFollowing() {
        isFollowing = false;
    }

    public void SetPivot(Vector3 newPivot) {
        transform.position += newPivot - pivotPoint;
        pivotPoint = newPivot;
        distanceToPivot = Vector3.Distance(transform.position, pivotPoint);
    }
}
