using System;
using System.Collections.Generic;
using NUnit.Framework;
using UnityEngine;
using UnityEngine.EventSystems;
using UnityEngine.InputSystem;

public class InfoButton : MonoBehaviour, IPointerDownHandler
{
    private List<Action> leftClickListeners = new();
    private List<Action> rightClickListeners = new();

    void IPointerDownHandler.OnPointerDown(PointerEventData eventData)
    {
        if (eventData.button == PointerEventData.InputButton.Left)
        {
            leftClickListeners.ForEach(a => a());
        }
        else if (eventData.button == PointerEventData.InputButton.Right)
        {
            rightClickListeners.ForEach(a => a());
        }
    }

    public void AddListener(MouseButtons button, Action action) {
        if(button == MouseButtons.LEFT) {
            leftClickListeners.Add(action);
        }
        else {
            rightClickListeners.Add(action);
        }
    }


}
public enum MouseButtons {
    RIGHT,
    LEFT
}