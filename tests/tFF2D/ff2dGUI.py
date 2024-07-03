#!/usr/bin/python
import sys
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.widgets import Button,Slider

import numpy as np

sys.path.append("../../")
from pyBall import FF2D as ff

# Initialize the figure and axis
fig, ax = plt.subplots()
ax.set_title('LMB:select,RMB:add(atom,connect by bond), LMB+Ctrl:Delete')
ax.set_xlim(0, 10)
ax.set_ylim(0, 10)

# List to store points and lines
points = []
lines  = []

# Variable to store the index of the point being dragged
dragging_point_index = None

icur_atom = -1

# Function to update the plot
def update_plot( bLabels=True):
    global icur_atom
    ax.clear()
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 10)
    ax.set_title('Interactive Plot: Add, Delete, Connect, and Move Points')

    #print( "----- py.update_plot() " )
    types, apos, neighs = ff.getAtoms()  #;print("apos ", apos) #;print("types", types) ;print("neighs", ne
    bonds               = ff.getBonds()  #;print("bonds ", bonds)

    # Plot points
    #x_values, y_values = zip(*points)
    #ax.plot(x_values, y_values, 'bo', markersize=5)  # Plot points
    ax.plot( apos[:,0], apos[:,1], 'bo', markersize=5 )
    if icur_atom >= 0:
        ax.plot( [apos[icur_atom,0]], [apos[icur_atom,1]], 'ro', markersize=5 )

    if bLabels:
        for i, (x, y) in enumerate(apos):
            ax.text(x, y, str(i), fontsize=10, ha='right')

    # Plot lines using LineCollection
    lines = [ (apos[i,:], apos[j,:]) for (i,j) in bonds ] 
    line_segments = LineCollection(lines, colors='b')
    ax.add_collection(line_segments)
    ax.set_aspect('equal', adjustable='datalim')
    fig.canvas.draw()

# Event handler function for mouse button press
def on_click(event):
    global icur_atom
    global dragging_point_index
    bUpdate = False
    if event.inaxes is not None:
        x = event.xdata
        y = event.ydata
        ia = ff.findAtomAt( x, y )
        if event.button == 1:  # Left mouse button  - Picking atoms 
            if event.key == 'control':  # Ctrl + Left Click to delete point
                #print( "CTRL+LMB: remove atom() ", ia )
                ff.removeAtom(ia)
                bUpdate = True
            else:  # Left Click to add point
                if icur_atom != ia: 
                    icur_atom = ia
                    bUpdate = True
        elif event.button == 3:  # Right mouse button to connect points
            ja = ia
            if ia<0:   # new atom
                ja = ff.addAtom( x, y )
                bUpdate   = True
            if icur_atom >= 0:   # bond to icur_atom
                ff.addBond(icur_atom, ja)
                bUpdate = True
            icur_atom = ja
        #elif event.button == 2:  # Middle mouse button to start dragging
        #    for (i, (x, y)) in enumerate(points):
        #        if abs(x - event.xdata) < 0.1 and abs(y - event.ydata) < 0.1:
        #            dragging_point_index = i
        #            break
        #ff.print_atoms()
        if bUpdate:
            update_plot()


key_actions = {
    #'r':         lambda: print("Key 'a' pressed"),
    #'d':         lambda: print("Key 'd' pressed"),
    'enter':     lambda: ( print("run"),                     ff.run(n=1000, dt=0.1, damp=0.1, Fconv=1e-2, bCleanForce=True ), update_plot() ),
    'delete':    lambda: ( print("delete atom ", icur_atom), ff.removeAtom(icur_atom),                                        update_plot() ),
    #'backspace': lambda: print("Backspace key pressed"),
    #'home':      lambda: print("Home key pressed"),
    # Add more key events as needed
}
def on_key(event):
    action = key_actions.get(event.key)
    if action:
        action()

# def on_key(event):
#     if event.key == 'r':
#         print("Key 'a' pressed => Relax")
#         ff.run(n=1000, dt=0.1, damp=0.1, Fconv=1e-2, bCleanForce=True )
#         update_plot()
#     elif event.key == 'delete':
#       print("Key 'd' pressed")

# Event handler function for mouse motion
def on_motion(event):
    global dragging_point_index
    if dragging_point_index is not None and event.inaxes is not None:
        points[dragging_point_index] = (event.xdata, event.ydata)
        lines[:] = [(points[i], points[i + 1]) for i in range(len(points) - 1)]
        update_plot()

# Event handler function for mouse button release
def on_release(event):
    global dragging_point_index
    dragging_point_index = None

# Connect the event handlers to the figure
fig.canvas.mpl_connect('button_press_event',   on_click   )
fig.canvas.mpl_connect('motion_notify_event',  on_motion  )
fig.canvas.mpl_connect('button_release_event', on_release )
fig.canvas.mpl_connect('key_press_event', on_key)



plt.subplots_adjust(bottom=0.2)  # Adjust the subplot to make room for widgets

# ===========  Adding Button
def on_button_clicked(event):
    print("Button clicked!")

# Add a button
ax_button = plt.axes([0.81, 0.05, 0.1, 0.075])  # Position [left, bottom, width, height]
button = Button(ax_button, 'Click Me')
button.on_clicked(on_button_clicked)

# ===========  Adding Slider

def on_slider_changed(val):
    print(f"Slider value: {val}")

# Add a slider
ax_slider = plt.axes([0.1, 0.1, 0.8, 0.05])  # Position [left, bottom, width, height]
slider = Slider(ax_slider, 'Value', 0, 100, valinit=50)
slider.on_changed(on_slider_changed)

# Initial plot update
update_plot()

plt.show()