#!/usr/bin/python
import sys
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
import numpy as np

sys.path.append("../../")
from pyBall import FF2D as ff

# Initialize the figure and axis
fig, ax = plt.subplots()
ax.set_title('Interactive Plot: Add, Delete, Connect, and Move Points')
ax.set_xlim(0, 10)
ax.set_ylim(0, 10)

# List to store points and lines
points = []
lines  = []

# Variable to store the index of the point being dragged
dragging_point_index = None

icur_atom = -1

# Function to update the plot
def update_plot():
    global icur_atom
    ax.clear()
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 10)
    ax.set_title('Interactive Plot: Add, Delete, Connect, and Move Points')

    print( "----- py.update_plot() " )
    types, apos, neighs = ff.getAtoms()  #;print("apos ", apos) #;print("types", types) ;print("neighs", ne
    bonds               = ff.getBonds()  #;print("bonds ", bonds)

    # Plot points
    #x_values, y_values = zip(*points)
    #ax.plot(x_values, y_values, 'bo', markersize=5)  # Plot points
    ax.plot( apos[:,0], apos[:,1], 'bo', markersize=5 )
    if icur_atom >= 0:
        ax.plot( [apos[icur_atom,0]], [apos[icur_atom,1]], 'ro', markersize=5 )

    # Plot lines using LineCollection
    lines = [ (apos[i,:], apos[j,:]) for (i,j) in bonds ] 
    line_segments = LineCollection(lines, colors='b')
    ax.add_collection(line_segments)
    
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
                print( "CTRL+LMB: remove atom() ", ia )
                ff.removeAtom(ia)
                bUpdate = True
            else:  # Left Click to add point
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
        if bUpdate:
            update_plot()

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

# Initial plot update
update_plot()

plt.show()