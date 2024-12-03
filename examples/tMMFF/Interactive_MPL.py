#!/usr/bin/python

import matplotlib.pyplot as plt

# Initialize the figure and axis
fig, ax = plt.subplots()
ax.set_title('Interactive Plot: Add, Delete, Connect, and Move Points')
ax.set_xlim(0, 10)
ax.set_ylim(0, 10)

# List to store points
points = []

# Variables to store the index of the point being dragged
dragging_point_index = None

# Function to update the plot
def update_plot():
    ax.clear()
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 10)
    ax.set_title('Interactive Plot: Add, Delete, Connect, and Move Points')
    if points:
        x_values, y_values = zip(*points)
        ax.plot(x_values, y_values, 'bo-', markersize=5)  # Plot points and lines
    fig.canvas.draw()

# Event handler function for mouse button press
def on_click(event):
    R = 0.1**2
    global dragging_point_index
    if event.inaxes is not None:
        if event.button == 1:  # Left mouse button
            if event.key == 'control':  # Ctrl + Left Click to delete point
                for (i, (x, y)) in enumerate(points):
                    r2 = (x - event.xdata)**2 + (y - event.ydata)**2
                    if r2 < R2:
                        points.pop(i)
                        break
            else:  # Left Click to add point
                points.append((event.xdata, event.ydata))
        elif event.button == 3:  # Right mouse button to connect points
            if len(points) > 1:
                x_values, y_values = zip(*points)
                ax.plot(x_values, y_values, 'b-')  # Draw connecting line
        elif event.button == 2:  # Middle mouse button to start dragging
            for (i, (x, y)) in enumerate(points):
                if abs(x - event.xdata) < 0.1 and abs(y - event.ydata) < 0.1:
                    dragging_point_index = i
                    break
        update_plot()

# Event handler function for mouse motion
def on_motion(event):
    global dragging_point_index
    if dragging_point_index is not None and event.inaxes is not None:
        points[dragging_point_index] = (event.xdata, event.ydata)
        update_plot()

# Event handler function for mouse button release
def on_release(event):
    global dragging_point_index
    dragging_point_index = None

# Connect the event handlers to the figure
fig.canvas.mpl_connect('button_press_event', on_click)
fig.canvas.mpl_connect('motion_notify_event', on_motion)
fig.canvas.mpl_connect('button_release_event', on_release)

# Initial plot update
update_plot()

plt.show()