#!/usr/bin/pythonimport bpy
import numpy as np
import sys
import os
import argparse

# Assuming 'AtomiSystem.py', 'elements.py', and 'atomicUtils.py' are in the Python path
from AtomicSystem import AtomicSystem
import elements

def setup_render_engine():
    """
    Initializes Blender's Cycles render engine and configures GPU rendering.
    """
    bpy.context.scene.render.engine = 'CYCLES'
    
    # Attempt to set the device to GPU
    bpy.context.scene.cycles.device = 'GPU'

    # Get the system's Cycles preferences and enable the first available GPU device
    prefs = bpy.context.preferences.addons['cycles'].preferences
    prefs.get_devices()
    
    # Set the compute device type (e.g., CUDA, OPTIX, HIP, METAL)
    # You might want to make this an argument if you use different hardware
    prefs.compute_device_type = 'CUDA' 
    
    print("Available devices:")
    for device in prefs.devices:
        print(f"- {device.name} ({device.type})")

    # Enable all compatible GPUs
    for device in prefs.devices:
        if device.type == prefs.compute_device_type:
            device.use = True
            print(f"Using device: {device.name}")

def clear_scene():
    """
    Removes all objects from the current scene to start fresh.
    """
    # Ensure we are in Object Mode
    if bpy.context.object and bpy.context.object.mode != 'OBJECT':
        bpy.ops.object.mode_set(mode='OBJECT')
        
    # Select all objects
    bpy.ops.object.select_all(action='SELECT')
    
    # Delete selected objects
    bpy.ops.object.delete()

def render_molecule(system, output_path):
    """
    Generates and renders the molecular visualization in Blender.

    Args:
        system (AtomiSystem): An instance of the AtomiSystem class containing the molecular data.
        output_path (str): The file path for the rendered output image.
    """
    clear_scene()
    setup_render_engine()

    # --- Create a collection for the molecule for better organization ---
    molecule_collection = bpy.data.collections.new("Molecule")
    bpy.context.scene.collection.children.link(molecule_collection)

    # --- Create a NURBS sphere for each atom ---
    for i, pos in enumerate(system.apos):
        ename = system.enames[i]
        
        # Get atomic properties from the 'elements' module
        radius = elements.ATOMIC_RADII.get(ename, 0.5) # Default radius if not found
        color = elements.CPK_COLORS.get(ename, elements.CPK_COLORS['DEFAULT']) # Default color
        
        # Create the sphere object
        bpy.ops.surface.primitive_nurbs_surface_sphere_add(
            radius=radius,
            location=tuple(pos),
            enter_editmode=False,
            align='WORLD'
        )
        atom_obj = bpy.context.object
        atom_obj.name = f"{ename}_{i}"

        # --- Create and assign a material for the atom ---
        mat = bpy.data.materials.new(name=f"{ename}_mat")
        mat.use_nodes = True
        bsdf = mat.node_tree.nodes.get("Principled BSDF")
        if bsdf:
            bsdf.inputs['Base Color'].default_value = color
            bsdf.inputs['Roughness'].default_value = 0.5
            bsdf.inputs['Sheen'].default_value = 0.5
        atom_obj.data.materials.append(mat)
        
        # --- Link the new atom to the molecule collection ---
        # Unlink from the scene's master collection first
        bpy.context.scene.collection.objects.unlink(atom_obj)
        # Link to our new molecule collection
        molecule_collection.objects.link(atom_obj)

    # --- Setup Camera ---
    bpy.ops.object.camera_add(location=(0, -15, 0)) # Position camera
    camera = bpy.context.object
    bpy.context.scene.camera = camera # Set as active camera

    # Point camera to the center of the molecule
    center_of_mass = np.mean(system.apos, axis=0)
    cam_location = camera.location
    direction = tuple(center_of_mass) - cam_location
    rot_quat = direction.to_track_quat('-Z', 'Y')
    camera.rotation_euler = rot_quat.to_euler()

    # --- Setup Lighting ---
    bpy.ops.object.light_add(type='SUN', align='WORLD', location=(10, -10, 10))
    light = bpy.context.object
    light.data.energy = 4.0 # Adjust brightness

    # --- Configure Render Settings ---
    scene = bpy.context.scene
    scene.render.image_settings.file_format = 'PNG'
    scene.render.filepath = output_path
    scene.render.resolution_x = 1920
    scene.render.resolution_y = 1080
    scene.render.resolution_percentage = 100
    scene.cycles.samples = 128 # Adjust for quality vs. speed

    # --- Render the final image ---
    print(f"Rendering scene to {output_path}...")
    bpy.ops.render.render(write_still=True)
    print("Rendering complete.")

if __name__ == "__main__":
    """
    Main function to parse arguments and run the rendering process.
    """
    # Blender's Python API requires parsing args passed after '--'
    argv = sys.argv
    if "--" not in argv:
        argv = []
    else:
        argv = argv[argv.index("--") + 1:]

    parser = argparse.ArgumentParser(description="Render a molecule from an .xyz file using Blender.")
    parser.add_argument("-i", "--input", required=True, help="Path to the input .xyz file.")
    parser.add_argument("-o", "--output", required=True, help="Path for the output render image.")
    args = parser.parse_args(argv)

    # Load the molecule using your AtomiSystem class
    print(f"Loading molecule from: {args.input}")
    molecule_system = AtomiSystem(fname=args.input)

    # Render the molecule
    render_molecule(molecule_system, args.output)
