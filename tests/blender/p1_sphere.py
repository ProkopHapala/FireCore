# blender --background --python 01_cube.py --render-frame 1 -- </path/to/output/image> <resolution_percentage> <num_samples>
# using https://github.com/yuki-koyama/blender-cli-rendering


import bpy
import sys
import math
import os

working_dir_path = os.path.dirname(os.path.abspath(__file__))
sys.path.append(working_dir_path)

import utils

'''
def set_scene_objects() -> bpy.types.Object:
    num_suzannes = 15
    for index in range(num_suzannes):
        utils.create_smooth_monkey(location=((index - (num_suzannes - 1) / 2) * 3.0, 0.0, 0.0),
                                   name="Suzanne" + str(index))
    return bpy.data.objects["Suzanne" + str(int((num_suzannes - 1) / 2))]
'''

def make_spheres( ps, rs=None, R=1.0, subdivision_level = 1, names=None ) -> bpy.types.Object:
    n = len(ps)
    if rs is None:
        rs = [ R for i in range(n) ]
    if names is None:
        names = [ ("sphere%03i" %i) for i in range(n) ]
    for i in range( n ):
        #utils.create_smooth_monkey(location=((index - (num_suzannes - 1) / 2) * 3.0, 0.0, 0.0),  name="Suzanne" + str(index))
        utils.create_smooth_sphere(location=( ps[i][0],ps[i][1],ps[i][2] ),   radius = rs[i],  subdivision_level= subdivision_level,  name = names[i] )
    return bpy.data.objects[ names[0] ]


# Args
output_file_path = bpy.path.relpath(str(sys.argv[sys.argv.index('--') + 1]))
resolution_percentage = int(sys.argv[sys.argv.index('--') + 2])
num_samples = int(sys.argv[sys.argv.index('--') + 3])

# Scene Building

## Reset
utils.clean_objects()

ps = [ (-5.0,0.0,0.0), (0.0,0.0,0.0), (+5.0,0.0,0.0) ]
objects = make_spheres( ps, rs=[1. , 1.5, 2.0] )

## Camera
camera_object = utils.create_camera(location=(10.0, -7.0, 0.0))

utils.add_track_to_constraint(camera_object, objects)
utils.set_camera_params(camera_object.data, objects, lens=50.0)

## Lights
utils.create_sun_light(rotation=(0.0, math.pi * 0.5, -math.pi * 0.1))

# Render Setting
scene = bpy.data.scenes["Scene"]
utils.set_output_properties(scene, resolution_percentage, output_file_path)
utils.set_cycles_renderer(scene, camera_object, num_samples)
