# blender --background --python 01_cube.py --render-frame 1 -- </path/to/output/image> <resolution_percentage> <num_samples>
# using https://github.com/yuki-koyama/blender-cli-rendering

import bpy
import sys
import math
import os
import codecs

from pyBall import atomicUtils as au
from pyBall import elements

working_dir_path = os.path.dirname(os.path.abspath(__file__))
sys.path.append(working_dir_path)

import utils

atoms = au.AtomicSystem('AT.xyz')

def colorFromStr( s ):
    r = int( s[1:3], 16 )
    g = int( s[3:5], 16 )
    b = int( s[5:7], 16 )
    return (r,g,b)

def make_materials( es ):
    elems = set(es)
    mats={}
    for e in elems:
        mat = bpy.data.materials.new("mat_"+e )
        epars = elements.ELEMENT_DICT[e]
        c = colorFromStr( epars[8] )
        print( e, epars[8], c )
        mat.diffuse_color = ( c[0]/256.0, c[1]/256.0, c[2]/256.0,  1.0 )
        mats[e]=mat
    return mats

def make_atoms( ps, es, mats=None, rs=None, subdivision_level = 1, names=None, R0=-1.0, Rsc=1.0 ) -> bpy.types.Object:
    n = len(ps)
    #if rs is None:     rs    = [ R for i in range(n) ]
    if names is None:  names = [ ("sphere%03i" %i) for i in range(n) ]
    if mats is None:   mats  = make_materials( atoms.enames )
    for i in range( n ):
        epars = elements.ELEMENT_DICT[ es[i] ]
        Ri   = (epars[7] + R0)*Rsc  
        #utils.create_smooth_monkey(location=((index - (num_suzannes - 1) / 2) * 3.0, 0.0, 0.0),  name="Suzanne" + str(index))
        o = utils.create_smooth_sphere(location=( ps[i][0],ps[i][1],ps[i][2] ),   radius = Ri,  subdivision_level= subdivision_level,  name = names[i] )
        o.active_material = mats[ es[i] ]

    return bpy.data.objects[ names[0] ]


# Args
output_file_path = bpy.path.relpath(str(sys.argv[sys.argv.index('--') + 1]))
resolution_percentage = int(sys.argv[sys.argv.index('--') + 2])
num_samples = int(sys.argv[sys.argv.index('--') + 3])

# Scene Building

## Reset
utils.clean_objects()

objects = make_atoms( atoms.apos, atoms.enames )

## Camera
camera_object = utils.create_camera(location=(0.0, 0.0, 100.0))

utils.add_track_to_constraint(camera_object, objects)
utils.set_camera_params(camera_object.data, objects, lens=50.0)

## Lights
#utils.create_sun_light(  rotation=(0.0, math.pi * 0.5, -math.pi * 0.1) )

utils.create_sun_light(  rotation=( 0.0, 0.0, -1.0 ) )

# Render Setting
scene = bpy.data.scenes["Scene"]
utils.set_output_properties(scene, resolution_percentage, output_file_path)
utils.set_cycles_renderer(scene, camera_object, num_samples)
