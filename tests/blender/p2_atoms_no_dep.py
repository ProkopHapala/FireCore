# blender --background --python 01_cube.py --render-frame 1 -- </path/to/output/image> <resolution_percentage> <num_samples>
# using https://github.com/yuki-koyama/blender-cli-rendering

import bpy
import sys
import math
import os
import codecs

from pyBall import atomicUtils as au
from pyBall import elements

from pyBall import blender_cli as blender

mol = au.AtomicSystem('AT.xyz')
mol.findBonds()
blender.elements = elements

# ------ Scene Building
blender.clean_objects()


cam_focus = blender.make_atoms( mol.apos, mol.enames )
cam_focus = blender.make_bonds( mol.apos, mol.bonds, r=0.1, mat=None )

# ----- Floor
floor_mat, floor_node = blender.make_pbr_mat("Floor_mat")
blender.set_principled_node_as_ceramic(floor_node)
current_object = blender.create_plane(size=1000.0, name="Floor", location=(0.0,0.0,-3.0))
current_object.data.materials.append(floor_mat)


blender.finishScene( cam_focus, cam_pos=(0.0, 0.0, 100.0), sun_dir=(-0.5,-0.5,-1.0), light_size=2.0, lens=50.0, fname_render="./out/render.png", fname_blend="tmp.blend",  resolution_percentage=100, num_samples=128 )
