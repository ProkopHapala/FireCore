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

#build_environment_texture_background(world: bpy.types.World, hdri_path: str, rotation: float = 0.0)

#blender.build_environment_texture_background( bpy.data.scenes["Scene"].world, hdri_path = "/home/prokop/git_SW/blender-cli-rendering/assets/HDRIs/royal_esplanade_2k.hdr",  rotation = 0.0  )
#blender.build_environment_texture_background( bpy.data.scenes["Scene"].world, hdri_path = "/home/prokop/git_SW/blender-cli-rendering/assets/HDRIs/green_point_park_2k.hdr",  rotation = 0.0  )
#blender.build_environment_texture_background( bpy.data.scenes["Scene"].world, hdri_path = "/home/prokop/git_SW/blender-cli-rendering/assets/HDRIs/syferfontein_6d_clear_2k.hdr",  rotation = 0.0  )
#blender.build_environment_texture_background( bpy.data.scenes["Scene"].world, hdri_path = "/home/prokop/git_SW/blender-cli-rendering/assets/HDRIs/old_depot_2k.hdr",  rotation = 0.0        )
#blender.finishScene( cam_focus, cam_pos=(10.0, 10.0, 100.0), bLight=False, lens=50.0, fname_render="./out/render.png", fname_blend="tmp.blend",  resolution_percentage=100, num_samples=64 )

blender.finishScene( cam_focus, cam_pos=(10.0, 10.0, 100.0), sun_dir=(-0.5,-0.5,-1.0), sun_energy=1000.0, light_size=20.0, light_pos=(50.0,50.0,100.0), lens=50.0, fname_render="./out/render.png", fname_blend="tmp.blend", resolution_percentage=100, num_samples=64, background_color=(0.05,0.05,0.05,1.0) )
#blender.finishScene( cam_focus, cam_pos=(10.0, 10.0, 100.0), sun_dir=(-0.5,-0.5,-1.0), sun_energy=1000.0, light_size=20.0, light_pos=(50.0,50.0,100.0), lens=50.0, fname_render="./out/render.png", fname_blend="tmp.blend", resolution_percentage=100, num_samples=64, background_color=(0.05,0.05,0.05,1.0), bCycles=False )
