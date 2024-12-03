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

def make_materials( es, metallic=0.0,  specular=0.5,   roughness=0.0, clearcoat=0.5, clearcoat_roughness=0.030,  ior=1.45,  transmission=0.98   ):
    elems = set(es)
    mats={}
    for e in elems:
        '''
        mat = bpy.data.materials.new("mat_"+e )
        epars = elements.ELEMENT_DICT[e]
        c = colorFromStr( epars[8] )
        print( e, epars[8], c )
        mat.diffuse_color = ( c[0]/256.0, c[1]/256.0, c[2]/256.0,  1.0 )
        mats[e]=mat
        '''
        

        '''
        mat = utils.add_material("Material_Left", use_nodes=True, make_node_tree_empty=True)
        nodes = mat.node_tree.nodes
        links = mat.node_tree.links
        output_node = nodes.new(type='ShaderNodeOutputMaterial')
        principled_node = nodes.new(type='ShaderNodeBsdfPrincipled')
        if( e=='C' ):
            set_principled_node_as_glass(principled_node)
        elif( e=='O' ):
            set_principled_node_as_gold(principled_node)
        elif( e=='N' ):
            set_principled_node_as_rough_blue(principled_node)
        elif( e=='H' ):    
            set_principled_node_as_ceramic(principled_node)
        links.new(principled_node.outputs['BSDF'], output_node.inputs['Surface'])
        #left_object.data.materials.append(mat)
        '''

        epars = elements.ELEMENT_DICT[e]
        c = colorFromStr( epars[8] )
        mat = utils.add_material("Material_Left", use_nodes=True, make_node_tree_empty=True)
        nodes = mat.node_tree.nodes
        links = mat.node_tree.links
        output_node = nodes.new(type='ShaderNodeOutputMaterial')
        principled_node = nodes.new(type='ShaderNodeBsdfPrincipled')
        utils.set_principled_node(principled_node=principled_node,
                              base_color= ( c[0]/256.0, c[1]/256.0, c[2]/256.0,  1.0 ),
                              metallic=metallic,
                              specular=specular,
                              roughness=roughness,
                              clearcoat=clearcoat,
                              clearcoat_roughness=clearcoat_roughness,
                              ior=ior,
                              transmission=transmission )
        links.new(principled_node.outputs['BSDF'], output_node.inputs['Surface'])
        #left_object.data.materials.append(mat)

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
        mat = mats[ es[i] ]
        #o.active_material = mat
        o.data.materials.append(mat)

    return bpy.data.objects[ names[0] ]


# Args
output_file_path = bpy.path.relpath(str(sys.argv[sys.argv.index('--') + 1]))
resolution_percentage = int(sys.argv[sys.argv.index('--') + 2])
num_samples = int(sys.argv[sys.argv.index('--') + 3])

hdri_path = os.path.join(working_dir_path, "assets/HDRIs/green_point_park_2k.hdr")
scene = bpy.data.scenes["Scene"]
world = scene.world

## Reset
utils.clean_objects()

objects = make_atoms( atoms.apos, atoms.enames )


# ----- Floor
floor_mat, floor_node = utils.make_pbr_mat("Floor_mat")
utils.set_principled_node_as_ceramic(floor_node)
current_object = utils.create_plane(size=1000.0, name="Floor", location=(0.0,0.0,-3.0))
current_object.data.materials.append(floor_mat)


## Camera
camera_object = utils.create_camera(location=(0.0, 0.2, 30.0))

utils.add_track_to_constraint(camera_object, objects)
utils.set_camera_params(camera_object.data, objects, lens=50.0)

## Lights
#utils.create_sun_light(  rotation=(0.0, math.pi * 0.5, -math.pi * 0.1) )
#utils.create_sun_light(  rotation=( 0.0, 0.0, -1.0 ) )

# Render Setting
utils.build_environment_texture_background(world, hdri_path)
utils.set_output_properties(scene, resolution_percentage, output_file_path)
utils.set_cycles_renderer(scene, camera_object, num_samples)
