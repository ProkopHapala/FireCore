import bpy    # bpy is the Blender Python API
import math   # math is the Python math module
from typing import Optional, Tuple, Iterable, Sequence
#from utils.node import arrange_nodes
#from utils.modifier import add_subdivision_surface_modifier


elements=None

# =============================================================================
#               Camera
# =============================================================================

def create_camera(location: Tuple[float, float, float]) -> bpy.types.Object:
    bpy.ops.object.camera_add(location=location)

    return bpy.context.object


def set_camera_params(camera: bpy.types.Camera,
                      focus_target_object: bpy.types.Object,
                      lens: float = 85.0,
                      fstop: float = 1.4) -> None:
    # Simulate Sony's FE 85mm F1.4 GM
    camera.sensor_fit = 'HORIZONTAL'
    camera.sensor_width = 36.0
    camera.sensor_height = 24.0
    camera.lens = lens
    camera.dof.use_dof = True
    camera.dof.focus_object = focus_target_object
    camera.dof.aperture_fstop = fstop
    camera.dof.aperture_blades = 11

# =============================================================================
#               Light
# =============================================================================

def create_area_light(location: Tuple[float, float, float] = (0.0, 0.0, 5.0),
                      rotation: Tuple[float, float, float] = (0.0, 0.0, 0.0),
                      size: float = 5.0,
                      color: Tuple[float, float, float, float] = (1.00, 0.90, 0.80, 1.00),
                      energy: float = 1000.0,
                      hide_render: bool = True,
                      name: Optional[str] = None) -> bpy.types.Object:
    if bpy.app.version >= (2, 80, 0):
        bpy.ops.object.light_add(type='AREA', location=location, rotation=rotation)
    else:
        bpy.ops.object.lamp_add(type='AREA', location=location, rotation=rotation)

    obj = bpy.context.object
    if name is not None:
        obj.name = name
    #obj.hide_render = hide_render
    #obj.cycles_visibility.camera = False
    #obj.use_ray_visibility.camera = False
    obj.visible_camera = False
    light = bpy.context.object.data
    #print( light.cycles_visibility )
    light.size      = size
    light.use_nodes = True
    light.node_tree.nodes["Emission"].inputs["Color"].default_value = color
    light.energy    = energy
    return bpy.context.object


def create_sun_light(location: Tuple[float, float, float] = (0.0, 0.0, 5.0),
                     rotation: Tuple[float, float, float] = (0.0, 0.0, 0.0),
                     energy=1.0,
                     name: Optional[str] = None) -> bpy.types.Object:
    bpy.ops.object.light_add(type='SUN', location=location, rotation=rotation )
    obj = bpy.context.object
    obj.data.energy = energy
    if name is not None: obj.name = name
    return obj


def build_environment_texture_background(world: bpy.types.World, hdri_path: str, rotation: float = 0.0) -> None:
    world.use_nodes = True
    node_tree = world.node_tree

    environment_texture_node = node_tree.nodes.new(type="ShaderNodeTexEnvironment")
    environment_texture_node.image = bpy.data.images.load(hdri_path)

    mapping_node = node_tree.nodes.new(type="ShaderNodeMapping")
    if bpy.app.version >= (2, 81, 0):
        mapping_node.inputs["Rotation"].default_value = (0.0, 0.0, rotation)
    else:
        mapping_node.rotation[2] = rotation

    tex_coord_node = node_tree.nodes.new(type="ShaderNodeTexCoord")

    node_tree.links.new(tex_coord_node.outputs["Generated"], mapping_node.inputs["Vector"])
    node_tree.links.new(mapping_node.outputs["Vector"], environment_texture_node.inputs["Vector"])
    node_tree.links.new(environment_texture_node.outputs["Color"], node_tree.nodes["Background"].inputs["Color"])

    #arrange_nodes(node_tree)

# =============================================================================
# Constraints
# =============================================================================

def add_track_to_constraint(camera_object: bpy.types.Object, track_to_target_object: bpy.types.Object) -> None:
    constraint = camera_object.constraints.new(type='TRACK_TO')
    constraint.target = track_to_target_object
    constraint.track_axis = 'TRACK_NEGATIVE_Z'
    constraint.up_axis = 'UP_Y'


def add_copy_location_constraint(copy_to_object: bpy.types.Object,
                                 copy_from_object: bpy.types.Object,
                                 use_x: bool,
                                 use_y: bool,
                                 use_z: bool,
                                 bone_name: str = '') -> None:
    '''
    https://docs.blender.org/api/current/bpy.types.CopyLocationConstraint.html
    '''
    constraint = copy_to_object.constraints.new(type='COPY_LOCATION')
    constraint.target = copy_from_object
    constraint.use_x = use_x
    constraint.use_y = use_y
    constraint.use_z = use_z
    if bone_name:
        constraint.subtarget = bone_name

# =============================================================================
#               Modifiers
# =============================================================================

def add_boolean_modifier(mesh_object: bpy.types.Object,
                         another_mesh_object: bpy.types.Object,
                         operation: str = "DIFFERENCE") -> None:
    '''
    https://docs.blender.org/api/current/bpy.types.BooleanModifier.html
    '''

    modifier: bpy.types.SubsurfModifier = mesh_object.modifiers.new(name="Boolean", type='BOOLEAN')

    modifier.object = another_mesh_object
    modifier.operation = operation


def add_subdivision_surface_modifier(mesh_object: bpy.types.Object, level: int, is_simple: bool = False) -> None:
    '''
    https://docs.blender.org/api/current/bpy.types.SubsurfModifier.html
    '''

    modifier: bpy.types.SubsurfModifier = mesh_object.modifiers.new(name="Subsurf", type='SUBSURF')

    modifier.levels = level
    modifier.render_levels = level
    modifier.subdivision_type = 'SIMPLE' if is_simple else 'CATMULL_CLARK'


def add_solidify_modifier(mesh_object: bpy.types.Object,
                          thickness: float = 0.01,
                          flip_normal: bool = False,
                          fill_rim: bool = True,
                          material_index_offset: int = 0,
                          shell_vertex_group: str = "",
                          rim_vertex_group: str = "") -> None:
    '''
    https://docs.blender.org/api/current/bpy.types.SolidifyModifier.html
    '''

    modifier: bpy.types.SolidifyModifier = mesh_object.modifiers.new(name="Solidify", type='SOLIDIFY')

    modifier.material_offset = material_index_offset
    modifier.thickness = thickness
    modifier.use_flip_normals = flip_normal
    modifier.use_rim = fill_rim

    # TODO: Check whether shell_vertex_group is either empty or defined
    # TODO: Check whether rim_vertex_group is either empty or defined

    modifier.shell_vertex_group = shell_vertex_group
    modifier.rim_vertex_group = rim_vertex_group


def add_displace_modifier(mesh_object: bpy.types.Object,
                          texture_name: str,
                          vertex_group: str = "",
                          mid_level: float = 0.5,
                          strength: float = 1.0) -> None:
    '''
    https://docs.blender.org/api/current/bpy.types.DisplaceModifier.html
    '''

    modifier = mesh_object.modifiers.new(name="Displace", type='DISPLACE')

    modifier.mid_level = mid_level
    modifier.strength = strength

    # TODO: Check whether texture_name is properly defined
    modifier.texture = bpy.data.textures[texture_name]

    # TODO: Check whether vertex_group is either empty or defined
    modifier.vertex_group = vertex_group

# =============================================================================
#               NURBS
# =============================================================================

# =============================================================================
#               Mesh
# =============================================================================

def set_smooth_shading(mesh: bpy.types.Mesh) -> None:
    for polygon in mesh.polygons:
        polygon.use_smooth = True

def create_mesh_from_pydata(scene: bpy.types.Scene,
                            vertices: Iterable[Iterable[float]],
                            faces: Iterable[Iterable[int]],
                            mesh_name: str,
                            object_name: str,
                            use_smooth: bool = True) -> bpy.types.Object:
    '''
    https://docs.blender.org/api/current/bpy.types.Mesh.html
    Add a new mesh and set vertices and faces
    In this case, it does not require to set edges
    After manipulating mesh data, update() needs to be called
    '''
    new_mesh: bpy.types.Mesh = bpy.data.meshes.new(mesh_name)
    new_mesh.from_pydata(vertices, [], faces)
    new_mesh.update()
    if use_smooth:
        set_smooth_shading(new_mesh)
    new_object: bpy.types.Object = bpy.data.objects.new(object_name, new_mesh)
    scene.collection.objects.link(new_object)
    return new_object

def create_cached_mesh_from_alembic(file_path: str, name: str) -> bpy.types.Object:
    bpy.ops.wm.alembic_import(filepath=file_path, as_background_job=False)
    bpy.context.active_object.name = name
    return bpy.context.active_object

def create_plane(location: Tuple[float, float, float] = (0.0, 0.0, 0.0),
                 rotation: Tuple[float, float, float] = (0.0, 0.0, 0.0),
                 size: float = 2.0,
                 name: Optional[str] = None) -> bpy.types.Object:
    '''
    https://docs.blender.org/api/current/bpy.ops.mesh.html#bpy.ops.mesh.primitive_plane_add
    '''
    bpy.ops.mesh.primitive_plane_add(size=size, location=location, rotation=rotation)
    obj = bpy.context.object
    if name is not None:
        obj.name = name
    return obj

def create_smooth_sphere(location: Tuple[float, float, float] = (0.0, 0.0, 0.0),
                         radius: float = 1.0,
                         subdivision_level: int = 1,
                         name: Optional[str] = None) -> bpy.types.Object:
    '''
    https://docs.blender.org/api/current/bpy.ops.mesh.html#bpy.ops.mesh.primitive_uv_sphere_add
    '''
    bpy.ops.mesh.primitive_uv_sphere_add(radius=radius, location=location, calc_uvs=True)
    obj = bpy.context.object
    if name is not None:
        obj.name = name
    set_smooth_shading(obj.data)
    add_subdivision_surface_modifier(obj, subdivision_level)
    return obj

def create_nurbs_sphere(location: Tuple[float, float, float] = (0.0, 0.0, 0.0),
                         radius: float = 1.0,
                         name: Optional[str] = None) -> bpy.types.Object:
    '''
    https://docs.blender.org/api/current/bpy.ops.mesh.html#bpy.ops.mesh.primitive_uv_sphere_add
    '''
    bpy.ops.surface.primitive_nurbs_surface_sphere_add(radius=radius, location=location)
    obj = bpy.context.object
    if name is not None: obj.name = name
    #set_smooth_shading(obj.data)
    #add_subdivision_surface_modifier(obj, subdivision_level)
    return obj

def create_cylinder_between( p1, p2, r, name=None):
    x2, y2, z2 = p2
    x1, y1, z1 = p1
    dx = x2 - x1
    dy = y2 - y1
    dz = z2 - z1    
    dist  = math.sqrt(dx**2 + dy**2 + dz**2)
    '''
    https://docs.blender.org/api/current/bpy.ops.mesh.html#bpy.ops.mesh.primitive_cylinder_add
    https://blender.stackexchange.com/questions/12445/apply-rotation-to-an-object-and-its-children
    '''
    bpy.ops.mesh.primitive_cylinder_add(      radius = r,       depth = dist,      location = (dx/2 + x1, dy/2 + y1, dz/2 + z1)  ) 
    obj = bpy.context.object
    phi   = math.atan2(dy, dx) 
    theta = math.acos(dz/dist) 
    obj.rotation_euler[1] = theta 
    obj.rotation_euler[2] = phi 
    if name is not None: obj.name = name
    return obj

# =============================================================================
#               Set Render
# =============================================================================

def set_output_properties(scene: bpy.types.Scene,
                          resolution_percentage: int = 100,
                          output_file_path: str = "",
                          res_x: int = 1920,
                          res_y: int = 1080) -> None:
    scene.render.resolution_percentage = resolution_percentage
    scene.render.resolution_x = res_x
    scene.render.resolution_y = res_y
    if output_file_path:
        scene.render.filepath = output_file_path

'''
def set_cycles_renderer(scene: bpy.types.Scene,
                        camera_object: bpy.types.Object,
                        num_samples: int,
                        use_denoising: bool = True,
                        use_motion_blur: bool = False,
                        use_transparent_bg: bool = False,
                        prefer_cuda_use: bool = True,
                        use_adaptive_sampling: bool = False) -> None:
    scene.camera = camera_object
    scene.render.image_settings.file_format = 'PNG'
    scene.render.engine = 'CYCLES'
    scene.render.use_motion_blur = use_motion_blur
    scene.render.film_transparent = use_transparent_bg
    scene.view_layers[0].cycles.use_denoising = use_denoising
    scene.cycles.use_adaptive_sampling = use_adaptive_sampling
    scene.cycles.samples = num_samples
    # Enable GPU acceleration
    # Source - https://blender.stackexchange.com/a/196702
    if prefer_cuda_use:
        bpy.context.scene.cycles.device = "GPU"
        # Change the preference setting
        bpy.context.preferences.addons["cycles"].preferences.compute_device_type = "CUDA"
    # Call get_devices() to let Blender detects GPU device (if any)
    bpy.context.preferences.addons["cycles"].preferences.get_devices()
    # Let Blender use all available devices, include GPU and CPU
    for d in bpy.context.preferences.addons["cycles"].preferences.devices:
        d["use"] = 1
    # Display the devices to be used for rendering
    print("----")
    print("The following devices will be used for path tracing:")
    for d in bpy.context.preferences.addons["cycles"].preferences.devices:
        print("- {}".format(d["name"]))
    print("----")
'''

def set_renderer(scene: bpy.types.Scene,
                        camera_object: bpy.types.Object,
                        num_samples: int,
                        use_denoising: bool = True,
                        use_motion_blur: bool = False,
                        use_transparent_bg: bool = False,
                        prefer_cuda_use: bool = True,
                        use_adaptive_sampling: bool = True,
                        engine : str ='CYCLES'
                        ) -> None:
    scene.camera = camera_object
    scene.render.image_settings.file_format = 'PNG'
    scene.render.engine = engine

    scene.render.use_motion_blur  = use_motion_blur
    scene.render.film_transparent = use_transparent_bg
    scene.view_layers[0].cycles.use_denoising = use_denoising

    if engine == 'CYCLES':
        scene.cycles.use_adaptive_sampling = use_adaptive_sampling
        scene.cycles.samples = num_samples
        # Enable GPU acceleration
        # Source - https://blender.stackexchange.com/a/196702
        if prefer_cuda_use:
            bpy.context.scene.cycles.device = "GPU"
            # Change the preference setting
            bpy.context.preferences.addons["cycles"].preferences.compute_device_type = "CUDA"
        # Call get_devices() to let Blender detects GPU device (if any)
        bpy.context.preferences.addons["cycles"].preferences.get_devices()
        # Let Blender use all available devices, include GPU and CPU
        for d in bpy.context.preferences.addons["cycles"].preferences.devices:
            d["use"] = 1
        # Display the devices to be used for rendering
        print("The following devices will be used for path tracing:")
        for d in bpy.context.preferences.addons["cycles"].preferences.devices:
            print("- {}".format(d["name"]))
    elif engine == 'BLENDER_EEVEE':
        eevee = bpy.context.scene.eevee
        eevee.use_soft_shadows = True
        eevee.use_ssr = True
        eevee.use_ssr_refraction = True
        eevee.use_gtao = True
        eevee.gtao_distance = 1
        eevee.use_volumetric_shadows = True
        eevee.volumetric_tile_size = '2'



def clean_objects() -> None:
    for item in bpy.data.objects:
        bpy.data.objects.remove(item)

# =============================================================================
#               Material
# =============================================================================

def clean_nodes(nodes: bpy.types.Nodes) -> None:
    for node in nodes:
        nodes.remove(node)

def set_principled_node(principled_node: bpy.types.Node,
                        base_color: Tuple[float, float, float, float] = (0.6, 0.6, 0.6, 1.0),
                        subsurface: float = 0.0,
                        subsurface_color: Tuple[float, float, float, float] = (0.8, 0.8, 0.8, 1.0),
                        subsurface_radius: Tuple[float, float, float] = (1.0, 0.2, 0.1),
                        metallic: float = 0.0,
                        specular: float = 0.5,
                        specular_tint: float = 0.0,
                        roughness: float = 0.5,
                        anisotropic: float = 0.0,
                        anisotropic_rotation: float = 0.0,
                        sheen: float = 0.0,
                        sheen_tint: float = 0.5,
                        clearcoat: float = 0.0,
                        clearcoat_roughness: float = 0.03,
                        ior: float = 1.45,
                        transmission: float = 0.0,
                        transmission_roughness: float = 0.0) -> None:
    principled_node.inputs['Base Color'].default_value = base_color
    principled_node.inputs['Subsurface'].default_value = subsurface
    principled_node.inputs['Subsurface Color'].default_value = subsurface_color
    principled_node.inputs['Subsurface Radius'].default_value = subsurface_radius
    principled_node.inputs['Metallic'].default_value = metallic
    principled_node.inputs['Specular'].default_value = specular
    principled_node.inputs['Specular Tint'].default_value = specular_tint
    principled_node.inputs['Roughness'].default_value = roughness
    principled_node.inputs['Anisotropic'].default_value = anisotropic
    principled_node.inputs['Anisotropic Rotation'].default_value = anisotropic_rotation
    principled_node.inputs['Sheen'].default_value = sheen
    principled_node.inputs['Sheen Tint'].default_value = sheen_tint
    principled_node.inputs['Clearcoat'].default_value = clearcoat
    principled_node.inputs['Clearcoat Roughness'].default_value = clearcoat_roughness
    principled_node.inputs['IOR'].default_value = ior
    principled_node.inputs['Transmission'].default_value = transmission
    principled_node.inputs['Transmission Roughness'].default_value = transmission_roughness

def build_pbr_nodes(node_tree: bpy.types.NodeTree,
                    base_color: Tuple[float, float, float, float] = (0.6, 0.6, 0.6, 1.0),
                    metallic: float = 0.0,
                    specular: float = 0.5,
                    roughness: float = 0.5,
                    sheen: float = 0.0) -> None:
    output_node = node_tree.nodes.new(type='ShaderNodeOutputMaterial')
    principled_node = node_tree.nodes.new(type='ShaderNodeBsdfPrincipled')
    node_tree.links.new(principled_node.outputs['BSDF'], output_node.inputs['Surface'])

    set_principled_node(principled_node=principled_node,
                        base_color=base_color,
                        metallic=metallic,
                        specular=specular,
                        roughness=roughness,
                        sheen=sheen)
    #arrange_nodes(node_tree)

def build_checker_board_nodes(node_tree: bpy.types.NodeTree, size: float) -> None:
    output_node = node_tree.nodes.new(type='ShaderNodeOutputMaterial')
    principled_node = node_tree.nodes.new(type='ShaderNodeBsdfPrincipled')
    checker_texture_node = node_tree.nodes.new(type='ShaderNodeTexChecker')

    set_principled_node(principled_node=principled_node)
    checker_texture_node.inputs['Scale'].default_value = size

    node_tree.links.new(checker_texture_node.outputs['Color'], principled_node.inputs['Base Color'])
    node_tree.links.new(principled_node.outputs['BSDF'], output_node.inputs['Surface'])
    #arrange_nodes(node_tree)

def build_emission_nodes(node_tree: bpy.types.NodeTree,
                         color: Tuple[float, float, float] = (0.0, 0.0, 0.0),
                         strength: float = 1.0) -> None:
    '''
    https://docs.blender.org/api/current/bpy.types.ShaderNodeEmission.html
    '''
    output_node = node_tree.nodes.new(type='ShaderNodeOutputMaterial')
    emission_node = node_tree.nodes.new(type='ShaderNodeEmission')

    emission_node.inputs["Color"].default_value = color + (1.0, )
    emission_node.inputs["Strength"].default_value = strength

    node_tree.links.new(emission_node.outputs['Emission'], output_node.inputs['Surface'])

    #arrange_nodes(node_tree)

def add_material(name: str = "Material",
                 use_nodes: bool = False,
                 make_node_tree_empty: bool = False) -> bpy.types.Material:
    '''
    https://docs.blender.org/api/current/bpy.types.BlendDataMaterials.html
    https://docs.blender.org/api/current/bpy.types.Material.html
    '''

    # TODO: Check whether the name is already used or not

    material = bpy.data.materials.new(name)
    material.use_nodes = use_nodes

    if use_nodes and make_node_tree_empty:
        clean_nodes(material.node_tree.nodes)

    return material

def make_pbr_mat(name):
    mat = add_material(name, use_nodes=True, make_node_tree_empty=True)
    nodes = mat.node_tree.nodes
    links = mat.node_tree.links
    output_node     = nodes.new(type='ShaderNodeOutputMaterial')
    principled_node = nodes.new(type='ShaderNodeBsdfPrincipled')
    links.new(principled_node.outputs['BSDF'], output_node.inputs['Surface'])
    return mat,principled_node

def set_principled_node_as_rough_blue(principled_node: bpy.types.Node) -> None:
    set_principled_node(
        principled_node=principled_node,
        base_color=(0.1, 0.2, 0.6, 1.0),
        metallic=0.5,
        specular=0.5,
        roughness=0.9,
    )

def set_principled_node_as_ceramic(principled_node: bpy.types.Node) -> None:
    set_principled_node(
        principled_node=principled_node,
        base_color=(0.8, 0.8, 0.8, 1.0),
        subsurface=0.1,
        subsurface_color=(0.9, 0.9, 0.9, 1.0),
        subsurface_radius=(1.0, 1.0, 1.0),
        metallic=0.2,
        specular=0.5,
        roughness=0.0,
    )

def set_principled_node_as_gold(principled_node: bpy.types.Node) -> None:
    set_principled_node(
        principled_node=principled_node,
        base_color=(1.00, 0.71, 0.22, 1.0),
        metallic=1.0,
        specular=0.5,
        roughness=0.1,
    )

def set_principled_node_as_glass(principled_node: bpy.types.Node) -> None:
    set_principled_node(principled_node=principled_node,
                              base_color=(0.95, 0.95, 0.95, 1.0),
                              metallic=0.0,
                              specular=0.5,
                              roughness=0.0,
                              clearcoat=0.5,
                              clearcoat_roughness=0.030,
                              ior=1.45,
                              transmission=0.98)

# =============================================================================
#               Macros
# =============================================================================

def colorFromStr( s ):
    r = int( s[1:3], 16 )
    g = int( s[3:5], 16 )
    b = int( s[5:7], 16 )
    return (r,g,b)

def make_atom_materials( es, metallic=0.5, roughness=0.001, specular_intensity=1.0 ):
    elems = set(es)
    mats={}
    for e in elems:
        mat   = bpy.data.materials.new("mat_"+e )
        epars = elements.ELEMENT_DICT[e]
        c = colorFromStr( epars[8] )
        print( e, epars[8], c )
        mat.diffuse_color = ( c[0]/256.0, c[1]/256.0, c[2]/256.0,  1.0 )
        mat.roughness = roughness
        mat.specular_intensity = specular_intensity
        mat.metallic = metallic
        mats[e]=mat
    return mats

def make_atom_materials_pbr( es, metallic=0.0, specular=1.0, roughness=0.05, clearcoat=0.0, clearcoat_roughness=0.030, ior=1.45, transmission=0.0 ):
    elems = set(es)
    mats={}
    for e in elems:
        epars = elements.ELEMENT_DICT[e]
        c = colorFromStr( epars[8] )
        mat = add_material("Material_Left", use_nodes=True, make_node_tree_empty=True)
        nodes = mat.node_tree.nodes
        links = mat.node_tree.links
        output_node = nodes.new(type='ShaderNodeOutputMaterial')
        principled_node = nodes.new(type='ShaderNodeBsdfPrincipled')
        set_principled_node(principled_node=principled_node,
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

def make_atoms( ps, es, mats=None, rs=None, subdivision_level=1, names=None, R0=-0.8, Rsc=0.5, bPBR=True ) -> bpy.types.Object:
    n = len(ps)
    #if rs is None:     rs    = [ R for i in range(n) ]
    if names is None:  names = [ ("atom_%s_%03i" %(es[i],i)) for i in range(n) ]
    
    if mats is None:
        if bPBR:
            mats = make_atom_materials_pbr( es )
        else:
            mats = make_atom_materials( es )
    for i in range( n ):
        epars = elements.ELEMENT_DICT[ es[i] ]
        Ri   = (epars[7] + R0)*Rsc  
        p = ps[i]
        #utils.create_smooth_monkey(location=((index - (num_suzannes - 1) / 2) * 3.0, 0.0, 0.0),  name="Suzanne" + str(index))
        #o = blender.create_smooth_sphere(location=( p[0],p[1],p[2] ),   radius = Ri,  subdivision_level= subdivision_level,  name = names[i] )
        o = create_nurbs_sphere( location=( p[0],p[1],p[2] ),   radius = Ri, name = names[i] )
        #bpy.ops.surface.primitive_nurbs_surface_sphere_add(radius=Ri, location=p ); o = bpy.context.object
        o.active_material = mats[ es[i] ]
    return bpy.data.objects[ names[0] ]

def make_bonds( ps, bonds, r=1., mat=None, c=(128,128,128), names=None ):
    if mat is None:
        mat = bpy.data.materials.new("mat_bond" )
        mat.diffuse_color = ( c[0]/256.0, c[1]/256.0, c[2]/256.0,  1.0 )
    for i,b in enumerate(bonds):
        p1 = ps[b[0]]
        p2 = ps[b[1]]
        o = create_cylinder_between( p1, p2, r, name = "bond_%i" %i )
        o.active_material = mat 
    return bpy.data.objects[ "bond_0" ]

def set_background_color( color ):
    world = bpy.context.scene.world
    world.use_nodes = True
    shader_node = world.node_tree.nodes.get('Background')
    if shader_node:
        shader_node.inputs['Color'].default_value = color
    else:
        print("Background shader node not found.")

def finishScene( cam_target, cam_pos=(0.0, 0.0, 100.0), sun_dir=(0.0,0.0,-1.0), sun_energy=10.0, light_size=None, light_pos=(50.0,50.0,100.0), lens=50.0, fname_render="./out/render.png", fname_blend="tmp.blend",  resolution_percentage=100, num_samples=128, bLight=True, background_color=(0,0,0,1), bCycles=True ):
    cam = create_camera(location=cam_pos )
    add_track_to_constraint(cam,      cam_target )
    set_camera_params      (cam.data, cam_target, lens=lens )
    set_background_color( background_color )
    if bLight:
        if light_size is not None:
            create_area_light(location=light_pos, rotation=sun_dir, size=light_size, color=(1.0,1.0,1.0,1.0),energy=sun_energy*100.0 )
        elif sun_dir is not None: 
            create_sun_light( rotation=sun_dir, energy=sun_energy )

    scene = bpy.data.scenes["Scene"]
    bpy.ops.wm.save_as_mainfile(filepath="tmp.blend")
    set_output_properties(scene, resolution_percentage, fname_render)
    if bCycles:
        set_renderer(scene, cam, num_samples,engine='CYCLES' )
    else:
        set_renderer(scene, cam, num_samples,engine='BLENDER_EEVEE' )
