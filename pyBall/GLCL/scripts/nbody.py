import numpy as np

config = {
    "simulation_name":  "NBody Simulation",
    "description":      "NBody simulation using OpenCL for computation and OpenGL for rendering.",
    "parameters": {
        #                  value, type, step
        "particle_count": (2048,  "int" , 1      ),
        "dt":             (0.001, "float", 0.0001 )
    },
    "opencl_source": ["nbody_sim.cl"],
    "kernels": {
        #               local_size, global_size     buffers                    parameters
        "nbody_sim" : ( (32,), ("particle_count"), ["positions", "velocities"], ["dt"] ) 
    },
    "opengl_shaders": {
        #               vertex,            fragment           uniforms
        "nbody_render" : ("points.glslv", "monocolor.glslf", ["positions"])
    },
    #                  shader,        element_count,    vertex_buffer,   index_buffer
    "draw_calls":   ( "nbody_render", "particle_count", "positions",     None )
}

def init():
    particle_count = config["particle_count"]
    positions      = np.random.rand(particle_count, 4) * 2 - 1
    velocities     = np.random.rand(particle_count, 4) * 0.1
    return {
        "positions":  positions.astype(np.float32),
        "velocities": velocities.astype(np.float32)
    }