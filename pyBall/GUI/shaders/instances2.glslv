#version 330 core
layout (location = 0) in vec3 aPos;
layout (location = 1) in vec3 aNormal;

// Per-instance attributes
layout (location = 2) in vec3 instancePosition_model;
layout (location = 3) in float instanceActualSphereRadius; // Used for spheres, can be repurposed
layout (location = 4) in vec4 instanceColor;
layout (location = 5) in mat4 instanceMatrix; // NEW: Per-instance transformation (rotation, scale)

// Uniforms
uniform mat4 projection;
uniform mat4 view;
uniform mat4 model; // Global model matrix (from trackball)

// Outputs to Fragment Shader
out vec3 fNormal_world;
out vec3 fpos_world;
out vec4 fColor;

// For sphere ray-tracing
out vec4 sphere_obj_world;
out vec3 viewPos_world;

void main() {
    // Transform normal and position by the instance's unique matrix first, then by the global model matrix
    mat4 model_instance = model * instanceMatrix;
    vec4 world_pos_vec4 = model_instance * vec4(aPos, 1.0) + vec4(instancePosition_model, 0.0);
    gl_Position = projection * view * world_pos_vec4;

    // Data for fragment shader
    fpos_world = world_pos_vec4.xyz;
    fNormal_world = mat3(transpose(inverse(model_instance))) * aNormal;
    fColor = instanceColor;

    // Data for sphere fragment shader
    vec4 sphere_center_world = model * vec4(instancePosition_model, 1.0);
    float sphere_radius_world = instanceActualSphereRadius; // Assuming uniform scaling for simplicity
    sphere_obj_world = vec4(sphere_center_world.xyz, sphere_radius_world);
    viewPos_world = (inverse(view) * vec4(0.0, 0.0, 0.0, 1.0)).xyz;
}