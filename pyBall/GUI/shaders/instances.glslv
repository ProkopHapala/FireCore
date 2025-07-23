#version 330 core

layout (location = 0) in vec3 aPos;      // Vertex position of the base sphere mesh
layout (location = 1) in vec3 aNormal;   // Vertex normal of the base sphere mesh (not used by raytracer directly)

// Per-instance attributes
layout (location = 2) in vec4 instanceMatrix_row0;
layout (location = 3) in vec4 instanceMatrix_row1;
layout (location = 4) in vec4 instanceMatrix_row2;
layout (location = 5) in vec4 instanceMatrix_row3;
layout (location = 6) in float instanceActualSphereRadius;  // Actual radius of the sphere to ray-trace
layout (location = 7) in vec4 instanceColor;    // Color (RGBA) of the sphere instance

// Outputs to Fragment Shader
out vec3 fpos_world;        // Fragment position on the bounding mesh in world space
out vec4 sphere_obj_world;  // Sphere center (xyz) and radius (w) in world space
out vec4 fColor;             // Pass through atom color

uniform mat4 projection;
uniform mat4 view;
uniform mat4 model; // Overall model orientation (trackball)

void main()
{
    mat4 instanceMatrix = mat4(
        instanceMatrix_row0,
        instanceMatrix_row1,
        instanceMatrix_row2,
        instanceMatrix_row3
    );

    // Apply instance matrix to the vertex position
    vec3 instanced_pos = vec3(instanceMatrix * vec4(aPos, 1.0));

    // World space center and radius for the sphere to be ray-traced by the fragment shader.
    vec3 actualSphereCenter_world_space = vec3(model * instanceMatrix[3]); // Use the translation part of the instance matrix
    sphere_obj_world = vec4(actualSphereCenter_world_space, instanceActualSphereRadius);

    // Scale the base mesh vertices by the instance's actual radius
    vec3 scaled_aPos = instanced_pos * instanceActualSphereRadius;

    fpos_world = vec3(model * vec4(scaled_aPos, 1.0));
    fColor = instanceColor;
    gl_Position = projection * view * model * vec4(scaled_aPos, 1.0);
}