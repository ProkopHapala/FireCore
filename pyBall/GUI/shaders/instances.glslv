#version 330 core

layout (location = 0) in vec3 aPos;      // Vertex position of the base sphere mesh
layout (location = 1) in vec3 aNormal;   // Vertex normal of the base sphere mesh (not used by raytracer directly)

// Per-instance attributes
layout (location = 2) in vec3  ins_pos; // Center of the sphere instance (model space)
layout (location = 3) in float ins_rad;  // Actual radius of the sphere to ray-trace
layout (location = 4) in vec4  ins_col;    // Color (RGBA) of the sphere instance

// Outputs to Fragment Shader
out vec3 fpos;        // Fragment position on the bounding mesh in world space
out vec4 opos;  // Sphere center (xyz) and radius (w) in world space
out vec4 col;     // Pass through atom color

uniform mat4 projection;
uniform mat4 view;
uniform mat4 model; // Overall model orientation (trackball)

void main()
{
    // World space center and radius for the sphere to be ray-traced by the fragment shader.
    // ins_rad is the small radius (e.g., 0.005 for an atom).
    vec3 actualSphereCenter_world_space = vec3(model * vec4(ins_pos, 1.0));
    opos = vec4(actualSphereCenter_world_space, ins_rad);

    // Create a model matrix for the bounding box mesh (whose vertices are in aPos).
    // aPos comes from a mesh already defined with its bounding box radius (e.g., 1.5 in octahedron_sphere_mesh).
    // This matrix only translates the bounding box to the instancePosition_model.
    // The overall 'model' matrix (from trackball) will handle rotation.
    mat4 boundingBoxTransformMatrix = mat4(1.0);
    boundingBoxTransformMatrix[3] = vec4(ins_pos, 1.0); // Translate
    mat4 finalBoundingBoxModelMatrix = model * boundingBoxTransformMatrix;

    // Scale the base mesh vertices by the instance's actual radius
    vec3 scaled_aPos = aPos * ins_rad;

    fpos = vec3(finalBoundingBoxModelMatrix * vec4(scaled_aPos, 1.0));
    col  = ins_col;
    gl_Position = projection * view * finalBoundingBoxModelMatrix * vec4(scaled_aPos, 1.0);
}