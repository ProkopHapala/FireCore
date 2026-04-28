#version 330 core

layout (location = 0) in vec3 aPos3D;       // 3D base position of the character quad (per instance)
layout (location = 1) in vec2 aLocalOffset; // Local 2D offset within the character quad (-0.5 to 0.5)
layout (location = 2) in vec2 aTexCoord;    // UV coordinates for texture sampling

uniform mat4  projection;
uniform mat4  view;
uniform mat4  model;
uniform float labelScale; // Global scale for all labels
uniform vec2  offset; 

out vec2 v_texCoord;

void main() {
    v_texCoord = aTexCoord;

    // To create a billboard, we transform the local offset by the inverse of the view matrix's rotation part.
    // This aligns the quad with the camera's view plane.
    vec3 cameraRight = vec3(view[0][0], view[1][0], view[2][0])*0.1;
    vec3 cameraUp    = vec3(view[0][1], view[1][1], view[2][1])*0.1;

    vec2 off = (offset + aLocalOffset)*labelScale;
    
    // Calculate the final world position for the vertex
    vec3 finalWorldPos = (model * vec4(aPos3D, 1.0)).xyz + (cameraRight * off.x) + (cameraUp * off.y);

    gl_Position = projection * view * vec4(finalWorldPos, 1.0);

    //gl_Position = vec4(aLocalOffset*10.0, 0.0, 1.0);

    //gl_Position = vec4( aPos3D + vec3( aLocalOffset.xyx ), 1.0 );
}
