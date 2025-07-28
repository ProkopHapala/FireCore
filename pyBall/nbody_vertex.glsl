#version 330 core
layout (location = 0) in vec4 position;
void main() { gl_Position = vec4(position.xyz, 1.0); gl_PointSize = 2.0; }
