precision highp float;
uniform sampler2D tPosOld;
uniform sampler2D tPosNew;
uniform float dt;

layout(location = 0) out vec4 outVel;

void main() {
    ivec2 coord = ivec2(gl_FragCoord.xy);
    vec3 p_old = texelFetch(tPosOld, coord, 0).rgb;
    vec3 p_new = texelFetch(tPosNew, coord, 0).rgb;
    vec3 v = (p_new - p_old) / dt;
    // Damping
    v *= 0.99; 
    outVel = vec4(v, 0.0);
}