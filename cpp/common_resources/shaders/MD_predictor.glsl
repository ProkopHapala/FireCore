precision highp float;
uniform sampler2D tPos;
uniform sampler2D tVel;
uniform float dt;
uniform int texWidth;

layout(location = 0) out vec4 outSn;

void main() {
    ivec2 coord = ivec2(gl_FragCoord.xy);
    // 1D Texture Logic: Y is 0
    vec4 pData = texelFetch(tPos, ivec2(coord.x, 0), 0);
    vec3 p = pData.rgb;
    float m = pData.a;
    
    vec3 v = texelFetch(tVel, ivec2(coord.x, 0), 0).rgb;
    
    // Simple gravity
    vec3 gravity = vec3(0.0, -9.81, 0.0);
    
    // s_n = p + v*dt + F_ext/m * dt^2
    vec3 sn = p + v * dt;
    if (m > 0.0) {
       sn += gravity * dt * dt;
    }
    
    outSn = vec4(sn, m);
}