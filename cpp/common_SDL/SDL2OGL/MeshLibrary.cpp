#include "MeshLibrary.h"
#include "Vec3.h"

static GLMesh<MPOS> makePointCross() {
    GLMesh<MPOS> mesh(GL_LINES);
    mesh.addVertex({-1, 0, 0});
    mesh.addVertex({ 1, 0, 0});
    mesh.addVertex({ 0,-1, 0});
    mesh.addVertex({ 0, 1, 0});
    mesh.addVertex({ 0, 0,-1});
    mesh.addVertex({ 0, 0, 1});
    return mesh;
}

static GLMesh<MPOS> makeLine() {
    GLMesh<MPOS> mesh(GL_LINES);
    mesh.addVertex({0, 0, 0});
    mesh.addVertex({1, 1, 1});
    return mesh;
}

static GLMesh<MPOS> makeLine2D() {
    GLMesh<MPOS> mesh(GL_LINES);
    mesh.addVertex({0, 0, 0});
    mesh.addVertex({1, 1, 0});
    return mesh;
}

static GLMesh<MPOS> makePoint() {
    GLMesh<MPOS> m = GLMesh<MPOS>(GL_POINTS);
    m.addVertex({0, 0, 0});
    return m;
}

static GLMesh<MPOS> makeWireCube() {
    GLMesh<MPOS> m = GLMesh<MPOS>(GL_LINES);
    // Bottom face
    m.addVertex({0, 0, 0}); m.addVertex({1, 0, 0});
    m.addVertex({1, 0, 0}); m.addVertex({1, 0, 1});
    m.addVertex({1, 0, 1}); m.addVertex({0, 0, 1});
    m.addVertex({0, 0, 1}); m.addVertex({0, 0, 0});
    
    // Top face
    m.addVertex({0, 1, 0}); m.addVertex({1, 1, 0});
    m.addVertex({1, 1, 0}); m.addVertex({1, 1, 1});
    m.addVertex({1, 1, 1}); m.addVertex({0, 1, 1});
    m.addVertex({0, 1, 1}); m.addVertex({0, 1, 0});
    
    // Connecting vertical lines
    m.addVertex({0, 0, 0}); m.addVertex({0, 1, 0});
    m.addVertex({1, 0, 0}); m.addVertex({1, 1, 0});
    m.addVertex({1, 0, 1}); m.addVertex({1, 1, 1});
    m.addVertex({0, 0, 1}); m.addVertex({0, 1, 1});
    
    return m;
}

static GLMesh<MPOS,MNORMAL> makeCubeWithNormals() {
    GLMesh<MPOS,MNORMAL> m = GLMesh<MPOS,MNORMAL>(GL_TRIANGLES);

    // positive x face
    m.addVertex({1, 0, 0}, {1, 0, 0});
    m.addVertex({1, 1, 0}, {1, 0, 0});
    m.addVertex({1, 1, 1}, {1, 0, 0});
    m.addVertex({1, 0, 0}, {1, 0, 0});
    m.addVertex({1, 1, 1}, {1, 0, 0});
    m.addVertex({1, 0, 1}, {1, 0, 0});

    // negative x face
    m.addVertex({0, 0, 0}, {-1, 0, 0});
    m.addVertex({0, 1, 0}, {-1, 0, 0});
    m.addVertex({0, 1, 1}, {-1, 0, 0});
    m.addVertex({0, 0, 0}, {-1, 0, 0});
    m.addVertex({0, 1, 1}, {-1, 0, 0});
    m.addVertex({0, 0, 1}, {-1, 0, 0});

    // positive y face
    m.addVertex({0, 1, 0}, {0, 1, 0});
    m.addVertex({1, 1, 0}, {0, 1, 0});
    m.addVertex({1, 1, 1}, {0, 1, 0});
    m.addVertex({0, 1, 0}, {0, 1, 0});
    m.addVertex({1, 1, 1}, {0, 1, 0});
    m.addVertex({0, 1, 1}, {0, 1, 0});

    // negative y face
    m.addVertex({0, 0, 0}, {0, -1, 0});
    m.addVertex({1, 0, 0}, {0, -1, 0});
    m.addVertex({1, 0, 1}, {0, -1, 0});
    m.addVertex({0, 0, 0}, {0, -1, 0});
    m.addVertex({1, 0, 1}, {0, -1, 0});
    m.addVertex({0, 0, 1}, {0, -1, 0});

    // positive z face
    m.addVertex({0, 0, 1}, {0, 0, 1});
    m.addVertex({1, 0, 1}, {0, 0, 1});
    m.addVertex({1, 1, 1}, {0, 0, 1});
    m.addVertex({0, 0, 1}, {0, 0, 1});
    m.addVertex({1, 1, 1}, {0, 0, 1});
    m.addVertex({0, 1, 1}, {0, 0, 1});

    // negative z face
    m.addVertex({0, 0, 0}, {0, 0, -1});
    m.addVertex({1, 0, 0}, {0, 0, -1});
    m.addVertex({1, 1, 0}, {0, 0, -1});
    m.addVertex({0, 0, 0}, {0, 0, -1});
    m.addVertex({1, 1, 0}, {0, 0, -1});
    m.addVertex({0, 1, 0}, {0, 0, -1});

    return m;
}

static GLMesh<MPOS> makeRectMesh(){
    GLMesh<MPOS> m = GLMesh<MPOS>(GL_TRIANGLE_FAN);
    m.addVertex( {0, 0, 0} );
    m.addVertex( {0, 1, 0} );
    m.addVertex( {1, 1, 0} );
    m.addVertex( {1, 0, 0} );
    return m;
}

static const GLMesh<MPOS> makeCircleMesh(){
	GLMesh<MPOS> m = GLMesh<MPOS>(GL_LINE_LOOP);
	const int n = 64;
	float dphi =  6.28318530718f / n;
	Vec2f drot; drot.fromAngle( dphi );
	Vec2f v = {1, 0};
	for ( int i=0; i<n; i++ ){
		m.addVertex( {v.x, v.y, 0} );
		v.mul_cmplx( drot );
	}
	return m;
}

static const char* vertexSphere = R"(
#version 300 es

in highp vec3 vPosition;

out highp vec3 fPos;
out highp vec3 fDir;

uniform mat4 uMVPinv;
uniform mat4 uMVPMatrix;
uniform vec3 uPos;
uniform float uRadius;

vec4 model2clip(vec4 pos){
    return uMVPMatrix * (pos*uRadius + vec4(uPos, 0.0)*pos.w);
}

vec4 clip2model(vec4 pos){
    return ((uMVPinv * pos) - vec4(uPos, 0.0)*pos.w) / uRadius;
}

void main(){
    mediump vec3 sphPos = model2clip(vec4(0.0, 0.0, 0.0, 1.0)).xyz;
    mediump vec3 vecRight = normalize(clip2model(vec4(1.0, 0.0, 0.0, 0.0)).xyz);
    mediump vec3 vecUp    = normalize(clip2model(vec4(0.0, 1.0, 0.0, 0.0)).xyz);
    
    vecRight = model2clip(vec4(vecRight, 0.0)).xyz;
    vecUp    = model2clip(vec4(vecUp   , 0.0)).xyz;

    mediump vec3 newPos = vecRight*vPosition.x + vecUp*vPosition.y;
    gl_Position = vec4(newPos + sphPos, 1.0);

    fPos = clip2model(vec4(gl_Position.xy, 0.0, 1.0)).xyz;
    fDir = clip2model(vec4(0.0, 0.0, 2.0, 0.0)).xyz;
}
)";

static const char* fragSphere = R"(
#version 300 es

in highp vec3 fPos;
in highp vec3 fDir;

layout(location=0) out mediump vec4 FragColor;

uniform mediump vec3 uColor;

void main(){
    FragColor = vec4(uColor, 1.0);

    highp float DirLen = length(fDir);
    highp vec3 Dir = fDir / DirLen;
    highp float Tc = -dot(fPos, Dir);
    highp float dsqr = dot(fPos, fPos) - Tc*Tc;
    if (dsqr >= 1.0) discard;
    highp float T1c = sqrt(1.0 - dsqr);

    highp float T1 = Tc - T1c;

    highp vec3 normal = fPos + T1*Dir;
    gl_FragDepth = clamp((T1 / DirLen) + .5, 0.0, 1.0);

    mediump float light = dot(normal, vec3(1.0, -1.0, 1.0));
    light = (light+1.0)/2.0;
    light = .3 + light*.6;
    FragColor = FragColor*vec4(light, light, light, 1.0);
}
)";

static GLMeshBase<MPOS> makeSphere(){
    GLMeshBase<MPOS> m = GLMeshBase<MPOS>(GL_TRIANGLES, GL_STATIC_DRAW, new Shader<MPOS>(vertexSphere, fragSphere));
    m.addVertex({-1, -1, -1});
    m.addVertex({ 1, -1, -1});
    m.addVertex({-1,  1, -1});

    m.addVertex({-1,  1, -1});
    m.addVertex({ 1, -1, -1});
    m.addVertex({ 1,  1, -1});
    return m;
}

static GLMesh<MPOS> makeCross() {
    GLMesh<MPOS> m = GLMesh<MPOS>(GL_LINES);
    float sz = 0.5f;
    m.addVertex({-sz, 0, 0}); m.addVertex({sz, 0, 0});
    m.addVertex({0, -sz, 0}); m.addVertex({0, sz, 0});
    return m;
}

static GLMesh<MPOS> makeXMark() {
    GLMesh<MPOS> m = GLMesh<MPOS>(GL_LINES);
    float sz = 0.5f;
    m.addVertex({-sz, -sz, 0}); m.addVertex({sz, sz, 0});
    m.addVertex({-sz, sz, 0}); m.addVertex({sz, -sz, 0});
    return m;
}

namespace MeshLibrary {
    GLMesh<MPOS> pointCross = makePointCross();
    GLMesh<MPOS> line = makeLine();
    GLMesh<MPOS> line2D = makeLine2D();
    GLMesh<MPOS> point = makePoint();
    GLMesh<MPOS> wireCube = makeWireCube();
    GLMesh<MPOS,MNORMAL> cubeWithNormals = makeCubeWithNormals();
    GLMesh<MPOS> rect = makeRectMesh();
    GLMesh<MPOS> circle = makeCircleMesh();
    GLMeshBase<MPOS> sphere = makeSphere();
    GLMesh<MPOS> cross = makeCross();
    GLMesh<MPOS> xmark = makeXMark();
}
