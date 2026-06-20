#!/usr/bin/env python3
"""Generate a standalone WebGL 2.0 HTML page for interactive single-crystal XRD.

Reads an XYZ file and embeds atom data + GLSL shaders into a single HTML file
that can be opened directly in a browser (no server required).

Usage:
    python scripts/generate_xrd_webgl.py \
        --xyz tests/tSiNCs/fixtures/vibration_parallel/structures/diamond_nc_R6_relaxed.xyz \
        --out OUT_XRD/diamond_nc_R6_single_crystal.html
"""
import argparse
import json
import os
import numpy as np


GLSL_VERT = """#version 300 es
in vec2 a_position;
void main() {
    gl_Position = vec4(a_position, 0.0, 1.0);
}
"""

GLSL_FRAG_TEMPLATE = """#version 300 es
precision highp float;
precision highp int;

out vec4 outColor;

uniform vec4 u_atoms[{MAX_ATOMS}];   // xyz + Z (as float)
uniform int  u_nAtoms;
uniform float u_wavelength;          // Å
uniform vec3  u_beamDir;             // unit incident direction
uniform mat3  u_rot;                 // crystal rotation
uniform float u_detectorL;         // detector distance (Å)
uniform float u_pixelSize;           // Å per pixel
uniform vec2  u_detectorSize;      // width, height in pixels
uniform float u_exposure;          // brightness scale
uniform float u_bloom;             // Gaussian bloom radius in pixels

float getFF(int Z, float Q) {{
    // Simplified Cromer-Mann-ish: f ≈ Z * exp(-b * s²) with s = Q/(4π)
    float s2 = Q * Q * 0.00633257; // (1/(4π))²
    float f = 0.0;
    if (Z == 6) {{      // C
        f = 2.31*exp(-20.8439*s2) + 1.02*exp(-10.2075*s2)
          + 1.5886*exp(-0.5687*s2) + 0.865*exp(-51.6512*s2) + 0.2156;
    }} else if (Z == 1) {{ // H
        f = 0.489918*exp(-20.6593*s2) + 0.262003*exp(-7.74039*s2)
          + 0.196767*exp(-49.5519*s2) + 0.049879*exp(-2.20159*s2) + 0.001305;
    }} else if (Z == 14) {{ // Si
        f = 6.2925*exp(-2.4386*s2) + 3.0353*exp(-32.3337*s2)
          + 1.9891*exp(-0.6785*s2) + 1.541*exp(-81.6937*s2) + 1.1407;
    }}
    return f;
}}

void main() {{
    vec2 pix = gl_FragCoord.xy - u_detectorSize * 0.5;
    float x = pix.x * u_pixelSize;
    float y = pix.y * u_pixelSize;

    // Build detector basis u,v perpendicular to beam
    vec3 u;
    if (abs(u_beamDir.z) < 0.9) {{
        u = normalize(cross(vec3(0.0, 0.0, 1.0), u_beamDir));
    }} else {{
        u = normalize(cross(vec3(0.0, 1.0, 0.0), u_beamDir));
    }}
    vec3 v = cross(u_beamDir, u);

    float k = 6.28318530718 / u_wavelength;
    vec3 kOutDir = normalize(u_beamDir * u_detectorL + u * x + v * y);
    vec3 kOut = kOutDir * k;
    vec3 kIn  = u_beamDir * k;
    vec3 Qvec = kOut - kIn;
    float Qmag = length(Qvec);

    vec2 F = vec2(0.0, 0.0);
    for (int i = 0; i < {MAX_ATOMS}; i++) {{
        if (i >= u_nAtoms) break;
        vec4 a = u_atoms[i];
        vec3 r = u_rot * a.xyz;
        int Z = int(a.w + 0.5);
        float f = getFF(Z, Qmag);
        float phase = dot(Qvec, r);
        F.x += f * cos(phase);
        F.y += f * sin(phase);
    }}

    float I = F.x*F.x + F.y*F.y;

    // Optional bloom (simple box blur approximation)
    float intensity = I * u_exposure;
    float logI = log(1.0 + intensity);
    float bright = logI / (1.0 + logI); // soft saturate
    outColor = vec4(vec3(bright), 1.0);
}}
"""

HTML_TEMPLATE = """<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>XRD Single Crystal — {title}</title>
<style>
  body {{ margin: 0; background: #050505; color: #ddd; font-family: sans-serif; overflow: hidden; }}
  #glcanvas {{ display: block; width: 100vw; height: 100vh; }}
  #crystalview {{
    position: absolute; top: 10px; right: 10px;
    width: 260px; height: 260px;
    background: rgba(10,10,15,0.9);
    border: 1px solid #333; border-radius: 6px;
    z-index: 10;
  }}
  #info {{
    position: absolute; top: 10px; left: 10px;
    background: rgba(0,0,0,0.8); padding: 12px 16px;
    border-radius: 6px; font-size: 13px; line-height: 1.5;
    pointer-events: none; max-width: 320px;
  }}
  #info h3 {{ margin: 0 0 6px 0; font-size: 15px; color: #8cf; }}
  #info p {{ margin: 2px 0; }}
  #controls {{
    position: absolute; bottom: 10px; left: 10px;
    background: rgba(0,0,0,0.8); padding: 10px 14px;
    border-radius: 6px; font-size: 12px;
  }}
  #controls label {{ display: inline-block; min-width: 110px; }}
  #controls input {{ vertical-align: middle; }}
</style>
</head>
<body>
<canvas id="glcanvas"></canvas>
<canvas id="crystalview" width="260" height="260"></canvas>
<div id="info">
  <h3>{title}</h3>
  <p><b>natoms:</b> {n_atoms}  |  <b>bonds:</b> {n_bonds}</p>
  <p><b>λ:</b> {wavelength} Å  |  <b>L:</b> {detector_L} Å</p>
  <p><b>Drag</b> to rotate crystal</p>
  <p><b>Scroll</b> to zoom exposure</p>
  <p><b>R</b> to reset rotation</p>
</div>
<div id="controls">
  <label>Exposure</label><input type="range" id="exposure" min="-6" max="2" step="0.1" value="{exposure_log10}"><br>
  <label>Detector L</label><input type="range" id="detL" min="2" max="5" step="0.05" value="{det_L_log10}"><br>
  <label>Wavelength</label><input type="range" id="wavelength" min="0.5" max="3.0" step="0.01" value="{wavelength}"><br>
</div>
<script>
'use strict';

// ── embedded atom data ──
const ATOMS = {atoms_json};
const N_ATOMS = ATOMS.length;
const MAX_ATOMS = {max_atoms};

// ── WebGL 2.0 setup ──
const canvas = document.getElementById('glcanvas');
const gl = canvas.getContext('webgl2');
if (!gl) {{ alert('WebGL 2.0 not supported'); throw new Error('No WebGL2'); }}

function compile(type, src) {{
    const s = gl.createShader(type);
    gl.shaderSource(s, src);
    gl.compileShader(s);
    if (!gl.getShaderParameter(s, gl.COMPILE_STATUS)) {{
        console.error(gl.getShaderInfoLog(s));
        throw new Error('Shader compile failed');
    }}
    return s;
}}

const prog = gl.createProgram();
gl.attachShader(prog, compile(gl.VERTEX_SHADER, {vert_json}));
gl.attachShader(prog, compile(gl.FRAGMENT_SHADER, {frag_json}));
gl.linkProgram(prog);
if (!gl.getProgramParameter(prog, gl.LINK_STATUS)) {{
    console.error(gl.getProgramInfoLog(prog));
    throw new Error('Program link failed');
}}
gl.useProgram(prog);

// Full-screen quad
const buf = gl.createBuffer();
gl.bindBuffer(gl.ARRAY_BUFFER, buf);
gl.bufferData(gl.ARRAY_BUFFER, new Float32Array([-1,-1, 1,-1, -1,1, 1,1]), gl.STATIC_DRAW);
const aPos = gl.getAttribLocation(prog, 'a_position');
gl.enableVertexAttribArray(aPos);
gl.vertexAttribPointer(aPos, 2, gl.FLOAT, false, 0, 0);

// Uniform locations
const u_atoms      = gl.getUniformLocation(prog, 'u_atoms');
const u_nAtoms     = gl.getUniformLocation(prog, 'u_nAtoms');
const u_wavelength = gl.getUniformLocation(prog, 'u_wavelength');
const u_beamDir    = gl.getUniformLocation(prog, 'u_beamDir');
const u_rot        = gl.getUniformLocation(prog, 'u_rot');
const u_detectorL  = gl.getUniformLocation(prog, 'u_detectorL');
const u_pixelSize  = gl.getUniformLocation(prog, 'u_pixelSize');
const u_detectorSize = gl.getUniformLocation(prog, 'u_detectorSize');
const u_exposure   = gl.getUniformLocation(prog, 'u_exposure');

// Upload atoms (pad to MAX_ATOMS with zeros)
const atomArr = new Float32Array(MAX_ATOMS * 4);
for (let i = 0; i < N_ATOMS; i++) {{
    atomArr[i*4+0] = ATOMS[i][0];
    atomArr[i*4+1] = ATOMS[i][1];
    atomArr[i*4+2] = ATOMS[i][2];
    atomArr[i*4+3] = ATOMS[i][3];
}}
gl.uniform4fv(u_atoms, atomArr);
gl.uniform1i(u_nAtoms, N_ATOMS);

// ── state ──
let wavelength = {wavelength};
let detectorL  = {detector_L};
let pixelSize  = {pixel_size};
let exposure   = {exposure};
let beamDir    = [1.0, 0.0, 0.0];

// Rotation quaternion
let rotQuat = [0,0,0,1];

function quatToMat3(q) {{
    const [x,y,z,w] = q;
    const x2=x+x, y2=y+y, z2=z+z;
    const xx=x*x2, xy=x*y2, xz=x*z2;
    const yy=y*y2, yz=y*z2, zz=z*z2;
    const wx=w*x2, wy=w*y2, wz=w*z2;
    return new Float32Array([
        1.0-(yy+zz), xy+wz,       xz-wy,
        xy-wz,       1.0-(xx+zz), yz+wx,
        xz+wy,       yz-wx,       1.0-(xx+yy)
    ]);
}}

function normalize(v) {{
    const len = Math.sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
    return [v[0]/len, v[1]/len, v[2]/len];
}}

function quatMul(a, b) {{
    return [
        a[3]*b[0] + a[0]*b[3] + a[1]*b[2] - a[2]*b[1],
        a[3]*b[1] - a[0]*b[2] + a[1]*b[3] + a[2]*b[0],
        a[3]*b[2] + a[0]*b[1] - a[1]*b[0] + a[2]*b[3],
        a[3]*b[3] - a[0]*b[0] - a[1]*b[1] - a[2]*b[2]
    ];
}}

function axisAngleQuat(axis, angle) {{
    const half = angle * 0.5;
    const s = Math.sin(half);
    return [axis[0]*s, axis[1]*s, axis[2]*s, Math.cos(half)];
}}

function resize() {{
    const w = canvas.clientWidth;
    const h = canvas.clientHeight;
    if (canvas.width !== w || canvas.height !== h) {{
        canvas.width = w;
        canvas.height = h;
        gl.viewport(0, 0, w, h);
    }}
}}

function draw() {{
    resize();
    gl.uniform1f(u_wavelength, wavelength);
    gl.uniform3f(u_beamDir, beamDir[0], beamDir[1], beamDir[2]);
    gl.uniformMatrix3fv(u_rot, false, quatToMat3(rotQuat));
    gl.uniform1f(u_detectorL, detectorL);
    gl.uniform1f(u_pixelSize, pixelSize);
    gl.uniform2f(u_detectorSize, canvas.width, canvas.height);
    gl.uniform1f(u_exposure, exposure);
    gl.drawArrays(gl.TRIANGLE_STRIP, 0, 4);
    drawCrystalView();
}}

// ── Crystal view renderer (small overlay) ──
const cvCanvas = document.getElementById('crystalview');
const cvGL = cvCanvas.getContext('webgl2', {{ alpha: false, antialias: true }});
if (!cvGL) {{ console.warn('Crystal view: no WebGL2'); }}

const BONDS = {bonds_json};

const CV_ATOM_VERT = `#version 300 es
in vec3 a_pos;
uniform mat4 u_mvp;
void main() {{
    gl_Position = u_mvp * vec4(a_pos, 1.0);
    gl_PointSize = 6.0;
}}`;
const CV_ATOM_FRAG = `#version 300 es
precision highp float;
out vec4 outColor;
uniform vec3 u_color;
void main() {{
    float r = length(gl_PointCoord - 0.5);
    if (r > 0.5) discard;
    float alpha = 1.0 - smoothstep(0.3, 0.5, r);
    outColor = vec4(u_color, alpha);
}}`;
const CV_BOND_VERT = `#version 300 es
in vec3 a_pos;
uniform mat4 u_mvp;
void main() {{
    gl_Position = u_mvp * vec4(a_pos, 1.0);
}}`;
const CV_BOND_FRAG = `#version 300 es
precision highp float;
out vec4 outColor;
void main() {{
    outColor = vec4(0.5, 0.5, 0.5, 1.0);
}}`;

function cvCompile(gl, type, src) {{
    const s = gl.createShader(type);
    gl.shaderSource(s, src);
    gl.compileShader(s);
    if (!gl.getShaderParameter(s, gl.COMPILE_STATUS)) {{
        console.error(gl.getShaderInfoLog(s));
        throw new Error('CV shader compile failed');
    }}
    return s;
}}

let atomProg = null, bondProg = null, atomBuf = null, bondBuf = null;
if (cvGL) {{
    atomProg = cvGL.createProgram();
    cvGL.attachShader(atomProg, cvCompile(cvGL, cvGL.VERTEX_SHADER, CV_ATOM_VERT));
    cvGL.attachShader(atomProg, cvCompile(cvGL, cvGL.FRAGMENT_SHADER, CV_ATOM_FRAG));
    cvGL.linkProgram(atomProg);

    bondProg = cvGL.createProgram();
    cvGL.attachShader(bondProg, cvCompile(cvGL, cvGL.VERTEX_SHADER, CV_BOND_VERT));
    cvGL.attachShader(bondProg, cvCompile(cvGL, cvGL.FRAGMENT_SHADER, CV_BOND_FRAG));
    cvGL.linkProgram(bondProg);

    // Atom positions buffer (interleaved: x,y,z for each atom)
    const atomPos = new Float32Array(N_ATOMS * 3);
    for (let i = 0; i < N_ATOMS; i++) {{
        atomPos[i*3+0] = ATOMS[i][0];
        atomPos[i*3+1] = ATOMS[i][1];
        atomPos[i*3+2] = ATOMS[i][2];
    }}
    atomBuf = cvGL.createBuffer();
    cvGL.bindBuffer(cvGL.ARRAY_BUFFER, atomBuf);
    cvGL.bufferData(cvGL.ARRAY_BUFFER, atomPos, cvGL.STATIC_DRAW);

    // Bond buffer (line segments: start,end for each bond)
    if (BONDS.length > 0) {{
        const bondPos = new Float32Array(BONDS.length * 6);
        for (let b = 0; b < BONDS.length; b++) {{
            const i = BONDS[b][0];
            const j = BONDS[b][1];
            bondPos[b*6+0] = ATOMS[i][0];
            bondPos[b*6+1] = ATOMS[i][1];
            bondPos[b*6+2] = ATOMS[i][2];
            bondPos[b*6+3] = ATOMS[j][0];
            bondPos[b*6+4] = ATOMS[j][1];
            bondPos[b*6+5] = ATOMS[j][2];
        }}
        bondBuf = cvGL.createBuffer();
        cvGL.bindBuffer(cvGL.ARRAY_BUFFER, bondBuf);
        cvGL.bufferData(cvGL.ARRAY_BUFFER, bondPos, cvGL.STATIC_DRAW);
    }}
}}

function drawCrystalView() {{
    if (!cvGL) return;
    cvGL.viewport(0, 0, cvCanvas.width, cvCanvas.height);
    cvGL.clearColor(0.04, 0.04, 0.06, 1.0);
    cvGL.clear(cvGL.COLOR_BUFFER_BIT);
    cvGL.disable(cvGL.DEPTH_TEST);

    // Compute bounding box to set scale
    let xmin = 1e9, xmax = -1e9, ymin = 1e9, ymax = -1e9, zmin = 1e9, zmax = -1e9;
    for (let i = 0; i < N_ATOMS; i++) {{
        xmin = Math.min(xmin, ATOMS[i][0]); xmax = Math.max(xmax, ATOMS[i][0]);
        ymin = Math.min(ymin, ATOMS[i][1]); ymax = Math.max(ymax, ATOMS[i][1]);
        zmin = Math.min(zmin, ATOMS[i][2]); zmax = Math.max(zmax, ATOMS[i][2]);
    }}
    const cx = (xmin + xmax) * 0.5;
    const cy = (ymin + ymax) * 0.5;
    const cz = (zmin + zmax) * 0.5;
    const maxR = Math.max(xmax - xmin, ymax - ymin, zmax - zmin) * 0.5;
    const scale = maxR > 0 ? 0.85 / maxR : 1.0;

    // Build MVP = S * R * T(-center) in column-major order for WebGL
    // quatToMat3 returns column-major data, so R[0..2]=col0, R[3..5]=col1, R[6..8]=col2
    const R = quatToMat3(rotQuat);
    const tx = -(R[0]*cx + R[3]*cy + R[6]*cz) * scale;
    const ty = -(R[1]*cx + R[4]*cy + R[7]*cz) * scale;
    const tz = -(R[2]*cx + R[5]*cy + R[8]*cz) * scale;
    const mvpRot = new Float32Array([
        R[0]*scale, R[1]*scale, R[2]*scale, 0,   // col 0
        R[3]*scale, R[4]*scale, R[5]*scale, 0,   // col 1
        R[6]*scale, R[7]*scale, R[8]*scale, 0,   // col 2
        tx,         ty,         tz,         1.0  // col 3 (translation)
    ]);

    // Draw bonds first
    if (bondBuf && BONDS.length > 0) {{
        cvGL.useProgram(bondProg);
        cvGL.uniformMatrix4fv(cvGL.getUniformLocation(bondProg, 'u_mvp'), false, mvpRot);
        cvGL.bindBuffer(cvGL.ARRAY_BUFFER, bondBuf);
        const aPosBond = cvGL.getAttribLocation(bondProg, 'a_pos');
        cvGL.enableVertexAttribArray(aPosBond);
        cvGL.vertexAttribPointer(aPosBond, 3, cvGL.FLOAT, false, 0, 0);
        cvGL.drawArrays(cvGL.LINES, 0, BONDS.length * 2);
    }}

    // Draw atoms — group by element for different colors
    cvGL.useProgram(atomProg);
    cvGL.uniformMatrix4fv(cvGL.getUniformLocation(atomProg, 'u_mvp'), false, mvpRot);
    cvGL.bindBuffer(cvGL.ARRAY_BUFFER, atomBuf);
    const aPosAtom = cvGL.getAttribLocation(atomProg, 'a_pos');
    cvGL.enableVertexAttribArray(aPosAtom);
    cvGL.vertexAttribPointer(aPosAtom, 3, cvGL.FLOAT, false, 0, 0);
    // All atoms same color for now (green for carbon)
    const ZtoColor = {{6: [0.3, 0.8, 0.3], 1: [0.8, 0.3, 0.3], 14: [0.5, 0.5, 0.9]}};
    const col = ZtoColor[ATOMS[0][3]] || [0.7, 0.7, 0.7];
    cvGL.uniform3f(cvGL.getUniformLocation(atomProg, 'u_color'), col[0], col[1], col[2]);
    cvGL.drawArrays(cvGL.POINTS, 0, N_ATOMS);
}}

// ── mouse controls ──
let dragging = false;
let lastX = 0, lastY = 0;

canvas.addEventListener('mousedown', e => {{ dragging = true; lastX = e.clientX; lastY = e.clientY; }});
window.addEventListener('mouseup', () => {{ dragging = false; }});
window.addEventListener('mousemove', e => {{
    if (!dragging) return;
    const dx = e.clientX - lastX;
    const dy = e.clientY - lastY;
    lastX = e.clientX;
    lastY = e.clientY;
    const sens = 0.005;
    const qx = axisAngleQuat([0,1,0], dx * sens);
    const qy = axisAngleQuat([1,0,0], dy * sens);
    rotQuat = quatMul(quatMul(qx, qy), rotQuat);
    draw();
}});

canvas.addEventListener('wheel', e => {{
    e.preventDefault();
    const delta = Math.exp(-e.deltaY * 0.001);
    exposure *= delta;
    document.getElementById('exposure').value = Math.log10(exposure);
    draw();
}}, {{passive: false}});

window.addEventListener('keydown', e => {{
    if (e.key === 'r' || e.key === 'R') {{
        rotQuat = [0,0,0,1];
        draw();
    }}
}});

// ── UI controls ──
document.getElementById('exposure').addEventListener('input', e => {{
    exposure = Math.pow(10, parseFloat(e.target.value));
    draw();
}});
document.getElementById('detL').addEventListener('input', e => {{
    detectorL = Math.pow(10, parseFloat(e.target.value));
    draw();
}});
document.getElementById('wavelength').addEventListener('input', e => {{
    wavelength = parseFloat(e.target.value);
    draw();
}});

// ── animate ──
function loop() {{
    draw();
    requestAnimationFrame(loop);
}}
requestAnimationFrame(loop);

</script>
</body>
</html>
"""


def load_xyz(path: str):
    with open(path, 'r') as f:
        lines = f.readlines()
    n = int(lines[0].strip())
    atoms = []
    symbol_to_Z = {'H': 1, 'C': 6, 'Si': 14}
    for line in lines[2:2 + n]:
        parts = line.split()
        sym = parts[0]
        Z = symbol_to_Z.get(sym, 0)
        if Z == 0:
            raise ValueError(f"Unknown element symbol: {sym}")
        x, y, z = float(parts[1]), float(parts[2]), float(parts[3])
        atoms.append([x, y, z, Z])
    return atoms


def build_bonds(atoms, max_bond={'CC': 1.8, 'CH': 1.2, 'HH': 1.0, 'SiSi': 2.4, 'SiC': 2.0, 'SiH': 1.6}):
    """Return list of bond indices [i,j] for pairs within covalent distance."""
    n = len(atoms)
    bonds = []
    elems = {1: 'H', 6: 'C', 14: 'Si'}
    for i in range(n):
        zi = int(atoms[i][3])
        ei = elems.get(zi, 'X')
        for j in range(i + 1, n):
            zj = int(atoms[j][3])
            ej = elems.get(zj, 'X')
            key = ''.join(sorted([ei, ej]))
            thr = max_bond.get(key, 2.0)
            dx = atoms[i][0] - atoms[j][0]
            dy = atoms[i][1] - atoms[j][1]
            dz = atoms[i][2] - atoms[j][2]
            r = np.sqrt(dx * dx + dy * dy + dz * dz)
            if r < thr:
                bonds.append([i, j])
    return bonds


def main():
    ap = argparse.ArgumentParser(description='Generate WebGL 2.0 single-crystal XRD viewer')
    ap.add_argument('--xyz', required=True, help='Input XYZ file')
    ap.add_argument('--out', required=True, help='Output HTML file')
    ap.add_argument('--wavelength', type=float, default=1.54, help='X-ray wavelength (Å)')
    ap.add_argument('--detector-L', type=float, default=5000.0, help='Detector distance (Å)')
    ap.add_argument('--pixel-size', type=float, default=20.0, help='Detector pixel size (Å)')
    ap.add_argument('--exposure', type=float, default=1e-3, help='Intensity scale factor')
    ap.add_argument('--title', default=None, help='Page title')
    args = ap.parse_args()

    atoms = load_xyz(args.xyz)
    bonds = build_bonds(atoms)
    n_atoms = len(atoms)
    max_atoms = 1
    while max_atoms < n_atoms:
        max_atoms *= 2
    if max_atoms < 256:
        max_atoms = 256

    title = args.title or os.path.splitext(os.path.basename(args.xyz))[0]

    vert_src = GLSL_VERT
    frag_src = GLSL_FRAG_TEMPLATE.format(MAX_ATOMS=max_atoms)

    html = HTML_TEMPLATE.format(
        title=title,
        n_atoms=n_atoms,
        n_bonds=len(bonds),
        wavelength=args.wavelength,
        detector_L=args.detector_L,
        pixel_size=args.pixel_size,
        exposure=args.exposure,
        exposure_log10=round(np.log10(args.exposure), 2) if args.exposure > 0 else -3,
        det_L_log10=round(np.log10(args.detector_L), 2),
        max_atoms=max_atoms,
        atoms_json=json.dumps(atoms),
        bonds_json=json.dumps(bonds),
        vert_json=json.dumps(vert_src),
        frag_json=json.dumps(frag_src),
    )

    os.makedirs(os.path.dirname(args.out) or '.', exist_ok=True)
    with open(args.out, 'w') as f:
        f.write(html)
    print(f"Wrote {args.out}  ({n_atoms} atoms, max_atoms={max_atoms})")


if __name__ == '__main__':
    main()
