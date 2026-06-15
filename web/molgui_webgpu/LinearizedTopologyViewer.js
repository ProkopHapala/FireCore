/// Three.js HTML debug viewers for linearized MMFFL stick topology (not topology math).
import { exportTopologySpringList } from './MMFFLTopology.js';

const STICK_VIEWER_SCENE_JS = `
const TAG_BY_CLASS = { 1: 'bond', 2: 'angle', 3: 'dihedral' };
const COLOR_MAP = { bond: 0xff3333, angle: 0x33cc33, dihedral: 0x3388ff };
const ATOM_COLORS = { C: 0x909090, Si: 0xdaa520, H: 0xffffff, N: 0x3050f8, O: 0xff0d0d };
const ATOM_RAD = { H: 0.07, C: 0.15, Si: 0.18, N: 0.14, O: 0.13 };
function atomRadius(sym) { return ATOM_RAD[sym] || 0.12; }
function symFromZ(z) { return ({1:'H',6:'C',7:'N',8:'O',14:'Si'})[z|0] || 'X'; }

function buildScene(DATA) {
  const scene = new THREE.Scene();
  scene.background = new THREE.Color(0xf0f0f0);
  const stickGroups = { bond: new THREE.Group(), angle: new THREE.Group(), dihedral: new THREE.Group() };
  scene.add(stickGroups.bond, stickGroups.angle, stickGroups.dihedral);
  let cx = 0, cy = 0, cz = 0;
  for (const a of DATA.atoms) { cx += a.x; cy += a.y; cz += a.z; }
  cx /= DATA.atoms.length; cy /= DATA.atoms.length; cz /= DATA.atoms.length;
  for (const a of DATA.atoms) {
    const r = atomRadius(a.sym);
    const geo = new THREE.SphereGeometry(r, 14, 14);
    const m = new THREE.Mesh(geo, new THREE.MeshPhongMaterial({ color: ATOM_COLORS[a.sym] || 0x888888 }));
    m.position.set(a.x - cx, a.y - cy, a.z - cz);
    scene.add(m);
  }
  for (const s of DATA.springs) {
    const a = DATA.atoms[s.i], b = DATA.atoms[s.j];
    const p1 = new THREE.Vector3(a.x - cx, a.y - cy, a.z - cz);
    const p2 = new THREE.Vector3(b.x - cx, b.y - cy, b.z - cz);
    const geo = new THREE.BufferGeometry().setFromPoints([p1, p2]);
    const col = COLOR_MAP[s.tag] || 0x888888;
    const line = new THREE.Line(geo, new THREE.LineBasicMaterial({ color: col }));
    (stickGroups[s.tag] || stickGroups.bond).add(line);
  }
  scene.add(new THREE.AmbientLight(0xffffff, 0.6));
  const dl = new THREE.DirectionalLight(0xffffff, 0.8);
  dl.position.set(5, 8, 6);
  scene.add(dl);
  const box = new THREE.Box3();
  scene.traverse(o => { if (o.isMesh || o.isLine) box.expandByObject(o); });
  const size = box.getSize(new THREE.Vector3()).length() || 10;
  return { scene, stickGroups, size };
}

function wireLegend(stickGroups) {
  const bind = (id, tag) => {
    const el = document.getElementById(id);
    if (!el) return;
    const apply = () => { stickGroups[tag].visible = el.checked; };
    el.addEventListener('change', apply);
    apply();
  };
  bind('tog-bond', 'bond');
  bind('tog-angle', 'angle');
  bind('tog-dihedral', 'dihedral');
}

function springsFromNpz(npz) {
  const n = npz.natoms | 0, m = npz.max_neighbors | 0;
  const pos = npz.pos, Z = npz.Z, neigh = npz.neigh_idx, cls = npz.stick_class;
  const atoms = [];
  for (let i = 0; i < n; i++) atoms.push({ x: pos[i*3], y: pos[i*3+1], z: pos[i*3+2], sym: symFromZ(Z[i]) });
  const seen = new Set();
  const springs = [];
  for (let i = 0; i < n; i++) {
    for (let k = 0; k < m; k++) {
      const j = neigh[i*m+k] | 0;
      if (j < 0) break;
      const c = cls[i*m+k] | 0;
      const tag = TAG_BY_CLASS[c];
      if (!tag) continue;
      const a = i < j ? i : j, b = i < j ? j : i;
      const pk = a + ',' + b + ',' + tag;
      if (seen.has(pk)) continue;
      seen.add(pk);
      springs.push({ i: a, j: b, tag });
    }
  }
  return { atoms, springs };
}
`;

function legendHtml(title) {
    return `<div id="legend"><b>${title}</b><br>
<label class="tog"><input type="checkbox" id="tog-bond" checked><span class="sw bond"></span> K₁₂ bond</label>
<label class="tog"><input type="checkbox" id="tog-angle" checked><span class="sw angle"></span> K₁₃ angle</label>
<label class="tog"><input type="checkbox" id="tog-dihedral" checked><span class="sw dihedral"></span> K₁₄ dihedral</label>
</div>`;
}

const STICK_VIEWER_CSS = `body{margin:0;overflow:hidden;font-family:sans-serif}
#legend{position:absolute;top:8px;left:8px;background:rgba(255,255,255,.9);padding:8px 12px;border-radius:6px;font-size:13px;line-height:1.7;user-select:none}
.tog{display:block;cursor:pointer;margin:2px 0}
.tog input{margin-right:6px;vertical-align:middle;cursor:pointer}
.sw{display:inline-block;width:12px;height:12px;margin-right:4px;vertical-align:middle;border-radius:2px}
.sw.bond{background:#c33}.sw.angle{background:#3a3}.sw.dihedral{background:#38f}
#status{position:absolute;bottom:8px;left:8px;background:rgba(255,255,255,.85);padding:4px 8px;border-radius:4px;font-size:12px}
#fileRow{position:absolute;top:8px;right:8px;background:rgba(255,255,255,.9);padding:6px 10px;border-radius:6px;font-size:12px}`;

function topoViewerPayload(topo, title) {
    const springs = exportTopologySpringList(topo);
    const atoms = [];
    for (let i = 0; i < (topo.n_real | 0); i++) {
        const p = topo.apos[i];
        const sym = (topo.type_names && topo.type_names[i]) ? String(topo.type_names[i]).split('_')[0] : 'X';
        atoms.push({ x: +p[0], y: +p[1], z: +p[2], sym });
    }
    return JSON.stringify({ title, atoms, springs });
}

/// Self-contained HTML viewer with embedded topology JSON.
export function exportStickViewerHTML(topo, title = 'MMFFL topology') {
    const payload = topoViewerPayload(topo, title);
    return `<!DOCTYPE html>
<html><head><meta charset="utf-8"><title>${title}</title>
<style>${STICK_VIEWER_CSS}</style>
<script type="importmap">{"imports":{"three":"https://cdn.jsdelivr.net/npm/three@0.160.0/build/three.module.js","three/addons/":"https://cdn.jsdelivr.net/npm/three@0.160.0/examples/jsm/"}}</script>
</head><body>
${legendHtml(title)}
<script type="module">
import * as THREE from 'three';
import { OrbitControls } from 'three/addons/controls/OrbitControls.js';
${STICK_VIEWER_SCENE_JS}
const DATA = ${payload};
const camera = new THREE.PerspectiveCamera(55, innerWidth/innerHeight, 0.01, 500);
const renderer = new THREE.WebGLRenderer({antialias:true}); renderer.setSize(innerWidth, innerHeight); document.body.appendChild(renderer.domElement);
const controls = new OrbitControls(camera, renderer.domElement); controls.enableDamping = true;
const { scene, stickGroups, size } = buildScene(DATA);
wireLegend(stickGroups);
camera.position.set(size*0.6, size*0.5, size*0.8); controls.target.set(0,0,0); controls.update();
function anim(){requestAnimationFrame(anim); controls.update(); renderer.render(scene,camera);} anim();
addEventListener('resize',()=>{camera.aspect=innerWidth/innerHeight; camera.updateProjectionMatrix(); renderer.setSize(innerWidth,innerHeight);});
</script></body></html>`;
}

/// HTML viewer that loads *_topology.npz (file picker or HTTP fetch).
export function exportStickViewerNpzLoaderHTML(title = 'MMFFL topology', npzBasename = 'topology') {
    const safeTitle = String(title).replace(/</g, '&lt;');
    const npzFile = `${npzBasename}_topology.npz`;
    return `<!DOCTYPE html>
<html><head><meta charset="utf-8"><title>${safeTitle} (NPZ)</title>
<style>${STICK_VIEWER_CSS} #fileRow input{margin-left:6px}</style>
<script type="importmap">{"imports":{"three":"https://cdn.jsdelivr.net/npm/three@0.160.0/build/three.module.js","three/addons/":"https://cdn.jsdelivr.net/npm/three@0.160.0/examples/jsm/","jszip":"https://cdn.jsdelivr.net/npm/jszip@3.10.1/+esm"}}</script>
</head><body>
${legendHtml(safeTitle)}
<div id="fileRow"><label>NPZ <input type="file" id="npzFile" accept=".npz"></label></div>
<div id="status">Load ${npzFile} or pick file</div>
<script type="module">
import * as THREE from 'three';
import { OrbitControls } from 'three/addons/controls/OrbitControls.js';
import JSZip from 'jszip';
${STICK_VIEWER_SCENE_JS}

function parseNpy(buf) {
  const u8 = new Uint8Array(buf);
  if (String.fromCharCode(...u8.slice(0, 6)) !== '\\x93NUMPY') throw new Error('Not NPY');
  const verMajor = u8[6];
  let headerLen, headerOff, dataOff;
  if (verMajor === 1) { headerLen = u8[8] | (u8[9] << 8); headerOff = 10; dataOff = 10 + headerLen; }
  else { headerLen = new DataView(buf).getUint32(8, true); headerOff = 12; dataOff = 12 + headerLen; }
  const header = new TextDecoder().decode(u8.slice(headerOff, dataOff));
  const descrM = header.match(/'descr': '([^']+)'/);
  if (!descrM) throw new Error('NPY header missing descr');
  const descr = descrM[1];
  if (descr === '|O' || descr.startsWith('|S')) return null;
  const shapeM = header.match(/'shape': \\(([^)]*)\\)/);
  const shape = shapeM ? shapeM[1].split(',').map(s => parseInt(s.trim(), 10)).filter(n => !isNaN(n)) : [];
  const nElem = shape.length ? shape.reduce((a, b) => a * b, 1) : 1;
  const raw = u8.slice(dataOff);
  let data;
  if (descr === '<f4' || descr === '>f4') data = new Float32Array(raw.buffer, raw.byteOffset, nElem);
  else if (descr === '<f8' || descr === '>f8') data = new Float64Array(raw.buffer, raw.byteOffset, nElem);
  else if (descr === '<i4' || descr === '>i4') data = new Int32Array(raw.buffer, raw.byteOffset, nElem);
  else if (descr === '|u1') data = new Uint8Array(raw.buffer, raw.byteOffset, nElem);
  else return null;
  return { data, shape };
}

const NPZ_VIEWER_KEYS = new Set(['pos', 'Z', 'neigh_idx', 'stick_class', 'natoms', 'max_neighbors']);

async function loadNpz(buf) {
  const zip = await JSZip.loadAsync(buf);
  const out = {};
  for (const [name, entry] of Object.entries(zip.files)) {
    if (entry.dir || !name.endsWith('.npy')) continue;
    const key = name.replace(/\\.npy$/, '');
    if (!NPZ_VIEWER_KEYS.has(key)) continue;
    const parsed = parseNpy(await entry.async('arraybuffer'));
    if (!parsed) continue;
    out[key] = parsed;
  }
  const flat = (k) => out[k] ? Array.from(out[k].data) : null;
  const scalar = (k) => out[k] ? Number(out[k].data[0]) : 0;
  if (!out.pos || !out.Z || !out.neigh_idx || !out.stick_class) throw new Error('NPZ missing pos/Z/neigh_idx/stick_class');
  return {
    natoms: scalar('natoms') || Math.floor(flat('pos').length / 3),
    max_neighbors: scalar('max_neighbors') || Math.floor(flat('neigh_idx').length / (flat('pos').length / 3)),
    pos: flat('pos'), Z: flat('Z'), neigh_idx: flat('neigh_idx'), stick_class: flat('stick_class'),
  };
}

let camera, renderer, controls, sceneRoot;
function mountScene(DATA) {
  if (renderer) { renderer.domElement.remove(); renderer.dispose(); }
  camera = new THREE.PerspectiveCamera(55, innerWidth/innerHeight, 0.01, 500);
  renderer = new THREE.WebGLRenderer({ antialias: true });
  renderer.setSize(innerWidth, innerHeight);
  document.body.appendChild(renderer.domElement);
  controls = new OrbitControls(camera, renderer.domElement);
  controls.enableDamping = true;
  const { scene, stickGroups, size } = buildScene(DATA);
  sceneRoot = scene;
  wireLegend(stickGroups);
  camera.position.set(size * 0.6, size * 0.5, size * 0.8);
  controls.target.set(0, 0, 0);
  controls.update();
  function anim() { requestAnimationFrame(anim); controls.update(); renderer.render(scene, camera); }
  anim();
}

async function onNpzBuffer(buf, label) {
  try {
    const npz = await loadNpz(buf);
    const DATA = springsFromNpz(npz);
    DATA.title = label;
    mountScene(DATA);
    document.getElementById('status').textContent = label + ' — atoms=' + DATA.atoms.length + ' springs=' + DATA.springs.length;
  } catch (e) {
    document.getElementById('status').textContent = 'Error: ' + e.message;
    console.error(e);
  }
}

document.getElementById('npzFile').addEventListener('change', async (ev) => {
  const f = ev.target.files[0];
  if (!f) return;
  onNpzBuffer(await f.arrayBuffer(), f.name);
});

const defaultNpz = '${npzFile}';
if (location.protocol === 'file:') {
  document.getElementById('status').textContent = 'Use file picker to load ' + defaultNpz + ' (file:// cannot auto-fetch)';
} else if (location.protocol.startsWith('http')) {
  fetch(defaultNpz).then(r => r.ok ? r.arrayBuffer() : Promise.reject()).then(b => onNpzBuffer(b, defaultNpz)).catch(() => {
    document.getElementById('status').textContent = 'Pick ' + defaultNpz + ' or place it beside this HTML';
  });
}
addEventListener('resize', () => {
  if (!camera || !renderer) return;
  camera.aspect = innerWidth / innerHeight;
  camera.updateProjectionMatrix();
  renderer.setSize(innerWidth, innerHeight);
});
</script></body></html>`;
}
