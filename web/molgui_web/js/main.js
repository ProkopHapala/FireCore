import { logger } from '../../common_js/Logger.js';
import { EditableMolecule } from './EditableMolecule.js';
import { MMParams } from './MMParams.js';
import { MoleculeRenderer, PackedMolecule } from './MoleculeRenderer.js';
import { GUI } from './GUI.js';
import { Editor } from './Editor.js';
import { ShortcutManager } from './ShortcutManager.js';
import { Vec3 } from '../../common_js/Vec3.js';
import { buildWireframeCellVerts, buildWireframeAABBVerts } from '../../common_js/Buckets.js';
import { ScriptRunner } from './ScriptRunner.js';

export class MolGUIApp {
    constructor() {
        this.container = document.getElementById('canvas-container');
        this.scene = null;
        this.camera = null;
        this.renderer = null;
        this.controls = null;
        this._rafPending = false;
        this._rafLoopActive = false;
        this.continuousRender = false;
    }

    requestRender() {
        if (this._rafPending) return;
        this._rafPending = true;
        requestAnimationFrame(() => {
            this._rafPending = false;
            if (this.renderers) {
                Object.values(this.renderers).forEach(r => r.update());
            }
            if (this.renderer && this.scene && this.camera) this.renderer.render(this.scene, this.camera);
        });
    }

    setContinuousRender(enabled) {
        enabled = !!enabled;
        if (this.continuousRender === enabled) return;
        this.continuousRender = enabled;
        if (enabled) {
            this._startLoop();
        } else {
            this._stopLoop();
            this.requestRender();
        }
    }

    _startLoop() {
        if (this._rafLoopActive) return;
        this._rafLoopActive = true;
        this.animate();
    }

    _stopLoop() {
        this._rafLoopActive = false;
    }

    async init() {
        logger.info("Initializing MolGUI...");

        // Load Shaders first
        try {
            const vPromise  = fetch('../common_resources/shaders/atom.glslv').then(r => r.text());
            const fPromise  = fetch('../common_resources/shaders/atom.glslf').then(r => r.text());
            const bvPromise = fetch('../common_resources/shaders/bond.glslv').then(r => r.text());
            const svPromise = fetch('../common_resources/shaders/selection.glslv').then(r => r.text());
            const bfPromise = fetch('../common_resources/shaders/bond_color.glslf').then(r => r.text());
            const cfPromise = fetch('../common_resources/shaders/color.glslf').then(r => r.text());
            const lvPromise = fetch('../common_resources/shaders/label.glslv').then(r => r.text());
            const lfPromise = fetch('../common_resources/shaders/label.glslf').then(r => r.text());

            const [vertex, fragment, bVertex, sVertex, bondFrag, colorFrag, lVertex, lFragment] = await Promise.all([
                vPromise, fPromise, bvPromise, svPromise, bfPromise, cfPromise, lvPromise, lfPromise
            ]);

            this.shaders = {
                atom:      { vertex,      fragment },
                bond:      { vertex: bVertex, fragment: bondFrag },   // bond_color.glslf (uses vColor)
                selection: { vertex: sVertex, fragment: colorFrag },  // color.glslf (uses uColor)
                label:     { vertex: lVertex, fragment: lFragment }
            };
            window.logger.info("Shaders loaded.");
        } catch (e) {
            window.logger.error("Failed to load shaders: " + e);
            return;
        }

        // 1. Scene
        this.scene = new THREE.Scene();
        this.scene.background = new THREE.Color(0x222222);

        // 2. Camera
        const width = this.container.clientWidth;
        const height = this.container.clientHeight;
        const aspect = width / height;
        this.frustumSize = 20;
        this.camera = new THREE.OrthographicCamera(
            this.frustumSize * aspect / -2,
            this.frustumSize * aspect / 2,
            this.frustumSize / 2,
            this.frustumSize / -2,
            -1000.,
             1000.
        );
        this.camera.position.set(10, 10, 10);
        this.camera.lookAt(0, 0, 0);

        // 3. Renderer
        this.renderer = new THREE.WebGLRenderer({
            antialias: true,
            powerPreference: "high-performance"
        });
        this.renderer.setSize(width, height);
        this.renderer.setPixelRatio(window.devicePixelRatio);
        // Make canvas focusable so keyboard shortcuts work when user clicks viewport
        this.renderer.domElement.tabIndex = 0;
        this.renderer.domElement.style.outline = 'none';
        this.renderer.domElement.addEventListener('pointerdown', () => {
            this.renderer.domElement.focus();
        });
        this.container.appendChild(this.renderer.domElement);

        // 4. Controls
        if (THREE.OrbitControls) {
            this.controls = new THREE.OrbitControls(this.camera, this.renderer.domElement);
            this.controls.enableDamping = false; // Disable inertia for raw performance feel
            this.controls.dampingFactor = 0.05;

            // Render only when something changes
            this.controls.addEventListener('change', () => this.requestRender());

            // Unity Style Controls:
            // LMB: Selection (Handled by Editor) -> Set to NULL here
            // MMB: Zoom/Dolly (Standard)
            // RMB: Rotate (Standard)
            // Shift + RMB: Pan (Custom)

            this.controls.mouseButtons = {
                LEFT: null, // Disable Left Click for Camera
                MIDDLE: THREE.MOUSE.DOLLY,
                RIGHT: THREE.MOUSE.ROTATE
            };

            // Listen for Shift Key to toggle Pan Mode on RMB
            window.addEventListener('keydown', (e) => {
                if (e.key === 'Shift') {
                    this.controls.mouseButtons.RIGHT = THREE.MOUSE.PAN;
                }
            });

            window.addEventListener('keyup', (e) => {
                if (e.key === 'Shift') {
                    this.controls.mouseButtons.RIGHT = THREE.MOUSE.ROTATE;
                }
            });
        } else {
            window.logger.error("OrbitControls not loaded!");
        }

        // 5. Molecule Systems (authoritative editable models)
        this.systems = {
            molecule:  new EditableMolecule(),
            substrate: new EditableMolecule()
        };
        this.system = this.systems.molecule; // default active system
        
        this.packedSystems = {
            molecule:  new PackedMolecule(),
            substrate: new PackedMolecule()
        };

        // 5. MMParams (Load resources)
        this.mmParams = new MMParams();
        await this.mmParams.loadResources('../common_resources/ElementTypes.dat', '../common_resources/AtomTypes.dat', '../common_resources/BondTypes.dat');

        // 6. Molecule Renderers
        this.renderers = {
            molecule:  new MoleculeRenderer(this.scene, this.packedSystems.molecule,  this.shaders, this.mmParams, this.systems.molecule),
            substrate: new MoleculeRenderer(this.scene, this.packedSystems.substrate, this.shaders, this.mmParams, this.systems.substrate)
        };
        this.molRenderer = this.renderers.molecule; // default active renderer

        // Default: show axes (GUI default is checked)
        Object.values(this.renderers).forEach(r => {
            if (r && r.toggleAxes) r.toggleAxes(true);
        });

        // 7. Script Runner
        this.scriptRunner = new ScriptRunner(this);

        // --- Replica / Lattice State ---
        this.lattices = new Map();
        this.getLattice = (name = 'default') => {
            if (!this.lattices.has(name)) {
                const lat = {
                    lvec: [new Vec3(10, 0, 0), new Vec3(0, 10, 0), new Vec3(0, 0, 10)],
                    nrep: { x: 1, y: 1, z: 1 },
                    show: false,
                    showBox: false,
                    group: new THREE.Group(),
                    box: null,
                    filter: null // optional function (atom) => boolean
                };
                this.scene.add(lat.group);
                this.lattices.set(name, lat);
            }
            return this.lattices.get(name);
        };
        this.lattice = this.getLattice('default'); // Default/Legacy lattice

        this.updateLatticeBox = (name = 'default') => {
            const lat = this.getLattice(name);
            if (lat.box) { lat.group.remove(lat.box); lat.box = null; }
            if (!lat.showBox) return;
            const [a, b, c] = lat.lvec;
            const { x: nx, y: ny, z: nz } = lat.nrep;
            const superA = a.clone().mulScalar(nx * 2 + 1);
            const superB = b.clone().mulScalar(ny * 2 + 1);
            const superC = c.clone().mulScalar(nz * 2 + 1);
            const origin = a.clone().mulScalar(-nx).addMul(b, -ny).addMul(c, -nz);
            const verts = buildWireframeCellVerts(superA, superB, superC, origin);
            const geom = new THREE.BufferGeometry();
            geom.setAttribute('position', new THREE.BufferAttribute(verts, 3));
            const mat = new THREE.LineBasicMaterial({ color: (name === 'substrate' ? 0xffff00 : 0x00ffff), opacity: 0.5, transparent: true });
            lat.box = new THREE.LineSegments(geom, mat);
            lat.group.add(lat.box);
        };

        this.updateReplicas = (name = 'default') => {
            const lat = this.getLattice(name);
            lat.group.clear();
            this.updateLatticeBox(name);
            if (!lat.show) { this.requestRender(); return; }
            const { x: nx, y: ny, z: nz } = lat.nrep;
            const [a, b, c] = lat.lvec;
            
            // Pick correct renderer based on lattice name
            let targetRenderer = this.molRenderer;
            if (name === 'substrate') targetRenderer = this.renderers.substrate;
            else if (name === 'molecule') targetRenderer = this.renderers.molecule;

            const baseMeshes = [
                targetRenderer.atomMesh, 
                targetRenderer.bondLines, 
                targetRenderer.labelMesh
            ].filter(m => m !== null);

            for (let ix = -nx; ix <= nx; ix++) {
                for (let iy = -ny; iy <= ny; iy++) {
                    for (let iz = -nz; iz <= nz; iz++) {
                        if (ix === 0 && iy === 0 && iz === 0) continue;
                        const shift = new Vec3();
                        shift.addMul(a, ix);
                        shift.addMul(b, iy);
                        shift.addMul(c, iz);
                        const rep = new THREE.Group();
                        rep.position.set(shift.x, shift.y, shift.z);
                        baseMeshes.forEach(mesh => {
                            let clone;
                            if (mesh instanceof THREE.LineSegments) {
                                clone = new THREE.LineSegments(mesh.geometry, mesh.material);
                            } else if (mesh instanceof THREE.InstancedMesh) {
                                clone = new THREE.InstancedMesh(mesh.geometry, mesh.material, mesh.count);
                                for (const key in mesh.geometry.attributes) {
                                    clone.geometry.setAttribute(key, mesh.geometry.attributes[key]);
                                }
                                clone.count = mesh.count;
                            } else {
                                clone = new THREE.Mesh(mesh.geometry, mesh.material);
                            }
                            clone.raycast = () => { }; // non-pickable
                            rep.add(clone);
                        });
                        lat.group.add(rep);
                    }
                }
            }
            this.requestRender();
        };

        this.bakeReplicas = (name = 'default') => {
            const lat = this.getLattice(name);
            const { x: nx, y: ny, z: nz } = lat.nrep;
            
            let targetSystem = this.system;
            if (name === 'substrate') targetSystem = this.systems.substrate;
            else if (name === 'molecule') targetSystem = this.systems.molecule;

            targetSystem.replicate([nx * 2 + 1, ny * 2 + 1, nz * 2 + 1], lat.lvec);
            
            // Trigger UI and Renderer update
            if (this.gui) this.gui.updateSelectionUI();
            this.requestRender();
            window.logger.info(`Baked replicas for '${name}': system now has ${targetSystem.atoms.length} atoms.`);
        };

        // --- Bucket overlay (debug visualization) ---
        this.bucketOverlay = null;
        {
            const geom = new THREE.BufferGeometry();
            geom.setAttribute('position', new THREE.BufferAttribute(new Float32Array(0), 3));
            const mat = new THREE.LineBasicMaterial({ color: 0xff8844, transparent: true, opacity: 0.6 });
            this.bucketOverlay = new THREE.LineSegments(geom, mat);
            this.bucketOverlay.renderOrder = 12;
            this.bucketOverlay.visible = false;
            this.scene.add(this.bucketOverlay);
        }

        // --- Bucket atom->center lines (debug visualization) ---
        this.bucketAtomLines = null;
        {
            const geom = new THREE.BufferGeometry();
            geom.setAttribute('position', new THREE.BufferAttribute(new Float32Array(0), 3));
            const mat = new THREE.LineBasicMaterial({ color: 0x55ccff, transparent: true, opacity: 0.6 });
            this.bucketAtomLines = new THREE.LineSegments(geom, mat);
            this.bucketAtomLines.renderOrder = 13;
            this.bucketAtomLines.visible = false;
            this.scene.add(this.bucketAtomLines);
        }

        // 7. GUI
        this.gui = new GUI(this.system, this.molRenderer);

        // 7. Editor (Selection, Gizmo)
        this.editor = new Editor(this.scene, this.camera, this.renderer, this.system, this.molRenderer);

        // defaults / shared state
        this.bondRecalcMode = this.bondRecalcMode ? String(this.bondRecalcMode) : 'brute';
        this.showBucketBoxes = !!this.showBucketBoxes;
        this.autoUpdateBuckets = (this.autoUpdateBuckets !== undefined) ? !!this.autoUpdateBuckets : true;
        this.showBucketAtomLines = !!this.showBucketAtomLines;

        this.refreshBucketDebug = () => {
            const bg = this.lastBucketGraph;
            if (!bg) return;
            if (typeof bg.toInds === 'function') bg.toInds(this.system);
            if (typeof bg.pruneEmptyBuckets === 'function') bg.pruneEmptyBuckets();
            if (typeof bg.recalcBounds === 'function') bg.recalcBounds(this.system);
            if (typeof this.updateBucketOverlay === 'function') this.updateBucketOverlay();
            if (!this.bucketAtomLines) return;
            const showLines = !!this.showBucketAtomLines;
            if (!showLines || !bg.buckets || bg.buckets.length === 0) {
                this.bucketAtomLines.visible = false;
                this.bucketAtomLines.geometry.setAttribute('position', new THREE.BufferAttribute(new Float32Array(0), 3));
                this.bucketAtomLines.geometry.computeBoundingSphere();
                this.requestRender();
                return;
            }
            let nPairs = 0;
            for (let ib = 0; ib < bg.buckets.length; ib++) nPairs += (bg.buckets[ib].atoms.length | 0);
            const verts = new Float32Array(nPairs * 2 * 3);
            let k = 0;
            const c = new Vec3();
            for (let ib = 0; ib < bg.buckets.length; ib++) {
                bg.getBucketCenterFromBounds(ib, c);
                const bs = bg.buckets[ib].atoms;
                for (let i = 0; i < bs.length; i++) {
                    const ia = bs[i] | 0;
                    const a = this.system.atoms[ia];
                    if (!a) continue;
                    const p = a.pos;
                    verts[k++] = p.x; verts[k++] = p.y; verts[k++] = p.z;
                    verts[k++] = c.x; verts[k++] = c.y; verts[k++] = c.z;
                }
            }
            const v2 = (k === (verts.length | 0)) ? verts : verts.subarray(0, k);
            this.bucketAtomLines.visible = true;
            this.bucketAtomLines.geometry.setAttribute('position', new THREE.BufferAttribute(v2, 3));
            this.bucketAtomLines.geometry.computeBoundingSphere();
            this.requestRender();
        };

        this.updateBucketOverlay = () => {
            if (!this.bucketOverlay) return;
            const show = !!this.showBucketBoxes;
            const bg = this.lastBucketGraph;
            if (!show || !bg || !bg.buckets || bg.buckets.length === 0) {
                this.bucketOverlay.visible = false;
                this.bucketOverlay.geometry.setAttribute('position', new THREE.BufferAttribute(new Float32Array(0), 3));
                this.bucketOverlay.geometry.computeBoundingSphere();
                this.requestRender();
                return;
            }
            const bs = bg.buckets;
            let verts = new Float32Array(0);
            let n3 = 0;
            if (bg.meta && bg.meta.kind === 'crystal_cells' && bg.meta.lvec && bg.meta.origin) {
                const A = bg.meta.lvec[0], B = bg.meta.lvec[1], C = bg.meta.lvec[2];
                const o0 = bg.meta.origin;
                const per = 24 * 3;
                verts = new Float32Array(bs.length * per);
                const O = new Vec3();
                for (let ib = 0; ib < bs.length; ib++) {
                    const m = bs[ib].meta;
                    O.setV(o0);
                    if (m) {
                        if (m.ix) O.addMul(A, m.ix);
                        if (m.iy) O.addMul(B, m.iy);
                        if (m.iz) O.addMul(C, m.iz);
                    }
                    buildWireframeCellVerts(A, B, C, O, verts, n3);
                    n3 += per;
                }
            } else {
                const per = 12 * 2 * 3;
                verts = new Float32Array(bs.length * per);
                for (let ib = 0; ib < bs.length; ib++) {
                    const b = bs[ib];
                    buildWireframeAABBVerts(b.pmin, b.pmax, verts, n3);
                    n3 += per;
                }
            }
            const v2 = (n3 === (verts.length | 0)) ? verts : verts.subarray(0, n3);
            this.bucketOverlay.visible = true;
            this.bucketOverlay.geometry.setAttribute('position', new THREE.BufferAttribute(v2, 3));
            this.bucketOverlay.geometry.computeBoundingSphere();
            this.requestRender();
        };

        // 8. Selection Rendering (Centralized in MoleculeRenderer)
        // No extra code needed here, MoleculeRenderer handles it.

        this.editor.onSelectionChange = () => {
            this.gui.updateSelectionUI();
            this.molRenderer.updateSelection();
            this.requestRender();
        };

        // 9. Shortcut Manager
        this.shortcuts = new ShortcutManager(this.editor);

        // Hook into Editor or System updates?
        // Ideally System should dispatch events, but we can just patch the update method or call it manually.
        // For now, let's monkey-patch the renderer update or just call it in animate loop?
        // Better: Pass it to Editor? Or make Editor call a global update?
        // Let's make Editor call it.
        // this.editor.onSelectionChange is already set above.

        // Also hook GUI input back to renderer
        // Also hook GUI input back to renderer
        this.gui.onSelectionChanged = () => {
            this.molRenderer.updateSelection();
        };

        // No default molecules added; scripts/builders should populate systems as needed.

        this.molRenderer.update();

        // Lights
        const ambientLight = new THREE.AmbientLight(0xffffff, 0.5);
        this.scene.add(ambientLight);
        const dirLight = new THREE.DirectionalLight(0xffffff, 1);
        dirLight.position.set(5, 10, 7);
        this.scene.add(dirLight);

        // Events
        window.addEventListener('resize', this.onWindowResize.bind(this));

        // Initial render
        this.requestRender();

        window.logger.info("Initialization Complete.");
    }

    onWindowResize() {
        const width = this.container.clientWidth;
        const height = this.container.clientHeight;
        const aspect = width / height;

        this.camera.left = -this.frustumSize * aspect / 2;
        this.camera.right = this.frustumSize * aspect / 2;
        this.camera.top = this.frustumSize / 2;
        this.camera.bottom = -this.frustumSize / 2;

        this.camera.updateProjectionMatrix();
        this.renderer.setSize(width, height);
        this.requestRender();
    }

    animate() {
        if (!this._rafLoopActive) return;
        requestAnimationFrame(this.animate.bind(this));
        if (!this.continuousRender) return;
        if (this.controls) this.controls.update();
        if (this.renderers) {
            Object.values(this.renderers).forEach(r => r.update());
        }
        if (this.renderer && this.scene && this.camera) this.renderer.render(this.scene, this.camera);
    }
}

// Start
window.onload = () => {
    const app = new MolGUIApp();
    window.app = app; // Expose for debugging
    app.init();
};
