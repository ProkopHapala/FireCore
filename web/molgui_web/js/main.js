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
import { PDSimulation } from './ProjectiveDynamics.js';

export class MolGUIApp {
    constructor() {
        this.container = document.getElementById('canvas-container');
        if (!this.container) throw new Error("Could not find #canvas-container");
        this.scene = null;
        this.camera = null;
        this.renderer = null;
        this.system = new EditableMolecule();
        this.controls = null;
        this._rafPending = false;
        this._rafLoopActive = false;
        this.continuousRender = false;
        this.pdDebugTexture = null;
    }

    requestRender() {
        if (this._rafPending) return;
        this._rafPending = true;
        requestAnimationFrame(() => {
            this._rafPending = false;
            if (this.renderers) {
                Object.values(this.renderers).forEach(r => r.update());
            }
            if (this.renderer && this.scene && this.camera) {
                this.renderer.render(this.scene, this.camera);
                if (this.pdSimulation && this.pdSimulation.debugTextureName) {
                    this.pdSimulation.renderDebugOverlay(this.renderer);
                }
            }
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
        this.camera.position.set(0, 0, 10);
        this.camera.lookAt(0, 0, 0);

        // 3. Renderer
        const canvas = document.createElement('canvas');
        const gl = canvas.getContext('webgl2', {
            antialias: true,
            powerPreference: "high-performance",
            preserveDrawingBuffer: false,
        });
        if (!gl) {
            throw new Error('WebGL2 context not available. Projective Dynamics requires WebGL2/GLSL 300 ES.');
        }
        this.renderer = new THREE.WebGLRenderer({
            canvas,
            context: gl
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
        if (!this.renderer.capabilities.isWebGL2) {
            throw new Error('Three.js fell back to WebGL1. Projective Dynamics shaders need WebGL2.');
        }

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
        await this.mmParams.loadResources('../common_resources/ElementTypes.dat', '../common_resources/AtomTypes.dat', '../common_resources/BondTypes.dat', '../common_resources/AngleTypes.dat');

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
        this.requestRender = this.requestRender.bind(this);

        // --- Replica / Lattice State ---
        this.lattices = new Map();
        this.getLattice = (name = 'default') => {
            const key = name || 'default';
            if (!this.lattices.has(key)) {
                this.lattices.set(key, {
                    lvec: [new Vec3(10, 0, 0), new Vec3(0, 10, 0), new Vec3(0, 0, 10)],
                    nrep: { x: 0, y: 0, z: 0 },
                    show: false,
                    showBox: false,
                    filter: null
                });
            }
            return this.lattices.get(key);
        };
        this.lattice = this.getLattice('default');

        this.updateReplicas = (name = 'default') => {
            const lat = this.getLattice(name);
            if (!lat) return;
            const renderer = this.renderers[name] || this.renderers.molecule || this.molRenderer;
            const logPrefix = `[updateReplicas:${name}]`;
            const logMsg = `${logPrefix} show=${lat.show} showBox=${lat.showBox} nrep=(${lat.nrep.x},${lat.nrep.y},${lat.nrep.z}) lvec=${lat.lvec.map(v => v ? `(${v.x.toFixed(2)},${v.y.toFixed(2)},${v.z.toFixed(2)})` : '(null)').join(' ')}`;
            if (window.logger && window.logger.info) window.logger.info(logMsg);
            else console.log(logMsg);
            if (renderer) {
                renderer.setReplicas({
                    lvec: lat.lvec,
                    nrep: lat.nrep,
                    show: lat.show,
                    showBox: lat.showBox
                });
                renderer.update();
            }
            this.requestRender();
            if (this.gui && typeof this.gui.refreshLatticeControls === 'function') {
                this.gui.refreshLatticeControls(name);
            }
        };

        // Alias for GUI checkbox; MoleculeRenderer handles box via showBox flag
        this.updateLatticeBox = (name = 'default') => {
            this.updateReplicas(name);
        };

        this.bakeReplicas = (name = 'default') => {
            const lat = this.getLattice(name);
            const { x: nx, y: ny, z: nz } = lat.nrep;
            let targetSystem = this.system;
            if (name === 'substrate') targetSystem = this.systems.substrate;
            else if (name === 'molecule') targetSystem = this.systems.molecule;

            targetSystem.replicate([nx * 2 + 1, ny * 2 + 1, nz * 2 + 1], lat.lvec);

            const renderer = name === 'substrate' ? this.renderers.substrate : (name === 'molecule' ? this.renderers.molecule : this.molRenderer);
            if (renderer) {
                lat.nrep = { x: 0, y: 0, z: 0 };
                lat.show = false;
                lat.showBox = false;
                renderer.setReplicas({ nrep: lat.nrep, show: false, showBox: false });
                renderer.update();
            }
            if (this.gui) this.gui.updateSelectionUI();
            this.requestRender();
            window.logger.info(`Baked replicas for '${name}': system now has ${targetSystem.atoms.length} atoms.`);
        };

        // --- Bucket overlay (debug visualization) ---
        this.refreshBucketDebug = () => {
            const bg = this.lastBucketGraph;
            if (!bg || !this.molRenderer.bucketRenderer) return;
            if (typeof bg.toInds === 'function') bg.toInds(this.system);
            if (typeof bg.pruneEmptyBuckets === 'function') bg.pruneEmptyBuckets();
            if (typeof bg.recalcBounds === 'function') bg.recalcBounds(this.system);

            this.molRenderer.bucketRenderer.updateBuckets(
                bg,
                this.system,
                !!this.showBucketBoxes,
                !!this.showBucketAtomLines
            );
            this.requestRender();
        };

        this.updateBucketOverlay = () => {
            this.refreshBucketDebug();
        };

        // 7. GUI
        this.gui = new GUI(this.system, this.molRenderer);

        // 7. Editor (Selection, Gizmo)
        this.editor = new Editor(this.scene, this.camera, this.renderer, this.system, this.molRenderer);

        // defaults / shared state
        this.bondRecalcMode = this.bondRecalcMode ? String(this.bondRecalcMode) : 'brute';
        this.showBucketBoxes = !!this.showBucketBoxes;

        // Default: show axes (GUI default is checked)
        Object.values(this.renderers).forEach(r => {
            if (r && r.toggleAxes) r.toggleAxes(true);
        });

        // 7. Script Runner
        this.scriptRunner = new ScriptRunner(this);

        // 8. Projective Dynamics Simulation (optional)
        this.pdSimulation = null;
        this.pdEnabled = false;
        this.pdDt = 0.016;

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

    async enablePD(force = false) {
        const nAtoms = this.system.atoms.length || 0;
        if (nAtoms === 0) {
            window.logger.warn("Cannot enable PD: no atoms in system");
            return;
        }
        if (this.pdSimulation && (force || this.pdSimulation.N !== nAtoms)) {
            this.pdSimulation.resize(nAtoms);
            this.pdEnabled = true;
            return;
        }
        if (!force && this.pdEnabled && this.pdSimulation && this.pdSimulation.N === nAtoms) return;
        this.pdSimulation = new PDSimulation(this.renderer, nAtoms, 16);
        await this.pdSimulation.init();
        this.pdEnabled = true;
        window.logger.info("Projective Dynamics simulation ready");
    }

    async runRelax({ dt = 0.1, iterations = 40, rebuild = true } = {}) {
        const nAtoms = this.system.atoms.length || 0;
        if (nAtoms === 0) {
            window.logger.warn("runRelax: no atoms to relax");
            return;
        }
        await this.enablePD();
        if (this.pdSimulation && this.pdSimulation.N !== nAtoms) {
            this.pdSimulation.resize(nAtoms);
        }
        if (!this.pdSimulation) {
            window.logger.error("runRelax: PD simulation not initialized");
            return;
        }

        // Ensure packed buffers reflect current geometry
        this.molRenderer.update();

        if (rebuild) {
            this.system.updateNeighborList();
            this.system.printSizes('beforePD');
            const bondCount = this.system.bonds ? this.system.bonds.length : 0;
            console.log(`[PD] setTopology rebuild: atoms=${nAtoms} bonds=${bondCount}`);
            this.pdSimulation.setTopology(this.system, this.mmParams);
        }
        this.pdSimulation.setPositions(this.packedSystems.molecule.pos);
        this.pdSimulation.setVelocities(null);
        this.pdSimulation.uploadInitialData();

        console.log('--- PD BEFORE STEP ---');
        this.pdSimulation.debugReadTexture('Initial Positions', this.pdSimulation.texPosInitial, 8);
        this.pdSimulation.debugReadTexture('Initial Velocities', this.pdSimulation.texVelInitial, 8);
        this.pdSimulation.debugReadTexture('texBonds', this.pdSimulation.texBonds, {
            count: Math.min(8, nAtoms),
            mode: 'bonds',
            width: this.pdSimulation.maxBonds,
            height: this.pdSimulation.N
        });

        const evenIterations = (iterations % 2 === 0) ? iterations : iterations + 1;
        const prevIter = this.pdSimulation.iterationCount;
        this.pdSimulation.iterationCount = evenIterations;
        this.pdSimulation.step(dt);
        this.pdSimulation.iterationCount = prevIter;

        console.log('--- PD AFTER STEP ---');
        this.pdSimulation.debugReadTexture('Final Positions', this.pdSimulation.getOutputTexture(), 8);
        this.pdSimulation.debugReadTexture('Velocities', this.pdSimulation.fboVel, 8);
        this.pdSimulation.renderDebugView();

        const buffer = this.pdSimulation.readPositions();
        if (!buffer) {
            window.logger.warn("runRelax: failed to read positions from PD");
            return;
        }

        const packedPos = this.packedSystems.molecule.pos;
        for (let i = 0; i < nAtoms; i++) {
            const a = this.system.atoms[i];
            if (!a) continue;
            const x = buffer[i * 4];
            const y = buffer[i * 4 + 1];
            const z = buffer[i * 4 + 2];
            a.pos.set(x, y, z);
            const base = i * 3;
            packedPos[base] = x;
            packedPos[base + 1] = y;
            packedPos[base + 2] = z;
        }
        this.system.dirtyGeom = true;
        this.system.dirtyExport = true;
        this.molRenderer.update();
        this.requestRender();
        window.logger.info(`runRelax: completed ${evenIterations} iterations`);
    }
}

// Start
window.onload = () => {
    const app = new MolGUIApp();
    window.app = app; // Expose for debugging
    app.init();
};
