import * as THREE from 'three';
import { WebGPURenderer } from 'three/webgpu';
import { OrbitControls } from 'three/addons/controls/OrbitControls.js';
import { TransformControls } from 'three/addons/controls/TransformControls.js';

import { logger } from '../common_js/Logger.js';
import { EditableMolecule } from './EditableMolecule.js';
import { installMoleculeIOMethods } from './MoleculeIO.js';
import { MMParams } from './MMParams.js';
import { MoleculeRenderer, PackedMolecule } from './MoleculeRenderer.js';
import { RawWebGPUAtomsRenderer } from './RawWebGPUAtomsRenderer.js';
import { GUI } from './GUI.js';
import { Editor } from './Editor.js';
import { ShortcutManager } from './ShortcutManager.js';
import { Vec3 } from '../common_js/Vec3.js';
import { ScriptRunner } from './ScriptRunner.js';
import { XPDB_WebGPU } from './XPDB_WebGPU.js';
import { buildXPDBTopology, getMaxRadius } from './XPDBTopology.js';

installMoleculeIOMethods(EditableMolecule);

export class MolGUIApp {
    constructor() {
        this.container = document.getElementById('canvas-container');
        if (!this.container) throw new Error("Could not find #canvas-container");
        this.scene = null;
        this.camera = null;
        this.renderer = null;
        this.raw = null;
        this.overlayRenderer = null;
        this.overlayCanvas = null;
        this.overlayScene = null;
        this.useRawWebGPU = !!window.USE_RAW_WEBGPU;
        this._rawLastLabelMode = 'none';
        this.system = new EditableMolecule();
        this.controls = null;
        this._rafPending = false;
        this.continuousRender = false;

        // XPDB WebGPU simulation state
        this.xpdb = null;
        this.xpdbEnabled = false;
        this.xpdbUseAngles = true;
        this.xpdbPaused = true;
        this.xpdbDisableBonds = false;
        this.xpdbDisableCollisions = false;
        this.xpdbParams = {
            dt: 0.01,
            iterations: 4,
            k_coll: 500.0,
            omega: 1.0,
            momentum_beta: 0.0
        };
        this.xpdbDirty = true; // Flag to re-upload topology
    }

    requestRender() {
        if (this._rafPending) return;
        this._rafPending = true;
        requestAnimationFrame(async () => {
            this._rafPending = false;

            // Step XPDB simulation if enabled
            if (this.xpdbEnabled) {
                await this.stepXPDB();
            }

            if (this.renderers) {
                Object.values(this.renderers).forEach(r => r.update());
            }

            if (this.useRawWebGPU && this.raw) {
                if (this.controls && this.controls.update) this.controls.update();
                this.camera.updateMatrixWorld(true);
                this.camera.matrixWorldInverse.copy(this.camera.matrixWorld).invert();
                this.camera.updateProjectionMatrix();

                const view = new Float32Array(this.camera.matrixWorldInverse.elements);
                const proj = new Float32Array(this.camera.projectionMatrix.elements);
                this.raw.updateCamera(view, proj);

                const sys = this.packedSystems ? this.packedSystems.molecule : null;
                if (sys && sys.nAtoms > 0) {
                    const n = sys.nAtoms | 0;
                    const pos = sys.pos;
                    const col = new Float32Array(n * 3);
                    const rad = new Float32Array(n);
                    const sel = this.systems.molecule.selection;
                    for (let i = 0; i < n; i++) {
                        const t = sys.types[i] | 0;
                        const c = this.mmParams ? this.mmParams.getColor(t) : [1, 0, 1];
                        col[i * 3] = c[0];
                        col[i * 3 + 1] = c[1];
                        col[i * 3 + 2] = c[2];
                        const r = this.mmParams ? (this.mmParams.getRadius(t) * 0.4) : 0.5;
                        rad[i] = r;
                    }
                    this.raw.updateAtoms(pos, col, rad, n);

                    // --- Bonds ---
                    if (sys.bonds && sys.bonds.length) {
                        this.raw.updateBonds(sys.bonds);
                    } else {
                        this.raw.updateBonds([]);
                    }

                    // --- Selection halo ---
                    const selIdx = sel ? Array.from(sel) : [];
                    this.raw.updateSelection(selIdx, n, 1.3);

                    // --- Labels ---
                    const mode = (this.molRenderer && this.molRenderer.labelMode) ? this.molRenderer.labelMode : 'none';
                    const wantLabels = (mode !== 'none');
                    console.log('[MolGUIApp/raw] labels', { mode, wantLabels, nAtoms: n });
                    this.raw.setLabelsVisible(wantLabels);
                    if (wantLabels) {
                        const strings = new Array(n);
                        for (let i = 0; i < n; i++) {
                            if (mode === 'id') {
                                strings[i] = i.toString();
                            } else if (mode === 'element') {
                                const type = sys.types[i] | 0;
                                if (this.mmParams && this.mmParams.byAtomicNumber && this.mmParams.byAtomicNumber[type]) {
                                    strings[i] = this.mmParams.byAtomicNumber[type].name;
                                } else {
                                    strings[i] = type.toString();
                                }
                            } else if (mode === 'type') {
                                strings[i] = '?';
                            } else {
                                strings[i] = '';
                            }
                        }
                        console.log('[MolGUIApp/raw] updateLabels sample', { mode, nAtoms: n, s0: strings[0], s1: strings[1], s2: strings[2], s3: strings[3] });
                        this.raw.updateLabels(strings, n);
                    } else {
                        this.raw.updateLabels(null, 0);
                    }
                } else {
                    this.raw.nAtoms = 0;
                    this.raw.updateBonds([]);
                    this.raw.updateLabels(null, 0);
                }

                this.raw.render();

                // Pure WebGPU mode: no WebGL overlay rendering.
                return;
            }

            if (this.renderer && this.scene && this.camera) {
                if (this.renderer.renderAsync) await this.renderer.renderAsync(this.scene, this.camera);
                else this.renderer.render(this.scene, this.camera);
            }
        });
    }

    setContinuousRender(enabled) {
        enabled = !!enabled;
        if (this.continuousRender === enabled) return;
        this.continuousRender = enabled;
        if (enabled) {
            this.renderer.setAnimationLoop(() => {
                this.requestRender();
            });
        } else {
            this.renderer.setAnimationLoop(null);
        }
    }

    async ensureXPDBReady(needTopology = false) {
        if ((this.system?.nAtoms | 0) <= 0) {
            console.warn('[ensureXPDBReady] No atoms loaded');
            return false;
        }
        if (!this.xpdb) {
            await this.initXPDB();
        }
        if (!this.xpdb) return false;
        if (needTopology || this.xpdbDirty) {
            this.xpdbEnabled = true;
            await this.updateXPDBTopology();
        }
        return true;
    }

    async init() {
        logger.info("Initializing MolGUI (WebGPU/TSL)...");

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

        // 3. WebGPU Hardware Diagnostic
        logger.info("--- WebGPU Hardware Diagnostic ---");
        if (!navigator.gpu) {
            logger.error("navigator.gpu is UNDEFINED. Are you in a secure context (HTTPS/localhost)?");
        } else {
            try {
                const adapter = await navigator.gpu.requestAdapter();
                if (adapter) {
                    logger.info(`Adapter: ${adapter.name || 'unknown'}`);
                    const device = await adapter.requestDevice();
                    if (device) {
                        logger.info("Device: Successfully created.");
                        // Raw compute test
                        const shaderCode = `@group(0) @binding(0) var<storage, read_write> data: array<f32>; @compute @workgroup_size(1) fn main(@builtin(global_invocation_id) id: vec3<u32>) { data[id.x] = data[id.x] * 2.0; }`;
                        const mod = device.createShaderModule({ code: shaderCode });
                        const pipe = device.createComputePipeline({ layout: 'auto', compute: { module: mod, entryPoint: 'main' } });
                        const buf = device.createBuffer({ size: 4, usage: GPUBufferUsage.STORAGE | GPUBufferUsage.COPY_SRC | GPUBufferUsage.COPY_DST });
                        device.queue.writeBuffer(buf, 0, new Float32Array([21.0]));
                        const bg = device.createBindGroup({ layout: pipe.getBindGroupLayout(0), entries: [{ binding: 0, resource: { buffer: buf } }] });
                        const enc = device.createCommandEncoder();
                        const p = enc.beginComputePass(); p.setPipeline(pipe); p.setBindGroup(0, bg); p.dispatchWorkgroups(1); p.end();
                        const readBuf = device.createBuffer({ size: 4, usage: GPUBufferUsage.COPY_DST | GPUBufferUsage.MAP_READ });
                        enc.copyBufferToBuffer(buf, 0, readBuf, 0, 4);
                        device.queue.submit([enc.finish()]);
                        await readBuf.mapAsync(GPUMapMode.READ);
                        const res = new Float32Array(readBuf.getMappedRange())[0];
                        readBuf.unmap();
                        logger.info(`Compute Test Result: ${res === 42 ? "PASS (42)" : "FAIL (" + res + ")"}`);
                        device.destroy(); // Cleanup diagnostic device
                    } else {
                        logger.error("Device: FAILED to create.");
                    }
                } else {
                    logger.error("Adapter: FAILED to find any adapter.");
                }
            } catch (e) {
                logger.error(`Diagnostic Error: ${e.message}`);
                console.error(e);
            }
        }
        logger.info("----------------------------------");

        // 4. Renderer
        if (this.useRawWebGPU) {
            const canvas = document.createElement('canvas');
            canvas.style.width = '100%';
            canvas.style.height = '100%';
            canvas.tabIndex = 0;
            canvas.style.outline = 'none';
            canvas.style.position = 'absolute';
            canvas.style.top = '0';
            canvas.style.left = '0';
            canvas.style.zIndex = '0';
            canvas.addEventListener('pointerdown', () => canvas.focus());
            this.container.style.position = 'relative';
            this.container.appendChild(canvas);

            this.raw = new RawWebGPUAtomsRenderer(canvas);
            await this.raw.init({ maxAtoms: 65536 });

            console.warn('[MolGUIApp] Raw WebGPU mode enabled. WebGL/Three overlay (gizmo) is DISABLED to keep renderer pure WebGPU.');

            // Provide minimal renderer-like object for Editor/controls
            this.renderer = {
                domElement: canvas,
                setSize: (w, h) => this.raw.resize(w, h),
                setPixelRatio: () => {},
                setAnimationLoop: () => {}
            };

            window.logger.info('Initialization Complete. Backend: WebGPU (Raw)');
        } else {
            // Three.js WebGPURenderer (may fallback)
            this.renderer = new WebGPURenderer({
                antialias: true,
                trackTimestamp: true,
                forceWebGPU: true,
                forceWebGL: false,
            });

            try {
                if (this.renderer.init) await this.renderer.init();
            } catch (e) {
                logger.error(`WebGPURenderer init failed: ${e?.message || e}`);
                console.error(e);
            }
            this.renderer.setSize(width, height);
            this.renderer.setPixelRatio(window.devicePixelRatio);
            this.renderer.domElement.tabIndex = 0;
            this.renderer.domElement.style.outline = 'none';
            this.renderer.domElement.addEventListener('pointerdown', () => this.renderer.domElement.focus());
            this.container.appendChild(this.renderer.domElement);

            const isWebGPU = this.renderer.backend && this.renderer.backend.isWebGPU;
            window.logger.info(`Initialization Complete. Backend: ${isWebGPU ? "WebGPU" : "WebGL (Fallback)"}`);
            if (!isWebGPU) {
                const backend = this.renderer.backend;
                window.logger.error(`WebGPU unavailable. backend: ${backend ? backend.type || 'unknown' : 'none'}`);
            }
        }
        if (OrbitControls) {
            this.controls = new OrbitControls(this.camera, this.renderer.domElement);
            this.controls.enableDamping = false;
            this.controls.dampingFactor = 0.05;
            this.controls.addEventListener('change', () => this.requestRender());
            this.controls.mouseButtons = {
                LEFT: null,
                MIDDLE: THREE.MOUSE.DOLLY,
                RIGHT: THREE.MOUSE.ROTATE
            };
            window.addEventListener('keydown', (e) => {
                if (e.key === 'Shift') this.controls.mouseButtons.RIGHT = THREE.MOUSE.PAN;
            });
            window.addEventListener('keyup', (e) => {
                if (e.key === 'Shift') this.controls.mouseButtons.RIGHT = THREE.MOUSE.ROTATE;
            });
        }

        // 5. Molecule Systems
        this.systems = {
            molecule: new EditableMolecule(),
            substrate: new EditableMolecule()
        };
        this.system = this.systems.molecule;

        this.packedSystems = {
            molecule: new PackedMolecule(),
            substrate: new PackedMolecule()
        };

        this.mmParams = new MMParams();
        await this.mmParams.loadResources('../common_resources/ElementTypes.dat', '../common_resources/AtomTypes.dat', '../common_resources/BondTypes.dat', '../common_resources/AngleTypes.dat');

        // 6. Molecule Renderers
        // Pass null for shaders, we generate them in MeshRenderer now
        this.renderers = {
            molecule: new MoleculeRenderer(this.scene, this.packedSystems.molecule, null, this.mmParams, this.systems.molecule),
            substrate: new MoleculeRenderer(this.scene, this.packedSystems.substrate, null, this.mmParams, this.systems.substrate)
        };
        this.molRenderer = this.renderers.molecule;

        await this._autoLoadDefaultMolecule();

        Object.values(this.renderers).forEach(r => {
            if (r && r.toggleAxes) r.toggleAxes(true);
        });

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

        // --- GUI & Editor ---
        this.gui = new GUI(this.system, this.molRenderer);
        this.editor = new Editor(this.scene, this.camera, this.renderer, this.system, this.molRenderer);
        if (this.useRawWebGPU && this.overlayScene && this.editor && this.editor.gizmo && this.editor.gizmo.getHelper) {
            this.overlayScene.add(this.editor.gizmo.getHelper());
        }
        this.scriptRunner = new ScriptRunner(this);
        this.shortcuts = new ShortcutManager(this.editor);

        this.editor.onSelectionChange = () => {
            this.gui.updateSelectionUI();
            this.molRenderer.updateSelection();
            this.requestRender();
        };

        this.gui.onSelectionChanged = () => {
            this.molRenderer.updateSelection();
        };

        this.molRenderer.update();

        const ambientLight = new THREE.AmbientLight(0xffffff, 0.5);
        this.scene.add(ambientLight);
        const dirLight = new THREE.DirectionalLight(0xffffff, 1);
        dirLight.position.set(5, 10, 7);
        this.scene.add(dirLight);

        window.addEventListener('resize', this.onWindowResize.bind(this));
        this.requestRender();
        const finalIsWebGPU = this.useRawWebGPU ? true : (this.renderer.backend && this.renderer.backend.isWebGPU);
        window.logger.info(`Initialization Complete. Backend: ${finalIsWebGPU ? "WebGPU" : "WebGL (Fallback)"}`);

        // Initialize XPDB WebGPU simulation
        this.initXPDB();
    }

    async initXPDB() {
        if (!navigator.gpu) {
            window.logger.warn("WebGPU not available, XPDB simulation disabled");
            return;
        }

        try {
            const nAtoms = this.system.nAtoms || 0;
            if (nAtoms === 0) {
                window.logger.info("XPDB: No atoms yet, will initialize when molecule is loaded");
                return;
            }

            this.xpdb = new XPDB_WebGPU(nAtoms);
            await this.xpdb.init();
            window.logger.info(`XPDB initialized with ${nAtoms} atoms`);
            this.xpdbDirty = true;
        } catch (e) {
            window.logger.error(`Failed to initialize XPDB: ${e.message}`);
            console.error(e);
        }
    }

    async updateXPDBTopology() {
        if (!this.xpdb || !this.xpdbEnabled) return;

        try {
            const nAtoms = this.system.nAtoms || 0;
            if (nAtoms === 0) return;

            // Re-initialize if atom count changed
            if (this.xpdb.numAtoms !== nAtoms) {
                this.xpdb = new XPDB_WebGPU(nAtoms);
                await this.xpdb.init();
            }

            // Upload atoms
            this.xpdb.uploadAtomsFromMolecule(this.system, this.mmParams);

            // Build topology
            const { bondsAdj } = buildXPDBTopology(this.system, this.mmParams, {
                includeAngleConstraints: !!this.xpdbUseAngles,
                maxBonds: 16
            });

            if ((window.VERBOSITY_LEVEL | 0) >= 3) {
                const counts = bondsAdj.map(bs => (bs ? bs.length : 0));
                console.log('[MolGUIApp.updateXPDBTopology][DEBUG] xpdbUseAngles=', this.xpdbUseAngles);
                console.log('[MolGUIApp.updateXPDBTopology][DEBUG] bondsAdj counts=', counts);
                console.log('[MolGUIApp.updateXPDBTopology][DEBUG] bondsAdj=', bondsAdj);
            }

            // Upload bonds
            this.xpdb.uploadBonds(bondsAdj, 0.4, 200, 120, false);

            this.xpdbDirty = false;
            window.logger.info(`XPDB topology updated: ${nAtoms} atoms`);
        } catch (e) {
            window.logger.error(`Failed to update XPDB topology: ${e.message}`);
            console.error(e);
        }
    }

    stepXPDB() {
        if (!this.xpdb || !this.xpdbEnabled || this.xpdbDirty || this.xpdbPaused) return;

        try {
            const maxRadius = getMaxRadius(this.system, this.mmParams, 1.0);
            this.xpdb.step(
                this.xpdbParams.dt,
                this.xpdbParams.iterations,
                this.xpdbParams.k_coll,
                this.xpdbParams.omega,
                this.xpdbParams.momentum_beta,
                null, // mousePos
                -1,   // pickedIdx
                maxRadius
            );

            // Update molecule positions from GPU
            this.syncXPDBPositions();
        } catch (e) {
            window.logger.error(`XPDB step failed: ${e.message}`);
            console.error(e);
        }
    }

    async dumpXPDBTopology() {
        if (!this.system || !this.mmParams) {
            console.warn('[dumpXPDBTopology] System or mmParams missing');
            return;
        }
        const ready = await this.ensureXPDBReady(true);
        if (!ready) return;
        console.log('[dumpXPDBTopology] === XPDB Topology Dump ===');
        console.log('[dumpXPDBTopology] nAtoms=', this.system.nAtoms);
        console.log('[dumpXPDBTopology] xpdbUseAngles=', this.xpdbUseAngles);
        try {
            const { bondsAdj, stats } = buildXPDBTopology(this.system, this.mmParams, {
                includeAngleConstraints: !!this.xpdbUseAngles,
                maxBonds: 16
            });
            console.log('[dumpXPDBTopology] stats=', stats);
            const counts = bondsAdj.map(bs => (bs ? bs.length : 0));
            const maxDeg = Math.max(...counts);
            const avgDeg = counts.reduce((a, b) => a + b, 0) / counts.length;
            console.log('[dumpXPDBTopology] bondsAdj counts (avg/max)=', avgDeg.toFixed(2), maxDeg);
            const lines = bondsAdj.map((bs, i) => {
                const parts = bs.map(e => `${e[0]}:${(+e[1]).toFixed(4)}/${(+e[2]).toFixed(1)}`);
                return `${i}: (${bs.length}) ${parts.join('  ')}`;
            });
            console.log('[dumpXPDBTopology] bondsAdj formatted:\n' + lines.join('\n'));
        } catch (e) {
            console.error('[dumpXPDBTopology] Failed to build topology:', e);
        }
    }

    async xpdbStepOnce() {
        const ready = await this.ensureXPDBReady(true);
        if (!ready) return;
        this.xpdbEnabled = true;
        this.xpdbPaused = true;
        this.xpdbDirty = false;
        try {
            const maxRadius = getMaxRadius(this.system, this.mmParams, 1.0);
            const kColl = this.xpdbDisableCollisions ? 0.0 : this.xpdbParams.k_coll;
            const bondScale = this.xpdbDisableBonds ? 0.0 : 1.0;
            if (this.xpdb.setBondStiffnessScale) this.xpdb.setBondStiffnessScale(bondScale);
            this.xpdb.step(
                this.xpdbParams.dt,
                this.xpdbParams.iterations,
                kColl,
                this.xpdbParams.omega,
                this.xpdbParams.momentum_beta,
                null, -1, maxRadius
            );
            if (this.xpdb.setBondStiffnessScale && this.xpdbDisableBonds) this.xpdb.setBondStiffnessScale(1.0);
            this.syncXPDBPositions();
            if (this.renderers?.molecule?.update) this.renderers.molecule.update();
            this.requestRender();
            this.logRealBondLengths('after step once');
            window.logger.info('XPDB single step completed');
        } catch (e) {
            window.logger.error(`XPDB step failed: ${e.message}`);
            console.error(e);
        }
    }

    async xpdbRelaxSteps(n = 10) {
        const ready = await this.ensureXPDBReady(true);
        if (!ready) return;
        this.xpdbEnabled = true;
        this.xpdbPaused = true;
        this.xpdbDirty = false;
        try {
            const maxRadius = getMaxRadius(this.system, this.mmParams, 1.0);
            const kColl = this.xpdbDisableCollisions ? 0.0 : this.xpdbParams.k_coll;
            const bondScale = this.xpdbDisableBonds ? 0.0 : 1.0;
            if (this.xpdb.setBondStiffnessScale) this.xpdb.setBondStiffnessScale(bondScale);
            for (let i = 0; i < n; i++) {
                this.xpdb.step(
                    this.xpdbParams.dt,
                    this.xpdbParams.iterations,
                    kColl,
                    this.xpdbParams.omega,
                    this.xpdbParams.momentum_beta,
                    null, -1, maxRadius
                );
                this.syncXPDBPositions();
            }
            if (this.xpdb.setBondStiffnessScale && this.xpdbDisableBonds) this.xpdb.setBondStiffnessScale(1.0);
            if (this.renderers?.molecule?.update) this.renderers.molecule.update();
            this.requestRender();
            this.logRealBondLengths(`after relax ${n} steps`);
            window.logger.info(`XPDB relaxed ${n} steps`);
        } catch (e) {
            window.logger.error(`XPDB relax failed: ${e.message}`);
            console.error(e);
        }
    }

    setXPDBPaused(paused) {
        this.xpdbPaused = !!paused;
    }

    setXPDBDisableBonds(disable) {
        this.xpdbDisableBonds = !!disable;
    }

    setXPDBDisableCollisions(disable) {
        this.xpdbDisableCollisions = !!disable;
    }

    logRealBondLengths(label = '') {
        if (!this.system || !this.system.bonds || this.system.bonds.length === 0) {
            console.warn('[logRealBondLengths] No bonds available');
            return;
        }
        const bnds = this.system.bonds;
        const atoms = this.system.atoms;
        const lens = new Array(bnds.length);
        let minL = 1e9, maxL = -1e9, sumL = 0;
        for (let i = 0; i < bnds.length; i++) {
            const b = bnds[i];
            if (typeof b.ensureIndices === 'function') b.ensureIndices(this.system);
            const ia = b.a ?? b.i;
            const ib = b.b ?? b.j;
            const A = atoms[ia];
            const B = atoms[ib];
            if (!A || !B) continue;
            const dx = A.pos.x - B.pos.x;
            const dy = A.pos.y - B.pos.y;
            const dz = A.pos.z - B.pos.z;
            const l = Math.sqrt(dx*dx + dy*dy + dz*dz);
            lens[i] = l;
            if (l < minL) minL = l;
            if (l > maxL) maxL = l;
            sumL += l;
        }
        const avgL = sumL / lens.length;
        const sample = lens.slice(0, Math.min(10, lens.length)).map(v => v.toFixed(4));
        console.log(`[logRealBondLengths] ${label} bonds=${lens.length} min=${minL.toFixed(4)} max=${maxL.toFixed(4)} avg=${avgL.toFixed(4)} sample=${sample.join(', ')}`);
    }

    async dumpXPDBBuffers() {
        if (!this.system) {
            console.warn('[dumpXPDBBuffers] System missing');
            return;
        }
        const ready = await this.ensureXPDBReady(true);
        if (!ready) return;
        console.log('[dumpXPDBBuffers] === XPDB Buffer Dump ===');
        console.log('[dumpXPDBBuffers] nAtoms=', this.system.nAtoms);
        console.log('[dumpXPDBBuffers] xpdbParams=', this.xpdbParams);
        try {
            const posData = await this.xpdb.readPositions();
            console.log('[dumpXPDBBuffers] pos0:');
            for (let i = 0; i < this.system.nAtoms; i++) {
                const x = posData[i * 3 + 0];
                const y = posData[i * 3 + 1];
                const z = posData[i * 3 + 2];
                console.log(`${x.toFixed(6)} ${y.toFixed(6)} ${z.toFixed(6)}`);
            }
        } catch (e) {
            console.error('[dumpXPDBBuffers] Failed to read positions:', e);
        }
    }

    async syncXPDBPositions() {
        if (!this.xpdb) return;

        try {
            const posData = await this.xpdb.readPositions();

            // Update molecule atom positions
            for (let i = 0; i < this.system.nAtoms; i++) {
                const atom = this.system.atoms[i];
                if (atom) {
                    atom.pos.x = posData[i * 3];
                    atom.pos.y = posData[i * 3 + 1];
                    atom.pos.z = posData[i * 3 + 2];
                }
            }

            // Mark geometry as dirty
            this.system._touchGeom();
        } catch (e) {
            window.logger.error(`Failed to sync XPDB positions: ${e.message}`);
            console.error(e);
        }
    }

    async toggleXPDB(enabled) {
        this.xpdbEnabled = !!enabled;
        if (this.xpdbEnabled) {
            window.logger.info("XPDB simulation enabled");
            if (!this.xpdb) {
                await this.initXPDB();
            }
            if (!this.xpdb) {
                window.logger.error("XPDB init failed; cannot enable simulation");
                this.xpdbEnabled = false;
                return;
            }
            await this.updateXPDBTopology();
        } else {
            window.logger.info("XPDB simulation disabled");
        }
    }

    setXPDBTopologyParams(params) {
        if (!params) return;
        if (params.useAngles !== undefined) this.xpdbUseAngles = !!params.useAngles;
        this.xpdbDirty = true;
    }

    setXPDBParams(params) {
        Object.assign(this.xpdbParams, params);
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
        if (this.overlayRenderer) this.overlayRenderer.setSize(width, height, false);
        this.requestRender();
    }

    async _autoLoadDefaultMolecule() {
        const url = '../common_resources/mol/H2O.mol2';
        try {
            console.log('[MolGUIApp] auto-loading default molecule', url);
            const resp = await fetch(url);
            if (!resp.ok) throw new Error(`HTTP ${resp.status}`);
            const text = await resp.text();
            const parsed = EditableMolecule.parseMol2(text);
            if (!parsed || !parsed.pos || parsed.pos.length <= 0) throw new Error('parsed molecule empty');
            if (typeof this.system.clear === 'function') this.system.clear();
            this.system.appendParsedSystem(parsed);
            this.molRenderer.setLabelMode('id');
            if (this.renderers?.molecule?.update) this.renderers.molecule.update();
            if (this.raw && typeof this.raw.setLabelsVisible === 'function') this.raw.setLabelsVisible(true);
            this.requestRender();
            logger.info('[MolGUIApp.init] Auto-loaded H2O with Atom ID labels');
        } catch (err) {
            logger.error(`[MolGUIApp.init] Failed to auto-load H2O: ${err?.message || err}`);
        }
    }
}

window.onload = () => {
    const app = new MolGUIApp();
    window.app = app;
    app.init();
};