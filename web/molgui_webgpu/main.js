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
import { XPDB_CPU } from './XPDB_CPU.js';
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
        this.debugRawLogs = false;
        this._rawLastLabelMode = 'none';
        this._rawLastLabelsVisible = false;
        this._rawLabelsDirty = true;
        this._rawAtomsDirty = true;
        this._rawBondsDirty = true;
        this._rawSelectionDirty = true;
        this.system = new EditableMolecule();
        this.controls = null;
        this._rafPending = false;
        this._needsExtraFrame = false;
        this.continuousRender = false;

        // XPDB WebGPU simulation state
        this.xpdb = null;
        this.xpdbCPU = null;
        this.xpdbEnabled = false;
        this.xpdbUseAngles = true;
        this.xpdbUseCPU = false;
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
        if (this._rafPending) {
            this._needsExtraFrame = true;
            return;
        }
        this._rafPending = true;
        requestAnimationFrame(async () => {
            this._rafPending = false;

            // Step XPDB simulation if enabled
            if (this.xpdbEnabled) {
                await this.stepXPDB();
            }

            if (!this.useRawWebGPU && this.renderers) {
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

                const activeEditable = this.systems ? this.systems.molecule : null;
                const sys = this.packedSystems ? this.packedSystems.molecule : null;

                if (activeEditable && sys && (activeEditable.dirtyExport || activeEditable.dirtyTopo || activeEditable.dirtyGeom)) {
                    console.log('[MolGUIApp/raw] exporting EditableMolecule -> PackedMolecule', {
                        dirtyExport: activeEditable.dirtyExport,
                        dirtyTopo: activeEditable.dirtyTopo,
                        dirtyGeom: activeEditable.dirtyGeom,
                        nAtoms: activeEditable.atoms.length
                    });
                    activeEditable.exportToMoleculeSystem(sys);
                    this._markRawAllDirty('exportEditable');
                }

                if (sys) {
                    if (this.debugRawLogs) {
                        console.log('[MolGUIApp/raw] packed snapshot', {
                            nAtoms: sys.nAtoms,
                            nBonds: sys.bonds ? sys.bonds.length : 0,
                            dirtyAtoms: this._rawAtomsDirty,
                            dirtyBonds: this._rawBondsDirty,
                            dirtyLabels: this._rawLabelsDirty,
                            dirtySel: this._rawSelectionDirty
                        });
                    }
                    const n = sys.nAtoms | 0;
                    const pos = sys.pos;
                    const sel = this.systems.molecule.selection;

                    // --- Labels visibility is controlled regardless of atom count ---
                    const mode = (this.molRenderer && this.molRenderer.labelMode) ? this.molRenderer.labelMode : 'none';
                    const wantLabels = (mode !== 'none');
                    if (wantLabels !== this._rawLastLabelsVisible) {
                        this.raw.setLabelsVisible(wantLabels);
                        this._rawLastLabelsVisible = wantLabels;
                    }

                    if (n > 0) {
                        if (this._rawAtomsDirty) {
                            const col = new Float32Array(n * 3);
                            const rad = new Float32Array(n);
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
                            console.log('[MolGUIApp/raw] updateAtoms', { nAtoms: n });
                            this._rawAtomsDirty = false;
                        }

                        if (this._rawBondsDirty) {
                            if (sys.bonds && sys.bonds.length) {
                                this.raw.updateBonds(sys.bonds);
                            } else {
                                this.raw.updateBonds([]);
                            }
                            console.log('[MolGUIApp/raw] updateBonds', { nBonds: sys.bonds ? sys.bonds.length : 0 });
                            this._rawBondsDirty = false;
                        }

                        if (this._rawSelectionDirty) {
                            const selIdx = sel ? Array.from(sel) : [];
                            this.raw.updateSelection(selIdx, n, 1.3);
                             console.log('[MolGUIApp/raw] updateSelection', { selectionSize: selIdx.length });
                            this._rawSelectionDirty = false;
                        }

                        if (wantLabels) {
                            if (this._rawLabelsDirty || mode !== this._rawLastLabelMode) {
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
                                this.raw.updateLabels(strings, n);
                                this._rawLabelsDirty = false;
                                this._rawLastLabelMode = mode;
                            }
                        } else if (this._rawLastLabelsVisible) {
                            this.raw.updateLabels(null, 0);
                            this._rawLabelsDirty = false;
                        }
                    } else {
                        if (this._rawAtomsDirty) {
                            this.raw.nAtoms = 0;
                            this._rawAtomsDirty = false;
                        }
                        if (this._rawBondsDirty) {
                            this.raw.updateBonds([]);
                            this._rawBondsDirty = false;
                        }
                        if (this._rawLabelsDirty || this._rawLastLabelsVisible) {
                            this.raw.updateLabels(null, 0);
                            this._rawLabelsDirty = false;
                        }
                        if (this._rawSelectionDirty) {
                            this.raw.updateSelection([], 0, 1.3);
                            this._rawSelectionDirty = false;
                        }
                    }
                }

                this.raw.render();
                this._maybeScheduleExtraFrame();

                // Pure WebGPU mode: no WebGL overlay rendering.
                return;
            }

            if (this.renderer && this.scene && this.camera) {
                if (this.renderer.renderAsync) await this.renderer.renderAsync(this.scene, this.camera);
                else this.renderer.render(this.scene, this.camera);
            }

            this._maybeScheduleExtraFrame();
        });
    }

    _maybeScheduleExtraFrame() {
        if (!this._needsExtraFrame) return;
        this._needsExtraFrame = false;
        this.requestRender();
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
            if (this.useRawWebGPU) this._markRawSelectionDirty('editorSelection');
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
            this.xpdbCPU = new XPDB_CPU(nAtoms);
            window.logger.info(`XPDB initialized with ${nAtoms} atoms`);
            this.xpdbDirty = true;
        } catch (e) {
            window.logger.error(`Failed to initialize XPDB: ${e.message}`);
            console.error(e);
        }
    }

    setXPDBUseCPU(useCPU) {
        this.xpdbUseCPU = !!useCPU;
        this.xpdbDirty = true;
        window.logger && window.logger.info(`[MolGUIApp.setXPDBUseCPU] xpdbUseCPU=${this.xpdbUseCPU}`);
        if (this.xpdbUseCPU) {
            if (this.xpdbCPU) {
                try {
                    this.xpdbCPU.uploadFromMolecule(this.system, this.mmParams);
                    if ((window.VERBOSITY_LEVEL | 0) >= 2) console.log('[XPDB_CPU] synced from system when enabling CPU mode');
                } catch (err) {
                    console.warn('[XPDB_CPU] failed to upload system when enabling CPU mode', err);
                }
            }
        } else {
            this.syncXPDBGPUFromSystem('switchToGPU');
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
                this.xpdbCPU = new XPDB_CPU(nAtoms);
            }

            // Upload atoms
            this.xpdb.uploadAtomsFromMolecule(this.system, this.mmParams);
            if (this.xpdbCPU) this.xpdbCPU.uploadFromMolecule(this.system, this.mmParams);

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

            // Keep CPU side topology for reference solver
            this._xpdbBondsAdj = bondsAdj;

            this.xpdbDirty = false;
            window.logger.info(`XPDB topology updated: ${nAtoms} atoms`);
        } catch (e) {
            window.logger.error(`Failed to update XPDB topology: ${e.message}`);
            console.error(e);
        }
    }

    async stepXPDB() {
        if (!this.xpdb || !this.xpdbEnabled || this.xpdbDirty || this.xpdbPaused) return;
        try {
            const maxRadius = getMaxRadius(this.system, this.mmParams, 1.0);
            if (this.xpdbUseCPU) {
                if (!this.xpdbCPU) throw new Error('XPDB CPU solver missing');
                if (!this._xpdbBondsAdj) throw new Error('XPDB CPU: missing bondsAdj (topology not built)');
                const kColl = this.xpdbDisableCollisions ? 0.0 : this.xpdbParams.k_coll;
                window.logger && window.logger.info('[XPDB][CPU] stepXPDB using CPU solver');
                this.xpdbCPU.step(this._xpdbBondsAdj, {
                    dt: this.xpdbParams.dt,
                    iterations: this.xpdbParams.iterations,
                    k_coll: kColl,
                    omega: this.xpdbParams.omega,
                    momentum_beta: this.xpdbParams.momentum_beta
                });
                if ((window.VERBOSITY_LEVEL | 0) >= 2) {
                    console.log(`[XPDB_CPU] stepXPDB iterations=${this.xpdbParams.iterations} dt=${this.xpdbParams.dt} k_coll=${kColl} lastDelta=${this.xpdbCPU.lastDelta.toExponential(4)}`);
                }
                this.xpdbCPU.writeToMolecule(this.system);
                this.system._touchGeom();
                this.debugLogCPUState('stepXPDB');
                this.syncXPDBGPUFromSystem('cpuStepXPDB');
            } else {
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
                await this.syncXPDBPositions();
            }
            this._markRawGeometryDirty('xpdbStep');
            this.requestRender();
        } catch (e) {
            window.logger.error(`XPDB step failed: ${e.message}`);
            console.error(e);
        }
    }

    async dumpXPDBBuffers() {
        if (!this.xpdb || !this.system) {
            console.warn('[dumpXPDBBuffers] XPDB or system not initialized');
            return;
        }
        const ready = await this.ensureXPDBReady(true);
        if (!ready) return;
        console.log('[dumpXPDBBuffers] === XPDB Buffer Dump ===');
        console.log('[dumpXPDBBuffers] nAtoms=', this.system.nAtoms);
        console.log('[dumpXPDBBuffers] xpdbParams=', this.xpdbParams);
        try {
            const { posData, ghostData, bondGlobalData, bondLocalData, bondLenStiffData } = await this.xpdb.readBuffersAsTyped();
            const nAtoms = this.system.nAtoms;
            const groupSize = this.xpdb.groupSize;
            const numGroups = this.xpdb.numGroups;
            const maxGhosts = this.xpdb.maxGhosts;
            const nMaxBonded = this.xpdb.nMaxBonded;

            console.log('[dumpXPDBBuffers] groupSize=', groupSize, 'numGroups=', numGroups, 'maxGhosts=', maxGhosts, 'nMaxBonded=', nMaxBonded);

            for (let g = 0; g < numGroups; g++) {
                const groupBase = g * groupSize;
                const ghostBase = g * (maxGhosts + 1);
                const ghostCount = ghostData[ghostBase + maxGhosts];
                console.log(`[dumpXPDBBuffers] === Group ${g} (global ${groupBase}-${groupBase + groupSize - 1}) ghostCount=${ghostCount} ===`);

                if (ghostCount > 0) {
                    console.log('[dumpXPDBBuffers] Ghost list:');
                    for (let k = 0; k < ghostCount; k++) {
                        const globalIdx = ghostData[ghostBase + k];
                        const px = posData[globalIdx * 4 + 0];
                        const py = posData[globalIdx * 4 + 1];
                        const pz = posData[globalIdx * 4 + 2];
                        const atom = this.system.atoms[globalIdx];
                        const typeName = atom ? (this.mmParams.getAtomTypeForAtom(atom)?.name || 'unknown') : 'missing';
                        console.log(`  ghost[${k}] -> global ${globalIdx} (${typeName}) pos=(${px.toFixed(3)}, ${py.toFixed(3)}, ${pz.toFixed(3)})`);
                    }
                }

                for (let lid = 0; lid < groupSize; lid++) {
                    const globalIdx = groupBase + lid;
                    if (globalIdx >= nAtoms) break;

                    const px = posData[globalIdx * 4 + 0];
                    const py = posData[globalIdx * 4 + 1];
                    const pz = posData[globalIdx * 4 + 2];
                    const atom = this.system.atoms[globalIdx];
                    const typeName = atom ? (this.mmParams.getAtomTypeForAtom(atom)?.name || 'unknown') : 'missing';
                    console.log(`[dumpXPDBBuffers] atom ${globalIdx} (lid ${lid}, group ${g}, type ${typeName}) pos=(${px.toFixed(3)}, ${py.toFixed(3)}, ${pz.toFixed(3)})`);

                    const bondBase = globalIdx * nMaxBonded;

                    console.log('  bond_indices_global:');
                    for (let slot = 0; slot < nMaxBonded; slot++) {
                        const targetGlobal = bondGlobalData[bondBase + slot];
                        if (targetGlobal === -1) break;
                        const L0 = bondLenStiffData[(bondBase + slot) * 2 + 0];
                        const K = bondLenStiffData[(bondBase + slot) * 2 + 1];
                        console.log(`    slot${slot}: global ${targetGlobal} L0=${L0.toFixed(4)} K=${K.toFixed(1)}`);
                    }

                    console.log('  bond_indices_local:');
                    for (let slot = 0; slot < nMaxBonded; slot++) {
                        const localIdx = bondLocalData[bondBase + slot];
                        if (localIdx === -1) break;
                        const targetGlobal = bondGlobalData[bondBase + slot];
                        const L0 = bondLenStiffData[(bondBase + slot) * 2 + 0];
                        const K = bondLenStiffData[(bondBase + slot) * 2 + 1];

                        let mappedGlobal = null;
                        let kind = 'unknown';
                        if (localIdx < groupSize) {
                            mappedGlobal = groupBase + localIdx;
                            kind = 'internal';
                        } else {
                            const ghostIdx = localIdx - groupSize;
                            if (ghostIdx < ghostCount) {
                                mappedGlobal = ghostData[ghostBase + ghostIdx];
                                kind = 'ghost';
                            } else {
                                kind = 'invalid';
                            }
                        }

                        const match = (mappedGlobal !== null && mappedGlobal === targetGlobal) ? 'MATCH' : 'MISMATCH';
                        console.log(`    slot${slot}: local ${localIdx} -> ${kind} ${mappedGlobal !== null ? mappedGlobal : '?'} (expected ${targetGlobal}) ${match} L0=${L0.toFixed(4)} K=${K.toFixed(1)}`);
                    }
                }
            }

            console.log('[dumpXPDBBuffers] === Real Bond Lengths (EditableMolecule.bonds) ===');
            if (this.system.bonds && this.system.bonds.length > 0) {
                const resolveIndices = (bond) => {
                    if (!bond) return [null, null];
                    if (Array.isArray(bond)) return bond;
                    if (typeof bond.ensureIndices === 'function') bond.ensureIndices(this.system);
                    const a = bond.a ?? bond.i ?? (bond.indices ? bond.indices[0] : null);
                    const b = bond.b ?? bond.j ?? (bond.indices ? bond.indices[1] : null);
                    return [a, b];
                };

                for (const bond of this.system.bonds) {
                    const [a, b] = resolveIndices(bond);
                    if (a === null || b === null) {
                        console.warn('[dumpXPDBBuffers] bond missing indices', bond);
                        continue;
                    }
                    if (a < 0 || a >= nAtoms || b < 0 || b >= nAtoms) {
                        console.warn(`[dumpXPDBBuffers] bond indices out of range: ${a}, ${b}`);
                        continue;
                    }
                    const ax = posData[a * 4 + 0], ay = posData[a * 4 + 1], az = posData[a * 4 + 2];
                    const bx = posData[b * 4 + 0], by = posData[b * 4 + 1], bz = posData[b * 4 + 2];
                    const dist = Math.sqrt((ax - bx) ** 2 + (ay - by) ** 2 + (az - bz) ** 2);
                    console.log(`[dumpXPDBBuffers] bond ${a}-${b}: dist=${dist.toFixed(4)} A=(${ax.toFixed(3)},${ay.toFixed(3)},${az.toFixed(3)}) B=(${bx.toFixed(3)},${by.toFixed(3)},${bz.toFixed(3)})`);
                }
            } else {
                console.log('[dumpXPDBBuffers] No bonds recorded in EditableMolecule');
            }

        } catch (e) {
            console.error('[dumpXPDBBuffers] Failed:', e);
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
            if (this.xpdbUseCPU) {
                if (!this.xpdbCPU) throw new Error('XPDB CPU solver missing');
                if (!this._xpdbBondsAdj) throw new Error('XPDB CPU: missing bondsAdj (topology not built)');
                window.logger && window.logger.info('[XPDB][CPU] xpdbStepOnce using CPU solver');
                this.xpdbCPU.step(this._xpdbBondsAdj, {
                    dt: this.xpdbParams.dt,
                    iterations: this.xpdbParams.iterations,
                    k_coll: kColl,
                    omega: this.xpdbParams.omega,
                    momentum_beta: this.xpdbParams.momentum_beta
                });
                console.log(`[XPDB_CPU] stepOnce iterations=${this.xpdbParams.iterations} dt=${this.xpdbParams.dt} k_coll=${kColl} lastDelta=${this.xpdbCPU.lastDelta.toExponential(4)}`);
                this.xpdbCPU.writeToMolecule(this.system);
                this.system._touchGeom();
                this.debugLogCPUState('stepOnce');
                this.syncXPDBGPUFromSystem('cpuStepOnce');
            } else {
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
                await this.syncXPDBPositions();
            }
            if (this.renderers?.molecule?.update) this.renderers.molecule.update();
            this._markRawGeometryDirty('xpdbStepOnce');
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
            if (this.xpdbUseCPU) {
                if (!this.xpdbCPU) throw new Error('XPDB CPU solver missing');
                if (!this._xpdbBondsAdj) throw new Error('XPDB CPU: missing bondsAdj (topology not built)');
                window.logger && window.logger.info(`[XPDB][CPU] xpdbRelaxSteps using CPU solver for ${n} outer steps`);
                for (let i = 0; i < n; i++) {
                    this.xpdbCPU.step(this._xpdbBondsAdj, {
                        dt: this.xpdbParams.dt,
                        iterations: this.xpdbParams.iterations,
                        k_coll: kColl,
                        omega: this.xpdbParams.omega,
                        momentum_beta: this.xpdbParams.momentum_beta
                    });
                    if ((window.VERBOSITY_LEVEL | 0) >= 3) {
                        console.log(`[XPDB_CPU] relax outer=${i} lastDelta=${this.xpdbCPU.lastDelta.toExponential(4)}`);
                    }
                }
                console.log(`[XPDB_CPU] relax completed lastDelta=${this.xpdbCPU.lastDelta.toExponential(4)}`);
                this.xpdbCPU.writeToMolecule(this.system);
                this.system._touchGeom();
                this.debugLogCPUState(`relax${n}`);
                this.syncXPDBGPUFromSystem('cpuRelax');
            } else {
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
                    await this.syncXPDBPositions();
                }
                if (this.xpdb.setBondStiffnessScale && this.xpdbDisableBonds) this.xpdb.setBondStiffnessScale(1.0);
            }
            if (this.renderers?.molecule?.update) this.renderers.molecule.update();
            this._markRawGeometryDirty('xpdbRelax');
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

    debugLogCPUState(label = 'cpu') {
        if (!this.xpdbCPU || !this.system) return;
        const verbosity = window.VERBOSITY_LEVEL | 0;
        if (verbosity < 1) return;
        const nAtoms = this.system.nAtoms || 0;
        console.log(`[XPDB_CPU][STATE] ${label} nAtoms=${nAtoms} lastDelta=${this.xpdbCPU.lastDelta.toExponential(6)}`);
        const maxAtoms = Math.min(4, nAtoms);
        for (let i = 0; i < maxAtoms; i++) {
            const atom = this.system.atoms[i];
            const px = this.xpdbCPU.pos[i * 3 + 0];
            const py = this.xpdbCPU.pos[i * 3 + 1];
            const pz = this.xpdbCPU.pos[i * 3 + 2];
            const labelName = atom?.label || atom?.element || atom?.name || 'atom';
            console.log(`  atom ${i} (${labelName}) -> (${px.toFixed(4)}, ${py.toFixed(4)}, ${pz.toFixed(4)})`);
        }
        const bonds = this.system.bonds || [];
        if (bonds.length > 0) {
            const maxBonds = Math.min(5, bonds.length);
            let logged = 0;
            const resolve = (bond) => {
                if (!bond) return [null, null];
                if (Array.isArray(bond)) return bond;
                if (typeof bond.ensureIndices === 'function') bond.ensureIndices(this.system);
                const a = bond.a ?? bond.i ?? (bond.indices ? bond.indices[0] : null);
                const b = bond.b ?? bond.j ?? (bond.indices ? bond.indices[1] : null);
                return [a, b];
            };
            for (const bond of bonds) {
                if (logged >= maxBonds) break;
                const [a, b] = resolve(bond);
                if (a === null || b === null) continue;
                const ax = this.xpdbCPU.pos[a * 3 + 0];
                const ay = this.xpdbCPU.pos[a * 3 + 1];
                const az = this.xpdbCPU.pos[a * 3 + 2];
                const bx = this.xpdbCPU.pos[b * 3 + 0];
                const by = this.xpdbCPU.pos[b * 3 + 1];
                const bz = this.xpdbCPU.pos[b * 3 + 2];
                const dx = ax - bx;
                const dy = ay - by;
                const dz = az - bz;
                const dist = Math.sqrt(dx * dx + dy * dy + dz * dz);
                console.log(`  bond ${a}-${b} length=${dist.toFixed(4)}`);
                logged++;
            }
        } else {
            console.log('  [XPDB_CPU][STATE] no bonds recorded');
        }

        if (verbosity >= 2) {
            this.logRealBondLengths(`[XPDB_CPU ${label}]`);
        }
    }

    syncXPDBGPUFromSystem(reason = 'manual') {
        if (!this.xpdb || !this.system || !this.mmParams) return;
        try {
            this.xpdb.uploadAtomsFromMolecule(this.system, this.mmParams);
            if ((window.VERBOSITY_LEVEL | 0) >= 2) {
                console.log(`[XPDB][syncGPU] reason=${reason} uploaded system positions to GPU buffers`);
            }
        } catch (err) {
            console.warn('[XPDB][syncGPU] Failed to upload atoms from system:', err);
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
            if (this.useRawWebGPU) this._markRawGeometryDirty('syncXPDBPositions');
            if (this.xpdbCPU && !this.xpdbUseCPU) {
                try {
                    this.xpdbCPU.uploadFromMolecule(this.system, this.mmParams);
                    if ((window.VERBOSITY_LEVEL | 0) >= 3) console.log('[XPDB_CPU] synced from system after GPU update');
                } catch (err) {
                    console.warn('[XPDB_CPU] failed to sync from GPU positions', err);
                }
            }
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
            this.markRawAllDirty('autoLoadDefault');
            if (this.raw && typeof this.raw.setLabelsVisible === 'function') this.raw.setLabelsVisible(true);
            this.requestRender();
            logger.info('[MolGUIApp.init] Auto-loaded H2O with Atom ID labels');
        } catch (err) {
            logger.error(`[MolGUIApp.init] Failed to auto-load H2O: ${err?.message || err}`);
        }
    }

    _markRawAllDirty(reason = 'manual') {
        if (!this.useRawWebGPU) return;
        this._rawAtomsDirty = true;
        this._rawBondsDirty = true;
        this._rawLabelsDirty = true;
        this._rawSelectionDirty = true;
        console.log('[MolGUIApp/raw] markAllDirty', { reason });
    }

    markRawAllDirty(reason = 'manualExternal') {
        this._markRawAllDirty(reason);
    }

    _markRawTopologyDirty(reason = 'topology') {
        this._markRawAllDirty(reason);
    }

    markRawTopologyDirty(reason = 'topologyExternal') {
        if (!this.useRawWebGPU) return;
        this._markRawTopologyDirty(reason);
    }

    _markRawGeometryDirty(reason = 'geometry') {
        if (!this.useRawWebGPU) return;
        this._rawAtomsDirty = true;
        this._rawLabelsDirty = true;
        console.log('[MolGUIApp/raw] markGeometryDirty', { reason });
    }

    markRawGeometryDirty(reason = 'geometryExternal') {
        this._markRawGeometryDirty(reason);
    }

    _markRawSelectionDirty(reason = 'manual') {
        if (!this.useRawWebGPU) return;
        this._rawSelectionDirty = true;
        this._rawLabelsDirty = true;
        console.log('[MolGUIApp/raw] markSelectionDirty', { reason });
    }

    _markRawLabelsDirty(reason = 'manual') {
        if (!this.useRawWebGPU) return;
        this._rawLabelsDirty = true;
        console.log('[MolGUIApp/raw] markLabelsDirty', { reason });
    }
}

window.onload = () => {
    const app = new MolGUIApp();
    window.app = app;
    app.init();
};