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
import { getMaxRadius } from './XPDBTopology.js';
import { buildXPDBInputsFromMol } from './MMFFLTopology.js';
import { XPDBTopologyRenderer } from './XPDBTopologyRenderer.js';

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

        // XPDB topology state (real vs dummies)
        this.xpdbNReal = 0;
        this.xpdbNAll = 0;

        // XPDB topology visualization overlay (separate from EditableMolecule renderer)
        this.xpdbTopoRenderer = null;
        this.xpdbTopoVisible = false;
        this.xpdbTopoVisibility = {
            atoms_real: true,
            atoms_dummy: true,
            primary: true,
            angle: true,
            pi: true,
            pi_align: true,
            epair: true,
            epair_pair: true,
        };
        this._xpdbTopoLast = null;
    }

    async _rebuildXPDBTopoCPU(reason = 'manual') {
        const verbosity = window.VERBOSITY_LEVEL | 0;
        const nReal = this.system.nAtoms || 0;
        if (!(nReal > 0)) {
            if (verbosity >= 1) console.warn(`[XPDBTopoCPU] skip rebuild reason=${reason} (no atoms)`);
            this._xpdbTopoLast = null;
            return null;
        }
        if (!this.mmParams) throw new Error(`[XPDBTopoCPU] mmParams missing (reason=${reason})`);
        try {
            const molXPDB = this._cloneSystemForXPDB();
            const packed = buildXPDBInputsFromMol(molXPDB, this.mmParams, {
                report_types: false,
                add_angle: !!this.xpdbUseAngles,
                add_pi: true,
                two_pi: true,
                add_pi_align: false,
                add_epair: true,
                add_epair_pairs: false,
                L_pi: 1.0,
                L_epair: 0.5,
                k_angle: 100.0,
                k_pi: 50.0,
                k_pi_orth: 30.0,
                k_pi_align: 15.0,
                k_ep: 40.0,
                k_ep_orth: 25.0,
                k_ep_pair: 10.0,
                defaultK: 200.0,
                nMaxBonded: 16,
            });
            this._xpdbTopoLast = packed.topo;
            if (verbosity >= 2) {
                console.log('[XPDBTopoCPU] rebuilt', {
                    reason,
                    n_real: packed.topo?.n_real | 0,
                    n_all: packed.topo?.n_all | 0,
                    bonds_primary: packed.topo?.bonds_primary ? packed.topo.bonds_primary.length : 0,
                });
            }
            return packed;
        } catch (e) {
            this._xpdbTopoLast = null;
            if (this.xpdbTopoRenderer) this.xpdbTopoRenderer.clear();
            throw e;
        }
    }

    async _refreshXPDBTopoOverlay(reason = 'manual') {
        const verbosity = window.VERBOSITY_LEVEL | 0;
        if (!this.xpdbTopoRenderer) return;

        if (!this.xpdbTopoVisible) {
            this.xpdbTopoRenderer.setEnabled(false);
            if (verbosity >= 3) console.log('[XPDBTopoOverlay] disabled', { reason });
            return;
        }

        // Ensure we have a topology even if user never enabled XPDB simulation.
        if (!this._xpdbTopoLast) {
            if (verbosity >= 2) console.log('[XPDBTopoOverlay] topo missing; rebuilding CPU topo', { reason });
            try {
                await this._rebuildXPDBTopoCPU(`overlay:${reason}`);
            } catch (e) {
                console.error('[XPDBTopoOverlay] rebuild failed', e);
                this.xpdbTopoRenderer.setEnabled(false);
                return;
            }
        }

        if (!this._xpdbTopoLast) {
            if (verbosity >= 1) console.warn('[XPDBTopoOverlay] still no topo after rebuild', { reason });
            this.xpdbTopoRenderer.setEnabled(false);
            return;
        }

        if (verbosity >= 2) {
            console.log('[XPDBTopoOverlay] update', {
                reason,
                visible: this.xpdbTopoVisible,
                n_all: this._xpdbTopoLast.n_all | 0,
                n_real: this._xpdbTopoLast.n_real | 0,
            });
        }
        this.xpdbTopoRenderer.updateFromTopo(this._xpdbTopoLast, { enabled: true, visibility: this.xpdbTopoVisibility });
    }

    _cloneSystemForXPDB() {
        // MMFFL topology builder mutates the molecule by appending dummy atoms.
        // Never run it on this.system directly.
        const m = new EditableMolecule();
        const parsed = this.system.exportAsParsed ? this.system.exportAsParsed() : null;
        if (!parsed) throw new Error('_cloneSystemForXPDB: system.exportAsParsed() missing');
        m.clear();
        m.appendParsedSystem(parsed);
        const v = window.VERBOSITY_LEVEL | 0;
        if (v >= 3) console.log('[MolGUIApp/_cloneSystemForXPDB] cloned system', { nAtoms: m.atoms.length, nBonds: m.bonds ? m.bonds.length : 0 });
        return m;
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

                // Optional WebGL overlay (topology debug, etc.) on top of raw WebGPU canvas.
                if (this.overlayRenderer && this.overlayScene && this.camera) {
                    try {
                        this.overlayRenderer.render(this.overlayScene, this.camera);
                    } catch (e) {
                        console.error('[MolGUIApp][overlayRenderer] render failed', e);
                    }
                }

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

            // WebGL overlay canvas for debug visuals (XPDB topology overlay etc.)
            // This does NOT affect the raw WebGPU renderer; it's a separate transparent layer.
            this.overlayCanvas = document.createElement('canvas');
            this.overlayCanvas.style.width = '100%';
            this.overlayCanvas.style.height = '100%';
            this.overlayCanvas.style.position = 'absolute';
            this.overlayCanvas.style.top = '0';
            this.overlayCanvas.style.left = '0';
            this.overlayCanvas.style.zIndex = '1';
            this.overlayCanvas.style.pointerEvents = 'none';
            this.container.appendChild(this.overlayCanvas);
            this.overlayRenderer = new THREE.WebGLRenderer({ canvas: this.overlayCanvas, alpha: true, antialias: true });
            this.overlayRenderer.setClearColor(0x000000, 0.0);
            this.overlayRenderer.setSize(width, height);
            this.overlayScene = new THREE.Scene();

            console.warn('[MolGUIApp] Raw WebGPU mode enabled. Using transparent WebGL overlay canvas for debug/topology visualization.');

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

        // 6b. XPDB topology overlay renderer (dummy atoms + constraint bonds)
        this.xpdbTopoRenderer = new XPDBTopologyRenderer(this.useRawWebGPU ? this.overlayScene : this.scene);
        this.xpdbTopoRenderer.setEnabled(!!this.xpdbTopoVisible);
        this.xpdbTopoRenderer.setVisibilityMap(this.xpdbTopoVisibility);

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

    setXPDBTopologyOverlayVisible(visible) {
        this.xpdbTopoVisible = !!visible;
        const verbosity = (window.VERBOSITY_LEVEL | 0);
        const dbg = !!window.DEBUG_XPDB_TOPO_VIZ;
        if (dbg) console.log('[MolGUIApp.setXPDBTopologyOverlayVisible]', { visible: this.xpdbTopoVisible, hasTopo: !!this._xpdbTopoLast });
        // async refresh (do not require checkbox handlers to be async)
        Promise.resolve().then(() => this._refreshXPDBTopoOverlay('setVisible')).catch((e) => {
            console.error('[XPDBTopoOverlay] refresh failed', e);
        });
        this.requestRender();
    }

    setXPDBTopologyOverlayVisibility(vis) {
        if (!vis) return;
        for (const k of Object.keys(vis)) this.xpdbTopoVisibility[k] = !!vis[k];
        const dbg = !!window.DEBUG_XPDB_TOPO_VIZ;
        if (dbg) console.log('[MolGUIApp.setXPDBTopologyOverlayVisibility]', { vis });
        // async refresh (may rebuild topo if missing)
        Promise.resolve().then(() => this._refreshXPDBTopoOverlay('setVisibility')).catch((e) => {
            console.error('[XPDBTopoOverlay] refresh failed', e);
        });
        this.requestRender();
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
            const nReal = this.system.nAtoms || 0;
            if (nReal === 0) return;

            const molXPDB = this._cloneSystemForXPDB();
            const packed = buildXPDBInputsFromMol(molXPDB, this.mmParams, {
                report_types: false,
                add_angle: !!this.xpdbUseAngles,
                add_pi: true,
                two_pi: true,
                add_pi_align: false,
                add_epair: true,
                add_epair_pairs: false,
                L_pi: 1.0,
                L_epair: 0.5,
                k_angle: 100.0,
                k_pi: 50.0,
                k_pi_orth: 30.0,
                k_pi_align: 15.0,
                k_ep: 40.0,
                k_ep_orth: 25.0,
                k_ep_pair: 10.0,
                defaultK: 200.0,
                nMaxBonded: 16,
            });

            const nAll = packed.nAtoms | 0;
            if (nAll !== (packed.topo?.n_all | 0)) {
                throw new Error(`updateXPDBTopology: packed.nAtoms mismatch topo.n_all packed=${nAll} topo.n_all=${packed.topo?.n_all}`);
            }
            if (nAll < nReal) throw new Error(`updateXPDBTopology: nAll=${nAll} < nReal=${nReal}`);

            // Re-initialize XPDB solvers if atom count changed (must match packed n_all)
            if (!this.xpdb || (this.xpdb.numAtoms !== nAll)) {
                this.xpdb = new XPDB_WebGPU(nAll);
                await this.xpdb.init();
            }
            if (!this.xpdbCPU || (this.xpdbCPU.nAtoms !== nAll)) {
                this.xpdbCPU = new XPDB_CPU(nAll);
            }

            // Upload packed atoms and topology to GPU
            this.xpdb.uploadPackedAtoms(packed);
            this.xpdb.uploadPackedTopology(packed);

            // Upload CPU positions/radii/mass from the packed molecule used for topology
            if (this.xpdbCPU) this.xpdbCPU.uploadFromMolecule(molXPDB, this.mmParams);

            this.xpdbNReal = nReal;
            this.xpdbNAll = nAll;
            this._xpdbBondsAdj = packed.bondsAdj;

            // Update topology visualization overlay (do not touch EditableMolecule rendering)
            this._xpdbTopoLast = packed.topo;
            if (this.xpdbTopoVisible) await this._refreshXPDBTopoOverlay('updateXPDBTopology');

            this.xpdbDirty = false;
            window.logger.info(`XPDB topology updated: n_real=${nReal} n_all=${nAll}`);
        } catch (e) {
            window.logger.error(`Failed to update XPDB topology: ${e.message}`);
            console.error(e);
        }
    }

    async dumpXPDBTopology() {
        const ready = await this.ensureXPDBReady(true);
        if (!ready) return;
        console.log('[dumpXPDBTopology] n_real=', this.xpdbNReal, 'n_all=', this.xpdbNAll, 'useAngles=', this.xpdbUseAngles);
        try {
            const molXPDB = this._cloneSystemForXPDB();
            const packed = buildXPDBInputsFromMol(molXPDB, this.mmParams, {
                report_types: false,
                add_angle: !!this.xpdbUseAngles,
                add_pi: true,
                two_pi: true,
                add_pi_align: false,
                add_epair: true,
                add_epair_pairs: false,
                L_pi: 1.0,
                L_epair: 0.5,
                k_angle: 100.0,
                k_pi: 50.0,
                k_pi_orth: 30.0,
                k_pi_align: 15.0,
                k_ep: 40.0,
                k_ep_orth: 25.0,
                k_ep_pair: 10.0,
                defaultK: 200.0,
                nMaxBonded: 16,
            });
            const topo = packed.topo;
            console.log('[dumpXPDBTopology] topo n_real=', topo.n_real, 'n_all=', topo.n_all);
            const logPairs = (label, arr) => {
                console.log(label + ' count=', arr?.length || 0);
                if (arr && arr.length) arr.forEach((p, idx) => console.log(`  ${label}[${idx}] = (${p[0]}, ${p[1]})`));
            };
            if (topo.type_names && topo.type_names.length) {
                console.log('[dumpXPDBTopology] type_names:');
                topo.type_names.forEach((name, idx) => console.log(`  type_names[${idx}] = ${name}`));
            }
            logPairs('[dumpXPDBTopology] bonds_primary', topo.bonds_primary);
            logPairs('[dumpXPDBTopology] bonds_angle', topo.bonds_angle);
            logPairs('[dumpXPDBTopology] bonds_pi', topo.bonds_pi);
            logPairs('[dumpXPDBTopology] bonds_epair', topo.bonds_epair);
        } catch (e) {
            console.error('[dumpXPDBTopology] failed', e);
        }
    }

    async dumpXPDBBuffers() {
        const ready = await this.ensureXPDBReady(true);
        if (!ready) return;
        console.log('[dumpXPDBBuffers] n_real=', this.xpdbNReal, 'n_all=', this.xpdbNAll, 'xpdbParams=', this.xpdbParams);
        try {
            await XPDB_WebGPU.dumpTypedState(this.xpdb, { label: 'XPDB_GUI', fixed: 6 });
        } catch (e) {
            console.error('[dumpXPDBBuffers] failed', e);
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
                this.xpdb.step({
                    dt: this.xpdbParams.dt,
                    iterations: this.xpdbParams.iterations,
                    k_coll: this.xpdbParams.k_coll,
                    omega: this.xpdbParams.omega,
                    momentum_beta: this.xpdbParams.momentum_beta,
                    maxRadius,
                    coll_scale: 2.0,
                    bbox_scale: 2.0,
                });

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
                this.xpdb.step({
                    dt: this.xpdbParams.dt,
                    iterations: this.xpdbParams.iterations,
                    k_coll: kColl,
                    omega: this.xpdbParams.omega,
                    momentum_beta: this.xpdbParams.momentum_beta,
                    maxRadius,
                    coll_scale: 2.0,
                    bbox_scale: 2.0,
                });
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
                    this.xpdb.step({
                        dt: this.xpdbParams.dt,
                        iterations: this.xpdbParams.iterations,
                        k_coll: kColl,
                        omega: this.xpdbParams.omega,
                        momentum_beta: this.xpdbParams.momentum_beta,
                        maxRadius,
                        coll_scale: 2.0,
                        bbox_scale: 2.0,
                    });
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

            const n = (this.xpdbNReal > 0) ? (this.xpdbNReal | 0) : (this.system.nAtoms | 0);
            for (let i = 0; i < n; i++) {
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
                    // Keep CPU solver in sync only for real atoms (avoid mismatch with dummy topology in GUI system)
                    // For now we disable this sync to prevent mol/cpu atom-count mismatch spam.
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
        // Any topology/atom-count change invalidates cached XPDB topo built on previous molecule.
        this.xpdbDirty = true;
        this._xpdbTopoLast = null;
        if (this.xpdbTopoRenderer) this.xpdbTopoRenderer.clear();
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