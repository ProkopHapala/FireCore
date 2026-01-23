import * as THREE from 'three';
import { WebGPURenderer } from 'three/webgpu';
import { OrbitControls } from 'three/addons/controls/OrbitControls.js';
import { TransformControls } from 'three/addons/controls/TransformControls.js';

import { logger } from '../common_js/Logger.js';
import { EditableMolecule } from './EditableMolecule.js';
import { MMParams } from './MMParams.js';
import { MoleculeRenderer, PackedMolecule } from './MoleculeRenderer.js';
import { RawWebGPUAtomsRenderer } from './RawWebGPUAtomsRenderer.js';
import { GUI } from './GUI.js';
import { Editor } from './Editor.js';
import { ShortcutManager } from './ShortcutManager.js';
import { Vec3 } from '../common_js/Vec3.js';
import { ScriptRunner } from './ScriptRunner.js';

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
        this.system = new EditableMolecule();
        this.controls = null;
        this._rafPending = false;
        this.continuousRender = false;
    }

    requestRender() {
        if (this._rafPending) return;
        this._rafPending = true;
        requestAnimationFrame(async () => {
            this._rafPending = false;
            if (this.renderers) {
                Object.values(this.renderers).forEach(r => r.update());
            }

            if (this.useRawWebGPU && this.raw) {
                if (this.controls && this.controls.update) this.controls.update();
                this.camera.updateMatrixWorld(true);
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
                        const isSel = sel && sel.has(i);
                        col[i * 3] = isSel ? 1.0 : c[0];
                        col[i * 3 + 1] = isSel ? 1.0 : c[1];
                        col[i * 3 + 2] = isSel ? 0.0 : c[2];
                        const r = this.mmParams ? (this.mmParams.getRadius(t) * 0.4) : 0.5;
                        rad[i] = isSel ? r * 1.3 : r;
                    }
                    this.raw.updateAtoms(pos, col, rad, n);
                } else {
                    this.raw.nAtoms = 0;
                }

                this.raw.render();

                // Overlay pass (gizmo/selection via Three WebGL)
                if (this.overlayRenderer) {
                    if (this.overlayScene) {
                        this.overlayRenderer.autoClear = true;
                        this.overlayRenderer.setClearColor(0x000000, 0.0);
                        this.overlayRenderer.clear(true, true, true);
                        this.overlayRenderer.render(this.overlayScene, this.camera);
                    }
                }
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

            // Provide minimal renderer-like object for Editor/controls
            this.renderer = {
                domElement: canvas,
                setSize: (w, h) => this.raw.resize(w, h),
                setPixelRatio: () => {},
                setAnimationLoop: () => {}
            };

            // Overlay canvas for gizmo/selection (Three WebGL)
            this.overlayCanvas = document.createElement('canvas');
            this.overlayCanvas.style.position = 'absolute';
            this.overlayCanvas.style.top = '0';
            this.overlayCanvas.style.left = '0';
            this.overlayCanvas.style.pointerEvents = 'none';
            this.overlayCanvas.style.zIndex = '1';
            this.overlayCanvas.style.backgroundColor = 'transparent';
            this.container.appendChild(this.overlayCanvas);

            this.overlayRenderer = new THREE.WebGLRenderer({ canvas: this.overlayCanvas, antialias: true, alpha: true });
            this.overlayRenderer.autoClear = true;
            this.overlayRenderer.setClearColor(0x000000, 0.0);
            this.overlayRenderer.setSize(width, height, false);

            // Dedicated overlay scene (do NOT render node-material meshes here)
            this.overlayScene = new THREE.Scene();

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
}

window.onload = () => {
    const app = new MolGUIApp();
    window.app = app;
    app.init();
};