import { Vec3 } from '../../common_js/Vec3.js';

const THREE = window.THREE || globalThis.THREE;
if (!THREE) {
    throw new Error('THREE.js not found on window. Include three.js before ProjectiveDynamics.js');
}

// --- Shader Loading Helper ---

async function loadShader(path) {
    try {
        const response = await fetch(new URL(path, import.meta.url));
        if (!response.ok) throw new Error(`Failed to load shader: ${path}`);
        let text = await response.text();
        text = text.replace(/#version\s+300\s+es\s*/gi, '');
        return text.trimStart();
    } catch (e) {
        console.error('loadShader error:', e);
        throw new Error(`loadShader failed for ${path}: ${e.message}`);
    }
}

// --- Shader Source Codes ---

const vertexShaderScreen = `
void main() {
    gl_Position = vec4(position, 1.0);
}
`;

// Pass 1: Prediction (Inertia)
let fragmentShaderPredict = null;
// Pass 2: Solver (Jacobi + Heavy Ball)
let fragmentShaderSolver = null;
// Pass 3: Update Velocities
let fragmentShaderUpdate = null;

export class PDSimulation {
    constructor(renderer, atomCount, maxBonds = 16) {
        this.renderer = renderer;
        this.N = atomCount;
        this.maxBonds = maxBonds;
        this.iterationCount = 10;
        this.omega = 0.85; // Heavy ball parameter
        this.initialized = false;

        // Debug helpers
        this.debugTextureName = null;
        this.debugScene = null;
        this.debugCamera = null;
        this.debugQuad = null;
        this.debugBufferSize = new THREE.Vector2();

        // Scratch targets
        this.readbackTargets = new Map();

        // Debug data for angle-derived constraints (transient, rebuilt on setTopology)
        this.debugAngleBonds = []; // [{a, c, center, l0, k, accumulated}]

        // Optional callback invoked whenever debugAngleBonds is rebuilt
        this.onAngleDebugUpdate = null;
    }

    disposeResources() {
        const maybeDispose = (obj) => { if (obj && typeof obj.dispose === 'function') obj.dispose(); };
        maybeDispose(this.texPosInitial);
        maybeDispose(this.texVelInitial);
        maybeDispose(this.texBonds);
        maybeDispose(this.fboVel);
        maybeDispose(this.fboReadback);
        maybeDispose(this.fboSolverA);
        maybeDispose(this.fboSolverB);
        this.readbackTargets.forEach(rt => maybeDispose(rt));
        this.readbackTargets.clear();
    }

    resize(atomCount, maxBonds = this.maxBonds) {
        console.log(`[PD] resize: N ${this.N} -> ${atomCount}, maxBonds ${this.maxBonds} -> ${maxBonds}`);
        this.disposeResources();
        this.N = atomCount;
        this.maxBonds = maxBonds;
        this.currentPosFBO = null;
        this.readbackTargets.clear();
        this.initTextures();
        this.initFBOs();
        this.initMaterials();
    }

    async init() {
        if (this.initialized) return;

        // Load shaders
        fragmentShaderPredict = await loadShader('../../common_resources/shaders/MD_predictor.glsl');
        fragmentShaderSolver = await loadShader('../../common_resources/shaders/ProjectiveDynamics_2D.glsl');
        fragmentShaderUpdate = await loadShader('../../common_resources/shaders/MD_move.glsl');

        // 1. Setup Camera for Fullscreen Quads
        this.camera = new THREE.OrthographicCamera(-1, 1, 1, -1, 0, 1);
        this.scene = new THREE.Scene();
        this.quadGeo = new THREE.PlaneGeometry(2, 2);

        // 2. Initialize Data Textures
        this.initTextures();

        // 3. Initialize Render Targets (FBOs)
        this.initFBOs();

        // 4. Initialize Materials
        this.initMaterials();

        this.initialized = true;
    }

    createFloatTexture(width, height, data) {
        const tex = new THREE.DataTexture(
            data, width, height, 
            THREE.RGBAFormat, THREE.FloatType
        );
        tex.minFilter = THREE.NearestFilter;
        tex.magFilter = THREE.NearestFilter;
        tex.needsUpdate = true;
        return tex;
    }

    createMRT(width, height, count = 2, options = {}) {
        const params = {
            minFilter: options.minFilter || THREE.NearestFilter,
            magFilter: options.magFilter || THREE.NearestFilter,
            type: options.type || THREE.FloatType,
            format: options.format || THREE.RGBAFormat,
            wrapS: options.wrapS || THREE.ClampToEdgeWrapping,
            wrapT: options.wrapT || THREE.ClampToEdgeWrapping,
            depthBuffer: options.depthBuffer || false,
            stencilBuffer: options.stencilBuffer || false,
            count: count
        };
        const rt = new THREE.WebGLRenderTarget(width, height, params);
        const attachments = rt.textures || rt.texture;
        if (Array.isArray(attachments)) {
            attachments.forEach((tex) => {
                if (!tex) return;
                tex.minFilter = params.minFilter;
                tex.magFilter = params.magFilter;
                tex.type = params.type;
                tex.format = params.format;
                tex.wrapS = params.wrapS;
                tex.wrapT = params.wrapT;
                tex.needsUpdate = true;
            });
        } else if (attachments) {
            attachments.minFilter = params.minFilter;
            attachments.magFilter = params.magFilter;
            attachments.type = params.type;
            attachments.format = params.format;
            attachments.wrapS = params.wrapS;
            attachments.wrapT = params.wrapT;
            attachments.needsUpdate = true;
        }
        return rt;
    }

    getRTTexture(rt, index = 0) {
        if (!rt) return null;
        if (rt.textures) return rt.textures[index] || null;
        if (Array.isArray(rt.texture)) return rt.texture[index] || null;
        return index === 0 ? (rt.texture || null) : null;
    }

    initTextures() {
        // Positions & Velocities (1D: Width=N, Height=1)
        const size1D = this.N * 4;
        this.dataPos = new Float32Array(size1D);
        this.dataVel = new Float32Array(size1D);
        this.dataBonds = new Float32Array(this.N * this.maxBonds * 4); // 2D: Width=16, Height=N

        // --- Fill with Dummy Data (Grid) ---
        // (User should overwrite these arrays before first render if needed)
        // ... see main example below for filling logic ...

        this.texPosInitial = this.createFloatTexture(this.N, 1, this.dataPos);
        this.texVelInitial = this.createFloatTexture(this.N, 1, this.dataVel);
        
        // Note: Bond texture width is 16, height is N
        this.texBonds = this.createFloatTexture(this.maxBonds, this.N, this.dataBonds);
    }

    initFBOs() {
        const options = {
            minFilter: THREE.NearestFilter,
            magFilter: THREE.NearestFilter,
            type: THREE.FloatType,
            format: THREE.RGBAFormat,
            stencilBuffer: false,
            depthBuffer: false,
        };

        // Velocity State
        this.fboVel = new THREE.WebGLRenderTarget(this.N, 1, options);

        // Readback target for CPU transfers
        this.fboReadback = new THREE.WebGLRenderTarget(this.N, 1, options);

        // Ping-Pong Buffers for Solver (Must be MRT)
        // [0] = Position, [1] = Momentum
        this.fboSolverA = this.createMRT(this.N, 1, 2, options);
        this.fboSolverB = this.createMRT(this.N, 1, 2, options);

        // Helper to initialize FBOs with data
        // We render the initial textures into the FBOs once at startup
        // (Implementation omitted for brevity, but essentially a simple copy pass)
    }

    initMaterials() {
        // 1. Prediction Material (kept for compatibility if needed)
        this.matPredict = new THREE.ShaderMaterial({
            vertexShader: vertexShaderScreen,
            glslVersion: THREE.GLSL3,
            depthTest: false,
            depthWrite: false,
            fragmentShader: fragmentShaderPredict,
            uniforms: {
                tPos: { value: null },
                tVel: { value: null },
                dt: { value: 0.016 },
                texWidth: { value: this.N }
            }
        });

        // 2. Solver Material
        this.matSolver = new THREE.ShaderMaterial({
            vertexShader: vertexShaderScreen,
            glslVersion: THREE.GLSL3,
            depthTest: false,
            depthWrite: false,
            fragmentShader: fragmentShaderSolver,
            uniforms: {
                tPos: { value: null },
                tMomentum: { value: null },
                tBonds: { value: this.texBonds },
                dt: { value: 0.016 },
                omega: { value: this.omega },
                texWidth: { value: this.N }
            }
        });

        // 3. Update Material
        this.matUpdate = new THREE.ShaderMaterial({
            vertexShader: vertexShaderScreen,
            glslVersion: THREE.GLSL3,
            depthTest: false,
            depthWrite: false,
            fragmentShader: fragmentShaderUpdate,
            uniforms: {
                tPosOld: { value: null },
                tPosNew: { value: null },
                dt: { value: 0.016 }
            }
        });
        
        // Copy pass helper (to init FBOs)
        this.matCopy = new THREE.ShaderMaterial({
            vertexShader: vertexShaderScreen,
            glslVersion: THREE.GLSL3,
            depthTest: false,
            depthWrite: false,
            fragmentShader: `
            uniform sampler2D tInput;
            layout(location=0) out vec4 outColor;
            void main() { outColor = texelFetch(tInput, ivec2(gl_FragCoord.xy), 0); }
            `,
            uniforms: { tInput: { value: null } }
        });

        // Debug view material (visualize any texture with NEAREST)
        this.matDebugView = new THREE.ShaderMaterial({
            vertexShader: vertexShaderScreen,
            glslVersion: THREE.GLSL3,
            depthTest: false,
            depthWrite: false,
            fragmentShader: `
            uniform sampler2D tTex;
            out vec4 outColor;
            void main() {
                ivec2 pix = ivec2(gl_FragCoord.xy);
                vec4 v = texelFetch(tTex, pix, 0);
                outColor = vec4(v.rgb, 1.0);
            }
            `,
            uniforms: { tTex: { value: null } }
        });
    }

    // Helper to run a full screen quad pass
    renderPass(material, target) {
        this.quadMesh = this.quadMesh || new THREE.Mesh(this.quadGeo, material);
        this.quadMesh.material = material;
        if(this.quadMesh.parent !== this.scene) this.scene.add(this.quadMesh);

        this.renderer.setRenderTarget(target);
        this.renderer.render(this.scene, this.camera);
        this.renderer.setRenderTarget(null);
    }
    
    // One-time init to copy JS arrays to GPU FBOs
    uploadInitialData() {
        const matInit = new THREE.ShaderMaterial({
            glslVersion: THREE.GLSL3,
            vertexShader: vertexShaderScreen,
            fragmentShader: `
            uniform sampler2D tInput;
            layout(location=0) out vec4 outColor0;
            layout(location=1) out vec4 outColor1;
            void main() {
                ivec2 uv = ivec2(gl_FragCoord.xy);
                outColor0 = texelFetch(tInput, uv, 0);
                outColor1 = vec4(0.0);
            }
            `,
            uniforms: { tInput: { value: this.texPosInitial } },
            depthTest: false,
            depthWrite: false
        });
        this.renderPass(matInit, this.fboSolverA);

        // Copy initial velocities to fboVel
        this.matCopy.uniforms.tInput.value = this.texVelInitial;
        this.renderPass(this.matCopy, this.fboVel);

        // Set current position FBO to A (contains initial positions)
        this.currentPosFBO = this.fboSolverA;
    }

    step(dt) {
        // Pointers
        let sourcePosTex = (this.currentPosFBO) ? this.getRTTexture(this.currentPosFBO, 0) : this.texPosInitial;
        let sourceMomTex = (this.currentPosFBO) ? this.getRTTexture(this.currentPosFBO, 1) : null; // Null on first frame
        const velTex     = this.fboVel.texture;

        // --- 1. SOLVER LOOP (start from predicted/current positions already in tPos) ---
        // Save the "Old" position (start of frame) for velocity update later
        const posAtStartOfFrame = sourcePosTex;

        // Current Input/Output Setup
        // On first frame, sourceMomTex is null, shader needs to handle or we init FBO to black.
        // (Assuming initFBOs cleared them to 0, which is standard).
        
        let readFBO  = (this.currentPosFBO === this.fboSolverA) ? this.fboSolverA : this.fboSolverB;
        let writeFBO = (this.currentPosFBO === this.fboSolverA) ? this.fboSolverB : this.fboSolverA;
        
        // If first frame, readFBO might be empty/null logic, let's force write to A
        if (!this.currentPosFBO) {
             readFBO = this.fboSolverB; // Assume B is empty zero
             writeFBO = this.fboSolverA;
             // Bind initial texture as "Previous Pos" for the very first solver step
             // Actually, the solver shader takes tPos.
        }

        this.matSolver.uniforms.dt.value = dt;
        this.matSolver.uniforms.omega.value = this.omega;
        this.matSolver.uniforms.texWidth.value = this.N;
        
        for (let i = 0; i < this.iterationCount; i++) {
            // Inputs
            // tPos: P_k (from read buffer)
            // tMomentum: D_k (from read buffer)
            if (i === 0) {
                 this.matSolver.uniforms.tPos.value = posAtStartOfFrame;
                 this.matSolver.uniforms.tMomentum.value = this.getRTTexture(readFBO, 1);
            } else {
                 this.matSolver.uniforms.tPos.value = this.getRTTexture(readFBO, 0);
                 this.matSolver.uniforms.tMomentum.value = this.getRTTexture(readFBO, 1);
            }

            // Render
            this.renderPass(this.matSolver, writeFBO);

            // Swap
            let temp = readFBO; readFBO = writeFBO; writeFBO = temp;
        }
        
        // After loop, 'readFBO' holds the final P_{new} and D_{new}
        this.currentPosFBO = readFBO;

        // --- 3. VELOCITY UPDATE ---
        this.matUpdate.uniforms.tPosOld.value = posAtStartOfFrame;
        this.matUpdate.uniforms.tPosNew.value = this.getRTTexture(this.currentPosFBO, 0);
        this.matUpdate.uniforms.dt.value = dt;
        this.renderPass(this.matUpdate, this.fboVel);
    }
    
    // Get the texture to display
    getOutputTexture() {
        if (this.currentPosFBO) return this.getRTTexture(this.currentPosFBO, 0);
        return this.texPosInitial;
    }

    setTopology(mol, mmParams, opts = {}) {
        if (!mol || !mmParams) throw new Error('setTopology: mol and mmParams required');

        const nAtoms = mol.nAtoms || mol.atoms?.length || 0;
        if (nAtoms !== this.N) {
            console.warn(`setTopology: atom count mismatch (mol=${nAtoms}, sim=${this.N})`);
        }

        const includeAngleConstraints = (opts.includeAngleConstraints !== undefined) ? !!opts.includeAngleConstraints : true;
        console.log(`[PD] setTopology: atoms=${nAtoms} maxBonds=${this.maxBonds} includeAngleConstraints=${includeAngleConstraints}`);

        // Clear bond texture
        this.dataBonds.fill(0);
        for (let i = 3; i < this.dataBonds.length; i += 4) {
            this.dataBonds[i] = -1.0; // mark slot inactive
        }

        // Helper to get atom type name
        const getAtomTypeName = (atom) => {
            if (atom.typeName) return atom.typeName;
            const at = mmParams.getAtomTypeForAtom(atom);
            return at ? at.name : '*';
        };

        // Helper to add bond entry to texture
        const addBondEntry = (i, j, l0, k) => {
            if (i < 0 || i >= this.N || j < 0 || j >= this.N) return -1;
            const row = i * this.maxBonds * 4;
            for (let slot = 0; slot < this.maxBonds; slot++) {
                const idx = row + slot * 4;
                if (this.dataBonds[idx + 3] < 0.0) {
                    this.dataBonds[idx] = j;
                    this.dataBonds[idx + 1] = l0;
                    this.dataBonds[idx + 2] = k;
                    this.dataBonds[idx + 3] = 1.0;
                    if (window.DEBUG_PD) {
                        console.debug(`[PD][CPU] bond atom=${i} slot=${slot} -> neighbor=${j} l0=${l0.toFixed(3)} k=${k.toFixed(3)}`);
                    }
                    return slot;
                }
            }
            console.warn(`[PD] Atom ${i} exceeds max bonds (${this.maxBonds})`);
            return -1;
        };

        // Helper to find existing bond entry
        const findBondEntry = (i, j) => {
            const row = i * this.maxBonds * 4;
            for (let slot = 0; slot < this.maxBonds; slot++) {
                const idx = row + slot * 4;
                if (this.dataBonds[idx + 3] > 0.0 && this.dataBonds[idx] === j) {
                    return idx;
                }
            }
            return -1;
        };

        // Build local neighbor list on the fly (transient, not stored on molecule)
        const localNeighbors = new Array(nAtoms);
        for (let i = 0; i < nAtoms; i++) {
            localNeighbors[i] = [];
        }

        // 1. Process real bonds
        const bonds = mol.bonds || [];
        let realConstraintsCount = 0;
        let atomsAtCapacity = 0;

        console.log(`[PD] Processing ${bonds.length} real bonds...`);
        for (const bond of bonds) {
            if (typeof bond.ensureIndices === 'function') {
                bond.ensureIndices(mol);
            }
            const i = (Number.isInteger(bond.a)) ? bond.a : bond.i;
            const j = (Number.isInteger(bond.b)) ? bond.b : bond.j;
            if (!Number.isInteger(i) || !Number.isInteger(j)) {
                throw new Error(`setTopology: bond ${bond.id ?? 'unknown'} missing atom indices (a=${bond.a}, b=${bond.b})`);
            }
            if (i < 0 || i >= this.N || j < 0 || j >= this.N) {
                throw new Error(`setTopology: bond ${bond.id ?? 'unknown'} out-of-range atoms i=${i} j=${j} N=${this.N}`);
            }

            const atomA = mol.atoms[i];
            const atomB = mol.atoms[j];
            if (!atomA || !atomB) continue;

            const zA = atomA.Z || 0;
            const zB = atomB.Z || 0;
            const bondData = mmParams.getBondL0(zA, zB);

            let l0 = bondData ? bondData.l0 : 1.5;
            let k = bondData ? bondData.k : 100.0;

            const slotA = addBondEntry(i, j, l0, k);
            const slotB = addBondEntry(j, i, l0, k);

            if (slotA >= 0 && slotB >= 0) {
                realConstraintsCount++;
                localNeighbors[i].push(j);
                localNeighbors[j].push(i);
            } else {
                atomsAtCapacity++;
            }
        }

        console.log(`[PD] Real bonds processed: ${realConstraintsCount} constraints, ${atomsAtCapacity} atoms at capacity`);

        // 2. Process angles to create virtual bonds (on-the-fly from localNeighbors)
        let virtualConstraintsCount = 0;
        let angleTriplesProcessed = 0;
        let angleParamsFound = 0;
        let angleParamsMissing = 0;
        let angleBondsAdded = 0;
        let angleBondsAccumulated = 0;

        // Clear debug data for angle-derived constraints
        this.debugAngleBonds.length = 0;

        if (includeAngleConstraints) {
            console.log(`[PD] Processing angle-derived constraints...`);
            
            for (let b = 0; b < nAtoms; b++) {
                const neighs = localNeighbors[b];
                if (!neighs || neighs.length < 2) continue;

                const atomB = mol.atoms[b];
                const nameB = getAtomTypeName(atomB);

                // Enumerate all pairs of neighbors to form angles A-B-C
                for (let ni = 0; ni < neighs.length; ni++) {
                    for (let nj = ni + 1; nj < neighs.length; nj++) {
                        const a = neighs[ni];
                        const c = neighs[nj];
                        if (a < 0 || a >= this.N || c < 0 || c >= this.N) continue;

                        const atomA = mol.atoms[a];
                        const atomC = mol.atoms[c];
                        if (!atomA || !atomC) continue;

                        const nameA = getAtomTypeName(atomA);
                        const nameC = getAtomTypeName(atomC);

                        angleTriplesProcessed++;

                        const angleParams = mmParams.getAngleParams(nameA, nameB, nameC);
                        if (!angleParams) {
                            angleParamsMissing++;
                            continue;
                        }

                        angleParamsFound++;

                        // Get bond lengths for AB and BC
                        let lab = 1.5, lbc = 1.5;
                        const idxAB = findBondEntry(a, b);
                        const idxBC = findBondEntry(b, c);
                        if (idxAB >= 0) lab = this.dataBonds[idxAB + 1];
                        if (idxBC >= 0) lbc = this.dataBonds[idxBC + 1];

                        // Convert angle to distance constraint
                        const virtBond = mmParams.convertAngleToDistance(
                            lab, lbc, angleParams.ang0, angleParams.k
                        );

                        // Add virtual bond A-C (accumulate if exists)
                        const idxAC = findBondEntry(a, c);
                        if (idxAC >= 0) {
                            this.dataBonds[idxAC + 2] += virtBond.stiffness;
                            angleBondsAccumulated++;
                            // Store debug data for accumulated bond
                            this.debugAngleBonds.push({ a, c, center: b, l0: virtBond.restLength, k: virtBond.stiffness, accumulated: true });
                            if (window.DEBUG_PD) {
                                console.debug(`[PD][Angle] accumulated A-C bond: a=${a} c=${c} center=${b} l0=${virtBond.restLength.toFixed(3)} k=${virtBond.stiffness.toFixed(3)}`);
                            }
                        } else {
                            const slotA = addBondEntry(a, c, virtBond.restLength, virtBond.stiffness);
                            const slotC = addBondEntry(c, a, virtBond.restLength, virtBond.stiffness);
                            if (slotA >= 0 && slotC >= 0) {
                                angleBondsAdded++;
                                virtualConstraintsCount++;
                                // Store debug data for newly added bond
                                this.debugAngleBonds.push({ a, c, center: b, l0: virtBond.restLength, k: virtBond.stiffness, accumulated: false });
                                if (window.DEBUG_PD) {
                                    console.debug(`[PD][Angle] added A-C bond: a=${a} c=${c} center=${b} l0=${virtBond.restLength.toFixed(3)} k=${virtBond.stiffness.toFixed(3)}`);
                                }
                            } else {
                                atomsAtCapacity++;
                            }
                        }
                    }
                }
            }

            console.log(`[PD] Angle processing complete: triples=${angleTriplesProcessed} paramsFound=${angleParamsFound} paramsMissing=${angleParamsMissing} bondsAdded=${angleBondsAdded} bondsAccumulated=${angleBondsAccumulated} virtualConstraints=${virtualConstraintsCount}`);
        } else {
            console.log(`[PD] Angle-derived constraints SKIPPED (includeAngleConstraints=false)`);
        }

        // 3. Finalize and report summary
        const totalConstraints = realConstraintsCount + virtualConstraintsCount;
        console.log(`[PD] setTopology summary: atoms=${nAtoms} realConstraints=${realConstraintsCount} virtualConstraints=${virtualConstraintsCount} totalConstraints=${totalConstraints} maxBonds=${this.maxBonds} atomsAtCapacity=${atomsAtCapacity}`);

        if (typeof this.onAngleDebugUpdate === 'function') {
            this.onAngleDebugUpdate(this.debugAngleBonds);
        }

        // Update bond texture
        this.texBonds.image.data.set(this.dataBonds);
        this.texBonds.needsUpdate = true;

        if (window.DEBUG_PD) {
            this.debugDumpBondsCPU(4);
        }
    }

    setPositions(posArray) {
        if (!posArray || posArray.length < this.N * 3) return;
        for (let i = 0; i < this.N; i++) {
            this.dataPos[i * 4] = posArray[i * 3];
            this.dataPos[i * 4 + 1] = posArray[i * 3 + 1];
            this.dataPos[i * 4 + 2] = posArray[i * 3 + 2];
            this.dataPos[i * 4 + 3] = 1.0;
        }
        this.texPosInitial.image.data.set(this.dataPos);
        this.texPosInitial.needsUpdate = true;
        this.currentPosFBO = null; // force re-seed on next upload/step
    }

    setVelocities(velArray) {
        if (!velArray || velArray.length < this.N * 3) {
            this.dataVel.fill(0);
        } else {
            for (let i = 0; i < this.N; i++) {
                this.dataVel[i * 4] = velArray[i * 3];
                this.dataVel[i * 4 + 1] = velArray[i * 3 + 1];
                this.dataVel[i * 4 + 2] = velArray[i * 3 + 2];
                this.dataVel[i * 4 + 3] = 0.0;
            }
        }
        this.texVelInitial.image.data.set(this.dataVel);
        this.texVelInitial.needsUpdate = true;
        this.currentPosFBO = null; // force re-seed on next upload/step
    }

    readPositions() {
        const result = this.readTextureBuffer(this.getOutputTexture());
        return result ? result.buffer : null;
    }

    getReadbackTarget(width, height) {
        const key = `${width}x${height}`;
        if (!this.readbackTargets.has(key)) {
            const rt = new THREE.WebGLRenderTarget(width, height, {
                minFilter: THREE.NearestFilter,
                magFilter: THREE.NearestFilter,
                format: THREE.RGBAFormat,
                type: THREE.FloatType,
                depthBuffer: false,
                stencilBuffer: false
            });
            this.readbackTargets.set(key, rt);
        }
        return this.readbackTargets.get(key);
    }

    readTextureBuffer(source, width, height) {
        if (!source) return null;

        if (source.isWebGLRenderTarget) {
            const rt = source;
            const w = width || rt.width;
            const h = height || rt.height;
            const buffer = new Float32Array(w * h * 4);
            this.renderer.readRenderTargetPixels(rt, 0, 0, w, h, buffer);
            return { buffer, width: w, height: h };
        }

        const tex = source.isTexture ? source : null;
        if (!tex) return null;

        const img = tex.image || tex.source?.data || {};
        const w = width || img.width || this.N;
        const h = height || img.height || 1;

        const target = this.getReadbackTarget(w, h);
        this.matCopy.uniforms.tInput.value = tex;
        this.renderPass(this.matCopy, target);

        const buffer = new Float32Array(w * h * 4);
        this.renderer.readRenderTargetPixels(target, 0, 0, w, h, buffer);
        return { buffer, width: w, height: h };
    }

    debugReadTexture(label, textureOrTarget, options = {}) {
        const opts = (typeof options === 'number') ? { count: options } : options;
        const count = opts.count ?? 5;
        const mode = opts.mode || 'vec4';
        const result = this.readTextureBuffer(textureOrTarget, opts.width, opts.height);
        if (!result) {
            console.warn(`debugReadTexture: ${label} unavailable`);
            return null;
        }

        const { buffer, width, height } = result;
        console.groupCollapsed(`GPU Debug: ${label} (${width}x${height})`);
        if (mode === 'bonds') {
            const atoms = Math.min(count, height);
            for (let atom = 0; atom < atoms; atom++) {
                let rowStr = `atom ${atom}: `;
                for (let slot = 0; slot < width; slot++) {
                    const base = (atom * width + slot) * 4;
                    const neighbor = buffer[base];
                    const l0 = buffer[base + 1];
                    const stiff = buffer[base + 2];
                    const active = buffer[base + 3];
                    if (active < 0.0) continue;
                    rowStr += `[slot ${slot} -> j=${neighbor.toFixed(0)} l0=${l0.toFixed(3)} k=${stiff.toFixed(3)} act=${active.toFixed(1)}] `;
                }
                console.log(rowStr);
            }
        } else {
            const rows = Math.min(count, width * height);
            for (let i = 0; i < rows; i++) {
                const base = i * 4;
                console.log(`[${i}] ${buffer[base].toFixed(4)} ${buffer[base+1].toFixed(4)} ${buffer[base+2].toFixed(4)} | w=${buffer[base+3].toFixed(4)}`);
            }
        }
        console.groupEnd();
        return buffer;
    }

    debugDumpBondsCPU(maxAtoms = 4) {
        const atoms = Math.min(maxAtoms, this.N);
        console.group('[PD][CPU] Bond slots');
        for (let atom = 0; atom < atoms; atom++) {
            let rowStr = `atom ${atom}: `;
            const row = atom * this.maxBonds * 4;
            for (let slot = 0; slot < this.maxBonds; slot++) {
                const base = row + slot * 4;
                if (this.dataBonds[base + 3] < 0.0) continue;
                rowStr += `[slot ${slot} -> j=${this.dataBonds[base]} l0=${this.dataBonds[base+1].toFixed(3)} k=${this.dataBonds[base+2].toFixed(3)}] `;
            }
            console.log(rowStr);
        }
        console.groupEnd();
    }

    debugDumpBondsGPU(maxAtoms = 4) {
        this.debugReadTexture('texBonds', this.texBonds, {
            count: maxAtoms,
            width: this.maxBonds,
            height: this.N,
            mode: 'bonds'
        });
    }

    setDebugTexture(name) {
        this.debugTextureName = name || null;
    }

    renderDebugOverlay(renderer = this.renderer) {
        if (!renderer || !this.debugTextureName) return;

        let tex = null;
        switch (this.debugTextureName) {
            case 'positions': tex = this.getOutputTexture(); break;
            case 'bonds': tex = this.texBonds; break;
            case 'velocities': tex = this.fboVel.texture; break;
            case 'prediction': tex = this.fboSn.texture; break;
            default: return;
        }
        if (!tex) return;

        if (!this.debugScene) {
            this.debugScene = new THREE.Scene();
            this.debugCamera = new THREE.OrthographicCamera(-1, 1, 1, -1, 0, 1);
            this.debugQuad = new THREE.Mesh(new THREE.PlaneGeometry(2, 2), this.matDebugView);
            this.debugScene.add(this.debugQuad);
        }

        this.matDebugView.uniforms.tTex.value = tex;

        const canvasWidth = renderer.domElement.width;
        const canvasHeight = renderer.domElement.height;
        const size = Math.min(256, Math.floor(canvasHeight * 0.35));
        const x = canvasWidth - size - 16;
        const y = 16;

        const prevScissor = renderer.getScissor(new THREE.Vector4());
        const prevViewport = renderer.getViewport(new THREE.Vector4());
        const prevScissorTest = renderer.getScissorTest();

        renderer.setScissorTest(true);
        renderer.setScissor(x, y, size, size);
        renderer.setViewport(x, y, size, size);
        renderer.render(this.debugScene, this.debugCamera);
        renderer.setScissorTest(prevScissorTest);
        renderer.setScissor(prevScissor);
        renderer.setViewport(prevViewport);
    }

    renderDebugView(renderer) {
        this.renderDebugOverlay(renderer);
    }

    logConstraint(prefix, atomIndex, slot, j, l0, l, dl, k) {
        console.log(`[CPU][${prefix}] atom ${atomIndex} slot ${slot} j=${j} l0=${l0.toFixed(6)} l=${l.toFixed(6)} dl=${dl.toExponential(6)} k=${k.toFixed(6)}`);
    }

    computeLinearRHS(pos, masses, dt, debugAtom = -1) {
        const N = this.N;
        const maxBonds = this.maxBonds;
        const rhs = new Float32Array(N * 3);
        const diag = new Float32Array(N);
        const readVec = (arr, i, out) => out.set(arr[i * 3], arr[i * 3 + 1], arr[i * 3 + 2]);
        const writeVec = (arr, i, v) => { arr[i * 3] = v.x; arr[i * 3 + 1] = v.y; arr[i * 3 + 2] = v.z; };
        const fmtVec = (v) => v.toString();
        const tmpPi = new Vec3();
        const tmpPj = new Vec3();
        const tmpDij = new Vec3();
        const tmpBi = new Vec3();
        const dt2 = dt * dt;
        const logConstraint = this.logConstraint.bind(this);

        console.log('[CPU][Jacobi] Precompute RHS/diag START');
        for (let i = 0; i < N; i++) {
            const shouldDebug = (debugAtom === -1 || debugAtom === i);
            const pi = readVec(pos, i, tmpPi);

            const Ii = masses[i] / dt2;
            let Aii = Ii;
            tmpBi.setV(pi).mulScalar(Ii);

            const row = i * maxBonds * 4;
            for (let slot = 0; slot < maxBonds; slot++) {
                const base = row + slot * 4;
                if (this.dataBonds[base + 3] < 0.0) continue;

                const j = this.dataBonds[base];
                const l0 = this.dataBonds[base + 1];
                const k = this.dataBonds[base + 2];

                if (j < 0 || j >= N) continue;

                const pj = readVec(pos, j, tmpPj);
                const dij = tmpDij.setSub(pi, pj);
                const dist = dij.norm();
                if (dist < 1e-12) continue;
                const dl = dist - l0;

                // Match C++/fly algebra for linear Jacobi solve:
                // p_new = ( Ii*pi + sum_j[ k*(pj + (l0/l)*(pi-pj)) ] ) / ( Ii + sum_j k )
                // so we store rhs_i = Ii*pi + sum_j[ k*(l0/l)*(pi-pj) ] and later add sum_j k*pj.
                tmpBi.addMul(dij, k * (l0 / dist));
                Aii += k;

                this.logConstraint('JacobiRHS', i, slot, j, l0, dist, dl, k);
            }

            if (Aii < 1e-12) {
                console.error(`[CPU][Jacobi] Atom ${i} diagonal too small: ${Aii}`);
                throw new Error(`runCpuJacobi: Atom ${i} has zero diagonal`);
            }

            diag[i] = Aii;
            writeVec(rhs, i, tmpBi);

            if (shouldDebug) {
                console.log(`[CPU][Jacobi] precompute atom ${i} summary`);
                console.log(`    pi: ${fmtVec(pi)}`);
                console.log(`    rhs: ${fmtVec(tmpBi)} diag=${Aii.toFixed(6)}`);
            }
        }
        console.log('[CPU][Jacobi] Precompute RHS/diag END');
        return { rhs, diag };
    }

    solveJacobiLinear(pos, rhs, diag, iterations, debugAtom = -1) {
        const N = this.N;
        const maxBonds = this.maxBonds;
        const readVec = (arr, i, out) => out.set(arr[i * 3], arr[i * 3 + 1], arr[i * 3 + 2]);
        const writeVec = (arr, i, v) => { arr[i * 3] = v.x; arr[i * 3 + 1] = v.y; arr[i * 3 + 2] = v.z; };
        const fmtVec = (v) => v.toString();
        const tmpPi = new Vec3();
        const tmpPj = new Vec3();
        const tmpSum = new Vec3();
        const tmpNumerator = new Vec3();
        const tmpBi = new Vec3();
        const tmpDelta = new Vec3();

        for (let iter = 0; iter < iterations; iter++) {
            console.log(`[CPU][Jacobi] Iteration ${iter} START`);

            const p_old = new Float32Array(pos);
            const p_new = new Float32Array(N * 3);

            for (let i = 0; i < N; i++) {
                const shouldDebug = (debugAtom === -1 || debugAtom === i);
                tmpSum.set(0, 0, 0);

                const row = i * maxBonds * 4;
                for (let slot = 0; slot < maxBonds; slot++) {
                    const base = row + slot * 4;
                    if (this.dataBonds[base + 3] < 0.0) continue;
                    const j = this.dataBonds[base];
                    const k = this.dataBonds[base + 2];
                    if (j < 0 || j >= N) continue;
                    const pj = readVec(p_old, j, tmpPj);
                    tmpSum.addMul(pj, k);
                    if (shouldDebug) {
                        console.log(`[CPU][Jacobi]   sum_j slot ${slot} j=${j} k=${k.toFixed(6)} contrib=${fmtVec(tmpSum)}`);
                    }
                }

                const rhs_i = readVec(rhs, i, tmpBi);
                tmpNumerator.setV(rhs_i).add(tmpSum);

                const Aii = diag[i];
                const p_new_i = tmpPj.setV(tmpNumerator).mulScalar(1.0 / Aii);
                writeVec(p_new, i, p_new_i);

                if (shouldDebug) {
                    const piOld = readVec(p_old, i, tmpPi);
                    console.log(`[CPU][Jacobi] Atom ${i}`);
                    console.log(`    rhs: ${fmtVec(rhs_i)} sum_j: ${fmtVec(tmpSum)}`);
                    console.log(`    numerator: ${fmtVec(tmpNumerator)} diag=${Aii.toFixed(6)}`);
                    console.log(`    p_new: ${fmtVec(p_new_i)}`);
                    const delta = tmpDelta.setSub(p_new_i, piOld);
                    console.log(`    delta: ${fmtVec(delta)} |delta|=${delta.norm().toFixed(6)}`);
                }
            }

            pos.set(p_new);
            console.log(`[CPU][Jacobi] Iteration ${iter} complete`);
            console.log(`[CPU][Jacobi] Iteration ${iter} END`);
        }
        return pos;
    }

    solveGaussSeidelLinear(pos, rhs, diag, iterations, debugAtom = -1) {
        const N = this.N;
        const maxBonds = this.maxBonds;
        const readVec = (arr, i, out) => out.set(arr[i * 3], arr[i * 3 + 1], arr[i * 3 + 2]);
        const writeVec = (arr, i, v) => { arr[i * 3] = v.x; arr[i * 3 + 1] = v.y; arr[i * 3 + 2] = v.z; };
        const fmtVec = (v) => v.toString();
        const tmpPi = new Vec3();
        const tmpPj = new Vec3();
        const tmpSum = new Vec3();
        const tmpBi = new Vec3();
        const tmpDelta = new Vec3();

        for (let iter = 0; iter < iterations; iter++) {
            console.log(`[CPU][GS] Iteration ${iter} START`);
            for (let i = 0; i < N; i++) {
                const shouldDebug = (debugAtom === -1 || debugAtom === i);
                tmpSum.set(0, 0, 0);

                const row = i * maxBonds * 4;
                for (let slot = 0; slot < maxBonds; slot++) {
                    const base = row + slot * 4;
                    if (this.dataBonds[base + 3] < 0.0) continue;
                    const j = this.dataBonds[base];
                    const k = this.dataBonds[base + 2];
                    if (j < 0 || j >= N) continue;
                    const pj = readVec(pos, j, tmpPj);
                    tmpSum.addMul(pj, k);
                }

                const rhs_i = readVec(rhs, i, tmpBi);
                tmpSum.add(rhs_i);
                const Aii = diag[i];
                const piOld = readVec(pos, i, tmpPi);
                const piNew = tmpPj.setV(tmpSum).mulScalar(1.0 / Aii);
                writeVec(pos, i, piNew);

                if (shouldDebug) {
                    const delta = tmpDelta.setSub(piNew, piOld);
                    console.log(`[CPU][GS] Atom ${i}`);
                    console.log(`    rhs: ${fmtVec(rhs_i)} sum_j: ${fmtVec(tmpSum)}`);
                    console.log(`    diag=${Aii.toFixed(6)} p_new=${fmtVec(piNew)} delta=${fmtVec(delta)}`);
                }
            }
            console.log(`[CPU][GS] Iteration ${iter} END`);
        }
        return pos;
    }

    solveJacobiFly(pos, masses, dt, iterations, debugAtom = -1, momentum = null) {
        const N = this.N;
        const maxBonds = this.maxBonds;
        const readVec = (arr, i, out) => out.set(arr[i * 3], arr[i * 3 + 1], arr[i * 3 + 2]);
        const writeVec = (arr, i, v) => { arr[i * 3] = v.x; arr[i * 3 + 1] = v.y; arr[i * 3 + 2] = v.z; };
        const fmtVec = (v) => v.toString();
        const tmpPi = new Vec3();
        const tmpPj = new Vec3();
        const tmpDij = new Vec3();
        const tmpSum = new Vec3();
        const tmpDelta = new Vec3();
        const tmpMom = new Vec3();
        const tmpAcc = new Vec3();
        const dt2 = dt * dt;
        const logConstraint = this.logConstraint.bind(this);
        const omega = this.omega || 0.0;
        let momPrev = momentum || new Float32Array(N * 3);
        let momNew = new Float32Array(N * 3);

        for (let iter = 0; iter < iterations; iter++) {
            console.log(`[CPU][JacobiFly] Iteration ${iter} START`);
            const p_old = new Float32Array(pos);
            const p_new = new Float32Array(N * 3);

            for (let i = 0; i < N; i++) {
                const shouldDebug = (debugAtom === -1 || debugAtom === i);
                const pi = readVec(p_old, i, tmpPi);
                let Aii = masses[i] / dt2;
                tmpSum.setV(pi).mulScalar(Aii);

                const row = i * maxBonds * 4;
                for (let slot = 0; slot < maxBonds; slot++) {
                    const base = row + slot * 4;
                    if (this.dataBonds[base + 3] < 0.0) continue;
                    const j = this.dataBonds[base];
                    const l0 = this.dataBonds[base + 1];
                    const k = this.dataBonds[base + 2];
                    if (j < 0 || j >= N) continue;
                    const pj = readVec(p_old, j, tmpPj);
                    const dij = tmpDij.setSub(pi, pj);
                    const dist = dij.norm();
                    if (dist < 1e-12) continue;
                    const dl = dist - l0;

                    tmpSum.addMul(dij, k * (l0 / dist));
                    tmpSum.addMul(pj, k);
                    Aii += k;

                    logConstraint('JacobiFly', i, slot, j, l0, dist, dl, k);
                }

                const p_new_i = tmpPj.setV(tmpSum).mulScalar(1.0 / Aii); // Jacobi solution

                const momPrev_i = readVec(momPrev, i, tmpMom);
                const p_acc = tmpAcc.setV(p_new_i).addMul(momPrev_i, omega);
                const delta = tmpDelta.setSub(p_acc, pi);

                writeVec(p_new, i, p_acc);
                writeVec(momNew, i, delta);

                if (shouldDebug) {
                    console.log(`[CPU][JacobiFly] atom ${i} p_new=${fmtVec(p_acc)} delta=${fmtVec(delta)} momPrev=${fmtVec(momPrev_i)}`);
                }
            }

            pos.set(p_new);
            momPrev = momNew;
            momNew = new Float32Array(N * 3);
            console.log(`[CPU][JacobiFly] Iteration ${iter} END`);
        }

        if (momentum) momentum.set(momPrev);
        return pos;
    }

    solveGaussSeidelFly(pos, masses, dt, iterations, debugAtom = -1, momentum = null) {
        const N = this.N;
        const maxBonds = this.maxBonds;
        const readVec = (arr, i, out) => out.set(arr[i * 3], arr[i * 3 + 1], arr[i * 3 + 2]);
        const writeVec = (arr, i, v) => { arr[i * 3] = v.x; arr[i * 3 + 1] = v.y; arr[i * 3 + 2] = v.z; };
        const fmtVec = (v) => v.toString();
        const tmpPi = new Vec3();
        const tmpPj = new Vec3();
        const tmpDij = new Vec3();
        const tmpSum = new Vec3();
        const tmpDelta = new Vec3();
        const tmpMom = new Vec3();
        const tmpAcc = new Vec3();
        const dt2 = dt * dt;
        const logConstraint = this.logConstraint.bind(this);
        const omega = this.omega || 0.0;
        let momPrev = momentum || new Float32Array(N * 3);

        for (let iter = 0; iter < iterations; iter++) {
            console.log(`[CPU][GSFly] Iteration ${iter} START`);
            for (let i = 0; i < N; i++) {
                const shouldDebug = (debugAtom === -1 || debugAtom === i);
                const pi = readVec(pos, i, tmpPi);
                let Aii = masses[i] / dt2;
                tmpSum.setV(pi).mulScalar(Aii);

                const row = i * maxBonds * 4;
                for (let slot = 0; slot < maxBonds; slot++) {
                    const base = row + slot * 4;
                    if (this.dataBonds[base + 3] < 0.0) continue;
                    const j = this.dataBonds[base];
                    const l0 = this.dataBonds[base + 1];
                    const k = this.dataBonds[base + 2];
                    if (j < 0 || j >= N) continue;
                    const pj = readVec(pos, j, tmpPj);
                    const dij = tmpDij.setSub(pi, pj);
                    const dist = dij.norm();
                    if (dist < 1e-12) continue;
                    const dl = dist - l0;

                    tmpSum.addMul(dij, k * (l0 / dist));
                    tmpSum.addMul(pj, k);
                    Aii += k;

                    logConstraint('GSFly', i, slot, j, l0, dist, dl, k);
                }

                const piNew = tmpPj.setV(tmpSum).mulScalar(1.0 / Aii); // GS solution
                const momPrev_i = readVec(momPrev, i, tmpMom);
                const piAcc = tmpDij.setV(piNew).addMul(momPrev_i, omega);
                const delta = tmpDelta.setSub(piAcc, pi);

                writeVec(pos, i, piAcc);
                writeVec(momPrev, i, delta);

                if (shouldDebug) {
                    console.log(`[CPU][GSFly] atom ${i} p_new=${fmtVec(piAcc)} delta=${fmtVec(delta)} momPrev=${fmtVec(momPrev_i)}`);
                }
            }
            console.log(`[CPU][GSFly] Iteration ${iter} END`);
        }

        return pos;
    }

    runCpuJacobi({ positions, dt = 0.1, iterations = 1, debugAtom = -1, masses = null, mode = 'jacobi_fly' } = {}) {
        if (!positions || positions.length < this.N * 3) {
            throw new Error('runCpuJacobi: positions array too small');
        }

        const N = this.N;
        const maxBonds = this.maxBonds;

        if (!masses) {
            masses = new Float32Array(N);
            masses.fill(1.0);
        }

        console.log('[CPU][Jacobi] Starting solver');
        console.log(`[CPU][Jacobi] N=${N} dt=${dt} iterations=${iterations} debugAtom=${debugAtom} mode=${mode}`);

        const pos = new Float32Array(positions);

        let result;
        const momBuf = new Float32Array(N * 3);
        if (mode === 'jacobi_lin' || mode === 'gs_lin') {
            const { rhs, diag } = this.computeLinearRHS(pos, masses, dt, debugAtom);
            if (mode === 'jacobi_lin') {
                result = this.solveJacobiLinear(pos, rhs, diag, iterations, debugAtom);
            } else {
                result = this.solveGaussSeidelLinear(pos, rhs, diag, iterations, debugAtom);
            }
        } else if (mode === 'jacobi_fly') {
            result = this.solveJacobiFly(pos, masses, dt, iterations, debugAtom, momBuf);
        } else if (mode === 'gs_fly') {
            result = this.solveGaussSeidelFly(pos, masses, dt, iterations, debugAtom, momBuf);
        } else {
            throw new Error(`runCpuJacobi: unknown mode "${mode}"`);
        }

        const readVecFinal = (arr, i, out) => out.set(arr[i * 3], arr[i * 3 + 1], arr[i * 3 + 2]);
        const tmpVec = new Vec3();
        console.log('[CPU][Jacobi] Final positions:');
        for (let i = 0; i < N; i++) {
            console.log(`[CPU][Jacobi]   atom ${i}: ${readVecFinal(result, i, tmpVec).toString()}`);
        }
        console.groupEnd();
        return result;
    }
}