export class XPDB_WebGPU {
    constructor(numAtoms, groupSize = 64) {
        this.numAtoms = numAtoms;
        this.groupSize = groupSize;
        this.numGroups = Math.ceil(numAtoms / groupSize);
        this.maxGhosts = 128;
        this.nMaxBonded = 16;
        
        this.device = null;
        this.pipelines = {};
        this.bindGroup = null;
        this.buffers = {};
    }

    static async dumpTypedState(sim, { label = 'XPDB', writeText = null, fixed = 6 } = {}) {
        const dbg = await import('./debugBuffers.js');
        const n = sim.numAtoms | 0;
        const ng = sim.numGroups | 0;
        const nb = sim.nMaxBonded | 0;
        const mg = sim.maxGhosts | 0;
        const state = await sim.readBuffersAsTyped();
        const sections = [];
        sections.push([`# ${label} numAtoms=${n} numGroups=${ng} nMaxBonded=${nb} maxGhosts=${mg}`]);
        sections.push(dbg.dumpVec4BufferLines('pos', state.posData, n, { stride: 4, cols: 3, fixed }));
        sections.push(dbg.dumpVec4BufferLines('pred_pos', state.predPosData, n, { stride: 4, cols: 3, fixed }));
        sections.push(dbg.dumpAtomParamsLines('atom_params', state.atomParamsData, n, { fixed }));
        sections.push(dbg.dumpVec4BufferLines('bboxes', state.bboxesData, ng * 2, { stride: 4, cols: 3, fixed }));
        sections.push(dbg.dumpGhostPackedLines('ghost_packed', state.ghostData, ng, mg));
        sections.push(dbg.dumpBondIndicesLines('bond_idx_global', state.bondGlobalData, n, nb));
        sections.push(dbg.dumpBondIndicesLines('bond_idx_local', state.bondLocalData, n, nb));
        sections.push(dbg.dumpBondLenStiffLines('bond_len_stiff', state.bondLenStiffData, n, nb, { fixed }));
        if (state.forceBondData) sections.push(dbg.dumpVec3Forces('force_bond', state.forceBondData, n, { fixed }));
        if (state.forceCollData) sections.push(dbg.dumpVec3Forces('force_coll', state.forceCollData, n, { fixed }));
        const text = dbg.joinSections(sections);
        if (typeof writeText === 'function') await writeText(text);
        else console.log(text);
        return { text, state };
    }

    async init() {
        if (!navigator.gpu) throw new Error("WebGPU not supported");
        const adapter = await navigator.gpu.requestAdapter();
        if (!adapter) throw new Error('XPDB_WebGPU.init: navigator.gpu.requestAdapter() returned null');
        const neededStorageBuffers = 10;
        const supported = adapter.limits?.maxStorageBuffersPerShaderStage;
        if (typeof supported === 'number' && supported < neededStorageBuffers) {
            throw new Error(`XPDB_WebGPU.init: adapter maxStorageBuffersPerShaderStage=${supported} < required ${neededStorageBuffers}`);
        }
        this.device = await adapter.requestDevice({
            requiredLimits: {
                maxStorageBuffersPerShaderStage: Math.max(neededStorageBuffers, supported || neededStorageBuffers)
            }
        });

        // Load Shader
        const shaderModule = this.device.createShaderModule({
            code: await fetch('./xpdb.wgsl').then(r => r.text())
        });

        // Explicit bind group layout with packed buffers to stay within storage limits
        this.bindGroupLayout = this.device.createBindGroupLayout({
            entries: [
                { binding: 0, visibility: GPUShaderStage.COMPUTE, buffer: { type: 'uniform' } },
                { binding: 1, visibility: GPUShaderStage.COMPUTE, buffer: { type: 'storage' } },
                { binding: 2, visibility: GPUShaderStage.COMPUTE, buffer: { type: 'storage' } },
                { binding: 3, visibility: GPUShaderStage.COMPUTE, buffer: { type: 'read-only-storage' } }, // atomParams
                { binding: 4, visibility: GPUShaderStage.COMPUTE, buffer: { type: 'storage' } }, // packed bbox min/max
                { binding: 5, visibility: GPUShaderStage.COMPUTE, buffer: { type: 'storage' } }, // packed ghosts (indices + counts)
                { binding: 6, visibility: GPUShaderStage.COMPUTE, buffer: { type: 'read-only-storage' } }, // bondIdxGlobal
                { binding: 7, visibility: GPUShaderStage.COMPUTE, buffer: { type: 'storage' } }, // bondIdxLocal
                { binding: 8, visibility: GPUShaderStage.COMPUTE, buffer: { type: 'read-only-storage' } }, // bondLenStiff packed
                { binding: 9, visibility: GPUShaderStage.COMPUTE, buffer: { type: 'storage' } }, // debug force bond
                { binding: 10, visibility: GPUShaderStage.COMPUTE, buffer: { type: 'storage' } }, // debug force coll
            ],
        });

        const pipelineLayout = this.device.createPipelineLayout({ bindGroupLayouts: [this.bindGroupLayout] });

        // Create Pipelines
        const createPipeline = (entryPoint) => {
            return this.device.createComputePipeline({
                layout: pipelineLayout,
                compute: { module: shaderModule, entryPoint }
            });
        };

        this.pipelines = {
            updateBBoxes: createPipeline('update_bboxes'),
            buildTopology: createPipeline('build_local_topology'),
            solve: createPipeline('solve_cluster_jacobi')
        };

        this._createBuffers();
    }

    _createBuffers() {
        const createBuf = (size, usage = GPUBufferUsage.STORAGE | GPUBufferUsage.COPY_DST | GPUBufferUsage.COPY_SRC) => {
            return this.device.createBuffer({ size, usage });
        };

        const f32 = 4;
        const vec4 = 16;
        const i32 = 4;
        const align = (n) => Math.ceil(n / 4) * 4; // Padding helper if needed

        // 0: Params Uniform
        // Struct size must be 16-byte aligned. 12 floats/ints = 48 bytes.
        this.buffers.params = this.device.createBuffer({
            size: 64, // Rounded up
            usage: GPUBufferUsage.UNIFORM | GPUBufferUsage.COPY_DST
        });

        // Atom Data
        this.buffers.pos = createBuf(this.numAtoms * vec4);
        this.buffers.predPos = createBuf(this.numAtoms * vec4);
        this.buffers.atomParams = createBuf(this.numAtoms * vec4);

        // BBoxes packed: [min,max] per group -> two vec4 per group
        this.buffers.bboxes = createBuf(this.numGroups * 2 * vec4);

        // Ghost Data packed: indices + count per group (count stored at offset maxGhosts)
        this.buffers.ghostPacked = createBuf(this.numGroups * (this.maxGhosts + 1) * i32);

        // Bonds
        this.buffers.bondIdxGlobal = createBuf(this.numAtoms * this.nMaxBonded * i32);
        this.buffers.bondIdxLocal = createBuf(this.numAtoms * this.nMaxBonded * i32);
        // Pack bondLen and bondStiff into vec2<f32>
        this.buffers.bondLenStiff = createBuf(this.numAtoms * this.nMaxBonded * (f32 * 2));

        this.buffers.forceBond = createBuf(this.numAtoms * vec4);
        this.buffers.forceColl = createBuf(this.numAtoms * vec4);

        // Bind Group
        const entries = [
            { binding: 0, resource: { buffer: this.buffers.params } },
            { binding: 1, resource: { buffer: this.buffers.pos } },
            { binding: 2, resource: { buffer: this.buffers.predPos } },
            { binding: 3, resource: { buffer: this.buffers.atomParams } },
            { binding: 4, resource: { buffer: this.buffers.bboxes } },
            { binding: 5, resource: { buffer: this.buffers.ghostPacked } },
            { binding: 6, resource: { buffer: this.buffers.bondIdxGlobal } },
            { binding: 7, resource: { buffer: this.buffers.bondIdxLocal } },
            { binding: 8, resource: { buffer: this.buffers.bondLenStiff } },
            { binding: 9, resource: { buffer: this.buffers.forceBond } },
            { binding: 10, resource: { buffer: this.buffers.forceColl } },
        ];

        this.bindGroup = this.device.createBindGroup({
            layout: this.bindGroupLayout,
            entries,
        });
    }

    uploadAtoms(posArray, radiusArray, massArray) {
        // posArray: Float32Array [x,y,z, x,y,z...]
        // Expand to vec4 [x,y,z,0]
        const posData = new Float32Array(this.numAtoms * 4);
        const paramsData = new Float32Array(this.numAtoms * 4);

        for (let i = 0; i < this.numAtoms; i++) {
            posData[i*4 + 0] = posArray[i*3 + 0];
            posData[i*4 + 1] = posArray[i*3 + 1];
            posData[i*4 + 2] = posArray[i*3 + 2];
            posData[i*4 + 3] = 0.0;

            paramsData[i*4 + 0] = radiusArray[i]; // radius
            paramsData[i*4 + 3] = massArray[i];   // mass
        }

        this.device.queue.writeBuffer(this.buffers.pos, 0, posData);
        this.device.queue.writeBuffer(this.buffers.predPos, 0, posData); // Init pred = curr
        this.device.queue.writeBuffer(this.buffers.atomParams, 0, paramsData);
    }

    uploadAtomsFromMolecule(mol, mmParams) {
        if (!mol || !mmParams) throw new Error('uploadAtomsFromMolecule: mol and mmParams required');
        const nAtoms = mol.nAtoms || mol.atoms?.length || 0;
        if (nAtoms !== this.numAtoms) throw new Error(`uploadAtomsFromMolecule: atom count mismatch (mol=${nAtoms}, sim=${this.numAtoms})`);

        const posData = new Float32Array(this.numAtoms * 4);
        const paramsData = new Float32Array(this.numAtoms * 4);

        for (let i = 0; i < this.numAtoms; i++) {
            const atom = mol.atoms[i];
            if (!atom) throw new Error(`uploadAtomsFromMolecule: missing atom at index ${i}`);

            posData[i*4 + 0] = atom.pos.x;
            posData[i*4 + 1] = atom.pos.y;
            posData[i*4 + 2] = atom.pos.z;
            posData[i*4 + 3] = 0.0;

            const at = mmParams.getAtomTypeForAtom(atom);
            const radius = (at && at.RvdW > 0) ? at.RvdW : 1.0;
            const mass = (atom.Z || 1) * 1.0;

            paramsData[i*4 + 0] = radius;
            paramsData[i*4 + 3] = mass;
        }

        this.device.queue.writeBuffer(this.buffers.pos, 0, posData);
        this.device.queue.writeBuffer(this.buffers.predPos, 0, posData);
        this.device.queue.writeBuffer(this.buffers.atomParams, 0, paramsData);
    }

    uploadPackedAtoms(packed) {
        if (!packed) throw new Error('uploadPackedAtoms: packed is null/undefined');
        if (!packed.pos4 || !packed.params4) throw new Error('uploadPackedAtoms: packed must provide pos4 and params4');
        if ((packed.pos4.length | 0) !== (this.numAtoms * 4)) throw new Error(`uploadPackedAtoms: pos4 length mismatch pos4=${packed.pos4.length} expected=${this.numAtoms * 4}`);
        if ((packed.params4.length | 0) !== (this.numAtoms * 4)) throw new Error(`uploadPackedAtoms: params4 length mismatch params4=${packed.params4.length} expected=${this.numAtoms * 4}`);
        this.device.queue.writeBuffer(this.buffers.pos, 0, packed.pos4);
        this.device.queue.writeBuffer(this.buffers.predPos, 0, packed.pos4);
        this.device.queue.writeBuffer(this.buffers.atomParams, 0, packed.params4);
    }

    uploadPackedTopology(packed) {
        if (!packed) throw new Error('uploadPackedTopology: packed is null/undefined');
        if (!packed.bondIndices || !packed.bondLenStiff) throw new Error('uploadPackedTopology: packed must provide bondIndices and bondLenStiff');
        this.device.queue.writeBuffer(this.buffers.bondIdxGlobal, 0, packed.bondIndices);
        this.device.queue.writeBuffer(this.buffers.bondLenStiff, 0, packed.bondLenStiff);
    }

    // DEPRECATED: this is legacy GUI-only path (bondsAdj + defaults). Use buildXPDBInputsFromMol/buildXPDBInputsFromXYZArgs + uploadPackedTopology().
    uploadBonds(bondsAdj, defaultL=0.4, defaultK=200, angleDeg=120, enableAngles=false) {
        const bondIndices = new Int32Array(this.numAtoms * this.nMaxBonded).fill(-1);
        const bondLenStiff = new Float32Array(this.numAtoms * this.nMaxBonded * 2).fill(0);

        const addBond = (i, j, L, K) => {
            // Find slot
            let slot = -1;
            const base = i * this.nMaxBonded;
            for(let k=0; k<this.nMaxBonded; k++) {
                if(bondIndices[base+k] === -1 || bondIndices[base+k] === j) {
                    slot = k; break;
                }
            }
            if(slot === -1) console.warn(`Atom ${i} max bonds exceeded`);
            else {
                bondIndices[base+slot] = j;
                bondLenStiff[(base+slot)*2 + 0] = L;
                bondLenStiff[(base+slot)*2 + 1] = K;
            }
        };

        // 1. Primary Bonds
        for (let i = 0; i < this.numAtoms; i++) {
            const neighs = bondsAdj[i] || [];
            for (let b of neighs) {
                // b: [neighbor_idx, length, stiffness]
                addBond(i, b[0], b[1] || defaultL, b[2] || defaultK);
            }
        }

        // 2. Angles (as distance constraints)
        // NOTE: Angle-derived constraints are expected to be built in `buildXPDBTopology()` using MMParams.
        // Doing it here as a second pass would double-add angles and can overflow nMaxBonded.
        //
        // if (enableAngles) {
        //     const rads = angleDeg * Math.PI / 180.0;
        //     for (let i = 0; i < this.numAtoms; i++) {
        //         const neighs = bondsAdj[i] || [];
        //         // Iterate pairs of neighbors
        //         for (let a = 0; a < neighs.length; a++) {
        //             for (let b = a + 1; b < neighs.length; b++) {
        //                 const j = neighs[a][0];
        //                 const k = neighs[b][0];
        //                 // Get lengths i-j and i-k
        //                 const Lij = neighs[a][1] || defaultL;
        //                 const Lik = neighs[b][1] || defaultL;
        //                 // Law of cosines
        //                 const Ljk = Math.sqrt(Math.max(0, Lij**2 + Lik**2 - 2*Lij*Lik*Math.cos(rads)));
        //
        //                 // Add constraints between j and k
        //                 addBond(j, k, Ljk, defaultK);
        //                 addBond(k, j, Ljk, defaultK);
        //             }
        //         }
        //     }
        // }

        if ((window.VERBOSITY_LEVEL | 0) >= 3) {
            const used = new Int32Array(this.numAtoms);
            for (let i = 0; i < this.numAtoms; i++) {
                let c = 0;
                const base = i * this.nMaxBonded;
                for (let k = 0; k < this.nMaxBonded; k++) if (bondIndices[base + k] >= 0) c++;
                used[i] = c;
            }
            console.log('[XPDB_WebGPU.uploadBonds][DEBUG] usedSlotsPerAtom=', Array.from(used));
            console.log('[XPDB_WebGPU.uploadBonds][DEBUG] bondIndices=', bondIndices);
            console.log('[XPDB_WebGPU.uploadBonds][DEBUG] bondLenStiff=', bondLenStiff);
        }

        this.uploadPackedTopology({ bondIndices, bondLenStiff });
    }

    async readBuffer(buffer, byteOffset = 0, byteLength = null) {
        if (!buffer) throw new Error('readBuffer: buffer is null/undefined');
        const size = byteLength !== null ? byteLength : buffer.size;
        const staging = this.device.createBuffer({
            size,
            usage: GPUBufferUsage.MAP_READ | GPUBufferUsage.COPY_DST
        });
        const commandEncoder = this.device.createCommandEncoder();
        commandEncoder.copyBufferToBuffer(buffer, byteOffset, staging, 0, size);
        this.device.queue.submit([commandEncoder.finish()]);
        await staging.mapAsync(GPUMapMode.READ);
        const mapped = staging.getMappedRange();
        const copy = new Uint8Array(mapped.byteLength);
        copy.set(new Uint8Array(mapped));
        staging.unmap();
        staging.destroy();
        return copy.buffer;
    }

    async readBuffersAsTyped() {
        const { numAtoms, numGroups, maxGhosts, nMaxBonded } = this;
        const vec4 = 16;
        const i32 = 4;
        const f32 = 4;

        const [posData, predPosData, atomParamsData, bboxesData, ghostData, bondGlobalData, bondLocalData, bondLenStiffData, forceBondData, forceCollData] = await Promise.all([
            this.readBuffer(this.buffers.pos).then(buf => new Float32Array(buf)),
            this.readBuffer(this.buffers.predPos).then(buf => new Float32Array(buf)),
            this.readBuffer(this.buffers.atomParams).then(buf => new Float32Array(buf)),
            this.readBuffer(this.buffers.bboxes).then(buf => new Float32Array(buf)),
            this.readBuffer(this.buffers.ghostPacked).then(buf => new Int32Array(buf)),
            this.readBuffer(this.buffers.bondIdxGlobal).then(buf => new Int32Array(buf)),
            this.readBuffer(this.buffers.bondIdxLocal).then(buf => new Int32Array(buf)),
            this.readBuffer(this.buffers.bondLenStiff).then(buf => new Float32Array(buf)),
            this.readBuffer(this.buffers.forceBond).then(buf => new Float32Array(buf)),
            this.readBuffer(this.buffers.forceColl).then(buf => new Float32Array(buf))
        ]);

        return { posData, predPosData, atomParamsData, bboxesData, ghostData, bondGlobalData, bondLocalData, bondLenStiffData, forceBondData, forceCollData };
    }

    step(arg0, iterations, k_coll, omega, momentum_beta, mousePos=null, pickedIdx=-1, maxRadius=null, coll_scale=2.0, bbox_scale=2.0) {
        const cfg = (arg0 !== null && typeof arg0 === 'object' && !Array.isArray(arg0))
            ? {
                dt: arg0.dt ?? 0.01,
                iterations: arg0.iterations ?? 0,
                k_coll: arg0.k_coll ?? 0.0,
                omega: arg0.omega ?? 0.0,
                momentum_beta: arg0.momentum_beta ?? 0.0,
                mousePos: arg0.mousePos ?? null,
                pickedIdx: arg0.pickedIdx ?? -1,
                maxRadius: arg0.maxRadius ?? null,
                coll_scale: arg0.coll_scale ?? 2.0,
                bbox_scale: arg0.bbox_scale ?? 2.0,
            }
            : {
                dt: arg0 ?? 0.01,
                iterations: iterations ?? 0,
                k_coll: k_coll ?? 0.0,
                omega: omega ?? 0.0,
                momentum_beta: momentum_beta ?? 0.0,
                mousePos,
                pickedIdx,
                maxRadius,
                coll_scale,
                bbox_scale,
            };

        const rMax = (cfg.maxRadius !== null && cfg.maxRadius > 0) ? cfg.maxRadius : 0.5;
        const bboxMargin = rMax * cfg.bbox_scale;
        const collMarginSq = (rMax * cfg.coll_scale) ** 2;

        const u32View = new Uint32Array(4);
        const f32View = new Float32Array(8);
        u32View[0] = this.numAtoms;
        u32View[1] = this.numGroups;
        u32View[2] = cfg.iterations;
        u32View[3] = 0;

        f32View[0] = cfg.dt;
        f32View[1] = cfg.k_coll;
        f32View[2] = cfg.omega;
        f32View[3] = cfg.momentum_beta;
        f32View[4] = collMarginSq;
        f32View[5] = bboxMargin;
        f32View[6] = 0;
        f32View[7] = 0;

        const paramBytes = new ArrayBuffer(64);
        new Uint32Array(paramBytes).set(u32View, 0);
        new Float32Array(paramBytes).set(f32View, 4);
        this.device.queue.writeBuffer(this.buffers.params, 0, paramBytes);

        const copyEncoder = this.device.createCommandEncoder();
        copyEncoder.copyBufferToBuffer(this.buffers.pos, 0, this.buffers.predPos, 0, this.numAtoms * 16);
        this.device.queue.submit([copyEncoder.finish()]);

        if (cfg.pickedIdx >= 0 && cfg.mousePos) {
            const mouseData = new Float32Array([cfg.mousePos[0], cfg.mousePos[1], 0.0, 0.0]);
            this.device.queue.writeBuffer(this.buffers.predPos, cfg.pickedIdx * 16, mouseData);
            this.device.queue.writeBuffer(this.buffers.pos, cfg.pickedIdx * 16, mouseData);
        }

        const encoder = this.device.createCommandEncoder();

        // Pass 1: BBoxes
        const pass1 = encoder.beginComputePass();
        pass1.setPipeline(this.pipelines.updateBBoxes);
        pass1.setBindGroup(0, this.bindGroup);
        pass1.dispatchWorkgroups(this.numGroups);
        pass1.end();

        // Pass 2: Topology
        const pass2 = encoder.beginComputePass();
        pass2.setPipeline(this.pipelines.buildTopology);
        pass2.setBindGroup(0, this.bindGroup);
        pass2.dispatchWorkgroups(this.numGroups);
        pass2.end();

        // Pass 3: Solve
        const pass3 = encoder.beginComputePass();
        pass3.setPipeline(this.pipelines.solve);
        pass3.setBindGroup(0, this.bindGroup);
        pass3.dispatchWorkgroups(this.numGroups);
        pass3.end();

        this.device.queue.submit([encoder.finish()]);
    }

    readPositions() {
        // Read positions from GPU buffer
        const bufferSize = this.numAtoms * 16; // vec4 per atom
        const readBuffer = this.device.createBuffer({
            size: bufferSize,
            usage: GPUBufferUsage.COPY_DST | GPUBufferUsage.MAP_READ
        });

        const encoder = this.device.createCommandEncoder();
        encoder.copyBufferToBuffer(this.buffers.pos, 0, readBuffer, 0, bufferSize);
        this.device.queue.submit([encoder.finish()]);

        // Map and read
        return readBuffer.mapAsync(GPUMapMode.READ).then(() => {
            const data = new Float32Array(readBuffer.getMappedRange());
            const pos3 = new Float32Array(this.numAtoms * 3);
            for (let i = 0; i < this.numAtoms; i++) {
                pos3[i * 3] = data[i * 4];
                pos3[i * 3 + 1] = data[i * 4 + 1];
                pos3[i * 3 + 2] = data[i * 4 + 2];
            }
            readBuffer.unmap();
            return pos3;
        });
    }

    async getPositions() {
        // Readback
        const size = this.numAtoms * 16;
        const readBuf = this.device.createBuffer({
            size,
            usage: GPUBufferUsage.COPY_DST | GPUBufferUsage.MAP_READ
        });

        const encoder = this.device.createCommandEncoder();
        encoder.copyBufferToBuffer(this.buffers.pos, 0, readBuf, 0, size);
        this.device.queue.submit([encoder.finish()]);

        await readBuf.mapAsync(GPUMapMode.READ);
        const data = new Float32Array(readBuf.getMappedRange());
        // Copy out because unmap invalidates
        const result = new Float32Array(data);
        readBuf.unmap();
        return result;
    }
}