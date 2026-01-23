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

    async init() {
        if (!navigator.gpu) throw new Error("WebGPU not supported");
        const adapter = await navigator.gpu.requestAdapter();
        this.device = await adapter.requestDevice();

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

    // Python 'build_bond_arrays_with_angles' equivalent
    uploadBonds(bondsAdj, defaultL=0.4, defaultK=200, angleDeg=120, enableAngles=true) {
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
        if (enableAngles) {
            const rads = angleDeg * Math.PI / 180.0;
            for (let i = 0; i < this.numAtoms; i++) {
                const neighs = bondsAdj[i] || [];
                // Iterate pairs of neighbors
                for (let a = 0; a < neighs.length; a++) {
                    for (let b = a + 1; b < neighs.length; b++) {
                        const j = neighs[a][0];
                        const k = neighs[b][0];
                        // Get lengths i-j and i-k
                        const Lij = neighs[a][1] || defaultL;
                        const Lik = neighs[b][1] || defaultL;
                        // Law of cosines
                        const Ljk = Math.sqrt(Math.max(0, Lij**2 + Lik**2 - 2*Lij*Lik*Math.cos(rads)));
                        
                        // Add constraints between j and k
                        addBond(j, k, Ljk, defaultK); // usually slightly weaker K for angles?
                        addBond(k, j, Ljk, defaultK);
                    }
                }
            }
        }

        this.device.queue.writeBuffer(this.buffers.bondIdxGlobal, 0, bondIndices);
        this.device.queue.writeBuffer(this.buffers.bondLenStiff, 0, bondLenStiff);
    }

    step(dt, iterations, k_coll, omega, momentum_beta, mousePos=null, pickedIdx=-1, maxRadius=null) {
        // 1. Update Params Uniforms
        // Struct: num_atoms(u32), num_groups, inner_iters, pad, dt(f32), k_coll, omega, beta, margin_sq, bbox_margin...
        const u32View = new Uint32Array(4);
        const f32View = new Float32Array(8);

        const rMax = (maxRadius !== null && maxRadius > 0) ? maxRadius : 0.5;
        const bboxScale = 2.0;
        const collScale = 2.0;
        const bboxMargin = rMax * bboxScale;
        const collMarginSq = (rMax * collScale) ** 2;

        u32View[0] = this.numAtoms;
        u32View[1] = this.numGroups;
        u32View[2] = iterations;
        u32View[3] = 0; // pad

        f32View[0] = dt;
        f32View[1] = k_coll;
        f32View[2] = omega;
        f32View[3] = momentum_beta;
        f32View[4] = collMarginSq;
        f32View[5] = bboxMargin;
        f32View[6] = 0; // pad
        f32View[7] = 0; // pad

        // Combine into one buffer write
        const paramBytes = new ArrayBuffer(64);
        new Uint32Array(paramBytes).set(u32View, 0);
        new Float32Array(paramBytes).set(f32View, 4); // Offset 4 u32s = 16 bytes
        this.device.queue.writeBuffer(this.buffers.params, 0, paramBytes);

        // 2. Handle Interaction (Mouse Drag)
        // In OpenCL code, pred_pos is copied from host. Here we can do a buffer copy or write.
        // Copy curr -> pred
        const copyEncoder = this.device.createCommandEncoder();
        copyEncoder.copyBufferToBuffer(this.buffers.pos, 0, this.buffers.predPos, 0, this.numAtoms * 16);
        this.device.queue.submit([copyEncoder.finish()]);

        if (pickedIdx >= 0 && mousePos) {
            // Write specific predicted position for picked atom
            const mouseData = new Float32Array([mousePos[0], mousePos[1], 0.0, 0.0]);
            this.device.queue.writeBuffer(this.buffers.predPos, pickedIdx * 16, mouseData);
            // Also force current pos to prevent spring-back lag
            this.device.queue.writeBuffer(this.buffers.pos, pickedIdx * 16, mouseData);
        }

        // 3. Dispatch Compute
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