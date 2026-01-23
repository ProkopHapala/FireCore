export class RawWebGPUAtomsRenderer {
    constructor(canvas) {
        this.canvas = canvas;
        this.context = null;
        this.adapter = null;
        this.device = null;
        this.format = null;
        this.depthTex = null;
        this.depthView = null;

        this.pipeline = null;
        this.bindGroup = null;

        this.ubuf = null;      // mat4 view, mat4 proj
        this.posBuf = null;    // vec4 xyzr per atom
        this.colBuf = null;    // vec4 rgb+pad per atom

        this.maxAtoms = 0;
        this.nAtoms = 0;

        this.size = { w: 0, h: 0 };
    }

    async init({ maxAtoms = 65536 } = {}) {
        if (!navigator.gpu) throw new Error('navigator.gpu is undefined');
        this.context = this.canvas.getContext('webgpu');
        if (!this.context) throw new Error('canvas.getContext(webgpu) failed');

        this.adapter = await navigator.gpu.requestAdapter();
        if (!this.adapter) throw new Error('navigator.gpu.requestAdapter() returned null');

        this.device = await this.adapter.requestDevice();
        if (!this.device) throw new Error('adapter.requestDevice() returned null');

        this.device.addEventListener('uncapturederror', (e) => {
            const msg = e && e.error ? (e.error.message || String(e.error)) : String(e);
            console.error('WebGPU uncapturederror:', e);
            throw new Error(`WebGPU uncapturederror: ${msg}`);
        });

        this.device.lost.then((info) => {
            console.error('WebGPU device lost:', info);
            throw new Error(`WebGPU device lost: ${info && info.message ? info.message : info}`);
        });

        this.format = navigator.gpu.getPreferredCanvasFormat();

        this.maxAtoms = maxAtoms | 0;
        this._allocBuffers();

        this._buildPipeline();

        this.resize(this.canvas.clientWidth || this.canvas.width, this.canvas.clientHeight || this.canvas.height);
    }

    _allocBuffers() {
        const dev = this.device;
        const n = this.maxAtoms;
        this.ubuf = dev.createBuffer({ size: 16 * 4 * 2, usage: GPUBufferUsage.UNIFORM | GPUBufferUsage.COPY_DST });
        this.posBuf = dev.createBuffer({ size: n * 16, usage: GPUBufferUsage.STORAGE | GPUBufferUsage.COPY_DST });
        this.colBuf = dev.createBuffer({ size: n * 16, usage: GPUBufferUsage.STORAGE | GPUBufferUsage.COPY_DST });
    }

    _buildPipeline() {
        const dev = this.device;

        const wgsl = `
struct UBO {
  view : mat4x4<f32>,
  proj : mat4x4<f32>,
};
@group(0) @binding(0) var<uniform> ubo : UBO;
@group(0) @binding(1) var<storage, read> posRad : array<vec4<f32>>;   // xyz,r
@group(0) @binding(2) var<storage, read> col    : array<vec4<f32>>;   // rgb,1

struct VSOut {
  @builtin(position) pos : vec4<f32>,
  @location(0) uv : vec2<f32>,
  @location(1) rgb : vec3<f32>,
};

@vertex
fn vs_main(@builtin(vertex_index) vid : u32, @builtin(instance_index) iid : u32) -> VSOut {
  var out : VSOut;

  // quad corners in [-1,+1]
  var q = array<vec2<f32>,4>(
    vec2<f32>(-1.0, -1.0),
    vec2<f32>( 1.0, -1.0),
    vec2<f32>(-1.0,  1.0),
    vec2<f32>( 1.0,  1.0)
  );

  let pr = posRad[iid];
  let centerW = vec4<f32>(pr.xyz, 1.0);
  let centerV = ubo.view * centerW;

  let r = pr.w;
  let offs = q[vid] * r;
  let v = vec4<f32>(centerV.x + offs.x, centerV.y + offs.y, centerV.z, 1.0);

  out.pos = ubo.proj * v;
  out.uv = (q[vid] * 0.5) + vec2<f32>(0.5, 0.5);
  out.rgb = col[iid].xyz;
  return out;
}

@fragment
fn fs_main(in : VSOut) -> @location(0) vec4<f32> {
  let d = (in.uv - vec2<f32>(0.5,0.5)) * 2.0;
  let r2 = dot(d,d);
  if (r2 > 1.0) { discard; }

  let z = sqrt(max(1.0 - r2, 0.0));
  let N = normalize(vec3<f32>(d.x, d.y, z));
  let L = normalize(vec3<f32>(0.4, 0.6, 1.0));
  let V = vec3<f32>(0.0, 0.0, 1.0);

  let diff = max(dot(N,L), 0.0);
  let R = reflect(-L, N);
  let spec = pow(max(dot(V, R), 0.0), 32.0);

  let c = in.rgb * (0.25 + 0.75*diff) + vec3<f32>(0.3) * spec;
  return vec4<f32>(c, 1.0);
}
`;

        const module = dev.createShaderModule({ code: wgsl });
        // Surface shader compilation messages
        module.getCompilationInfo().then((info) => {
            const msgs = info && info.messages ? info.messages : [];
            if (msgs.length > 0) {
                for (const m of msgs) {
                    const t = m.type || 'info';
                    const loc = (m.lineNum !== undefined) ? `${m.lineNum}:${m.linePos}` : '';
                    console[t === 'error' ? 'error' : 'warn'](`WGSL ${t} ${loc} ${m.message}`);
                }
            }
        });

        const bindGroupLayout = dev.createBindGroupLayout({
            entries: [
                { binding: 0, visibility: GPUShaderStage.VERTEX, buffer: { type: 'uniform' } },
                { binding: 1, visibility: GPUShaderStage.VERTEX, buffer: { type: 'read-only-storage' } },
                { binding: 2, visibility: GPUShaderStage.VERTEX, buffer: { type: 'read-only-storage' } },
            ]
        });

        const pipelineLayout = dev.createPipelineLayout({ bindGroupLayouts: [bindGroupLayout] });

        try {
            this.pipeline = dev.createRenderPipeline({
            layout: pipelineLayout,
            vertex: { module, entryPoint: 'vs_main' },
            fragment: { module, entryPoint: 'fs_main', targets: [{ format: this.format }] },
            primitive: { topology: 'triangle-strip', stripIndexFormat: undefined, cullMode: 'none' },
            depthStencil: { format: 'depth24plus', depthWriteEnabled: true, depthCompare: 'less' },
            });
        } catch (e) {
            console.error('createRenderPipeline failed', e);
            throw e;
        }

        this.bindGroup = dev.createBindGroup({
            layout: bindGroupLayout,
            entries: [
                { binding: 0, resource: { buffer: this.ubuf } },
                { binding: 1, resource: { buffer: this.posBuf } },
                { binding: 2, resource: { buffer: this.colBuf } },
            ]
        });
    }

    resize(w, h) {
        const dpr = Math.max(1, window.devicePixelRatio || 1);
        const W = Math.max(1, Math.round(w * dpr));
        const H = Math.max(1, Math.round(h * dpr));
        if (W === this.size.w && H === this.size.h) return;
        this.size.w = W;
        this.size.h = H;

        this.canvas.width = W;
        this.canvas.height = H;
        this.canvas.style.width = `${w}px`;
        this.canvas.style.height = `${h}px`;

        this.context.configure({ device: this.device, format: this.format, alphaMode: 'premultiplied' });

        this.depthTex = this.device.createTexture({ size: [W, H], format: 'depth24plus', usage: GPUTextureUsage.RENDER_ATTACHMENT });
        this.depthView = this.depthTex.createView();
    }

    // view, proj are Float32Array length 16 each
    updateCamera(view, proj) {
        if (!view || !proj) throw new Error('updateCamera expects view and proj matrices');
        this.device.queue.writeBuffer(this.ubuf, 0, view.buffer, view.byteOffset, 64);
        this.device.queue.writeBuffer(this.ubuf, 64, proj.buffer, proj.byteOffset, 64);
    }

    // pos is Float32Array [x,y,z]*n; col is Float32Array [r,g,b]*n; rad is Float32Array [r]*n
    updateAtoms(pos, col, rad, n) {
        n = n | 0;
        if (n > this.maxAtoms) throw new Error(`updateAtoms: n=${n} exceeds maxAtoms=${this.maxAtoms}`);
        this.nAtoms = n;

        const pr = new Float32Array(n * 4);
        const cc = new Float32Array(n * 4);
        for (let i = 0; i < n; i++) {
            pr[i * 4] = pos[i * 3];
            pr[i * 4 + 1] = pos[i * 3 + 1];
            pr[i * 4 + 2] = pos[i * 3 + 2];
            pr[i * 4 + 3] = rad ? rad[i] : 0.5;
            cc[i * 4] = col ? col[i * 3] : 1.0;
            cc[i * 4 + 1] = col ? col[i * 3 + 1] : 1.0;
            cc[i * 4 + 2] = col ? col[i * 3 + 2] : 1.0;
            cc[i * 4 + 3] = 1.0;
        }

        this.device.queue.writeBuffer(this.posBuf, 0, pr.buffer);
        this.device.queue.writeBuffer(this.colBuf, 0, cc.buffer);
    }

    render() {
        const dev = this.device;
        const texView = this.context.getCurrentTexture().createView();

        const encoder = dev.createCommandEncoder();
        const pass = encoder.beginRenderPass({
            colorAttachments: [{
                view: texView,
                clearValue: { r: 0.13, g: 0.13, b: 0.13, a: 1.0 },
                loadOp: 'clear',
                storeOp: 'store'
            }],
            depthStencilAttachment: {
                view: this.depthView,
                depthClearValue: 1.0,
                depthLoadOp: 'clear',
                depthStoreOp: 'store'
            }
        });

        pass.setPipeline(this.pipeline);
        pass.setBindGroup(0, this.bindGroup);
        if (this.nAtoms > 0) pass.draw(4, this.nAtoms, 0, 0);
        pass.end();

        dev.queue.submit([encoder.finish()]);
    }
}
