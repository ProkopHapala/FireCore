export class RawWebGPUAtomsRenderer {
    constructor(canvas) {
        this.canvas = canvas;
        this.context = null;
        this.adapter = null;
        this.device = null;
        this.format = null;
        this.depthTex = null;
        this.depthView = null;

        // --- Pipelines ---
        this.pipelineAtoms = null;
        this.pipelineHalo  = null;
        this.pipelineBonds = null;
        this.pipelineLabels = null;
        this.pipelineClear = null;
        this.pipelineAtlasDebug = null;

        // --- Bind groups ---
        this.bindGroupAtoms = null;
        this.bindGroupHalo  = null;
        this.bindGroupBonds = null;
        this.bindGroupLabels = null;

        this.ubuf = null;      // mat4 view, mat4 proj
        this.posBuf = null;    // vec4 xyzr per atom
        this.colBuf = null;    // vec4 rgb+pad per atom

        // selection halo scale per atom (f32). 0.0 -> not selected, >0 -> draw halo.
        this.selBuf = null;

        // bonds: per-bond indices
        this.bondBuf = null;
        this.maxBonds = 0;
        this.nBonds = 0;

        // labels: one instance per character quad
        this.fontTex = null;
        this.fontView = null;
        this.fontSampler = null;
        this.labelCharBuf = null;
        this.labelUBuf = null; // label settings
        this.maxLabelChars = 0;
        this.nLabelChars = 0;
        this.labelsVisible = false;
        this.debugShowLabelAtlas = true;

        this.maxAtoms = 0;
        this.nAtoms = 0;

        this.size = { w: 0, h: 0 };

        // Controls / settings
        this.bondsVisible = true;
        this.haloVisible  = true;
        this.labelScale = 0.5;
        this.labelColor = [1.0, 1.0, 1.0];
        this.bondWidth = 0.08; // in view-space units; tweak later
        this.debugForceSolidLabels = false;
        this.debugLabelMode = 0; // Default: normal rendering
    }

    async init({ maxAtoms = 65536 } = {}) {
        if (!navigator.gpu) throw new Error('navigator.gpu is undefined');
        this.context = this.canvas.getContext('webgpu');
        if (!this.context) throw new Error('canvas.getContext(webgpu) failed');

        this.adapter = await navigator.gpu.requestAdapter();
        if (!this.adapter) throw new Error('navigator.gpu.requestAdapter() returned null');

        this.device = await this.adapter.requestDevice();
        if (!this.device) throw new Error('adapter.requestDevice() returned null');

        this.lastUncapturedError = null;
        this.device.addEventListener('uncapturederror', (e) => {
            const msg = e && e.error ? (e.error.message || String(e.error)) : String(e);
            this.lastUncapturedError = msg;
            console.error('[RawWebGPUAtomsRenderer] WebGPU uncapturederror', { msg, error: e && e.error ? e.error : e });
            // IMPORTANT: do NOT throw here; throwing hides follow-up diagnostics and may break the app event loop.
            // We want a loud console + continued execution to locate the triggering call.
        });

        this.device.lost.then((info) => {
            console.error('[RawWebGPUAtomsRenderer] WebGPU device lost', info);
        });

        this.format = navigator.gpu.getPreferredCanvasFormat();

        // Configure swapchain immediately (will be reconfigured on resize and each render as well)
        this.device.pushErrorScope('validation');
        this._configureSwapChain(this.canvas.clientWidth || this.canvas.width || 1, this.canvas.clientHeight || this.canvas.height || 1);
        this.device.popErrorScope().then((err) => {
            if (err) console.error('[RawWebGPUAtomsRenderer] validation error during _configureSwapChain(init)', err);
            else console.log('[RawWebGPUAtomsRenderer] _configureSwapChain(init) ok');
        });

        this.maxAtoms = maxAtoms | 0;
        this.maxBonds = (this.maxAtoms * 4) | 0;
        this.maxLabelChars = (this.maxAtoms * 8) | 0;
        this._allocBuffers();

        await this._buildPipelines();
    }

    _allocBuffers() {
        const dev = this.device;
        const n = this.maxAtoms;
        this.ubuf = dev.createBuffer({ size: 16 * 4 * 2, usage: GPUBufferUsage.UNIFORM | GPUBufferUsage.COPY_DST });
        this.posBuf = dev.createBuffer({ size: n * 16, usage: GPUBufferUsage.STORAGE | GPUBufferUsage.COPY_DST });
        this.colBuf = dev.createBuffer({ size: n * 16, usage: GPUBufferUsage.STORAGE | GPUBufferUsage.COPY_DST });

        this.selBuf = dev.createBuffer({ size: n * 4, usage: GPUBufferUsage.STORAGE | GPUBufferUsage.COPY_DST });

        // bond buffer: u32 atomA, u32 atomB, u32 matID, u32 pad
        this.bondBuf = dev.createBuffer({ size: this.maxBonds * 16, usage: GPUBufferUsage.STORAGE | GPUBufferUsage.COPY_DST });

        // label settings UBO: vec4 color, vec4 params (scale, aspect, unused, unused)
        this.labelUBuf = dev.createBuffer({ size: 32, usage: GPUBufferUsage.UNIFORM | GPUBufferUsage.COPY_DST });
        // label char buffer: u32 atomID, u32 charCode, f32 offX, f32 offY
        this.labelCharBuf = dev.createBuffer({ size: this.maxLabelChars * 16, usage: GPUBufferUsage.STORAGE | GPUBufferUsage.COPY_DST });
    }

    async _ensureFontAtlas() {
        if (this.fontTex && this.fontView && this.fontSampler) return;
        console.debug('[RawWebGPUAtomsRenderer] building font atlas texture');
        const size = 512;
        const cols = 16;
        const rows = 6; // 96 printable ASCII chars
        const canvas = document.createElement('canvas');
        canvas.width = size;
        canvas.height = Math.ceil((96 / cols) * (size / cols));
        const ctx = canvas.getContext('2d', { willReadFrequently: true });
        // ctx.fillStyle = '#000000';
        // ctx.fillRect(0, 0, canvas.width, canvas.height);
        ctx.clearRect(0, 0, canvas.width, canvas.height);
        const grad = ctx.createLinearGradient(0, 0, canvas.width, canvas.height);
        grad.addColorStop(0.0, 'rgba(255,10,10,0.5)'); // Slightly more opaque
        grad.addColorStop(0.5, 'rgba(10,255,10,0.5)');
        grad.addColorStop(1.0, 'rgba(10,10,255,0.5)');
        ctx.fillStyle = grad;
        ctx.fillRect(0, 0, canvas.width, canvas.height);
        ctx.fillStyle = 'white'; // Standard white
        // Ensure text is opaque
        ctx.globalAlpha = 1.0;
        ctx.font = `${Math.floor(size / cols * 0.8)}px monospace`;
        ctx.textBaseline = 'middle';
        ctx.textAlign = 'center';
        const cw = canvas.width / cols;
        const ch = canvas.height / Math.ceil(96 / cols);
        console.log('[RawWebGPUAtomsRenderer] font atlas layout', { size, cols, rows, canvasW: canvas.width, canvasH: canvas.height, cw, ch });
        for (let i = 0; i < 96; i++) {
            const col = i % cols;
            const row = Math.floor(i / cols);
            const x = col * cw + cw * 0.5;
            const y = row * ch + ch * 0.5;
            ctx.fillText(String.fromCharCode(32 + i), x, y);
        }
        // Pixel probes: sample background and a few known glyph centers.
        // NOTE: We intentionally do NOT sample at (0,0) only; that's often empty.
        {
            const bg = ctx.getImageData(0, 0, 1, 1).data;
            console.log('[RawWebGPUAtomsRenderer] font atlas probe bg(0,0) rgba', Array.from(bg));

            const probeChar = (chAscii) => {
                const i = (chAscii - 32) | 0;
                const col = i % cols;
                const row = Math.floor(i / cols);
                const px = Math.floor(col * cw + cw * 0.5);
                const py = Math.floor(row * ch + ch * 0.5);
                const rgba = ctx.getImageData(px, py, 1, 1).data;
                console.log('[RawWebGPUAtomsRenderer] font atlas probe glyph center', { ch: String.fromCharCode(chAscii), chAscii, i, col, row, px, py, rgba: Array.from(rgba) });
            };
            probeChar(48);  // '0'
            probeChar(49);  // '1'
            probeChar(65);  // 'A'
            probeChar(88);  // 'X'
        }
        const bmp = await createImageBitmap(canvas);
        console.log('[RawWebGPUAtomsRenderer] font atlas bitmap created', { width: bmp.width, height: bmp.height });
        try {
            if (!this._fontAtlasDebugImg) {
                const img = document.createElement('img');
                img.style.position = 'fixed';
                img.style.bottom = '8px';
                img.style.left = '8px';
                img.style.width = '128px';
                img.style.border = '1px solid red';
                img.style.zIndex = '9999';
                img.style.pointerEvents = 'none';
                img.alt = 'FontAtlasDebug';
                document.body.appendChild(img);
                this._fontAtlasDebugImg = img;
            }
            if (this._fontAtlasDebugImg) this._fontAtlasDebugImg.src = canvas.toDataURL();
        } catch (e) {
            console.warn('[RawWebGPUAtomsRenderer] font atlas debug image failed', e);
        }
        this.fontTex = this.device.createTexture({
            size: [bmp.width, bmp.height],
            format: 'rgba8unorm',
            usage: GPUTextureUsage.TEXTURE_BINDING | GPUTextureUsage.COPY_DST | GPUTextureUsage.COPY_SRC | GPUTextureUsage.RENDER_ATTACHMENT
        });
        this.device.queue.copyExternalImageToTexture({ source: bmp }, { texture: this.fontTex }, [bmp.width, bmp.height]);
        console.log('[RawWebGPUAtomsRenderer] font atlas copied to GPU texture');
        this.fontView = this.fontTex.createView();
        this.fontSampler = this.device.createSampler({ magFilter: 'linear', minFilter: 'linear' });
        console.debug('[RawWebGPUAtomsRenderer] font atlas ready', { width: bmp.width, height: bmp.height });
        this._refreshLabelBindGroup('fontAtlasReady');
        const gpuProbePoints = [
            { label: 'bg00', x: 0, y: 0 },
            { label: 'glyph_0_center', x: Math.floor(0 * cw + cw * 0.5), y: Math.floor(1 * ch + ch * 0.5) },
            { label: 'glyph_1_center', x: Math.floor(1 * cw + cw * 0.5), y: Math.floor(1 * ch + ch * 0.5) },
            { label: 'glyph_A_center', x: Math.floor((65 - 32) % cols * cw + cw * 0.5), y: Math.floor(Math.floor((65 - 32) / cols) * ch + ch * 0.5) },
            { label: 'glyph_X_center', x: Math.floor((88 - 32) % cols * cw + cw * 0.5), y: Math.floor(Math.floor((88 - 32) / cols) * ch + ch * 0.5) },
        ];
        await this._probeFontAtlasGPU(gpuProbePoints);
    }

    async _probeFontAtlasGPU(points) {
        if (!this.device || !this.fontTex || !points || !points.length) return;
        if (!(this.fontTex.usage & GPUTextureUsage.COPY_SRC)) {
            console.warn('[RawWebGPUAtomsRenderer] font atlas GPU probe skipped (COPY_SRC unsupported)');
            return;
        }
        const bytesPerRow = 256; // WebGPU requires 256-byte alignment
        const bufSize = bytesPerRow * points.length;
        const probeBuf = this.device.createBuffer({ size: bufSize, usage: GPUBufferUsage.COPY_DST | GPUBufferUsage.MAP_READ });
        try {
            const encoder = this.device.createCommandEncoder();
            for (let i = 0; i < points.length; i++) {
                const pt = points[i];
                encoder.copyTextureToBuffer(
                    { texture: this.fontTex, origin: { x: Math.max(0, pt.x | 0), y: Math.max(0, pt.y | 0), z: 0 } },
                    { buffer: probeBuf, offset: i * bytesPerRow, bytesPerRow, rowsPerImage: 1 },
                    { width: 1, height: 1, depthOrArrayLayers: 1 }
                );
            }
            this.device.queue.submit([encoder.finish()]);
            await this.device.queue.onSubmittedWorkDone();
            const mapMode = (globalThis.GPUMapMode && GPUMapMode.READ) ? GPUMapMode.READ : 1;
            await probeBuf.mapAsync(mapMode);
            const data = new Uint8Array(probeBuf.getMappedRange());
            const results = points.map((pt, i) => {
                const off = i * bytesPerRow;
                return {
                    label: pt.label,
                    coord: { x: pt.x, y: pt.y },
                    rgba: Array.from(data.slice(off, off + 4))
                };
            });
            probeBuf.unmap();
            console.log('[RawWebGPUAtomsRenderer] font atlas GPU probes', results);
        } catch (err) {
            console.warn('[RawWebGPUAtomsRenderer] font atlas GPU probe failed', err);
        }
    }

    async _buildPipelines() {
        const dev = this.device;

        await this._ensureFontAtlas();

        // A minimal clear pipeline to avoid Dawn's internal clear path which can require COPY_DST on some backends.
        // We draw a fullscreen quad with depth=1 to effectively clear both color and depth.
        {
            const clearWGSL = `
struct VSOut { @builtin(position) pos: vec4<f32> };
@vertex
fn vs(@builtin(vertex_index) vid:u32) -> VSOut {
  var out:VSOut;
  var q = array<vec2<f32>,4>( vec2<f32>(-1.0,-1.0), vec2<f32>(1.0,-1.0), vec2<f32>(-1.0,1.0), vec2<f32>(1.0,1.0) );
  out.pos = vec4<f32>(q[vid], 1.0, 1.0);
  return out;
}
@fragment
fn fs() -> @location(0) vec4<f32> {
  return vec4<f32>(0.13, 0.13, 0.13, 1.0);
}
`;
            const m = dev.createShaderModule({ code: clearWGSL });
            m.getCompilationInfo().then((info) => {
                const msgs = info && info.messages ? info.messages : [];
                for (const mm of msgs) {
                    const t = mm.type || 'info';
                    const loc = (mm.lineNum !== undefined) ? `${mm.lineNum}:${mm.linePos}` : '';
                    console[t === 'error' ? 'error' : 'warn'](`[clear] WGSL ${t} ${loc} ${mm.message}`);
                }
            });
            this.pipelineClear = dev.createRenderPipeline({
                layout: 'auto',
                vertex: { module: m, entryPoint: 'vs' },
                fragment: { module: m, entryPoint: 'fs', targets: [{ format: this.format }] },
                primitive: { topology: 'triangle-strip', cullMode: 'none' },
                depthStencil: { format: 'depth24plus', depthWriteEnabled: true, depthCompare: 'always' },
            });
            console.log('[RawWebGPUAtomsRenderer] pipelineClear created', { ok: !!this.pipelineClear });
        }

        const commonWGSL = `
struct UBO { view : mat4x4<f32>, proj : mat4x4<f32>, };
@group(0) @binding(0) var<uniform> ubo : UBO;
@group(0) @binding(1) var<storage, read> posRad : array<vec4<f32>>;   // xyz,r
@group(0) @binding(2) var<storage, read> col    : array<vec4<f32>>;   // rgb,1
`;

        const atomsWGSL = commonWGSL + `
struct VSOut { @builtin(position) pos:vec4<f32>, @location(0) uv:vec2<f32>, @location(1) rgb:vec3<f32>, };

@vertex
fn vs_atoms(@builtin(vertex_index) vid:u32, @builtin(instance_index) iid:u32) -> VSOut {
  var out:VSOut;
  var q = array<vec2<f32>,4>( vec2<f32>(-1.0,-1.0), vec2<f32>(1.0,-1.0), vec2<f32>(-1.0,1.0), vec2<f32>(1.0,1.0) );
  let pr = posRad[iid];
  let centerV = ubo.view * vec4<f32>(pr.xyz, 1.0);
  let r = pr.w;
  let offs = q[vid] * r;
  let v = vec4<f32>(centerV.x + offs.x, centerV.y + offs.y, centerV.z, 1.0);
  out.pos = ubo.proj * v;
  out.uv = (q[vid] * 0.5) + vec2<f32>(0.5,0.5);
  out.rgb = col[iid].xyz;
  return out;
}

@fragment
fn fs_atoms(in:VSOut) -> @location(0) vec4<f32> {
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

        const haloWGSL = commonWGSL + `
@group(0) @binding(3) var<storage, read> selScale : array<f32>;
struct VSOut { @builtin(position) pos:vec4<f32>, @location(0) uv:vec2<f32>, @location(1) a:f32, };

@vertex
fn vs_halo(@builtin(vertex_index) vid:u32, @builtin(instance_index) iid:u32) -> VSOut {
  var out:VSOut;
  var q = array<vec2<f32>,4>( vec2<f32>(-1.0,-1.0), vec2<f32>(1.0,-1.0), vec2<f32>(-1.0,1.0), vec2<f32>(1.0,1.0) );
  let s = selScale[iid];
  if (s <= 0.0) { out.pos=vec4<f32>(2.0,2.0,2.0,1.0); out.uv=vec2<f32>(0.0); out.a=0.0; return out; }
  let pr = posRad[iid];
  let centerV = ubo.view * vec4<f32>(pr.xyz, 1.0);
  let r = pr.w * s;
  let offs = q[vid] * r;
  let v = vec4<f32>(centerV.x + offs.x, centerV.y + offs.y, centerV.z, 1.0);
  out.pos = ubo.proj * v;
  out.uv = (q[vid] * 0.5) + vec2<f32>(0.5,0.5);
  out.a = 0.22;
  return out;
}

@fragment
fn fs_halo(in:VSOut) -> @location(0) vec4<f32> {
  if (in.a <= 0.0) { discard; }
  let d = (in.uv - vec2<f32>(0.5,0.5)) * 2.0;
  let r2 = dot(d,d);
  if (r2 > 1.0) { discard; }
  return vec4<f32>(1.0, 1.0, 0.0, in.a);
}
`;

        const bondsWGSL = commonWGSL + `
struct Bond { a:u32, b:u32, mat:u32, pad:u32 };
@group(0) @binding(3) var<storage, read> bonds : array<Bond>;
struct VSOut { @builtin(position) pos:vec4<f32>, };

@vertex
fn vs_bonds(@builtin(vertex_index) vid:u32, @builtin(instance_index) iid:u32) -> VSOut {
  var out:VSOut;
  let bond = bonds[iid];
  let a = posRad[bond.a];
  let b = posRad[bond.b];
  let A = (ubo.view * vec4<f32>(a.xyz,1.0)).xyz;
  let B = (ubo.view * vec4<f32>(b.xyz,1.0)).xyz;
  let d = B - A;
  let L = max(length(d), 1e-6);
  let dir = d / L;
  let up = normalize(cross(dir, vec3<f32>(0.0, 0.0, -1.0)));
  let w = 0.08;
  let side = select(-1.0, 1.0, (vid & 1u) == 1u);
  let t = select(0.0, 1.0, vid >= 2u);
  let P = A + dir * (t * L) + up * (side * w);
  out.pos = ubo.proj * vec4<f32>(P, 1.0);
  return out;
}

@fragment
fn fs_bonds(in:VSOut) -> @location(0) vec4<f32> { return vec4<f32>(1.0,1.0,1.0,0.9); }
`;

        const labelsWGSL = commonWGSL + `
struct LabelChar { atom:u32, code:u32, offx:f32, offy:f32 }; // code is 0..95 (printable ASCII - 32)
@group(0) @binding(3) var<storage, read> labelChars : array<LabelChar>;
struct LabelUBO { color: vec4<f32>, params: vec4<f32> };
@group(0) @binding(4) var<uniform> labelU : LabelUBO;
@group(0) @binding(5) var fontTex : texture_2d<f32>;
@group(0) @binding(6) var fontSamp : sampler;

struct VSOut { @builtin(position) pos:vec4<f32>, @location(0) uv:vec2<f32>, @location(1) @interpolate(flat) code:f32 };

@vertex
fn vs_labels(@builtin(vertex_index) vid:u32, @builtin(instance_index) iid:u32) -> VSOut {
  var out:VSOut;
  var q = array<vec2<f32>,4>( vec2<f32>(-1.0,-1.0), vec2<f32>(1.0,-1.0), vec2<f32>(-1.0,1.0), vec2<f32>(1.0,1.0) );
  let ch = labelChars[iid];
  let pr = posRad[ch.atom];
  let centerV = (ubo.view * vec4<f32>(pr.xyz,1.0)).xyz;
  let scale = labelU.params.x;
  let aspect = labelU.params.y;
  let offs = vec2<f32>(ch.offx / aspect, ch.offy);
  let P = vec3<f32>(centerV.x + offs.x, centerV.y + offs.y, centerV.z);
  let corner = q[vid] * scale;
  out.pos = ubo.proj * vec4<f32>(P.x + corner.x, P.y + corner.y, P.z, 1.0);
  // Flip UV.y because canvas/texture top is row 0, but quad corner.y is positive at the top.
  out.uv = (q[vid] * 0.5) + vec2<f32>(0.5,0.5);
  out.uv.y = 1.0 - out.uv.y; 
  out.code = f32(ch.code);
  return out;
}

@fragment
fn fs_labels(in:VSOut) -> @location(0) vec4<f32> {
  let dbg = labelU.params.w;
  if (labelU.params.z > 0.5) {
    return vec4<f32>(1.0, 1.0, 0.0, 1.0);
  }
  let code = in.code;
  let cols = 16.0;
  let rows = 6.0;
  let colI = floor(code % cols);
  let rowI = floor(code / cols);
  let cell = vec2<f32>(1.0/cols, 1.0/rows);
  let base = vec2<f32>(colI*cell.x, rowI*cell.y);
  let fuv = base + in.uv*cell;
  
  let texel = textureSample(fontTex, fontSamp, fuv);
  
  if (dbg > 0.5 && dbg < 1.5) { // Mode 1: Debug UVs mapped to Atlas
     // If we see gradients here, it means we are sampling the background.
     // If we see black, it means we are sampling outside or texture is empty.
     return vec4<f32>(fuv.x, fuv.y, 0.0, 1.0);
  }
  if (dbg >= 1.5 && dbg < 2.5) { // Mode 2: Atlas Sampled Color
    return vec4<f32>(texel.rgb, 1.0);
  }
  if (dbg >= 2.5 && dbg < 3.5) { // Mode 3: Atlas Sampled Alpha
    return vec4<f32>(texel.aaa, 1.0);
  }
  if (dbg >= 3.5 && dbg < 4.5) { // Mode 4: Character Code Visualizer
    return vec4<f32>(code/96.0, fract(code/10.0), fract(code/16.0), 1.0);
  }

  // Normal mode: Alpha blending with white text
  // We expect alpha to be > 0.5 for glyphs if background is 0.5
  if (texel.a < 0.6) { discard; }
  return vec4<f32>(labelU.color.xyz, 1.0);
}
`;

        const mkModule = (code, tag) => {
            const module = dev.createShaderModule({ code });
            module.getCompilationInfo().then((info) => {
                const msgs = info && info.messages ? info.messages : [];
                for (const m of msgs) {
                    const t = m.type || 'info';
                    const loc = (m.lineNum !== undefined) ? `${m.lineNum}:${m.linePos}` : '';
                    console[t === 'error' ? 'error' : 'warn'](`[${tag}] WGSL ${t} ${loc} ${m.message}`);
                }
            });
            return module;
        };

        const atomsModule  = mkModule(atomsWGSL, 'atoms');
        const haloModule   = mkModule(haloWGSL, 'halo');
        const bondsModule  = mkModule(bondsWGSL, 'bonds');
        const labelsModule = mkModule(labelsWGSL, 'labels');

        const bglAtoms = dev.createBindGroupLayout({ entries: [
            { binding: 0, visibility: GPUShaderStage.VERTEX, buffer: { type: 'uniform' } },
            { binding: 1, visibility: GPUShaderStage.VERTEX, buffer: { type: 'read-only-storage' } },
            { binding: 2, visibility: GPUShaderStage.VERTEX, buffer: { type: 'read-only-storage' } },
        ]});
        const bglHalo = dev.createBindGroupLayout({ entries: [
            { binding: 0, visibility: GPUShaderStage.VERTEX, buffer: { type: 'uniform' } },
            { binding: 1, visibility: GPUShaderStage.VERTEX, buffer: { type: 'read-only-storage' } },
            { binding: 2, visibility: GPUShaderStage.VERTEX, buffer: { type: 'read-only-storage' } },
            { binding: 3, visibility: GPUShaderStage.VERTEX, buffer: { type: 'read-only-storage' } },
        ]});
        const bglBonds = dev.createBindGroupLayout({ entries: [
            { binding: 0, visibility: GPUShaderStage.VERTEX, buffer: { type: 'uniform' } },
            { binding: 1, visibility: GPUShaderStage.VERTEX, buffer: { type: 'read-only-storage' } },
            { binding: 2, visibility: GPUShaderStage.VERTEX, buffer: { type: 'read-only-storage' } },
            { binding: 3, visibility: GPUShaderStage.VERTEX, buffer: { type: 'read-only-storage' } },
        ]});
        const bglLabels = dev.createBindGroupLayout({ entries: [
            { binding: 0, visibility: GPUShaderStage.VERTEX, buffer: { type: 'uniform' } },
            { binding: 1, visibility: GPUShaderStage.VERTEX, buffer: { type: 'read-only-storage' } },
            { binding: 2, visibility: GPUShaderStage.VERTEX, buffer: { type: 'read-only-storage' } },
            { binding: 3, visibility: GPUShaderStage.VERTEX, buffer: { type: 'read-only-storage' } },
            { binding: 4, visibility: GPUShaderStage.VERTEX | GPUShaderStage.FRAGMENT, buffer: { type: 'uniform' } },
            { binding: 5, visibility: GPUShaderStage.FRAGMENT, texture: { sampleType: 'float' } },
            { binding: 6, visibility: GPUShaderStage.FRAGMENT, sampler: { type: 'filtering' } },
        ]});
        this._bglLabels = bglLabels;
        const bglAtlasDebug = dev.createBindGroupLayout({ entries: [
            { binding: 0, visibility: GPUShaderStage.FRAGMENT, texture: { sampleType: 'float' } },
            { binding: 1, visibility: GPUShaderStage.FRAGMENT, sampler: { type: 'filtering' } },
        ]});

        const mkPipe = (module, vs, fs, bgl, blend, disableDepth) => {
            const layout = dev.createPipelineLayout({ bindGroupLayouts: [bgl] });
            const targets = [{
                format: this.format,
                blend: blend ? {
                    color: { srcFactor: 'src-alpha', dstFactor: 'one-minus-src-alpha', operation: 'add' },
                    alpha: { srcFactor: 'one', dstFactor: 'one-minus-src-alpha', operation: 'add' }
                } : undefined
            }];
            const depthStencil = disableDepth
                ? { format: 'depth24plus', depthWriteEnabled: false, depthCompare: 'always' }
                : blend
                    ? { format: 'depth24plus', depthWriteEnabled: false, depthCompare: 'less' }
                    : { format: 'depth24plus', depthWriteEnabled: true, depthCompare: 'less' };
            return dev.createRenderPipeline({
                layout,
                vertex: { module, entryPoint: vs },
                fragment: { module, entryPoint: fs, targets },
                primitive: { topology: 'triangle-strip', cullMode: 'none' },
                depthStencil,
            });
        };

        try {
            this.pipelineAtoms  = mkPipe(atomsModule,  'vs_atoms',  'fs_atoms',  bglAtoms,  false, false);
            this.pipelineHalo   = mkPipe(haloModule,   'vs_halo',   'fs_halo',   bglHalo,   true, false);
            this.pipelineBonds  = mkPipe(bondsModule,  'vs_bonds',  'fs_bonds',  bglBonds,  true, false);
            this.pipelineLabels = mkPipe(labelsModule, 'vs_labels', 'fs_labels', bglLabels, true, true);
            console.debug('[RawWebGPUAtomsRenderer] pipelines created OK');
        } catch (e) {
            console.error('RawWebGPUAtomsRenderer: createRenderPipeline failed', e);
            throw e;
        }

        this.bindGroupAtoms = dev.createBindGroup({ layout: bglAtoms, entries: [
            { binding: 0, resource: { buffer: this.ubuf } },
            { binding: 1, resource: { buffer: this.posBuf } },
            { binding: 2, resource: { buffer: this.colBuf } },
        ]});
        this.bindGroupHalo = dev.createBindGroup({ layout: bglHalo, entries: [
            { binding: 0, resource: { buffer: this.ubuf } },
            { binding: 1, resource: { buffer: this.posBuf } },
            { binding: 2, resource: { buffer: this.colBuf } },
            { binding: 3, resource: { buffer: this.selBuf } },
        ]});
        this.bindGroupBonds = dev.createBindGroup({ layout: bglBonds, entries: [
            { binding: 0, resource: { buffer: this.ubuf } },
            { binding: 1, resource: { buffer: this.posBuf } },
            { binding: 2, resource: { buffer: this.colBuf } },
            { binding: 3, resource: { buffer: this.bondBuf } },
        ]});
        this._refreshLabelBindGroup('buildPipelines');

        // Debug atlas pipeline/bind group (optional overlay to visualize font texture)
        const atlasDebugWGSL = `
struct VSOut { @builtin(position) pos:vec4<f32>, @location(0) uv:vec2<f32> };
@vertex
fn vs(@builtin(vertex_index) vid:u32) -> VSOut {
  var positions = array<vec2<f32>,6>(
    vec2<f32>(-1.0, -1.0), vec2<f32>(-0.4, -1.0), vec2<f32>(-1.0, -0.4),
    vec2<f32>(-1.0, -0.4), vec2<f32>(-0.4, -1.0), vec2<f32>(-0.4, -0.4)
  );
  var uvs = array<vec2<f32>,6>(
    vec2<f32>(0.0, 1.0), vec2<f32>(1.0, 1.0), vec2<f32>(0.0, 0.0),
    vec2<f32>(0.0, 0.0), vec2<f32>(1.0, 1.0), vec2<f32>(1.0, 0.0)
  );
  var out:VSOut;
  out.pos = vec4<f32>(positions[vid], 0.0, 1.0);
  out.uv = uvs[vid];
  return out;
}
@group(0) @binding(0) var fontTex : texture_2d<f32>;
@group(0) @binding(1) var fontSamp : sampler;
@fragment
fn fs(in:VSOut) -> @location(0) vec4<f32> {
  let c = textureSample(fontTex, fontSamp, in.uv);
  return vec4<f32>(c.rgb, 1.0);
}
`;
        this.pipelineAtlasDebug = dev.createRenderPipeline({
            layout: dev.createPipelineLayout({ bindGroupLayouts: [bglAtlasDebug] }),
            vertex: { module: dev.createShaderModule({ code: atlasDebugWGSL }), entryPoint: 'vs' },
            fragment: { module: dev.createShaderModule({ code: atlasDebugWGSL }), entryPoint: 'fs', targets: [{ format: this.format }] },
            primitive: { topology: 'triangle-list', cullMode: 'none' },
            depthStencil: { format: 'depth24plus', depthWriteEnabled: false, depthCompare: 'always' },
        });
        this.bindGroupAtlasDebug = dev.createBindGroup({ layout: bglAtlasDebug, entries: [
            { binding: 0, resource: this.fontView },
            { binding: 1, resource: this.fontSampler },
        ]});

        // Initialize label UBO and selection buffer
        this._uploadLabelUBO();
        this.device.queue.writeBuffer(this.selBuf, 0, new Float32Array(this.maxAtoms));
    }

    resize(w, h) {
        this._configureSwapChain(w, h);
    }

    _configureSwapChain(widthPixels, heightPixels) {
        const dpr = Math.max(1, window.devicePixelRatio || 1);
        const W = Math.max(1, Math.round((widthPixels || 1) * dpr));
        const H = Math.max(1, Math.round((heightPixels || 1) * dpr));
        const needsConfig = (W !== this.size.w) || (H !== this.size.h) || !this.depthView;
        if (!needsConfig) return;

        console.log('[RawWebGPUAtomsRenderer] configure swapchain !!!!', widthPixels, heightPixels);
        console.log('[RawWebGPUAtomsRenderer] --- configure swapchain !!!!', widthPixels, heightPixels);

        this.size.w = W;
        this.size.h = H;

        this.canvas.width = W;
        this.canvas.height = H;
        this.canvas.style.width = `${widthPixels}px`;
        this.canvas.style.height = `${heightPixels}px`;

        try { this.context.unconfigure(); } catch (e) { /* ignore */ }

        const cfg = {
            device: this.device,
            format: this.format,
            alphaMode: 'premultiplied',
            // IMPORTANT: keep this minimal. On some Dawn/Vulkan configurations requesting COPY_DST still yields
            // a swapchain texture without COPY_DST but Dawn then tries to clear via a copy path and throws.
            usage: GPUTextureUsage.RENDER_ATTACHMENT | GPUTextureUsage.COPY_DST
        };
        console.debug('[RawWebGPUAtomsRenderer] configure swapchain', cfg, 'size', W, H);
        console.log('[RawWebGPUAtomsRenderer] _configureSwapChain stage=configure begin');
        this.context.configure(cfg);
        console.log('[RawWebGPUAtomsRenderer] _configureSwapChain stage=configure end');

        // NOTE: Do NOT call getCurrentTexture() here; it must be called exactly once per frame (inside render()).
        console.log('[RawWebGPUAtomsRenderer] _configureSwapChain stage=getCurrentTexture skipped (handled in render)');

        if (this.depthTex) {
            // Commented-out old pointer-only reassignment (kept for reference)
            // this.depthTex = this.device.createTexture({ size: [W, H], format: 'depth24plus', usage: GPUTextureUsage.RENDER_ATTACHMENT | GPUTextureUsage.COPY_SRC | GPUTextureUsage.COPY_DST });
            this.depthTex.destroy();
        }
        this.depthTex = this.device.createTexture({ size: [W, H], format: 'depth24plus', usage: GPUTextureUsage.RENDER_ATTACHMENT });
        this.depthView = this.depthTex.createView();

        this._uploadLabelUBO();
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

    updateSelection(selectedIndices, nAtoms, haloScale = 1.3) {
        nAtoms = (nAtoms | 0);
        if (nAtoms <= 0) {
            this.device.queue.writeBuffer(this.selBuf, 0, new Float32Array(this.maxAtoms));
            return;
        }
        const arr = new Float32Array(nAtoms);
        if (selectedIndices && selectedIndices.length) {
            for (let i = 0; i < selectedIndices.length; i++) {
                const ia = selectedIndices[i] | 0;
                if (ia >= 0 && ia < nAtoms) arr[ia] = haloScale;
            }
        }
        this.device.queue.writeBuffer(this.selBuf, 0, arr.buffer);
    }

    updateBonds(pairs) {
        const bonds = pairs || [];
        const n = bonds.length | 0;
        if (n > this.maxBonds) throw new Error(`updateBonds: nBonds=${n} exceeds maxBonds=${this.maxBonds}`);
        this.nBonds = n;
        const buf = new Uint32Array(n * 4);
        for (let i = 0; i < n; i++) {
            const e = bonds[i];
            buf[i * 4] = e[0] >>> 0;
            buf[i * 4 + 1] = e[1] >>> 0;
            buf[i * 4 + 2] = (e.length >= 3) ? (e[2] >>> 0) : 0;
            buf[i * 4 + 3] = 0;
        }
        this.device.queue.writeBuffer(this.bondBuf, 0, buf.buffer);
    }

    setLabelsVisible(visible) {
        console.debug('[RawWebGPUAtomsRenderer] setLabelsVisible', { visible, previous: this.labelsVisible });
        this.labelsVisible = visible;
    }

    setBondsVisible(f) { this.bondsVisible = !!f; }
    setHaloVisible(f) { this.haloVisible = !!f; }
    setLabelStyle(rgb, scale) {
        if (rgb) this.labelColor = [rgb[0], rgb[1], rgb[2]];
        if (scale !== undefined && scale !== null) this.labelScale = parseFloat(scale);
        this._uploadLabelUBO();
    }

    _uploadLabelUBO() {
        if (!this.device || !this.labelUBuf) return;
        const aspect = (this.size && this.size.h > 0) ? (this.size.w / this.size.h) : 1.0;
        const data = new Float32Array([
            this.labelColor[0], this.labelColor[1], this.labelColor[2], 1.0,
            this.labelScale, aspect, this.debugForceSolidLabels ? 1.0 : 0.0, this.debugLabelMode | 0
        ]);
        this.device.queue.writeBuffer(this.labelUBuf, 0, data.buffer);
        console.debug('[RawWebGPUAtomsRenderer] label UBO uploaded', { color: this.labelColor, scale: this.labelScale, aspect, debugForceSolid: this.debugForceSolidLabels, debugLabelMode: this.debugLabelMode });
    }
    setDebugForceSolidLabels(flag) {
        const prev = this.debugForceSolidLabels;
        this.debugForceSolidLabels = !!flag;
        console.log('[RawWebGPUAtomsRenderer] setDebugForceSolidLabels', { prev, next: this.debugForceSolidLabels });
        this._uploadLabelUBO();
    }
    setDebugLabelMode(mode) {
        const next = mode | 0;
        const prev = this.debugLabelMode | 0;
        this.debugLabelMode = next;
        console.log('[RawWebGPUAtomsRenderer] setDebugLabelMode', { prev, next });
        console.trace('[RawWebGPUAtomsRenderer] debugLabelMode stack');
        this._uploadLabelUBO();
    }

    logLabelBindingState() {
        console.log('[RawWebGPUAtomsRenderer] label binding state', {
            ubuf: !!this.ubuf,
            posBuf: !!this.posBuf,
            colBuf: !!this.colBuf,
            labelCharBuf: { exists: !!this.labelCharBuf, nLabelChars: this.nLabelChars },
            labelUBuf: !!this.labelUBuf,
            fontTex: !!this.fontTex,
            fontView: !!this.fontView,
            fontSampler: !!this.fontSampler,
            debugForceSolidLabels: this.debugForceSolidLabels,
            debugLabelMode: this.debugLabelMode,
        });
    }

    _refreshLabelBindGroup(reason = 'manual') {
        if (!this.device || !this._bglLabels) return;
        if (!this.ubuf || !this.posBuf || !this.colBuf || !this.labelCharBuf || !this.labelUBuf || !this.fontView || !this.fontSampler) {
            console.warn('[RawWebGPUAtomsRenderer] skip _refreshLabelBindGroup (resources missing)', { reason,
                ubuf: !!this.ubuf, posBuf: !!this.posBuf, colBuf: !!this.colBuf,
                labelCharBuf: !!this.labelCharBuf, labelUBuf: !!this.labelUBuf,
                fontView: !!this.fontView, fontSampler: !!this.fontSampler });
            return;
        }
        this.bindGroupLabels = this.device.createBindGroup({ layout: this._bglLabels, entries: [
            { binding: 0, resource: { buffer: this.ubuf } },
            { binding: 1, resource: { buffer: this.posBuf } },
            { binding: 2, resource: { buffer: this.colBuf } },
            { binding: 3, resource: { buffer: this.labelCharBuf } },
            { binding: 4, resource: { buffer: this.labelUBuf } },
            { binding: 5, resource: this.fontView },
            { binding: 6, resource: this.fontSampler },
        ]});
        console.log('[RawWebGPUAtomsRenderer] label bind group refreshed', { reason, bindGroupLabels: !!this.bindGroupLabels });
    }

    updateLabels(strings, nAtoms) {
        console.log('[RawWebGPUAtomsRenderer] updateLabels begin', {
            labelsVisible: this.labelsVisible,
            nAtoms,
            stringsType: Array.isArray(strings) ? 'array' : (strings === null ? 'null' : typeof strings),
            s0: (strings && strings[0]) ? strings[0] : null,
            s1: (strings && strings[1]) ? strings[1] : null,
            s2: (strings && strings[2]) ? strings[2] : null,
            s3: (strings && strings[3]) ? strings[3] : null,
        });
        if (!strings || nAtoms <= 0) {
            this.nLabelChars = 0;
            return;
        }
        if (!this.labelsVisible) { console.debug('[RawWebGPUAtomsRenderer] updateLabels skipped: labelsVisible=false'); return; }
        const maxChars = 8;
        const advance = 1.1; // in view-space units; multiplied by labelScale in shader
        let ptr = 0;
        const maxTotal = Math.min(this.maxLabelChars, nAtoms * maxChars);
        const buf = new ArrayBuffer(maxTotal * 16);
        const u32 = new Uint32Array(buf);
        const f32 = new Float32Array(buf);
        for (let ia = 0; ia < nAtoms; ia++) {
            const s = strings[ia] || '';
            const len = Math.min(s.length, maxChars);
            const half = (len - 1) * 0.5;
            for (let k = 0; k < maxChars; k++) {
                if (ptr >= maxTotal) break;
                let code = (k < len) ? (s.charCodeAt(k) | 0) : 0;
                if (code >= 32 && code <= 127) {
                    code = (code - 32) | 0; // atlas indices are 0..95 (printable ASCII offset by 32)
                } else {
                    code = 0; // blank cell
                }
                const offx = (k - half) * advance;
                const offy = 0.0;
                const base = ptr * 4;
                u32[base] = ia >>> 0;
                u32[base + 1] = code >>> 0;
                f32[base + 2] = offx;
                f32[base + 3] = offy;
                ptr++;
            }
        }
        this.nLabelChars = ptr;
        console.log('[RawWebGPUAtomsRenderer] updateLabels built', { nAtoms, maxTotal, nLabelChars: this.nLabelChars });
        {
            const nDump = Math.min(this.nLabelChars, 16);
            const dump = [];
            for (let i = 0; i < nDump; i++) {
                const b = i * 4;
                dump.push({ i, atom: u32[b], code: u32[b + 1], offx: f32[b + 2], offy: f32[b + 3] });
            }
            console.log('[RawWebGPUAtomsRenderer] labelChars dump (first)', dump);
            try { console.table ? console.table(dump) : null; } catch (_) {}
            for (const entry of dump) {
                const ascii = (entry.code >= 0 && entry.code <= 95) ? String.fromCharCode(entry.code + 32) : '?';
                const offx = Number.isFinite(entry.offx) ? entry.offx.toFixed(3) : entry.offx;
                const offy = Number.isFinite(entry.offy) ? entry.offy.toFixed(3) : entry.offy;
                console.log(`[RawWebGPUAtomsRenderer] labelChar[${entry.i}] atom=${entry.atom} code=${entry.code} (\"${ascii}\") offx=${offx} offy=${offy}`);
            }
        }
        if (this.nLabelChars > 0) {
            this.device.queue.writeBuffer(this.labelCharBuf, 0, buf);
            console.debug('[RawWebGPUAtomsRenderer] uploaded label chars', { nLabelChars: this.nLabelChars });
        } else {
            console.debug('[RawWebGPUAtomsRenderer] no label chars to upload');
        }
    }

    render() {
        // Ensure swap chain is configured with correct usage before rendering
        this._configureSwapChain(this.canvas.clientWidth || this.canvas.width || 1, this.canvas.clientHeight || this.canvas.height || 1);

        const dev = this.device;
        if (this.lastUncapturedError) {
            console.warn('[RawWebGPUAtomsRenderer] render continuing despite lastUncapturedError', this.lastUncapturedError);
        }

        // Put a validation scope around the whole render so we can see which call trips.
        dev.pushErrorScope('validation');
        console.log('[RawWebGPUAtomsRenderer] render stage=getCurrentTexture begin');
        const ct = this.context.getCurrentTexture();
        console.log('[RawWebGPUAtomsRenderer] render getCurrentTexture', { usage: ct.usage, size: this.size, nAtoms: this.nAtoms, nBonds: this.nBonds, nLabelChars: this.nLabelChars });
        console.log('[RawWebGPUAtomsRenderer] render stage=createView begin');
        const texView = ct.createView();
        console.log('[RawWebGPUAtomsRenderer] render stage=createView end');

        console.log('[RawWebGPUAtomsRenderer] render stage=beginRenderPass begin');
        const encoder = dev.createCommandEncoder();
        const pass = encoder.beginRenderPass({
            colorAttachments: [{
                view: texView,
                loadOp: 'load',
                storeOp: 'store'
            }],
            depthStencilAttachment: {
                view: this.depthView,
                depthLoadOp: 'load',
                depthStoreOp: 'store'
            }
        });
        console.log('[RawWebGPUAtomsRenderer] render stage=beginRenderPass end');

        // 0) Clear color+depth explicitly (avoids Dawn internal clear path requiring COPY_DST)
        if (!this.pipelineClear) throw new Error('pipelineClear is missing');
        pass.setPipeline(this.pipelineClear);
        pass.draw(4, 1, 0, 0);

        // 1) Bonds (behind atoms)
        if (this.bondsVisible && this.pipelineBonds && this.bindGroupBonds && this.nBonds > 0) {
            pass.setPipeline(this.pipelineBonds);
            pass.setBindGroup(0, this.bindGroupBonds);
            pass.draw(4, this.nBonds, 0, 0);
        }

        // 2) Atoms
        if (this.pipelineAtoms && this.bindGroupAtoms && this.nAtoms > 0) {
            pass.setPipeline(this.pipelineAtoms);
            pass.setBindGroup(0, this.bindGroupAtoms);
            pass.draw(4, this.nAtoms, 0, 0);
        }

        // 3) Selection halo (overlay)
        if (this.haloVisible && this.pipelineHalo && this.bindGroupHalo && this.nAtoms > 0) {
            pass.setPipeline(this.pipelineHalo);
            pass.setBindGroup(0, this.bindGroupHalo);
            pass.draw(4, this.nAtoms, 0, 0);
        }

        // 4) Labels (on top)
        if (this.labelsVisible && this.pipelineLabels && this.bindGroupLabels && this.nLabelChars > 0) {
            console.debug('[RawWebGPUAtomsRenderer] draw labels', { nLabelChars: this.nLabelChars, labelsVisible: this.labelsVisible });
            if (this._debugLogLabelBindGroupOnce !== false) {
                this._debugLogLabelBindGroupOnce = false;
                console.log('[RawWebGPUAtomsRenderer] label bind group snapshot', {
                    ubuf: !!this.ubuf,
                    posBuf: !!this.posBuf,
                    colBuf: !!this.colBuf,
                    labelCharBuf: !!this.labelCharBuf,
                    labelUBuf: !!this.labelUBuf,
                    fontView: !!this.fontView,
                    fontSampler: !!this.fontSampler,
                    debugLabelMode: this.debugLabelMode,
                    debugForceSolidLabels: this.debugForceSolidLabels,
                });
            }
            pass.setPipeline(this.pipelineLabels);
            pass.setBindGroup(0, this.bindGroupLabels);
            pass.draw(4, this.nLabelChars, 0, 0);
        } else {
            console.debug('[RawWebGPUAtomsRenderer] skip labels', { labelsVisible: this.labelsVisible, pipelineLabels: !!this.pipelineLabels, bindGroupLabels: !!this.bindGroupLabels, nLabelChars: this.nLabelChars });
        }

        if (this.debugShowLabelAtlas && this.pipelineAtlasDebug && this.bindGroupAtlasDebug) {
            pass.setPipeline(this.pipelineAtlasDebug);
            pass.setBindGroup(0, this.bindGroupAtlasDebug);
            pass.draw(6, 1, 0, 0);
        }

        pass.end();

        console.log('[RawWebGPUAtomsRenderer] render stage=submit begin');
        dev.queue.submit([encoder.finish()]);
        console.log('[RawWebGPUAtomsRenderer] render stage=submit end');

        dev.popErrorScope().then((err) => {
            if (err) {
                console.error('[RawWebGPUAtomsRenderer] validation error during render()', err);
                this.lastUncapturedError = err.message || String(err);
            }
        });
    }
}
