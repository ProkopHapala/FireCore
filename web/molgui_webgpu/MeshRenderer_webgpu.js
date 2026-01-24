import * as THREE from 'three';
import { Draw3D } from './Draw3D_webgpu.js';
import { MeshBasicNodeMaterial, LineBasicNodeMaterial } from 'three/webgpu';
// Import TSL Nodes via namespace to avoid missing named exports
import * as TSL from 'three/tsl';
const {
    float, vec2, vec3, vec4, color,
    positionLocal, uv,
    texture, uniform, attribute,
    instanceIndex,
    mix, normalize, dot, sqrt, max, pow, reflect,
    modelViewMatrix, cameraProjectionMatrix,
    varying
} = TSL;
const Fn = TSL.Fn || TSL.fn; // fallback for differing export names
const If = TSL.If || TSL.if;

export class MeshRenderer {
    constructor(scene, shaders, capacity) {
        this.scene = scene;
        // shaders param is ignored/deprecated, we build materials here
        this.capacity = capacity;

        this.posTexture = null;
        this.posData = null;
        this.texSize = 64;
        this.texHeight = 0;

        this.atomMesh = null;
        this.bondLines = null;
        this.labelMesh = null;
        this.selectionMesh = null;
        this.fontTexture = null;

        // Bonds CPU fallback data (for parity / debugging)
        this._bondPairs = null;
        this._bondCount = 0;

        this.init();
    }

    ensureCapacity(newCapacity) {
        const cap = newCapacity | 0;
        if (cap <= this.capacity) return;

        // Dispose existing
        [this.atomMesh, this.selectionMesh, this.bondLines, this.labelMesh].forEach(obj => {
            if (obj) {
                this.scene.remove(obj);
                if (obj.geometry) obj.geometry.dispose();
                if (obj.material) obj.material.dispose();
            }
        });

        this.posTexture.dispose();
        this.posData = null;

        this.capacity = cap;
        this.init();
    }

    init() {
        // --- 1. Position Texture ---
        this.texHeight = Math.ceil(this.capacity / this.texSize);
        const size = this.texSize * this.texHeight * 4;
        this.posData = new Float32Array(size);

        this.posTexture = new THREE.DataTexture(
            this.posData,
            this.texSize,
            this.texHeight,
            THREE.RGBAFormat,
            THREE.FloatType
        );
        this.posTexture.needsUpdate = true;

        // --- TSL Uniforms ---
        // We wrap Three.js objects in TSL uniforms
        const uPosTexNode = texture(this.posTexture);
        const uTexSizeNode = uniform(new THREE.Vector2(this.texSize, this.texHeight));

        // Helper function to read position from texture based on ID
        const getAtomPos = Fn(([id]) => {
            const tx = id.mod(uTexSizeNode.x).add(0.5).div(uTexSizeNode.x);
            const ty = id.div(uTexSizeNode.x).floor().add(0.5).div(uTexSizeNode.y);
            return uPosTexNode.sample(vec2(tx, ty)); // returns vec4(x,y,z,radius)
        });

        // --- 2. Atoms (Sprite / Impostor) ---
        {
            try {
            // Geometry: Simple Quad
            const atomGeo = new THREE.PlaneGeometry(1, 1);

            // Attributes
            const aAtomID = attribute('aAtomID', 'float');
            const instanceColor = attribute('instanceColor', 'vec3');
            const instanceScale = attribute('instanceScale', 'float');

            // --- Vertex Stage ---
            const atomPosData = getAtomPos(aAtomID);
            const centerPos = atomPosData.xyz;
            const baseRadius = atomPosData.w;

            // Global Uniforms
            const uPointScale = uniform(1.0);
            const uAlpha = uniform(1.0);

            // Billboarding & Scaling
            // In TSL, we modify 'positionLocal' or calculate 'positionWorld'
            // For billboarding instanced mesh:
            const finalRadius = instanceScale.mul(baseRadius).mul(uPointScale);

            // Create Material
            // We use NodeMaterial to construct the pipeline
            const material = new MeshBasicNodeMaterial();
            material.transparent = true;
            material.side = THREE.DoubleSide;

            // Vertex Node
            // Billboard logic: View-aligned offset added to center position
            // ProjectionMatrix * (ModelViewMatrix * Center + Offset)
            material.vertexNode = Fn(() => {
                const viewPos = modelViewMatrix.mul(vec4(centerPos, 1.0));
                viewPos.xy.addAssign(positionLocal.xy.mul(finalRadius));
                // Transform to Clip Space
                return cameraProjectionMatrix.mul(viewPos);
            })();

            // --- Fragment Stage (Impostor Sphere) ---
            const vUv = uv();

            material.colorNode = Fn(() => {
                const uvOffset = vUv.sub(0.5);
                const distSq = dot(uvOffset, uvOffset);

                // Outside circle -> fully transparent
                const mask = distSq.lessThan(0.25);

                // Normal calc (Impostor)
                const z = sqrt(max(float(0.25).sub(distSq), float(0.0)));
                const normal = normalize(vec3(uvOffset.x, uvOffset.y, z));

                // Simple Lighting
                const lightDir = normalize(vec3(0.5, 0.5, 1.0));
                const diffuse = max(dot(normal, lightDir), 0.0);
                const ambient = vec3(0.3);

                // Specular
                const viewDir = vec3(0.0, 0.0, 1.0);
                const reflectDir = reflect(lightDir.negate(), normal);
                const spec = pow(max(dot(viewDir, reflectDir), 0.0), 32.0);
                const specular = vec3(0.5).mul(spec);

                const finalColor = instanceColor.mul(ambient.add(diffuse)).add(specular);
                const alpha = uAlpha.mul(mask.select(1.0, 0.0));
                return vec4(finalColor, alpha);
            })();

            // Create Mesh
            this.atomMesh = new THREE.InstancedMesh(atomGeo, material, this.capacity);
            // Setup Attributes
            const ids = new Float32Array(this.capacity);
            const cols = new Float32Array(this.capacity * 3);
            const scales = new Float32Array(this.capacity);
            for (let i = 0; i < this.capacity; i++) { ids[i] = i; scales[i] = 0; cols[i * 3] = 1; cols[i * 3 + 1] = 1; cols[i * 3 + 2] = 1; }

            this.atomMesh.geometry.setAttribute('aAtomID', new THREE.InstancedBufferAttribute(ids, 1));
            this.atomMesh.geometry.setAttribute('instanceColor', new THREE.InstancedBufferAttribute(cols, 3));
            this.atomMesh.geometry.setAttribute('instanceScale', new THREE.InstancedBufferAttribute(scales, 1));
            this.atomMesh.count = 0;
            this.scene.add(this.atomMesh);

            // --- 2b. Selection Overlay ---
            // Same impostor as atoms, but slightly larger and translucent.
            {
                const selMat = new MeshBasicNodeMaterial();
                selMat.transparent = true;
                selMat.side = THREE.DoubleSide;
                selMat.depthWrite = false;
                selMat.depthTest = true;

                const uSelPointScale = uniform(1.08);
                const uSelAlpha = uniform(0.22);

                selMat.vertexNode = Fn(() => {
                    const atomPosData2 = getAtomPos(aAtomID);
                    const centerPos2 = atomPosData2.xyz;
                    const baseRadius2 = atomPosData2.w;
                    const finalRadius2 = instanceScale.mul(baseRadius2).mul(uSelPointScale);
                    const viewPos2 = modelViewMatrix.mul(vec4(centerPos2, 1.0));
                    viewPos2.xy.addAssign(positionLocal.xy.mul(finalRadius2));
                    return cameraProjectionMatrix.mul(viewPos2);
                })();

                selMat.colorNode = Fn(() => {
                    const vUv2 = uv();
                    const uvOffset2 = vUv2.sub(0.5);
                    const distSq2 = dot(uvOffset2, uvOffset2);
                    const mask2 = distSq2.lessThan(0.25);
                    const alpha2 = uSelAlpha.mul(mask2.select(1.0, 0.0));
                    return vec4(instanceColor, alpha2);
                })();

                this.selectionMesh = new THREE.InstancedMesh(atomGeo.clone(), selMat, this.capacity);
                // We reuse attribute arrays where it is safe (ids), but keep own colors/scales.
                this.selectionMesh.geometry.setAttribute('aAtomID', new THREE.InstancedBufferAttribute(ids, 1));
                this.selectionMesh.geometry.setAttribute('instanceColor', new THREE.InstancedBufferAttribute(new Float32Array(this.capacity * 3), 3));
                this.selectionMesh.geometry.setAttribute('instanceScale', new THREE.InstancedBufferAttribute(new Float32Array(this.capacity), 1));
                this.selectionMesh.count = 0;
                this.selectionMesh.renderOrder = 998;
                this.scene.add(this.selectionMesh);
            }
            } catch (e) {
                console.error('Atom material build failed', e);
                if (window.logger) window.logger.error(`[MeshRenderer] atom material failed: ${e?.message||e}`);
                throw e;
            }
        }

        // --- 3. Bonds (Line Segments) ---
        // CPU-updated line buffer as a safe fallback (works in WebGL and WebGPU backends).
        // Later we can optimize this back to an ID+texture driven shader.
        {
            const maxBonds = this.capacity * 4;
            const geom = new THREE.BufferGeometry();
            geom.setAttribute('position', new THREE.BufferAttribute(new Float32Array(maxBonds * 2 * 3), 3));
            // Keep these attributes for future shader parity (currently unused by LineBasicMaterial)
            geom.setAttribute('aAtomID', new THREE.BufferAttribute(new Float32Array(maxBonds * 2), 1));
            geom.setAttribute('aMatID', new THREE.BufferAttribute(new Float32Array(maxBonds * 2), 1));
            const mat = new THREE.LineBasicMaterial({ color: 0xffffff, transparent: true, opacity: 0.9 });
            this.bondLines = new THREE.LineSegments(geom, mat);
            this.bondLines.frustumCulled = false;
            this.bondLines.geometry.setDrawRange(0, 0);
            this.scene.add(this.bondLines);
        }

        // --- 4. Labels ---
        // Porting the font atlas logic to TSL.
        // We reuse Draw3D_webgpu geometry/attribute packing but provide a real NodeMaterial here.
        {
            this.fontTexture = Draw3D.createFontTexture();
            this.labelMesh = Draw3D.createLabelInstancedMesh(this.capacity, null, this.fontTexture, null);

            const uFontTexNode = texture(this.fontTexture);
            const uFontGrid = uniform(new THREE.Vector2(16, 16));
            const uLabelScale = uniform(0.5);
            const uLabelColor = uniform(new THREE.Vector3(1.0, 1.0, 1.0));
            const uAspect = uniform(1.0);

            // Expose for external styling / resize updates
            this._labelUniforms = { uLabelScale, uLabelColor, uAspect, uFontGrid };

            const aAtomID_L = attribute('aAtomID', 'float');
            const aLabel1 = attribute('aLabel1', 'vec4');
            const aLabel2 = attribute('aLabel2', 'vec4');
            const aStrLen = attribute('aStrLen', 'float');
            const aCharPos = attribute('aCharPos', 'float');

            // Choose packed ASCII code for this character quad (0..7)
            const getCharCode = Fn(([cpos]) => {
                // Nested selects (avoid try/catch; fail loudly if TSL changes semantics)
                const c0 = aLabel1.x;
                const c1 = aLabel1.y;
                const c2 = aLabel1.z;
                const c3 = aLabel1.w;
                const c4 = aLabel2.x;
                const c5 = aLabel2.y;
                const c6 = aLabel2.z;
                const c7 = aLabel2.w;
                let code = c7;
                code = cpos.lessThan(6.5).select(c6, code);
                code = cpos.lessThan(5.5).select(c5, code);
                code = cpos.lessThan(4.5).select(c4, code);
                code = cpos.lessThan(3.5).select(c3, code);
                code = cpos.lessThan(2.5).select(c2, code);
                code = cpos.lessThan(1.5).select(c1, code);
                code = cpos.lessThan(0.5).select(c0, code);
                return code;
            });

            const labelMat = new MeshBasicNodeMaterial();
            labelMat.transparent = true;
            labelMat.depthTest = false;
            labelMat.depthWrite = false;
            labelMat.side = THREE.DoubleSide;

            labelMat.vertexNode = Fn(() => {
                const atomPosDataL = getAtomPos(aAtomID_L);
                const centerPosL = atomPosDataL.xyz;

                // Center label by length
                const len = aStrLen;
                const halfSpan = len.sub(1.0).mul(0.5);
                const cx = aCharPos.sub(halfSpan);

                // per-character quad offset in view space (billboard)
                const charAdvance = float(1.1); // spacing
                const offs = vec2(cx.mul(charAdvance), float(0.0));
                const viewPosL = modelViewMatrix.mul(vec4(centerPosL, 1.0));
                // Apply aspect correction only in X so text does not stretch with viewport
                viewPosL.x.addAssign(offs.x.mul(uLabelScale).div(uAspect));
                viewPosL.y.addAssign(offs.y.mul(uLabelScale));
                viewPosL.xy.addAssign(positionLocal.xy.mul(uLabelScale));
                return cameraProjectionMatrix.mul(viewPosL);
            })();

            labelMat.colorNode = Fn(() => {
                const cpos = aCharPos;
                const len = aStrLen;
                const visible = cpos.lessThan(len);
                const code = getCharCode(cpos);

                // If string has no character here, make it transparent
                const active = visible.select(1.0, 0.0);

                // Map ASCII code into 16x16 atlas
                const cols = uFontGrid.x;
                const rows = uFontGrid.y;
                const col = code.mod(cols);
                const row = code.div(cols).floor();

                const cell = vec2(float(1.0).div(cols), float(1.0).div(rows));
                const base = vec2(col.mul(cell.x), row.mul(cell.y));

                const uvLocal = uv();
                const fuv = base.add(vec2(uvLocal.x.mul(cell.x), uvLocal.y.mul(cell.y)));
                const texel = uFontTexNode.sample(fuv);
                // Alpha from atlas (white on transparent)
                const a = texel.r.mul(active);
                return vec4(vec3(uLabelColor.x, uLabelColor.y, uLabelColor.z), a);
            })();

            this.labelMesh.material = labelMat;
            this.labelMesh.visible = false;
            this.scene.add(this.labelMesh);
        }
    }

    // ... Keep updatePositions, updateParticles, updateBonds, updateSelection, etc. ...
    // ... Copy them from previous MeshRenderer.js, they work on Data Buffers (CPU side) ...
    // ... Only change is remove references to `this.shaders` ...

    updatePositions(posArray, count) {
        if ((count | 0) > (this.capacity | 0)) throw new Error(`Capacity exceeded`);
        for (let i = 0; i < count; i++) {
            this.posData[i * 4] = posArray[i * 3];
            this.posData[i * 4 + 1] = posArray[i * 3 + 1];
            this.posData[i * 4 + 2] = posArray[i * 3 + 2];
            this.posData[i * 4 + 3] = 1.0;
        }
        this.posTexture.needsUpdate = true;

        // Bonds CPU fallback: update line vertex positions whenever atom positions change.
        if (this.bondLines && this._bondPairs && this._bondCount > 0) {
            this._updateBondLinePositions();
        }
    }

    updateParticles(count, colorGetter, scaleGetter) {
        if (!this.atomMesh) return;
        this.atomMesh.count = count;
        const colorAttr = this.atomMesh.geometry.getAttribute('instanceColor');
        const scaleAttr = this.atomMesh.geometry.getAttribute('instanceScale');
        for (let i = 0; i < count; i++) {
            const col = colorGetter(i);
            const scale = scaleGetter(i);
            colorAttr.setXYZ(i, col[0], col[1], col[2]);
            scaleAttr.setX(i, scale);
        }
        colorAttr.needsUpdate = true;
        scaleAttr.needsUpdate = true;
    }

    updateBonds(pairs) {
        if (!this.bondLines) return;

        this._bondPairs = pairs || [];
        this._bondCount = this._bondPairs.length | 0;

        const bondIDAttr = this.bondLines.geometry.getAttribute('aAtomID');
        const matIDAttr = this.bondLines.geometry.getAttribute('aMatID');
        let ptr = 0;
        for (let i = 0; i < this._bondCount; i++) {
            const entry = this._bondPairs[i];
            const matID = (entry.length >= 3) ? entry[2] : 0;
            bondIDAttr.setX(ptr, entry[0]); matIDAttr.setX(ptr, matID); ptr++;
            bondIDAttr.setX(ptr, entry[1]); matIDAttr.setX(ptr, matID); ptr++;
        }
        bondIDAttr.needsUpdate = true;
        matIDAttr.needsUpdate = true;

        this.bondLines.geometry.setDrawRange(0, this._bondCount * 2);
        this._updateBondLinePositions();

        if ((window.VERBOSITY_LEVEL | 0) >= 4) {
            console.debug(`[MeshRenderer_webgpu.updateBonds] nBonds=${this._bondCount}`);
        }
    }

    updateSelection(indices) {
        if (!this.selectionMesh) return;
        const maxCount = Math.min(indices.length, this.capacity);
        const idAttr = this.selectionMesh.geometry.getAttribute('aAtomID');
        const scaleAttr = this.selectionMesh.geometry.getAttribute('instanceScale');
        const colAttr = this.selectionMesh.geometry.getAttribute('instanceColor');

        // Hide all
        for (let i = 0; i < this.capacity; i++) scaleAttr.array[i] = 0.0;

        for (let i = 0; i < maxCount; i++) {
            idAttr.setX(i, indices[i]);
            scaleAttr.setX(i, 1.2);
            colAttr.setXYZ(i, 0.95, 0.95, 0.95);
        }
        idAttr.needsUpdate = true;
        scaleAttr.needsUpdate = true;
        colAttr.needsUpdate = true;
        this.selectionMesh.count = maxCount;

        if ((window.VERBOSITY_LEVEL | 0) >= 4) {
            console.debug(`[MeshRenderer_webgpu.updateSelection] nSelected=${maxCount}`);
        }
    }

    updateLabels(stringGetter, count) {
        if (!this.labelMesh || !this.labelMesh.visible) return;
        Draw3D.updateLabelBuffers(this.labelMesh, stringGetter, count);

        if ((window.VERBOSITY_LEVEL | 0) >= 4) {
            console.debug(`[MeshRenderer_webgpu.updateLabels] nLabels=${count}`);
        }
    }

    _updateBondLinePositions() {
        if (!this.bondLines) return;
        if (!this._bondPairs || (this._bondCount | 0) <= 0) return;

        const posAttr = this.bondLines.geometry.getAttribute('position');
        const arr = posAttr.array;
        let k = 0;
        for (let i = 0; i < this._bondCount; i++) {
            const e = this._bondPairs[i];
            const i0 = e[0] | 0;
            const i1 = e[1] | 0;
            const o0 = i0 * 4;
            const o1 = i1 * 4;
            // posData is RGBA-packed: xyzw
            arr[k++] = this.posData[o0];
            arr[k++] = this.posData[o0 + 1];
            arr[k++] = this.posData[o0 + 2];
            arr[k++] = this.posData[o1];
            arr[k++] = this.posData[o1 + 1];
            arr[k++] = this.posData[o1 + 2];
        }
        posAttr.needsUpdate = true;
        this.bondLines.geometry.computeBoundingSphere();
    }

    // Toggles
    setAtomsVisible(f) { if (this.atomMesh) this.atomMesh.visible = !!f; }
    setBondsVisible(f) { if (this.bondLines) this.bondLines.visible = !!f; }
    setLabelsVisible(f) { if (this.labelMesh) this.labelMesh.visible = !!f; }
    setSelectionVisible(f) { if (this.selectionMesh) this.selectionMesh.visible = !!f; }
}