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
            for (let i = 0; i < this.capacity; i++) { ids[i] = i; scales[i] = 0; cols[i * 3] = 1; }

            this.atomMesh.geometry.setAttribute('aAtomID', new THREE.InstancedBufferAttribute(ids, 1));
            this.atomMesh.geometry.setAttribute('instanceColor', new THREE.InstancedBufferAttribute(cols, 3));
            this.atomMesh.geometry.setAttribute('instanceScale', new THREE.InstancedBufferAttribute(scales, 1));
            this.atomMesh.count = 0;
            this.scene.add(this.atomMesh);
            } catch (e) {
                console.error('Atom material build failed', e);
                if (window.logger) window.logger.error(`[MeshRenderer] atom material failed: ${e?.message||e}`);
                throw e;
            }
        }

        // --- 3. Bonds (Line Segments) ---
        // TEMP disable bonds to isolate WebGPU errors
        this.bondLines = null;

        // --- 4. Labels ---
        // Porting the font atlas logic to TSL
        // TEMP disable labels to isolate WebGPU errors
        this.labelMesh = null;
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
        const bondIDAttr = this.bondLines.geometry.getAttribute('aAtomID');
        const matIDAttr = this.bondLines.geometry.getAttribute('aMatID');
        let ptr = 0;
        for (let i = 0; i < pairs.length; i++) {
            const entry = pairs[i];
            const matID = (entry.length >= 3) ? entry[2] : 0;
            bondIDAttr.setX(ptr, entry[0]); matIDAttr.setX(ptr, matID); ptr++;
            bondIDAttr.setX(ptr, entry[1]); matIDAttr.setX(ptr, matID); ptr++;
        }
        bondIDAttr.needsUpdate = true;
        matIDAttr.needsUpdate = true;
        this.bondLines.geometry.setDrawRange(0, pairs.length * 2);
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
    }

    updateLabels(stringGetter, count) {
        if (!this.labelMesh || !this.labelMesh.visible) return;
        Draw3D.updateLabelBuffers(this.labelMesh, stringGetter, count);
    }

    // Toggles
    setAtomsVisible(f) { if (this.atomMesh) this.atomMesh.visible = !!f; }
    setBondsVisible(f) { if (this.bondLines) this.bondLines.visible = !!f; }
    setLabelsVisible(f) { if (this.labelMesh) this.labelMesh.visible = !!f; }
    setSelectionVisible(f) { if (this.selectionMesh) this.selectionMesh.visible = !!f; }
}