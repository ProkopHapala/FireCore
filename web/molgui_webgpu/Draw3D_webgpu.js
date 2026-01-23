import * as THREE from 'three';

export class Draw3D {
    static drawSphereOctLines(n, R, pos, vertices) {
        const a = { x: 1, y: 0, z: 0 };
        const b = { x: 0, y: 1, z: 0 };
        const c = { x: 0, y: 0, z: 1 };

        this.drawCircleAxis(n, pos, a, c, R, vertices);
        this.drawCircleAxis(n, pos, b, a, R, vertices);
        this.drawCircleAxis(n, pos, c, b, R, vertices);
    }

    static drawCircleAxis(n, pos, v0, uaxis, R, vertices) {
        const dphi = 2 * Math.PI / n;
        const dca = Math.cos(dphi);
        const dsa = Math.sin(dphi);

        let x = v0.x;
        let y = v0.y;
        let z = v0.z;

        let px = pos.x + x * R;
        let py = pos.y + y * R;
        let pz = pos.z + z * R;

        for (let i = 0; i < n; i++) {
            const cx = uaxis.y * z - uaxis.z * y;
            const cy = uaxis.z * x - uaxis.x * z;
            const cz = uaxis.x * y - uaxis.y * x;

            const nx = x * dca + cx * dsa;
            const ny = y * dca + cy * dsa;
            const nz = z * dca + cz * dsa;

            const p2x = pos.x + nx * R;
            const p2y = pos.y + ny * R;
            const p2z = pos.z + nz * R;

            vertices.push(px, py, pz);
            vertices.push(p2x, p2y, p2z);

            x = nx;
            y = ny;
            z = nz;
            px = p2x;
            py = p2y;
            pz = p2z;
        }
    }

    static createOctSphereGeometry(n, R) {
        const vertices = [];
        this.drawSphereOctLines(n, R, { x: 0, y: 0, z: 0 }, vertices);
        const geometry = new THREE.BufferGeometry();
        geometry.setAttribute('position', new THREE.Float32BufferAttribute(vertices, 3));
        return geometry;
    }

    static createFontTexture() {
        const canvas = document.createElement('canvas');
        const size = 512;
        canvas.width = size;
        canvas.height = size;
        const ctx = canvas.getContext('2d');

        const cols = 16;
        const rows = 16;
        const charW = size / cols;
        const charH = size / rows;

        ctx.fillStyle = '#00000000';
        ctx.fillRect(0, 0, size, size);

        ctx.font = 'bold 24px monospace';
        ctx.fillStyle = '#FFFFFF';
        ctx.textAlign = 'center';
        ctx.textBaseline = 'middle';

        for (let i = 32; i < 127; i++) {
            const col = i % cols;
            const row = Math.floor(i / cols);

            const x = col * charW + charW / 2;
            const y = row * charH + charH / 2;

            ctx.fillText(String.fromCharCode(i), x, y + 2);
        }

        const texture = new THREE.CanvasTexture(canvas);
        texture.minFilter = THREE.LinearFilter;
        texture.magFilter = THREE.LinearFilter;
        texture.needsUpdate = true;
        return texture;
    }

    static createLabelInstancedMesh(capacity, shaders, fontTexture, uniforms) {
        const maxChars = 8;
        const baseGeo = new THREE.PlaneGeometry(1, 1);
        const basePos = baseGeo.attributes.position.array;
        const baseUv = baseGeo.attributes.uv.array;
        const baseIndex = baseGeo.index.array;

        const vertices = [];
        const uvs = [];
        const indices = [];
        const charPosAttr = [];

        for (let c = 0; c < maxChars; c++) {
            const vOffset = (vertices.length / 3);

            for (let i = 0; i < basePos.length; i += 3) {
                vertices.push(basePos[i], basePos[i + 1], basePos[i + 2]);
                charPosAttr.push(c);
            }

            for (let i = 0; i < baseUv.length; i += 2) {
                uvs.push(baseUv[i], baseUv[i + 1]);
            }

            for (let i = 0; i < baseIndex.length; i++) {
                indices.push(baseIndex[i] + vOffset);
            }
        }

        const geometry = new THREE.BufferGeometry();
        geometry.setAttribute('position', new THREE.Float32BufferAttribute(vertices, 3));
        geometry.setAttribute('uv', new THREE.Float32BufferAttribute(uvs, 2));
        geometry.setAttribute('aCharPos', new THREE.Float32BufferAttribute(charPosAttr, 1));
        geometry.setIndex(indices);

        const instID = new Float32Array(capacity);
        for (let i = 0; i < capacity; i++) instID[i] = i;

        const instLabel1 = new Float32Array(capacity * 4);
        const instLabel2 = new Float32Array(capacity * 4);
        const instStrLen = new Float32Array(capacity);

        const instGeo = new THREE.InstancedBufferGeometry();
        instGeo.copy(geometry);
        instGeo.setAttribute('aAtomID', new THREE.InstancedBufferAttribute(instID, 1));
        instGeo.setAttribute('aLabel1', new THREE.InstancedBufferAttribute(instLabel1, 4));
        instGeo.setAttribute('aLabel2', new THREE.InstancedBufferAttribute(instLabel2, 4));
        instGeo.setAttribute('aStrLen', new THREE.InstancedBufferAttribute(instStrLen, 1));

        const material = new THREE.MeshBasicMaterial({ visible: false });

        const mesh = new THREE.InstancedMesh(instGeo, material, capacity);
        mesh.frustumCulled = false;
        mesh.renderOrder = 999;
        mesh.count = 0;
        return mesh;
    }

    static updateLabelBuffers(mesh, stringGetter, count) {
        mesh.count = count;

        const attr1 = mesh.geometry.getAttribute('aLabel1');
        const attr2 = mesh.geometry.getAttribute('aLabel2');
        const attrLen = mesh.geometry.getAttribute('aStrLen');
        const arr1 = attr1.array;
        const arr2 = attr2.array;
        const arrLen = attrLen.array;

        for (let i = 0; i < count; i++) {
            const str = stringGetter(i) || "";
            const len = Math.min(str.length, 8);
            arrLen[i] = len;

            for (let k = 0; k < 8; k++) {
                if (k < 4) arr1[i * 4 + k] = 0;
                else arr2[i * 4 + (k - 4)] = 0;
            }

            for (let k = 0; k < len; k++) {
                const code = str.charCodeAt(k);
                if (k < 4) arr1[i * 4 + k] = code;
                else arr2[i * 4 + (k - 4)] = code;
            }
        }

        attr1.needsUpdate = true;
        attr2.needsUpdate = true;
        attrLen.needsUpdate = true;
    }
}