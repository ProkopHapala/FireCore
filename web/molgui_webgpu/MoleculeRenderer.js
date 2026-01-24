import * as THREE from 'three';
import { MeshRenderer } from './MeshRenderer_webgpu.js';
import { Draw3D } from './Draw3D_webgpu.js';
import { Logger } from '../common_js/Logger.js';
import { Vec3 } from '../common_js/Vec3.js';

export class PackedMolecule {
    constructor(capacity = 1024) {
        this.capacity = capacity | 0;
        this.nAtoms = 0;
        this.pos = new Float32Array(this.capacity * 3);
        this.types = new Uint8Array(this.capacity);
        this.bonds = [];
        this.selection = new Set();
        this.atomIds = new Int32Array(this.capacity);
        this.isDirty = true;
    }

    clear() {
        this.nAtoms = 0;
        this.bonds = [];
        this.selection = new Set();
        // keep buffers allocated; mark dirty so renderer uploads empty state
        this.isDirty = true;
    }

    resize(newCapacity) {
        newCapacity = newCapacity | 0;
        if (newCapacity <= (this.capacity | 0)) return;
        const newPos = new Float32Array(newCapacity * 3);
        const newTypes = new Uint8Array(newCapacity);
        const newAtomIds = new Int32Array(newCapacity);
        newPos.set(this.pos);
        newTypes.set(this.types);
        newAtomIds.set(this.atomIds);
        this.pos = newPos;
        this.types = newTypes;
        this.atomIds = newAtomIds;
        this.capacity = newCapacity;
        this.isDirty = true;
    }

    setBonds(bonds) {
        this.bonds = bonds || [];
        this.isDirty = true;
    }
}

export class BucketOverlayRenderer {
    constructor(scene) {
        this.scene = scene;
        this.group = new THREE.Group();
        this.scene.add(this.group);

        // --- Bucket overlay (boxes) ---
        const geomBox = new THREE.BufferGeometry();
        geomBox.setAttribute('position', new THREE.BufferAttribute(new Float32Array(0), 3));
        const matBox = new THREE.LineBasicMaterial({ color: 0xff8844, transparent: true, opacity: 0.6 });
        this.bucketOverlay = new THREE.LineSegments(geomBox, matBox);
        this.bucketOverlay.renderOrder = 12;
        this.bucketOverlay.visible = false;
        this.group.add(this.bucketOverlay);

        // --- Bucket atom->center lines ---
        const geomLines = new THREE.BufferGeometry();
        geomLines.setAttribute('position', new THREE.BufferAttribute(new Float32Array(0), 3));
        const matLines = new THREE.LineBasicMaterial({ color: 0x55ccff, transparent: true, opacity: 0.6 });
        this.bucketAtomLines = new THREE.LineSegments(geomLines, matLines);
        this.bucketAtomLines.renderOrder = 13;
        this.bucketAtomLines.visible = false;
        this.group.add(this.bucketAtomLines);
    }

    clear() {
        this.bucketOverlay.visible = false;
        this.bucketAtomLines.visible = false;
        this.bucketOverlay.geometry.setAttribute('position', new THREE.BufferAttribute(new Float32Array(0), 3));
        this.bucketAtomLines.geometry.setAttribute('position', new THREE.BufferAttribute(new Float32Array(0), 3));
    }

    updateBuckets(bg, mol, showBoxes, showLines) {
        if (!bg || !bg.buckets || bg.buckets.length === 0) {
            this.clear();
            return;
        }

        // 1. Update Boxes
        if (!showBoxes) {
            this.bucketOverlay.visible = false;
        } else {
            this.bucketOverlay.visible = true;
            let n3 = 0;
            let verts = new Float32Array(0);
            if (bg.meta && bg.meta.kind === 'crystal_cells' && bg.meta.lvec && bg.meta.origin) {
                const A = bg.meta.lvec[0], B = bg.meta.lvec[1], C = bg.meta.lvec[2];
                const o0 = bg.meta.origin;
                const n = bg.buckets.length;
                verts = new Float32Array(n * 12 * 2 * 3);
                for (let ib = 0; ib < n; ib++) {
                    const g = bg.buckets[ib].gridPos;
                    if (!g) continue;
                    const origin = o0.clone().addMul(A, g.x).addMul(B, g.y).addMul(C, g.z);
                    const bverts = this._buildWireframeCellVerts(A, B, C, origin);
                    verts.set(bverts, n3);
                    n3 += bverts.length;
                }
            } else {
                const n = bg.buckets.length;
                verts = new Float32Array(n * 12 * 2 * 3);
                for (let ib = 0; ib < n; ib++) {
                    const b = bg.buckets[ib];
                    const lo = b.min, hi = b.max;
                    if (!lo || !hi) continue;
                    const bverts = this._buildWireframeAABBVerts(lo, hi);
                    verts.set(bverts, n3);
                    n3 += bverts.length;
                }
            }
            const vFinal = (n3 === verts.length) ? verts : verts.subarray(0, n3);
            this.bucketOverlay.geometry.setAttribute('position', new THREE.BufferAttribute(vFinal, 3));
            this.bucketOverlay.geometry.computeBoundingSphere();
        }

        // 2. Update Atom Lines
        if (!showLines) {
            this.bucketAtomLines.visible = false;
        } else {
            this.bucketAtomLines.visible = true;
            let nPairs = 0;
            for (let ib = 0; ib < bg.buckets.length; ib++) nPairs += (bg.buckets[ib].atoms.length | 0);
            const verts = new Float32Array(nPairs * 2 * 3);
            let k = 0;
            const c = new Vec3();
            for (let ib = 0; ib < bg.buckets.length; ib++) {
                bg.getBucketCenterFromBounds(ib, c);
                const bs = bg.buckets[ib].atoms;
                for (let i = 0; i < bs.length; i++) {
                    const ia = bs[i] | 0;
                    const a = mol.atoms[ia];
                    if (!a) continue;
                    const p = a.pos;
                    verts[k++] = p.x; verts[k++] = p.y; verts[k++] = p.z;
                    verts[k++] = c.x; verts[k++] = c.y; verts[k++] = c.z;
                }
            }
            const vFinal = (k === verts.length) ? verts : verts.subarray(0, k);
            this.bucketAtomLines.geometry.setAttribute('position', new THREE.BufferAttribute(vFinal, 3));
            this.bucketAtomLines.geometry.computeBoundingSphere();
        }
    }

    _buildWireframeCellVerts(a, b, c, origin) {
        const v = new Float32Array(24 * 3);
        const corners = [
            origin.clone(), origin.clone().add(a), origin.clone().add(b), origin.clone().add(a).add(b),
            origin.clone().add(c), origin.clone().add(a).add(c), origin.clone().add(b).add(c), origin.clone().add(a).add(b).add(c)
        ];
        const edges = [[0, 1], [0, 2], [0, 4], [1, 3], [1, 5], [2, 3], [2, 6], [3, 7], [4, 5], [4, 6], [5, 7], [6, 7]];
        for (let i = 0; i < edges.length; i++) {
            const [i0, i1] = edges[i];
            const p0 = corners[i0], p1 = corners[i1];
            v[i * 6] = p0.x; v[i * 6 + 1] = p0.y; v[i * 6 + 2] = p0.z;
            v[i * 6 + 3] = p1.x; v[i * 6 + 4] = p1.y; v[i * 6 + 5] = p1.z;
        }
        return v;
    }

    _buildWireframeAABBVerts(lo, hi) {
        const v = new Float32Array(24 * 3);
        const c = [
            new Vec3(lo.x, lo.y, lo.z), new Vec3(hi.x, lo.y, lo.z), new Vec3(lo.x, hi.y, lo.z), new Vec3(hi.x, hi.y, lo.z),
            new Vec3(lo.x, lo.y, hi.z), new Vec3(hi.x, lo.y, hi.z), new Vec3(lo.x, hi.y, hi.z), new Vec3(hi.x, hi.y, hi.z)
        ];
        const edges = [[0, 1], [0, 2], [0, 4], [1, 3], [1, 5], [2, 3], [2, 6], [3, 7], [4, 5], [4, 6], [5, 7], [6, 7]];
        for (let i = 0; i < edges.length; i++) {
            const [i0, i1] = edges[i];
            const p0 = c[i0], p1 = c[i1];
            v[i * 6] = p0.x; v[i * 6 + 1] = p0.y; v[i * 6 + 2] = p0.z;
            v[i * 6 + 3] = p1.x; v[i * 6 + 4] = p1.y; v[i * 6 + 5] = p1.z;
        }
        return v;
    }
}

export class MoleculeRenderer extends MeshRenderer {
    /**
     * @param {THREE.Scene} scene
     * @param {PackedMolecule} system - packed buffer backing the renderer (positions/types/bonds/selection)
     * @param {*} shaders
     * @param {*} mmParams
     * @param {EditableMolecule|null} source - authoritative editable molecule; if provided, export on update()
     */
    constructor(scene, system, shaders, mmParams, source = null) {
        super(scene, shaders, system.capacity);
        this.system = system;
        this.source = source;
        this.mmParams = mmParams;

        this.axesHelper = null;
        this.labelMode = 'none';

        // selection is handled by MeshRenderer

        // --- Replicas / Lattice ---
        this.replicaGroup = new THREE.Group();
        this.scene.add(this.replicaGroup);
        this.replicaBox = null;
        this.replicaConfig = {
            lvec: [new Vec3(10, 0, 0), new Vec3(0, 10, 0), new Vec3(0, 0, 10)],
            nrep: { x: 0, y: 0, z: 0 },
            show: false,
            showBox: false,
            filter: null,
            dirty: false
        };
        this.replicaClones = [];

        // --- Buckets ---
        this.bucketRenderer = new BucketOverlayRenderer(this.scene);

        // --- Angle Constraint Debug Visualization ---
        this.angleConstraintLines = null;
        this.showAngleConstraints = false;
        this.pdSimulation = null; // Reference to PDSimulation for debug data
    }

    clear() {
        // Clear rendered buffers (PackedMolecule) and mark source dirty so export writes empty state.
        if (this.system && this.system.clear) this.system.clear(); // PackedMolecule
        if (this.source) { // EditableMolecule
            this.source.dirtyTopo = true;
            this.source.dirtyGeom = true;
            this.source.dirtyExport = true;
        }
        this.clearReplicas();
        if (this.bucketRenderer) this.bucketRenderer.clear();
        this.update();
    }

    setReplicas(cfg) {
        if (window.logger && window.logger.debug) {
            window.logger.debug(`[MoleculeRenderer.setReplicas] show=${cfg.show} showBox=${cfg.showBox} nrep=${cfg.nrep ? JSON.stringify(cfg.nrep) : 'null'}`);
        } else {
            console.debug('[MoleculeRenderer.setReplicas]', cfg);
        }
        // If lattice vectors not provided, fall back to source EditableMolecule.lvec; if still missing, disable replicas.
        if (cfg.lvec) {
            for (let i = 0; i < 3; i++) {
                if (cfg.lvec[i]) this.replicaConfig.lvec[i].setV(cfg.lvec[i]);
            }
        } else if (this.source && this.source.lvec) {
            for (let i = 0; i < 3; i++) {
                this.replicaConfig.lvec[i].setV(this.source.lvec[i]);
            }
        } else {
            this.replicaConfig.show = false;
        }
        if (cfg.nrep) {
            this.replicaConfig.nrep.x = cfg.nrep.x !== undefined ? cfg.nrep.x : this.replicaConfig.nrep.x;
            this.replicaConfig.nrep.y = cfg.nrep.y !== undefined ? cfg.nrep.y : this.replicaConfig.nrep.y;
            this.replicaConfig.nrep.z = cfg.nrep.z !== undefined ? cfg.nrep.z : this.replicaConfig.nrep.z;
        }
        if (cfg.show !== undefined) this.replicaConfig.show = cfg.show;
        if (cfg.showBox !== undefined) this.replicaConfig.showBox = cfg.showBox;
        if (cfg.filter !== undefined) this.replicaConfig.filter = cfg.filter;
        this.replicaConfig.dirty = true;
    }

    clearReplicas() {
        this.replicaGroup.clear();
        if (this.replicaBox) {
            this.replicaGroup.remove(this.replicaBox);
            if (this.replicaBox.geometry) this.replicaBox.geometry.dispose();
            if (this.replicaBox.material) this.replicaBox.material.dispose();
            this.replicaBox = null;
        }
        this.replicaClones = [];
        this.replicaConfig.dirty = false;
        this.replicaConfig.show = false;
        this.replicaConfig.showBox = false;
    }

    rebuildReplicas() {
        this.replicaGroup.clear();
        this.replicaBox = null; // updateReplicasBox will handle it
        this.updateReplicasBox();
        this.replicaClones = [];

        // Require valid lattice vectors; otherwise skip
        const lvecOK = this.replicaConfig.lvec && this.replicaConfig.lvec.length === 3;
        if (!this.replicaConfig.show || !lvecOK) {
            this.replicaConfig.dirty = false;
            return;
        }

        const { x: nx, y: ny, z: nz } = this.replicaConfig.nrep;
        const [a, b, c] = this.replicaConfig.lvec;

        const baseMeshes = [
            this.atomMesh,
            this.bondLines,
            this.labelMesh
        ].filter(m => m !== null);

        for (let ix = -nx; ix <= nx; ix++) {
            for (let iy = -ny; iy <= ny; iy++) {
                for (let iz = -nz; iz <= nz; iz++) {
                    if (ix === 0 && iy === 0 && iz === 0) continue;
                    const shift = new Vec3();
                    shift.addMul(a, ix);
                    shift.addMul(b, iy);
                    shift.addMul(c, iz);
                    const rep = new THREE.Group();
                    rep.position.set(shift.x, shift.y, shift.z);
                    baseMeshes.forEach(mesh => {
                        let clone;
                        if (mesh instanceof THREE.LineSegments) {
                            clone = new THREE.LineSegments(mesh.geometry, mesh.material);
                        } else if (mesh instanceof THREE.InstancedMesh) {
                            clone = new THREE.InstancedMesh(mesh.geometry, mesh.material, mesh.count);
                            clone.count = mesh.count;
                            this.replicaClones.push({ clone, source: mesh });
                        } else {
                            clone = new THREE.Mesh(mesh.geometry, mesh.material);
                        }
                        clone.raycast = () => { }; // non-pickable
                        rep.add(clone);
                    });
                    this.replicaGroup.add(rep);
                }
            }
        }
        this.replicaConfig.dirty = false;
        this._syncReplicaInstancedMeshes();
    }

    updateReplicasBox() {
        if (this.replicaBox) {
            this.replicaGroup.remove(this.replicaBox);
            if (this.replicaBox.geometry) this.replicaBox.geometry.dispose();
            if (this.replicaBox.material) this.replicaBox.material.dispose();
            this.replicaBox = null;
        }
        if (!this.replicaConfig.showBox) return;

        const [a, b, c] = this.replicaConfig.lvec;
        const origin = new Vec3(); // draw only the base unit cell, independent of nrep

        const verts = this._buildWireframeCellVerts(a, b, c, origin);
        const geom = new THREE.BufferGeometry();
        geom.setAttribute('position', new THREE.BufferAttribute(verts, 3));
        const mat = new THREE.LineBasicMaterial({ color: 0x00ffff, opacity: 0.5, transparent: true });
        this.replicaBox = new THREE.LineSegments(geom, mat);
        this.replicaGroup.add(this.replicaBox);
    }

    _buildWireframeCellVerts(a, b, c, origin) {
        const v = new Float32Array(24 * 3);
        const corners = [
            origin.clone(),
            origin.clone().add(a),
            origin.clone().add(b),
            origin.clone().add(a).add(b),
            origin.clone().add(c),
            origin.clone().add(a).add(c),
            origin.clone().add(b).add(c),
            origin.clone().add(a).add(b).add(c)
        ];
        const edges = [
            [0, 1], [0, 2], [0, 4], [1, 3], [1, 5], [2, 3], [2, 6], [3, 7], [4, 5], [4, 6], [5, 7], [6, 7]
        ];
        for (let i = 0; i < edges.length; i++) {
            const [i0, i1] = edges[i];
            const p0 = corners[i0], p1 = corners[i1];
            v[i * 6] = p0.x; v[i * 6 + 1] = p0.y; v[i * 6 + 2] = p0.z;
            v[i * 6 + 3] = p1.x; v[i * 6 + 4] = p1.y; v[i * 6 + 5] = p1.z;
        }
        return v;
    }

    update() {
        // Main update loop
        if (this.source && (this.source.dirtyExport || this.source.dirtyTopo || this.source.dirtyGeom)) {
            this.source.exportToMoleculeSystem(this.system);
        }
        if (this.system && (this.system.capacity | 0) > (this.capacity | 0)) {
            this.ensureCapacity(this.system.capacity);
        }
        if (this.system.isDirty) {
            this.updateStructure();
            this.updatePositions(); // Structure change implies position update too
            this.system.isDirty = false;
            this._syncReplicaInstancedMeshes();
        }

        if (this.replicaConfig.dirty) {
            this.rebuildReplicas();
        }

        // Update Labels (visibility check handled in parent)
        // But we need to update content if mode changes or structure changes.
        // Here we just ensure visibility logic if needed.
    }

    updatePositions() {
        super.updatePositions(this.system.pos, this.system.nAtoms);
    }

    updateStructure() {
        const nAtoms = this.system.nAtoms;

        // Update Particles
        this.updateParticles(nAtoms, (i) => {
            const type = this.system.types[i];
            // Use MMParams
            let col = [1, 0, 1]; // Default magenta
            if (this.mmParams) {
                col = this.mmParams.getColor(type);
            }
            return col;
        }, (i) => {
            const type = this.system.types[i];
            let radius = 1.0;
            if (this.mmParams) {
                radius = this.mmParams.getRadius(type) * 0.4;
            }
            return radius;
        });

        // Update Bonds
        this.updateBonds(this.system.bonds);

        this.updateSelection();
        this.updateLabelsContent(); // Update labels when structure changes

        if (window.VERBOSITY_LEVEL >= Logger.INFO) {
            // window.logger.info(`Renderer structure updated: ${ nAtoms } atoms.`);
        }
    }

    updateSelection() {
        if (this.source && this.source.dirtyExport) {
            this.source.exportToMoleculeSystem(this.system);
        }
        const selectedIDs = Array.from(this.system.selection);
        super.updateSelection(selectedIDs);
    }

    updateLabelsContent() {
        if (!this.labelMesh || this.labelMode === 'none') return;

        const nAtoms = this.system.nAtoms;

        this.updateLabels((i) => {
            if (this.labelMode === 'id') {
                return i.toString();
            } else if (this.labelMode === 'element') {
                const type = this.system.types[i];
                if (this.mmParams && this.mmParams.byAtomicNumber[type]) {
                    return this.mmParams.byAtomicNumber[type].name;
                }
                return type.toString();
            } else if (this.labelMode === 'type') {
                return "?"; // Placeholder as per previous logic
            }
            return "";
        }, nAtoms);
    }

    toggleAxes(visible) {
        if (!this.axesHelper) {
            this.axesHelper = new THREE.AxesHelper(5); // 5 units length
            this.scene.add(this.axesHelper);
        }
        this.axesHelper.visible = visible;
    }

    setLabelMode(mode) {
        if (this.labelMode !== mode) {
            this.labelMode = mode;
            if (window && window.app) {
                console.debug('[MoleculeRenderer] label mode changed', mode);
                if (typeof window.app.requestRender === 'function') window.app.requestRender();
            }
            if (this.labelMesh) {
                this.labelMesh.visible = (mode !== 'none');
                if (this.labelMesh.visible) {
                    this.updateLabelsContent();
                }
            }
        }
    }

    setLabelStyle(colorHex, scale) {
        if (!this.labelMesh) return;
        if (!this._labelUniforms) throw new Error('setLabelStyle: missing this._labelUniforms (MeshRenderer_webgpu init did not expose label uniforms)');

        if (colorHex !== undefined && colorHex !== null) {
            const c = new THREE.Color(colorHex);
            this._labelUniforms.uLabelColor.value.set(c.r, c.g, c.b);
        }
        if (scale !== undefined && scale !== null) {
            this._labelUniforms.uLabelScale.value = parseFloat(scale);
        }
    }

    setPDSimulation(pdSim) {
        if (this.pdSimulation && this._pdAngleCallback && this.pdSimulation.onAngleDebugUpdate === this._pdAngleCallback) {
            this.pdSimulation.onAngleDebugUpdate = null;
        }
        this.pdSimulation = pdSim || null;
        if (this.pdSimulation) {
            this._pdAngleCallback = () => this.updateAngleConstraintLines();
            this.pdSimulation.onAngleDebugUpdate = this._pdAngleCallback;
            this.updateAngleConstraintLines();
        } else {
            this._pdAngleCallback = null;
            if (this.angleConstraintLines) this.angleConstraintLines.visible = false;
        }
    }

    toggleAngleConstraints(visible) {
        this.showAngleConstraints = !!visible;
        if (!this.showAngleConstraints) {
            if (this.angleConstraintLines) {
                this.angleConstraintLines.visible = false;
            }
        } else {
            this.updateAngleConstraintLines();
        }
    }

    updateAngleConstraintLines() {
        if (!this.pdSimulation || !this.showAngleConstraints) {
            if (this.angleConstraintLines) {
                this.angleConstraintLines.visible = false;
            }
            return;
        }

        const bonds = this.pdSimulation.debugAngleBonds || [];
        if (bonds.length === 0) {
            if (this.angleConstraintLines) {
                this.angleConstraintLines.visible = false;
            }
            return;
        }

        const nBonds = bonds.length;
        const verts = new Float32Array(nBonds * 6);
        let floatCount = 0;

        for (let i = 0; i < nBonds; i++) {
            const bond = bonds[i];
            const a = bond.a | 0;
            const c = bond.c | 0;

            if (a < 0 || a >= this.system.nAtoms || c < 0 || c >= this.system.nAtoms) {
                continue;
            }

            const posA = this.system.pos;
            const idxA = a * 3;
            const idxC = c * 3;

            verts[floatCount++] = posA[idxA];
            verts[floatCount++] = posA[idxA + 1];
            verts[floatCount++] = posA[idxA + 2];
            verts[floatCount++] = posA[idxC];
            verts[floatCount++] = posA[idxC + 1];
            verts[floatCount++] = posA[idxC + 2];
        }

        if (floatCount === 0) {
            if (this.angleConstraintLines) this.angleConstraintLines.visible = false;
            return;
        }

        const usedVerts = (floatCount === verts.length) ? verts : verts.slice(0, floatCount);

        if (!this.angleConstraintLines) {
            const geometry = new THREE.BufferGeometry();
            geometry.setAttribute('position', new THREE.BufferAttribute(usedVerts, 3));
            geometry.computeBoundingSphere();

            const material = new THREE.LineBasicMaterial({
                color: 0x00ff00,
                linewidth: 1,
                transparent: true,
                opacity: 0.8
            });

            this.angleConstraintLines = new THREE.LineSegments(geometry, material);
            this.scene.add(this.angleConstraintLines);
        } else {
            const geometry = this.angleConstraintLines.geometry;
            let posAttr = geometry.getAttribute('position');
            if (!posAttr || posAttr.array.length !== usedVerts.length) {
                posAttr = new THREE.BufferAttribute(usedVerts, 3);
                geometry.setAttribute('position', posAttr);
            } else {
                posAttr.array.set(usedVerts);
                posAttr.needsUpdate = true;
            }
            geometry.computeBoundingSphere();
        }

        this.angleConstraintLines.geometry.setDrawRange(0, usedVerts.length / 3);
        this.angleConstraintLines.visible = true;
    }

    _syncReplicaInstancedMeshes() {
        if (!this.replicaClones || this.replicaClones.length === 0) return;
        for (let i = 0; i < this.replicaClones.length; i++) {
            const entry = this.replicaClones[i];
            if (!entry || !entry.clone || !entry.source) continue;
            entry.clone.count = entry.source.count;
            if (entry.clone.instanceMatrix && entry.source.instanceMatrix) {
                entry.clone.instanceMatrix.copy(entry.source.instanceMatrix);
                entry.clone.instanceMatrix.needsUpdate = true;
            }
        }
    }
}

