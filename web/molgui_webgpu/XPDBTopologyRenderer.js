import * as THREE from 'three';

function _hexToThreeColor(hex) {
    const c = new THREE.Color();
    c.setHex(hex >>> 0);
    return c;
}

function _pairsToLinePositions(pairs, pos, nAtoms, label) {
    if (!pairs || pairs.length === 0) return new Float32Array(0);
    const verts = new Float32Array(pairs.length * 2 * 3);
    let k = 0;
    for (let i = 0; i < pairs.length; i++) {
        const p = pairs[i];
        if (!p || p.length < 2) throw new Error(`${label}: invalid pair at i=${i}`);
        const a = p[0] | 0;
        const b = p[1] | 0;
        if (a < 0 || a >= nAtoms || b < 0 || b >= nAtoms) {
            throw new Error(`${label}: pair out of range i=${i} a=${a} b=${b} nAtoms=${nAtoms}`);
        }
        const ia = a * 3;
        const ib = b * 3;
        verts[k++] = pos[ia + 0];
        verts[k++] = pos[ia + 1];
        verts[k++] = pos[ia + 2];
        verts[k++] = pos[ib + 0];
        verts[k++] = pos[ib + 1];
        verts[k++] = pos[ib + 2];
    }
    return verts;
}

function _buildAtomPositionsFromTopo(topo) {
    if (!topo) throw new Error('_buildAtomPositionsFromTopo: topo is null');
    const nAll = topo.n_all | 0;
    if (!(nAll > 0)) return new Float32Array(0);
    if (!topo.apos || topo.apos.length !== nAll) {
        throw new Error(`_buildAtomPositionsFromTopo: topo.apos missing or length mismatch aposLen=${topo.apos ? topo.apos.length : 'null'} n_all=${nAll}`);
    }
    const pos = new Float32Array(nAll * 3);
    for (let i = 0; i < nAll; i++) {
        const p = topo.apos[i];
        if (!p || p.length < 3) throw new Error(`_buildAtomPositionsFromTopo: invalid apos[${i}]`);
        pos[i * 3 + 0] = +p[0];
        pos[i * 3 + 1] = +p[1];
        pos[i * 3 + 2] = +p[2];
    }
    return pos;
}

function _splitAtomPositions(topo, pos3) {
    const nAll = topo.n_all | 0;
    const nReal = topo.n_real | 0;
    if (pos3.length !== nAll * 3) throw new Error(`_splitAtomPositions: pos3 length mismatch pos3=${pos3.length} expected=${nAll * 3}`);
    const nDummy = Math.max(0, nAll - nReal);
    const posReal = new Float32Array(nReal * 3);
    const posDummy = new Float32Array(nDummy * 3);
    if (nReal > 0) posReal.set(pos3.subarray(0, nReal * 3));
    if (nDummy > 0) posDummy.set(pos3.subarray(nReal * 3, nAll * 3));
    return { posReal, posDummy, nReal, nDummy, nAll };
}

export class XPDBTopologyRenderer {
    constructor(scene) {
        if (!scene) throw new Error('XPDBTopologyRenderer: scene is required');
        this.scene = scene;
        this.group = new THREE.Group();
        this.group.name = 'XPDBTopologyOverlay';
        this.scene.add(this.group);

        this.enabled = false;

        this.atom = {
            real: null,
            dummy: null,
        };

        this.lines = new Map();

        this.style = {
            atomReal: { color: 0xffffff, size: 6.0, opacity: 0.9 },
            atomDummy: { color: 0xffcc55, size: 5.0, opacity: 0.8 },
            bonds: {
                primary: { color: 0xffffff, opacity: 0.65 },
                angle: { color: 0x55ccff, opacity: 0.5 },
                pi: { color: 0xff55ff, opacity: 0.55 },
                pi_align: { color: 0xaa55ff, opacity: 0.55 },
                epair: { color: 0xffaa00, opacity: 0.55 },
                epair_pair: { color: 0xff4400, opacity: 0.55 },
            },
        };

        this.visibility = {
            atoms_real: true,
            atoms_dummy: true,
            primary: true,
            angle: true,
            pi: true,
            pi_align: true,
            epair: true,
            epair_pair: true,
        };

        this.clear();
        this.setEnabled(false);
    }

    setEnabled(enabled) {
        this.enabled = !!enabled;
        this.group.visible = this.enabled;
    }

    setCategoryVisible(category, visible) {
        if (!(category in this.visibility)) throw new Error(`XPDBTopologyRenderer.setCategoryVisible: unknown category='${category}'`);
        this.visibility[category] = !!visible;
        if (category === 'atoms_real' && this.atom.real) this.atom.real.visible = this.visibility.atoms_real;
        if (category === 'atoms_dummy' && this.atom.dummy) this.atom.dummy.visible = this.visibility.atoms_dummy;
        const obj = this.lines.get(category);
        if (obj) obj.visible = !!visible;
    }

    setVisibilityMap(vis) {
        if (!vis) return;
        for (const k of Object.keys(vis)) this.setCategoryVisible(k, !!vis[k]);
    }

    clear() {
        if (this.atom.real) { this.group.remove(this.atom.real); this._disposeObject(this.atom.real); }
        if (this.atom.dummy) { this.group.remove(this.atom.dummy); this._disposeObject(this.atom.dummy); }
        this.atom.real = null;
        this.atom.dummy = null;
        for (const obj of this.lines.values()) {
            this.group.remove(obj);
            this._disposeObject(obj);
        }
        this.lines.clear();
    }

    updateFromTopo(topo, opts = {}) {
        if (!topo) throw new Error('XPDBTopologyRenderer.updateFromTopo: topo is null');
        const enabled = (opts.enabled !== undefined) ? !!opts.enabled : this.enabled;
        this.setEnabled(enabled);
        if (!this.enabled) return;

        const verbosity = (typeof window !== 'undefined') ? (window.VERBOSITY_LEVEL | 0) : 0;
        if (verbosity >= 2) {
            console.log('[XPDBTopologyRenderer.updateFromTopo]', {
                n_all: topo.n_all | 0,
                n_real: topo.n_real | 0,
                bonds_primary: topo.bonds_primary ? topo.bonds_primary.length : 0,
                bonds_angle: topo.bonds_angle ? topo.bonds_angle.length : 0,
                bonds_pi: topo.bonds_pi ? topo.bonds_pi.length : 0,
                bonds_pi_align: topo.bonds_pi_align ? topo.bonds_pi_align.length : 0,
                bonds_epair: topo.bonds_epair ? topo.bonds_epair.length : 0,
                bonds_epair_pair: topo.bonds_epair_pair ? topo.bonds_epair_pair.length : 0,
            });
        }

        const pos3 = _buildAtomPositionsFromTopo(topo);
        const { posReal, posDummy, nReal, nDummy, nAll } = _splitAtomPositions(topo, pos3);

        // Rebuild everything each update (topology can change, categories can be toggled)
        this.clear();

        if (nReal > 0) {
            this.atom.real = this._makePoints(posReal, this.style.atomReal);
            this.atom.real.renderOrder = 20;
            this.atom.real.visible = !!this.visibility.atoms_real;
            this.group.add(this.atom.real);
        }
        if (nDummy > 0) {
            this.atom.dummy = this._makePoints(posDummy, this.style.atomDummy);
            this.atom.dummy.renderOrder = 21;
            this.atom.dummy.visible = !!this.visibility.atoms_dummy;
            this.group.add(this.atom.dummy);
        }

        const addLines = (key, pairs, styleKey) => {
            const verts = _pairsToLinePositions(pairs, pos3, nAll, `XPDBTopologyRenderer:${key}`);
            const obj = this._makeLines(verts, this.style.bonds[styleKey]);
            obj.renderOrder = 10;
            obj.visible = !!this.visibility[key];
            this.lines.set(key, obj);
            this.group.add(obj);
        };

        addLines('primary', topo.bonds_primary, 'primary');
        addLines('angle', topo.bonds_angle, 'angle');
        addLines('pi', topo.bonds_pi, 'pi');
        addLines('pi_align', topo.bonds_pi_align, 'pi_align');
        addLines('epair', topo.bonds_epair, 'epair');
        addLines('epair_pair', topo.bonds_epair_pair, 'epair_pair');

        // apply any per-update overrides last
        if (opts.visibility) this.setVisibilityMap(opts.visibility);
    }

    _makePoints(pos3, st) {
        const geom = new THREE.BufferGeometry();
        geom.setAttribute('position', new THREE.BufferAttribute(pos3, 3));
        geom.computeBoundingSphere();
        const mat = new THREE.PointsMaterial({
            color: _hexToThreeColor(st.color),
            size: st.size,
            sizeAttenuation: true,
            transparent: true,
            opacity: st.opacity,
            depthTest: false,
            depthWrite: false,
        });
        return new THREE.Points(geom, mat);
    }

    _makeLines(pos3, st) {
        const geom = new THREE.BufferGeometry();
        geom.setAttribute('position', new THREE.BufferAttribute(pos3, 3));
        geom.computeBoundingSphere();
        const mat = new THREE.LineBasicMaterial({
            color: _hexToThreeColor(st.color),
            transparent: true,
            opacity: st.opacity,
            depthWrite: false,
        });
        return new THREE.LineSegments(geom, mat);
    }

    _disposeObject(obj) {
        if (!obj) return;
        if (obj.geometry && obj.geometry.dispose) obj.geometry.dispose();
        if (obj.material) {
            if (Array.isArray(obj.material)) obj.material.forEach(m => m && m.dispose && m.dispose());
            else if (obj.material.dispose) obj.material.dispose();
        }
    }
}
