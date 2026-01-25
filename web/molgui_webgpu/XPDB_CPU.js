export class XPDB_CPU {
    constructor(nAtoms) {
        this.nAtoms = nAtoms;
        this.pos = new Float32Array(nAtoms * 3);
        this.posPrev = new Float32Array(nAtoms * 3);
        this.posNew = new Float32Array(nAtoms * 3);
        this.pred = new Float32Array(nAtoms * 3);
        this.rad = new Float32Array(nAtoms);
        this.mass = new Float32Array(nAtoms);
        this._bondedSets = null;
        this.lastDelta = 0.0;
    }

    uploadFromMolecule(mol, mmParams) {
        const nAtoms = mol.nAtoms || mol.atoms?.length || 0;
        if (nAtoms !== this.nAtoms) throw new Error(`XPDB_CPU.uploadFromMolecule: atom count mismatch (mol=${nAtoms}, cpu=${this.nAtoms})`);
        for (let i = 0; i < this.nAtoms; i++) {
            const a = mol.atoms[i];
            if (!a) throw new Error(`XPDB_CPU.uploadFromMolecule: missing atom ${i}`);
            this.pos[i * 3 + 0] = a.pos.x;
            this.pos[i * 3 + 1] = a.pos.y;
            this.pos[i * 3 + 2] = a.pos.z;
            const at = mmParams.getAtomTypeForAtom(a);
            this.rad[i] = (at && at.RvdW > 0) ? at.RvdW : 1.0;
            this.mass[i] = (a.Z || 1) * 1.0;
        }
    }

    _ensureBondedSets(bondsAdj) {
        if (!bondsAdj) throw new Error('XPDB_CPU._ensureBondedSets: bondsAdj is null');
        const bonded = new Array(this.nAtoms);
        for (let i = 0; i < this.nAtoms; i++) bonded[i] = new Set();
        for (let i = 0; i < this.nAtoms; i++) {
            const bs = bondsAdj[i] || [];
            for (let k = 0; k < bs.length; k++) {
                const j = bs[k][0] | 0;
                if (j >= 0 && j < this.nAtoms) bonded[i].add(j);
            }
        }
        this._bondedSets = bonded;
    }

    step(bondsAdj, { dt, iterations, k_coll, omega, momentum_beta }) {
        if (!bondsAdj) throw new Error('XPDB_CPU.step: bondsAdj is null');
        if (!Number.isFinite(dt) || dt <= 0) throw new Error(`XPDB_CPU.step: bad dt=${dt}`);
        if (!Number.isFinite(iterations) || iterations < 1) throw new Error(`XPDB_CPU.step: bad iterations=${iterations}`);

        if (!this._bondedSets) this._ensureBondedSets(bondsAdj);

        this.pred.set(this.pos);
        this.posPrev.set(this.pos);

        const n = this.nAtoms;
        const pos = this.pos;
        const posPrev = this.posPrev;
        const posNew = this.posNew;
        const pred = this.pred;
        const rad = this.rad;
        const mass = this.mass;
        const bonded = this._bondedSets;

        const dbgV = (window.VERBOSITY_LEVEL | 0);
        const dbg = dbgV >= 4;
        const dbgAtoms = 4;
        const dbgIters = 1;
        let maxDeltaOverall = 0.0;

        for (let it = 0; it < iterations; it++) {
            let maxDeltaIter = 0.0;
            for (let i = 0; i < n; i++) {
                const ix = pos[i * 3 + 0];
                const iy = pos[i * 3 + 1];
                const iz = pos[i * 3 + 2];

                const alpha = mass[i] / (dt * dt);

                let fx = alpha * (pred[i * 3 + 0] - ix);
                let fy = alpha * (pred[i * 3 + 1] - iy);
                let fz = alpha * (pred[i * 3 + 2] - iz);
                let ksum = alpha;

                if (dbg && it < dbgIters && i < dbgAtoms) {
                    console.log(`[XPDB_CPU][DEBUG] it=${it} atom=${i} p=(${ix.toFixed(4)},${iy.toFixed(4)},${iz.toFixed(4)}) alpha=${alpha.toFixed(2)}`);
                }

                const bs = bondsAdj[i] || [];
                for (let s = 0; s < bs.length; s++) {
                    const j = bs[s][0] | 0;
                    const L = +bs[s][1];
                    const st = +bs[s][2];
                    if (!(j >= 0 && j < n)) continue;

                    const jx = pos[j * 3 + 0];
                    const jy = pos[j * 3 + 1];
                    const jz = pos[j * 3 + 2];

                    const dx = ix - jx;
                    const dy = iy - jy;
                    const dz = iz - jz;
                    const dist2 = dx * dx + dy * dy + dz * dz;
                    if (dist2 > 1e-12) {
                        const dist = Math.sqrt(dist2);
                        const inv = 1.0 / dist;
                        const C = dist - L;
                        const f = -st * C;
                        fx += f * dx * inv;
                        fy += f * dy * inv;
                        fz += f * dz * inv;
                        ksum += st;

                        if (dbg && it < dbgIters && i < dbgAtoms) {
                            console.log(`[XPDB_CPU][DEBUG]   bond s=${s} j=${j} L=${L.toFixed(4)} K=${st.toFixed(2)} dist=${dist.toFixed(4)} C=${C.toFixed(4)} f=${f.toFixed(4)}`);
                        }
                    }
                }

                if (k_coll !== 0.0) {
                    const ri = rad[i];
                    const bonded_i = bonded[i];
                    for (let j = 0; j < n; j++) {
                        if (j === i) continue;
                        if (bonded_i.has(j)) continue;

                        const jx = pos[j * 3 + 0];
                        const jy = pos[j * 3 + 1];
                        const jz = pos[j * 3 + 2];

                        const dx = ix - jx;
                        const dy = iy - jy;
                        const dz = iz - jz;
                        const dist2 = dx * dx + dy * dy + dz * dz;
                        const rsum = ri + rad[j];
                        if (dist2 < rsum * rsum && dist2 > 1e-12) {
                            const dist = Math.sqrt(dist2);
                            const inv = 1.0 / dist;
                            const overlap = dist - rsum;
                            const f = -k_coll * overlap;
                            fx += f * dx * inv;
                            fy += f * dy * inv;
                            fz += f * dz * inv;
                            ksum += k_coll;

                            if (dbg && it < dbgIters && i < dbgAtoms) {
                                console.log(`[XPDB_CPU][DEBUG]   coll j=${j} ri=${ri.toFixed(3)} rj=${rad[j].toFixed(3)} dist=${dist.toFixed(4)} overlap=${overlap.toFixed(4)} f=${f.toFixed(4)}`);
                            }
                        }
                    }
                }

                let nx = ix + (fx / ksum) * omega;
                let ny = iy + (fy / ksum) * omega;
                let nz = iz + (fz / ksum) * omega;
                if (momentum_beta !== 0.0) {
                    nx += (ix - posPrev[i * 3 + 0]) * momentum_beta;
                    ny += (iy - posPrev[i * 3 + 1]) * momentum_beta;
                    nz += (iz - posPrev[i * 3 + 2]) * momentum_beta;
                }

                posNew[i * 3 + 0] = nx;
                posNew[i * 3 + 1] = ny;
                posNew[i * 3 + 2] = nz;

                if (dbg && it < dbgIters && i < dbgAtoms) {
                    console.log(`[XPDB_CPU][DEBUG]   sum f=(${fx.toFixed(4)},${fy.toFixed(4)},${fz.toFixed(4)}) ksum=${ksum.toFixed(4)} pNew=(${nx.toFixed(4)},${ny.toFixed(4)},${nz.toFixed(4)})`);
                }

                const dnx = nx - ix;
                const dny = ny - iy;
                const dnz = nz - iz;
                const deltaLen = Math.sqrt(dnx * dnx + dny * dny + dnz * dnz);
                if (deltaLen > maxDeltaIter) maxDeltaIter = deltaLen;
            }

            posPrev.set(pos);
            pos.set(posNew);
            if (maxDeltaIter > maxDeltaOverall) maxDeltaOverall = maxDeltaIter;
        }

        this.lastDelta = maxDeltaOverall;
    }

    writeToMolecule(mol) {
        const nAtoms = mol.nAtoms || mol.atoms?.length || 0;
        if (nAtoms !== this.nAtoms) throw new Error(`XPDB_CPU.writeToMolecule: atom count mismatch (mol=${nAtoms}, cpu=${this.nAtoms})`);
        for (let i = 0; i < this.nAtoms; i++) {
            const a = mol.atoms[i];
            if (!a) continue;
            a.pos.x = this.pos[i * 3 + 0];
            a.pos.y = this.pos[i * 3 + 1];
            a.pos.z = this.pos[i * 3 + 2];
        }
    }
}
