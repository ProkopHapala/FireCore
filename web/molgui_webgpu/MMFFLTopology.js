import { Vec3 } from "../common_js/Vec3.js";
import { Z_TO_SYMBOL } from "./MoleculeIO.js";

function lawOfCosines(rab, rbc, cosTheta) {
    const t = rab * rab + rbc * rbc - 2.0 * rab * rbc * cosTheta;
    return Math.sqrt(Math.max(0.0, t));
}

export function packBondArrays(bondsAdj, nAtoms, nMaxBonded) {
    const bondIndices = new Int32Array(nAtoms * nMaxBonded);
    bondIndices.fill(-1);
    const bondLenStiff = new Float32Array(nAtoms * nMaxBonded * 2);
    bondLenStiff.fill(0);
    for (let i = 0; i < nAtoms; i++) {
        const neighs = bondsAdj[i] || [];
        if (neighs.length > nMaxBonded) {
            console.log(`[WARN] Atom ${i} has ${neighs.length} neighbors > n_max_bonded=${nMaxBonded}; truncating`);
        }
        const base = i * nMaxBonded;
        const m = Math.min(neighs.length, nMaxBonded);
        for (let k = 0; k < m; k++) {
            const b = neighs[k];
            bondIndices[base + k] = (b[0] | 0);
            bondLenStiff[(base + k) * 2 + 0] = +b[1];
            bondLenStiff[(base + k) * 2 + 1] = +b[2];
        }
    }
    return { bondIndices, bondLenStiff };
}

export function buildXPDBInputsFromMol(mol, mm, opts = {}) {
    if (!mol || !mm) throw new Error('buildXPDBInputsFromMol: mol and mm required');
    const nMaxBonded = (opts.nMaxBonded !== undefined) ? (opts.nMaxBonded | 0) : 16;
    if (!(nMaxBonded > 0)) throw new Error(`buildXPDBInputsFromMol: invalid nMaxBonded=${nMaxBonded}`);
    const topo = buildMMFFLTopology(mol, mm, opts);

    const nAtoms = topo.n_all | 0;
    const bondsAdj = new Array(nAtoms);
    for (let i = 0; i < nAtoms; i++) bondsAdj[i] = [];

    const addEdge = (i, j, L, K) => { bondsAdj[i].push([j | 0, +L, +K]); };
    const uniq = new Map();
    const bondKey = (a, b) => {
        const i = (a < b) ? a : b;
        const j = (a < b) ? b : a;
        return `${i},${j}`;
    };
    const derived = [];
    if (opts.add_angle) derived.push(...topo.bonds_angle);
    if (opts.add_pi) derived.push(...topo.bonds_pi);
    if (opts.add_pi && opts.add_pi_align) derived.push(...topo.bonds_pi_align);
    if (opts.add_epair) derived.push(...topo.bonds_epair);
    if (opts.add_epair && opts.add_epair_pairs) derived.push(...topo.bonds_epair_pair);
    const all0 = [...topo.bonds_primary, ...derived];
    for (const e of all0) {
        const a = e[0] | 0;
        const b = e[1] | 0;
        const i = (a < b) ? a : b;
        const j = (a < b) ? b : a;
        uniq.set(bondKey(a, b), [i, j]);
    }
    const bondsAll = Array.from(uniq.values());
    bondsAll.sort((p, q) => (p[0] - q[0]) || (p[1] - q[1]));

    const dist = (a, b) => {
        const pa = topo.apos[a];
        const pb = topo.apos[b];
        const dx = pa[0] - pb[0];
        const dy = pa[1] - pb[1];
        const dz = pa[2] - pb[2];
        return Math.sqrt(dx * dx + dy * dy + dz * dz);
    };
    const K0 = (opts.defaultK !== undefined) ? +opts.defaultK : 200.0;
    for (const [a, b] of bondsAll) {
        const L = dist(a, b);
        addEdge(a, b, L, K0);
        addEdge(b, a, L, K0);
    }

    const { bondIndices, bondLenStiff } = packBondArrays(bondsAdj, nAtoms, nMaxBonded);

    const pos4 = new Float32Array(nAtoms * 4);
    const params4 = new Float32Array(nAtoms * 4);
    const atom_rad = (opts.atom_rad !== undefined) ? +opts.atom_rad : 0.2;
    const atom_mass = (opts.atom_mass !== undefined) ? +opts.atom_mass : 1.0;
    for (let i = 0; i < nAtoms; i++) {
        const p = topo.apos[i];
        pos4[i * 4 + 0] = +p[0];
        pos4[i * 4 + 1] = +p[1];
        pos4[i * 4 + 2] = +p[2];
        pos4[i * 4 + 3] = 0.0;
        params4[i * 4 + 0] = atom_rad;
        params4[i * 4 + 3] = atom_mass;
    }

    return { topo, bondsAdj, nAtoms, nMaxBonded, bondIndices, bondLenStiff, pos4, params4 };
}

function piDomainDirForVSEPR(mol, neighs, hostPos) {
    const tmp = new Vec3();
    const v1 = new Vec3();
    const v2 = new Vec3();
    if (!neighs || neighs.length === 0) return null;
    if (neighs.length >= 2) {
        const j1 = neighs[0][0] | 0;
        const j2 = neighs[1][0] | 0;
        tmp.setV(mol.atoms[j1].pos); tmp.sub(hostPos); if (!(tmp.normalize() > 0)) return null; v1.setV(tmp);
        tmp.setV(mol.atoms[j2].pos); tmp.sub(hostPos); if (!(tmp.normalize() > 0)) return null; v2.setV(tmp);
        const cr = new Vec3().setCross(v1, v2);
        if (cr.norm() > 1e-6) { cr.mulScalar(1.0 / cr.norm()); return cr; }
        // fallthrough to nb=1-like if bonds are colinear
    }
    // neighs.length === 1 or colinear case
    const j = neighs[0][0] | 0;
    tmp.setV(mol.atoms[j].pos);
    tmp.sub(hostPos);
    if (!(tmp.normalize() > 0)) return null;
    const axis = tmp.clone();
    const u = new Vec3();
    const v = new Vec3();
    orthonormalBasisFromDir(axis, u, v);
    return u;
}

function addUndirectedPair(list, a, b) {
    const i = a < b ? a : b;
    const j = a < b ? b : a;
    list.push([i, j]);
}

function uniqPairs(pairs) {
    pairs.sort((p, q) => (p[0] - q[0]) || (p[1] - q[1]));
    const out = [];
    let pa = -1, pb = -1;
    for (let k = 0; k < pairs.length; k++) {
        const a = pairs[k][0] | 0;
        const b = pairs[k][1] | 0;
        if (a !== pa || b !== pb) out.push([a, b]);
        pa = a; pb = b;
    }
    return out;
}

function extendMolWithDummy(mol, Z, pos, atypeIndex = -1) {
    if (!mol || typeof mol.addAtom !== 'function' || typeof mol.getAtomIndex !== 'function') {
        throw new Error('extendMolWithDummy: mol must provide addAtom() and getAtomIndex()');
    }
    const id = mol.addAtom(pos.x, pos.y, pos.z, Z);
    const i = mol.getAtomIndex(id);
    if (i < 0) throw new Error('extendMolWithDummy: getAtomIndex() returned <0');
    mol.atoms[i].atype = atypeIndex;
    return { id, i };
}

function buildAngleBonds(mol, mmParams, bondsAdj1, opts, outLinear) {
    const n = mol.atoms.length;
    const kAng = (opts && (opts.k_angle !== undefined)) ? +opts.k_angle : 0.0;
    for (let ib = 0; ib < n; ib++) {
        const neighs = bondsAdj1[ib] || [];
        if (neighs.length < 2) continue;
        const symB = Z_TO_SYMBOL[mol.atoms[ib].Z] || 'X';
        const tBname = mmParams.resolveTypeNameTable(symB);
        const tB = mmParams.atomTypes[tBname];
        const thetaDeg = (tB && (tB.Ass !== undefined)) ? +tB.Ass : 120.0;
        const theta = thetaDeg * Math.PI / 180.0;
        const cosT = Math.cos(theta);
        const K = (kAng !== 0.0) ? kAng : ((tB && (tB.Kss !== undefined)) ? +tB.Kss : 0.0);
        if (!(K !== 0.0)) continue;

        for (let ia = 0; ia < neighs.length; ia++) {
            const ja = neighs[ia][0] | 0;
            const rab = +neighs[ia][1];
            if (!(rab > 0)) continue;
            for (let ic = ia + 1; ic < neighs.length; ic++) {
                const jc = neighs[ic][0] | 0;
                if (ja === jc) continue;
                const rbc = +neighs[ic][1];
                if (!(rbc > 0)) continue;
                const l0 = lawOfCosines(rab, rbc, cosT);
                if (!(l0 > 0)) continue;
                outLinear.push([ja, jc, l0, K, ['angle', [ja, ib, jc]]]);
            }
        }
    }
}

function countNeighborsFromAdj(adj) {
    let c = 0;
    for (let i = 0; i < adj.length; i++) if ((adj[i][0] | 0) >= 0) c++;
    return c;
}

function computeAtomiPiDirectionFromNeighs(mol, ia, nbSorted, tmpA, tmpB, tmpC) {
    // Python AtomicSystem.get_atomi_pi_direction():
    // vectors = normalize(rj-ri) for neighbors
    // dir = sum normalize(cross(v_k, v_{k+1})) with cyclic wrap
    const hostPos = mol.atoms[ia].pos;
    const vecs = [];
    for (let k = 0; k < nbSorted.length; k++) {
        const jb = nbSorted[k] | 0;
        if (jb < 0 || jb >= mol.atoms.length) continue;
        tmpA.setV(mol.atoms[jb].pos);
        tmpA.sub(hostPos);
        if (tmpA.normalize() > 1e-12) vecs.push(tmpA.clone());
    }
    tmpC.set(0, 0, 0);
    const m = vecs.length;
    if (m < 2) return 0.0;
    for (let k = 0; k < m; k++) {
        const a = vecs[k];
        const b = vecs[(k + 1) % m];
        tmpB.setCross(a, b);
        const n = tmpB.norm();
        if (n > 1e-12) tmpC.add(tmpB.mulScalar(1.0 / n));
    }
    return tmpC.norm();
}

function computePiOrientations(mol, bondsAdj1, npiList, alignPiVectors = false, debugTopo = null) {
    const n = mol.atoms.length;
    const pipos = new Array(n);
    const axisFallback = new Array(n);
    const hostPos = new Vec3();
    const ajPos = new Vec3();
    const dir = new Vec3();
    const acc = new Vec3();
    const cr = new Vec3();
    const tmpU = new Vec3();
    const tmpV = new Vec3();
    const tmpA = new Vec3();
    const tmpB = new Vec3();
    const tmpC = new Vec3();
    for (let ia = 0; ia < n; ia++) {
        const out = new Vec3(0, 0, 0);
        pipos[ia] = out;
        const npi = npiList[ia] | 0;
        if (npi <= 0) continue;
        const neighs = bondsAdj1[ia] || [];
        const nbSorted = [];
        for (let k = 0; k < neighs.length; k++) nbSorted.push(neighs[k][0] | 0);
        nbSorted.sort((a, b) => (a | 0) - (b | 0));

        // Python parity: initial pi dir from AtomicSystem.get_atomi_pi_direction()
        tmpC.set(0, 0, 0);
        const rawNorm = computeAtomiPiDirectionFromNeighs(mol, ia, nbSorted, tmpA, tmpB, tmpC);
        if (rawNorm > 1e-12) out.setV(tmpC.mulScalar(1.0 / rawNorm));
        if (nbSorted.length === 1) {
            // store axis fallback (used only if pi remains unresolved)
            hostPos.setV(mol.atoms[ia].pos);
            const jb = nbSorted[0] | 0;
            if (jb >= 0 && jb < n) {
                ajPos.setV(mol.atoms[jb].pos);
                dir.setV(ajPos).sub(hostPos);
                if (dir.normalize() > 0) axisFallback[ia] = dir.clone();
            }
        }
        if (debugTopo) debugTopo(`PI_INIT ia=${String(ia).padStart(4)} nb=${String(nbSorted.length).padStart(2)} npi=${String(npi).padStart(2)} neigh=${JSON.stringify(nbSorted)} raw=(${(rawNorm > 1e-12 ? out.x : 0).toFixed(6)},${(rawNorm > 1e-12 ? out.y : 0).toFixed(6)},${(rawNorm > 1e-12 ? out.z : 0).toFixed(6)}) pi=(${out.x.toFixed(6)},${out.y.toFixed(6)},${out.z.toFixed(6)})`);
    }

    // Python parity: propagate low-norm pi dirs from neighbors (min_norm=0.7, max_iter=4)
    const minNorm = 0.7;
    for (let it = 0; it < 4; it++) {
        let updated = 0;
        for (let ia = 0; ia < n; ia++) {
            if ((npiList[ia] | 0) <= 0) continue;
            const vHost = pipos[ia];
            const hostNorm = vHost.norm();
            if (hostNorm >= minNorm) continue;
            acc.set(0, 0, 0);
            const neighs = bondsAdj1[ia] || [];
            const nbSorted = [];
            for (let k = 0; k < neighs.length; k++) nbSorted.push(neighs[k][0] | 0);
            nbSorted.sort((a, b) => (a | 0) - (b | 0));
            for (let k = 0; k < nbSorted.length; k++) {
                const ja = nbSorted[k] | 0;
                if (ja < 0 || ja >= n) continue;
                if ((npiList[ja] | 0) <= 0) continue;
                const vj = pipos[ja];
                const nj = vj.norm();
                if (nj < 1e-6) continue;
                tmpU.setV(vj).mulScalar(1.0 / nj);
                if (hostNorm > 1e-6 && vHost.dot(tmpU) < 0) tmpU.mulScalar(-1);
                else if (acc.norm() > 1e-6 && acc.dot(tmpU) < 0) tmpU.mulScalar(-1);
                acc.add(tmpU);
            }
            const accNorm = acc.norm();
            if (accNorm >= minNorm) {
                vHost.setV(acc.mulScalar(1.0 / accNorm));
                if (debugTopo) debugTopo(`PI_PROP ia=${String(ia).padStart(4)} neigh=${JSON.stringify(nbSorted)} acc=(${acc.x.toFixed(6)},${acc.y.toFixed(6)},${acc.z.toFixed(6)}) acc_norm=${accNorm.toFixed(6)} pi=(${vHost.x.toFixed(6)},${vHost.y.toFixed(6)},${vHost.z.toFixed(6)})`);
                updated++;
            }
        }
        if (updated === 0) break;
    }

    if (alignPiVectors) {
        // Python parity: optional align signs (max_iter=6) against first available pi neighbor
        for (let it = 0; it < 6; it++) {
            let flippedAny = 0;
            for (let ia = 0; ia < n; ia++) {
                if ((npiList[ia] | 0) <= 0) continue;
                const vHost = pipos[ia];
                const hostNorm = vHost.norm();
                if (hostNorm < 1e-6) continue;
                const neighs = bondsAdj1[ia] || [];
                const nbSorted = [];
                for (let k = 0; k < neighs.length; k++) nbSorted.push(neighs[k][0] | 0);
                nbSorted.sort((a, b) => (a | 0) - (b | 0));
                for (let k = 0; k < nbSorted.length; k++) {
                    const ja = nbSorted[k] | 0;
                    if (ja < 0 || ja >= n) continue;
                    if ((npiList[ja] | 0) <= 0) continue;
                    const vj = pipos[ja];
                    const nj = vj.norm();
                    if (nj < 1e-6) continue;
                    tmpU.setV(vj).mulScalar(1.0 / nj);
                    const dot0 = vHost.dot(tmpU);
                    if (dot0 < 0) {
                        vHost.mulScalar(-1);
                        if (debugTopo) debugTopo(`PI_FLIP ia=${String(ia).padStart(4)} ja=${String(ja).padStart(4)} dot=${dot0.toFixed(6)} pi=(${vHost.x.toFixed(6)},${vHost.y.toFixed(6)},${vHost.z.toFixed(6)})`);
                        flippedAny++;
                    }
                    break;
                }
            }
            if (flippedAny === 0) break;
        }
    }

    // Last resort fallback for unresolved 1-neighbor cases: choose any perpendicular dir
    for (let ia = 0; ia < n; ia++) {
        if ((npiList[ia] | 0) <= 0) continue;
        if (pipos[ia].norm() > 1e-8) continue;
        const axis = axisFallback[ia];
        if (!axis) continue;
        orthonormalBasisFromDir(axis, tmpU, tmpV);
        pipos[ia].setV(tmpU);
    }

    for (let ia = 0; ia < n; ia++) {
        const npi = npiList[ia] | 0;
        if (npi <= 0) continue;
        const v = pipos[ia];
        const nn = v ? v.norm() : 0;
        if (debugTopo) debugTopo(`PI_FINAL ia=${String(ia).padStart(4)} npi=${String(npi).padStart(2)} pi=(${v.x.toFixed(6)},${v.y.toFixed(6)},${v.z.toFixed(6)}) norm=${nn.toFixed(6)}`);
    }
    return pipos;
}

function orthonormalBasisFromDir(dir, outU, outV) {
    const a = Math.abs(dir.x), b = Math.abs(dir.y), c = Math.abs(dir.z);
    const tmp = (a < 0.9) ? new Vec3(1, 0, 0) : ((b < 0.9) ? new Vec3(0, 1, 0) : new Vec3(0, 0, 1));
    outU.setV(tmp);
    outU.subMul(dir, outU.dot(dir));
    if (!(outU.normalize() > 0)) throw new Error('orthonormalBasisFromDir: failed');
    outV.setCross(dir, outU);
    if (!(outV.normalize() > 0)) throw new Error('orthonormalBasisFromDir: failed (v)');
}

function collectNeighborDirs(mol, ia, neighs, limit = 3) {
    const dirs = [];
    const hostPos = mol.atoms[ia].pos;
    const tmp = new Vec3();
    const nbSorted = [];
    for (let k = 0; k < neighs.length; k++) nbSorted.push(neighs[k][0] | 0);
    nbSorted.sort((a, b) => (a | 0) - (b | 0));
    for (let k = 0; k < nbSorted.length && dirs.length < limit; k++) {
        const jb = nbSorted[k] | 0;
        const atomB = mol.atoms[jb];
        if (!atomB) continue;
        tmp.setV(atomB.pos);
        tmp.sub(hostPos);
        if (tmp.normalize() > 1e-8) dirs.push(tmp.clone());
    }
    return dirs;
}

function epairDirsFromMMFF(mol, ia, neighs, pipos, npi, nep, debugTopo = null) {
    const dirs = collectNeighborDirs(mol, ia, neighs, 3);
    const nb = dirs.length;
    const out = [];
    if (nep <= 0) return out;
    const eps = 1e-8;

    if (npi === 0) {
        if (nb === 3) {
            const v10 = dirs[1].clone().sub(dirs[0]);
            const v20 = dirs[2].clone().sub(dirs[0]);
            const base = new Vec3().setCross(v10, v20);
            if (!(base.normalize() > eps)) base.set(0, 0, 1);
            const sum = dirs[0].clone().add(dirs[1]).add(dirs[2]);
            if (sum.dot(base) > 0) base.mulScalar(-1);
            out.push(base.clone());
        } else if (nb === 2) {
            const m_c = dirs[0].clone().add(dirs[1]);
            if (!(m_c.normalize() > eps)) m_c.set(0, 0, 1);
            const m_b = new Vec3().setCross(dirs[0], dirs[1]);
            if (!(m_b.normalize() > eps)) {
                const u = new Vec3();
                const v = new Vec3();
                orthonormalBasisFromDir(m_c, u, v);
                m_b.setV(u);
            }
            const cc = 0.57735026919;
            const cb = 0.81649658092;
            const d1 = m_c.clone().mulScalar(-cc).addMul(m_b, cb);
            const d2 = m_c.clone().mulScalar(-cc).addMul(m_b, -cb);
            if (d1.normalize() > eps) out.push(d1);
            if (d2.normalize() > eps) out.push(d2);
        }
    } else if (npi === 1) {
        if (nb === 2) {
            const m_c = dirs[0].clone().add(dirs[1]);
            if (!(m_c.normalize() > eps)) m_c.set(0, 0, 1);
            const neg = m_c.clone().mulScalar(-1);
            if (neg.normalize() > eps) out.push(neg);
        } else if (nb === 1) {
            const v1 = dirs[0].clone();
            const piDir = (pipos && pipos[ia]) ? pipos[ia] : null;
            const m_b = piDir ? piDir.clone() : new Vec3(0, 0, 1);
            if (!(m_b.normalize() > eps)) m_b.set(0, 0, 1);
            const m_c = new Vec3().setCross(v1, m_b);
            if (!(m_c.normalize() > eps)) {
                const u = new Vec3();
                const v = new Vec3();
                orthonormalBasisFromDir(m_b, u, v);
                m_c.setV(u);
            }
            const d1 = v1.clone().mulScalar(-0.5).addMul(m_c, +0.86602540378);
            const d2 = v1.clone().mulScalar(-0.5).addMul(m_c, -0.86602540378);
            if (d1.normalize() > eps) out.push(d1);
            if (d2.normalize() > eps) out.push(d2);
        }
    }

    if (out.length > nep) out.length = nep;
    if (debugTopo && out.length > 0) {
        const ss = out.map(v => `(${v.x.toFixed(6)},${v.y.toFixed(6)},${v.z.toFixed(6)})`).join(' ');
        debugTopo(`EP_DIR ia=${String(ia).padStart(4)} nb=${String(nb).padStart(2)} npi=${String(npi).padStart(2)} nep=${String(nep).padStart(2)} dirs=${ss}`);
    }
    return out;
}

function buildPiDummies(mol, mmParams, bondsAdj1, npiList, pipos, opts, outLinear, outPiDummies, debugTopo = null) {
    const L_pi = (opts.L_pi !== undefined) ? +opts.L_pi : 1.0;
    const K_pi_host = (opts.k_pi !== undefined) ? +opts.k_pi : 0.0;
    const K_pi_orth = (opts.k_pi_orth !== undefined) ? +opts.k_pi_orth : 0.0;
    const twoPi = !!opts.two_pi;

    const nBefore = mol.atoms.length;
    const tmp = new Vec3();
    for (let ia = 0; ia < nBefore; ia++) {
        const npi = npiList[ia] | 0;
        if (npi <= 0) continue;
        const dir = pipos[ia];
        if (!(dir.norm() > 1e-6)) continue;

        const host = mol.atoms[ia];
        const hostPos = host.pos;
        const signs = twoPi ? [1.0, -1.0] : [1.0];

        for (let si = 0; si < signs.length; si++) {
            const sign = signs[si];
            tmp.setV(dir);
            tmp.mulScalar(L_pi * sign);
            tmp.add(hostPos);

            const d = extendMolWithDummy(mol, 201, tmp, -1);
            outPiDummies.push({ index: d.i, host: ia, sign });
            if (debugTopo) debugTopo(`PI_DUMMY idx=${String(d.i).padStart(4)} host=${String(ia).padStart(4)} sign=${sign.toFixed(1)} pos=(${tmp.x.toFixed(6)},${tmp.y.toFixed(6)},${tmp.z.toFixed(6)})`);

            outLinear.push([ia, d.i, L_pi, K_pi_host, ['pi-host', [ia, d.i]]]);

            const neighs = bondsAdj1[ia] || [];
            for (let k = 0; k < neighs.length; k++) {
                const jb = neighs[k][0] | 0;
                const symH = Z_TO_SYMBOL[host.Z] || 'X';
                const symB = Z_TO_SYMBOL[mol.atoms[jb].Z] || 'X';
                const tH = mmParams.resolveTypeNameTable(symH);
                const tB = mmParams.resolveTypeNameTable(symB);
                const bp = mmParams.getBondParams(tH, tB);
                const rab = (bp && bp.l0 > 0) ? bp.l0 : 1.5;
                const l0 = Math.sqrt(L_pi * L_pi + rab * rab);
                outLinear.push([d.i, jb, l0, K_pi_orth, ['pi-orth', [d.i, jb]]]);
            }
        }
    }
}

function buildEpairDummies(mol, mmParams, bondsAdj1, npiList, nepList, pipos, opts, outLinear, outEpDummies, debugTopo = null) {
    const L_epair = (opts.L_epair !== undefined) ? +opts.L_epair : 0.5;
    const K_ep_host = (opts.k_ep !== undefined) ? +opts.k_ep : 0.0;
    const K_ep_orth = (opts.k_ep_orth !== undefined) ? +opts.k_ep_orth : 0.0;

    const tmp = new Vec3();
    const dirs = [];

    const nBefore = mol.atoms.length;
    for (let ia = 0; ia < nBefore; ia++) {
        const nep = nepList[ia] | 0;
        if (nep <= 0) continue;
        const host = mol.atoms[ia];
        const neighs = bondsAdj1[ia] || [];
        const npi = npiList[ia] | 0;

        const dirList = epairDirsFromMMFF(mol, ia, neighs, pipos, npi, nep, debugTopo);
        dirs.length = 0;
        for (let k = 0; k < dirList.length; k++) dirs.push(dirList[k]);
        for (let k = 0; k < dirs.length; k++) {
            const ddir = dirs[k];
            tmp.setV(ddir);
            tmp.mulScalar(L_epair);
            tmp.add(host.pos);
            const d = extendMolWithDummy(mol, 200, tmp, -1);
            outEpDummies.push({ index: d.i, host: ia, slot: k });
            if (debugTopo) debugTopo(`EP_DUMMY idx=${String(d.i).padStart(4)} host=${String(ia).padStart(4)} slot=${String(k).padStart(2)} pos=(${tmp.x.toFixed(6)},${tmp.y.toFixed(6)},${tmp.z.toFixed(6)})`);

            outLinear.push([ia, d.i, L_epair, K_ep_host, ['ep-host', [ia, d.i]]]);
            for (let kk = 0; kk < neighs.length; kk++) {
                const jb = neighs[kk][0] | 0;
                const symH = Z_TO_SYMBOL[host.Z] || 'X';
                const symB = Z_TO_SYMBOL[mol.atoms[jb].Z] || 'X';
                const tH = mmParams.resolveTypeNameTable(symH);
                const tB = mmParams.resolveTypeNameTable(symB);
                const bp = mmParams.getBondParams(tH, tB);
                const rab = (bp && bp.l0 > 0) ? bp.l0 : 1.5;
                const l0 = Math.sqrt(L_epair * L_epair + rab * rab);
                outLinear.push([d.i, jb, l0, K_ep_orth, ['ep-orth', [d.i, jb]]]);
            }
        }
    }
}

function buildEpairPairBonds(mol, epDummies, opts, outLinear) {
    const K_ep_pair = (opts.k_ep_pair !== undefined) ? +opts.k_ep_pair : 0.0;
    if (!(K_ep_pair !== 0.0)) return;
    const byHost = new Map();
    for (let i = 0; i < epDummies.length; i++) {
        const d = epDummies[i];
        const h = d.host | 0;
        if (!byHost.has(h)) byHost.set(h, []);
        byHost.get(h).push(d.index | 0);
    }
    for (const [host, idxs] of byHost.entries()) {
        if (idxs.length < 2) continue;
        for (let a = 0; a < idxs.length; a++) {
            for (let b = a + 1; b < idxs.length; b++) {
                const ia = idxs[a] | 0;
                const ib = idxs[b] | 0;
                const pa = mol.atoms[ia].pos;
                const pb = mol.atoms[ib].pos;
                const dx = pa.x - pb.x, dy = pa.y - pb.y, dz = pa.z - pb.z;
                const l0 = Math.sqrt(dx * dx + dy * dy + dz * dz);
                outLinear.push([ia, ib, l0, K_ep_pair, ['ep-pair', [host, ia, ib]]]);
            }
        }
    }
}

function buildPiAlignBonds(mol, piDummies, primaryPairs, opts, outLinear) {
    const K_pi_align = (opts.k_pi_align !== undefined) ? +opts.k_pi_align : 0.0;
    if (!(K_pi_align !== 0.0)) return;
    const byHost = new Map();
    for (let i = 0; i < piDummies.length; i++) {
        const d = piDummies[i];
        const h = d.host | 0;
        if (!byHost.has(h)) byHost.set(h, []);
        byHost.get(h).push(d.index | 0);
    }
    for (let k = 0; k < primaryPairs.length; k++) {
        const a = primaryPairs[k][0] | 0;
        const b = primaryPairs[k][1] | 0;
        const la = byHost.get(a);
        const lb = byHost.get(b);
        if (!la || !lb) continue;
        for (let ia = 0; ia < la.length; ia++) {
            for (let ib = 0; ib < lb.length; ib++) {
                const da = la[ia] | 0;
                const db = lb[ib] | 0;
                const pa = mol.atoms[da].pos;
                const pb = mol.atoms[db].pos;
                const dx = pa.x - pb.x, dy = pa.y - pb.y, dz = pa.z - pb.z;
                const l0 = Math.sqrt(dx * dx + dy * dy + dz * dz);
                outLinear.push([da, db, l0, K_pi_align, ['pi-align', [a, b, da, db]]]);
            }
        }
    }
}

export function buildMMFFLTopology(mol, mmParams, opts = {}) {
    if (!mol || !mmParams) throw new Error('buildMMFFLTopology: mol and mmParams required');
    if (!mol.atoms || mol.atoms.length === 0) throw new Error('buildMMFFLTopology: empty molecule');
    mmParams.ensureTypeIndex();

    // Ensure neighbor lists are correct.
    if (typeof mol.updateNeighborList === 'function') mol.updateNeighborList();

    const nReal = mol.atoms.length;

    // Assign type names from table for real atoms.
    const typeNamesReal = new Array(nReal);
    for (let i = 0; i < nReal; i++) {
        const sym = Z_TO_SYMBOL[mol.atoms[i].Z] || 'X';
        typeNamesReal[i] = mmParams.resolveTypeNameTable(sym);
    }

    // Primary bonds list (pairs)
    const primaryPairs = [];
    for (let i = 0; i < mol.bonds.length; i++) {
        const b = mol.bonds[i];
        b.ensureIndices(mol);
        addUndirectedPair(primaryPairs, b.a | 0, b.b | 0);
    }

    // Primary adjacency for real atoms
    const bondsAdj1 = new Array(nReal);
    for (let i = 0; i < nReal; i++) bondsAdj1[i] = [];
    for (let k = 0; k < primaryPairs.length; k++) {
        const a = primaryPairs[k][0] | 0;
        const b = primaryPairs[k][1] | 0;
        const ta = typeNamesReal[a];
        const tb = typeNamesReal[b];
        const bp = mmParams.getBondParams(ta, tb);
        const l0 = (bp && bp.l0 > 0) ? bp.l0 : ((opts.defaultL !== undefined) ? +opts.defaultL : 1.5);
        const K = (bp && bp.k > 0) ? bp.k : ((opts.defaultK !== undefined) ? +opts.defaultK : 100.0);
        bondsAdj1[a].push([b, l0, K]);
        bondsAdj1[b].push([a, l0, K]);
    }

    const addAngle = (opts.add_angle !== undefined) ? !!opts.add_angle : true;
    const addPi = (opts.add_pi !== undefined) ? !!opts.add_pi : false;
    const addPiAlign = (opts.add_pi_align !== undefined) ? !!opts.add_pi_align : false;
    const addEpair = (opts.add_epair !== undefined) ? !!opts.add_epair : false;
    const addEpairPairs = (opts.add_epair_pairs !== undefined) ? !!opts.add_epair_pairs : false;

    // npi/nep derived from AtomTypes + nbond (same as python report_types)
    const npiList = new Array(nReal);
    const nepList = new Array(nReal);
    for (let i = 0; i < nReal; i++) {
        const t = mmParams.atomTypes[typeNamesReal[i]];
        const nbond = countNeighborsFromAdj(bondsAdj1[i]);
        const nepair = (t && (t.nepair !== undefined)) ? (t.nepair | 0) : 0;
        let npi = (t && (t.valence !== undefined)) ? ((t.valence | 0) - nbond) : 0;
        if (npi < 0) npi = 0;
        npiList[i] = npi;
        nepList[i] = nepair;
    }

    if (!!opts.report_types) {
        for (let i = 0; i < nReal; i++) {
            const sym = Z_TO_SYMBOL[mol.atoms[i].Z] || 'X';
            const nbond = countNeighborsFromAdj(bondsAdj1[i]);
            const npi = npiList[i] | 0;
            const nepair = nepList[i] | 0;
            console.log(`TYPE ia=${String(i).padStart(4)} ename=${sym} type=${String(typeNamesReal[i]).padEnd(8)} nbond=${nbond} npi=${npi} nepair=${nepair}`);
        }
    }

    const debugTopo = (opts && (typeof opts.debugTopo === 'function')) ? opts.debugTopo : null;
    const alignPiVectors = !!(opts && opts.align_pi_vectors);
    const pipos = computePiOrientations(mol, bondsAdj1, npiList, alignPiVectors, debugTopo);

    const linearBonds = [];
    const piDummies = [];
    const epDummies = [];

    if (addAngle) buildAngleBonds(mol, mmParams, bondsAdj1, opts, linearBonds);
    if (addPi) buildPiDummies(mol, mmParams, bondsAdj1, npiList, pipos, opts, linearBonds, piDummies, debugTopo);
    if (addEpair) buildEpairDummies(mol, mmParams, bondsAdj1, npiList, nepList, pipos, opts, linearBonds, epDummies, debugTopo);
    if (addEpair && addEpairPairs) buildEpairPairBonds(mol, epDummies, opts, linearBonds);
    if (addPi && addPiAlign) buildPiAlignBonds(mol, piDummies, uniqPairs(primaryPairs), opts, linearBonds);

    const nAll = mol.atoms.length;
    const isDummy = new Array(nAll);
    const typeNamesAll = new Array(nAll);
    for (let i = 0; i < nAll; i++) {
        isDummy[i] = (i >= nReal);
        if (i < nReal) typeNamesAll[i] = typeNamesReal[i];
        else {
            const Z = mol.atoms[i].Z | 0;
            typeNamesAll[i] = (Z === 200) ? 'E' : 'Pi';
        }
    }

    const bondsPrimary = uniqPairs(primaryPairs);
    const bondsAngle = [];
    const bondsPi = [];
    const bondsPiAlign = [];
    const bondsEp = [];
    const bondsEpPair = [];

    for (let i = 0; i < linearBonds.length; i++) {
        const tag = linearBonds[i][4];
        const kind = Array.isArray(tag) ? String(tag[0]) : String(tag);
        const a = linearBonds[i][0] | 0;
        const b = linearBonds[i][1] | 0;
        if (kind === 'angle') addUndirectedPair(bondsAngle, a, b);
        else if (kind.startsWith('pi-align')) addUndirectedPair(bondsPiAlign, a, b);
        else if (kind.startsWith('pi')) addUndirectedPair(bondsPi, a, b);
        else if (kind.startsWith('ep-pair')) addUndirectedPair(bondsEpPair, a, b);
        else if (kind.startsWith('ep')) addUndirectedPair(bondsEp, a, b);
    }

    return {
        apos: mol.atoms.map(a => [a.pos.x, a.pos.y, a.pos.z]),
        type_names: typeNamesAll,
        is_dummy: isDummy,
        bonds_primary: uniqPairs(bondsPrimary),
        bonds_angle: uniqPairs(bondsAngle),
        bonds_pi: uniqPairs(bondsPi),
        bonds_pi_align: uniqPairs(bondsPiAlign),
        bonds_epair: uniqPairs(bondsEp),
        bonds_epair_pair: uniqPairs(bondsEpPair),
        bonds_linear: linearBonds,
        n_real: nReal,
        n_all: nAll,
        pi_dummies: piDummies,
        ep_dummies: epDummies,
    };
}
