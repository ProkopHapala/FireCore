import { Vec3 } from "../common_js/Vec3.js";
import { Z_TO_SYMBOL } from "./MoleculeIO.js";

function lawOfCosines(rab, rbc, cosTheta) {
    const t = rab * rab + rbc * rbc - 2.0 * rab * rbc * cosTheta;
    return Math.sqrt(Math.max(0.0, t));
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

function countNeighborsFromAdj(adj) {
    let c = 0;
    for (let i = 0; i < adj.length; i++) if ((adj[i][0] | 0) >= 0) c++;
    return c;
}

function computePiOrientations(mol, bondsAdj1, npiList) {
    const n = mol.atoms.length;
    const pipos = new Array(n);
    const hostPos = new Vec3();
    const ajPos = new Vec3();
    const dir = new Vec3();
    const acc = new Vec3();
    const cr = new Vec3();
    const tmpU = new Vec3();
    const tmpV = new Vec3();
    for (let ia = 0; ia < n; ia++) {
        const out = new Vec3(0, 0, 0);
        pipos[ia] = out;
        const npi = npiList[ia] | 0;
        if (npi <= 0) continue;
        const neighs = bondsAdj1[ia] || [];
        hostPos.setV(mol.atoms[ia].pos);
        const dirs = [];
        for (let k = 0; k < neighs.length; k++) {
            const jb = neighs[k][0] | 0;
            const atomB = mol.atoms[jb];
            if (!atomB) continue;
            ajPos.setV(atomB.pos);
            dir.setV(ajPos);
            dir.sub(hostPos);
            if (dir.normalize() > 0) dirs.push(dir.clone());
        }
        if (dirs.length >= 2) {
            acc.set(0, 0, 0);
            for (let a = 0; a < dirs.length; a++) {
                for (let b = a + 1; b < dirs.length; b++) {
                    cr.setCross(dirs[a], dirs[b]);
                    const norm = cr.norm();
                    if (norm > 1e-8) acc.add(cr.clone().mulScalar(1.0 / norm));
                }
            }
            const accNorm = acc.norm();
            if (accNorm > 1e-8) {
                out.setV(acc.mulScalar(1.0 / accNorm));
                continue;
            }
        }
        if (dirs.length === 1) {
            const axis = dirs[0];
            orthonormalBasisFromDir(axis, tmpU, tmpV);
            out.setV(tmpU);
            continue;
        }
        out.set(0, 0, 0);
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
    for (let k = 0; k < neighs.length && dirs.length < limit; k++) {
        const jb = neighs[k][0] | 0;
        const atomB = mol.atoms[jb];
        if (!atomB) continue;
        tmp.setV(atomB.pos);
        tmp.sub(hostPos);
        if (tmp.normalize() > 1e-8) dirs.push(tmp.clone());
    }
    return dirs;
}

function epairDirsFromMMFF(mol, ia, neighs, pipos, npi, nep) {
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
                throw new Error(`epairDirsFromMMFF: cannot build orthonormal basis for host ${ia} (nb=1, npi=1)`);
            }
            const c1 = v1.clone().mulScalar(-0.5).addMul(m_c, 0.86602540378);
            const c2 = v1.clone().mulScalar(-0.5).addMul(m_c, -0.86602540378);
            if (c1.normalize() > eps) out.push(c1);
            if (c2.normalize() > eps) out.push(c2);
        }
    }

    if (out.length > nep) out.length = nep;
    return out;
}

function missingDirsVSEPR(vs, nMissing, totalDomains, outDirs) {
    outDirs.length = 0;
    const nb = vs.length | 0;
    if (nMissing <= 0) return outDirs;
    if (nb === 0) {
        if (totalDomains === 4) {
            const a = 0.57735026919;
            outDirs.push(new Vec3( a, a, a));
            outDirs.push(new Vec3( a,-a,-a));
            outDirs.push(new Vec3(-a, a,-a));
            outDirs.push(new Vec3(-a,-a, a));
        } else if (totalDomains === 3) {
            outDirs.push(new Vec3(1, 0, 0));
            outDirs.push(new Vec3(-0.5, 0.86602540378, 0));
            outDirs.push(new Vec3(-0.5, -0.86602540378, 0));
        } else if (totalDomains === 2) {
            outDirs.push(new Vec3(1, 0, 0));
            outDirs.push(new Vec3(-1, 0, 0));
        } else {
            throw new Error(`missingDirsVSEPR: unsupported totalDomains=${totalDomains} for nb=0`);
        }
        while (outDirs.length > nMissing) outDirs.pop();
        return outDirs;
    }

    if (totalDomains === 2) {
        if (nb !== 1 || nMissing !== 1) throw new Error(`missingDirsVSEPR: linear expects nb=1,nMissing=1 got nb=${nb} nMissing=${nMissing}`);
        outDirs.push(vs[0].clone().mulScalar(-1));
        return outDirs;
    }

    if (totalDomains === 3) {
        if (nb === 2 && nMissing === 1) {
            const m = vs[0].clone().add(vs[1]);
            if (!(m.normalize() > 0)) throw new Error('missingDirsVSEPR: nb=2 planar but v1+v2 is zero');
            outDirs.push(m.mulScalar(-1));
            return outDirs;
        }
        if (nb === 1 && nMissing === 2) {
            const axis = vs[0].clone().mulScalar(-1);
            if (!(axis.normalize() > 0)) throw new Error('missingDirsVSEPR: nb=1 planar axis zero');
            const u = new Vec3();
            const v = new Vec3();
            orthonormalBasisFromDir(axis, u, v);
            const ca = -0.5;
            const sa = 0.86602540378;
            outDirs.push(axis.clone().mulScalar(ca).addMul(u, sa));
            outDirs.push(axis.clone().mulScalar(ca).addMul(u, -sa));
            return outDirs;
        }
        throw new Error(`missingDirsVSEPR: unsupported planar nb=${nb} nMissing=${nMissing}`);
    }

    if (totalDomains === 4) {
        if (nb === 3 && nMissing === 1) {
            const m = vs[0].clone().add(vs[1]).add(vs[2]);
            if (!(m.normalize() > 0)) throw new Error('missingDirsVSEPR: nb=3 tetra but sum is zero');
            outDirs.push(m.mulScalar(-1));
            return outDirs;
        }
        if (nb === 2 && nMissing === 2) {
            const m_c = vs[0].clone().add(vs[1]);
            if (!(m_c.normalize() > 0)) throw new Error('missingDirsVSEPR: nb=2 tetra but v1+v2 is zero');
            const m_b = new Vec3().setCross(vs[0], vs[1]);
            if (!(m_b.normalize() > 0)) {
                const u = new Vec3();
                const v = new Vec3();
                orthonormalBasisFromDir(m_c, u, v);
                m_b.setV(u);
            }
            const cc = 0.57735026919;
            const cb = 0.81649658092;
            const d1 = m_c.clone().mulScalar(-cc).addMul(m_b, cb);
            if (!(d1.normalize() > 0)) throw new Error('missingDirsVSEPR: failed normalize tetra dir1');
            const d2 = m_c.clone().mulScalar(-cc).addMul(m_b, -cb);
            if (!(d2.normalize() > 0)) throw new Error('missingDirsVSEPR: failed normalize tetra dir2');
            outDirs.push(d1);
            outDirs.push(d2);
            return outDirs;
        }
        if (nb === 1 && nMissing === 3) {
            const v1 = vs[0];
            const u = new Vec3();
            const v = new Vec3();
            orthonormalBasisFromDir(v1, u, v);
            const a = -1.0 / 3.0;
            const b = Math.sqrt(8.0 / 9.0);
            const c120 = -0.5;
            const s120 = 0.86602540378;
            const u2 = u.clone().mulScalar(c120).addMul(v, s120);
            const u3 = u.clone().mulScalar(c120).addMul(v, -s120);
            const d1 = v1.clone().mulScalar(a).addMul(u, b);
            const d2 = v1.clone().mulScalar(a).addMul(u2, b);
            const d3 = v1.clone().mulScalar(a).addMul(u3, b);
            if (!(d1.normalize() > 0 && d2.normalize() > 0 && d3.normalize() > 0)) throw new Error('missingDirsVSEPR: failed normalize tetra nb=1');
            outDirs.push(d1);
            outDirs.push(d2);
            outDirs.push(d3);
            return outDirs;
        }
        throw new Error(`missingDirsVSEPR: unsupported tetra nb=${nb} nMissing=${nMissing}`);
    }

    throw new Error(`missingDirsVSEPR: unsupported totalDomains=${totalDomains}`);
}

function buildAngleBonds(mol, mmParams, bondsAdj1, opts, outLinear) {
    const nAtoms = mol.atoms.length;
    const defaultL = (opts.defaultL !== undefined) ? +opts.defaultL : 1.5;
    const kAngle = (opts.k_angle !== undefined) ? +opts.k_angle : 0.0;

    for (let b = 0; b < nAtoms; b++) {
        const neighs = bondsAdj1[b] || [];
        if (neighs.length < 2) continue;
        const symB = Z_TO_SYMBOL[mol.atoms[b].Z] || 'X';
        const typeB = mmParams.resolveTypeNameTable(symB);
        const atB = mmParams.atomTypes[typeB];
        if (!atB) continue;
        const thetaDeg = (atB.Ass > 0) ? atB.Ass : 109.5;
        const theta = thetaDeg * Math.PI / 180.0;
        const cosT = Math.cos(theta);
        const K = (kAngle !== 0.0) ? kAngle : +atB.Kss;

        for (let ni = 0; ni < neighs.length; ni++) {
            for (let nj = ni + 1; nj < neighs.length; nj++) {
                let a = neighs[ni][0] | 0;
                let c = neighs[nj][0] | 0;
                if (a === c) continue;
                if (c < a) { const t = a; a = c; c = t; }

                const symA = Z_TO_SYMBOL[mol.atoms[a].Z] || 'X';
                const symC = Z_TO_SYMBOL[mol.atoms[c].Z] || 'X';
                const typeA = mmParams.resolveTypeNameTable(symA);
                const typeC = mmParams.resolveTypeNameTable(symC);
                const btAB = mmParams.getBondParams(typeB, typeA);
                const btBC = mmParams.getBondParams(typeB, typeC);
                const rab = (btAB && btAB.l0 > 0) ? btAB.l0 : defaultL;
                const rbc = (btBC && btBC.l0 > 0) ? btBC.l0 : defaultL;
                const l0 = lawOfCosines(rab, rbc, cosT);
                outLinear.push([a, c, l0, K, ['angle', [a, b, c]]]);
            }
        }
    }
}

function extendMolWithDummy(mol, Z, pos, atypeIndex = -1) {
    const id = mol.addAtom(pos.x, pos.y, pos.z, Z);
    const i = mol.getAtomIndex(id);
    mol.atoms[i].atype = atypeIndex;
    return { id, i };
}

function buildPiDummies(mol, mmParams, bondsAdj1, npiList, pipos, opts, outLinear, outPiDummies) {
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

        for (let s = 0; s < signs.length; s++) {
            const sign = signs[s];
            tmp.setV(dir);
            tmp.mulScalar(L_pi * sign);
            tmp.add(hostPos);

            const d = extendMolWithDummy(mol, 201, tmp, -1);
            outPiDummies.push({ index: d.i, host: ia, sign });

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

function buildEpairDummies(mol, mmParams, bondsAdj1, npiList, nepList, pipos, opts, outLinear, outEpDummies) {
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

        const dirList = epairDirsFromMMFF(mol, ia, neighs, pipos, npi, nep);
        dirs.length = 0;
        for (let k = 0; k < dirList.length; k++) dirs.push(dirList[k]);
        for (let k = 0; k < dirs.length; k++) {
            const ddir = dirs[k];
            tmp.setV(ddir);
            tmp.mulScalar(L_epair);
            tmp.add(host.pos);
            const d = extendMolWithDummy(mol, 200, tmp, -1);
            outEpDummies.push({ index: d.i, host: ia, slot: k });

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

    const pipos = computePiOrientations(mol, bondsAdj1, npiList);

    const linearBonds = [];
    const piDummies = [];
    const epDummies = [];

    if (addAngle) buildAngleBonds(mol, mmParams, bondsAdj1, opts, linearBonds);
    if (addPi) buildPiDummies(mol, mmParams, bondsAdj1, npiList, pipos, opts, linearBonds, piDummies);
    if (addEpair) buildEpairDummies(mol, mmParams, bondsAdj1, npiList, nepList, pipos, opts, linearBonds, epDummies);
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
