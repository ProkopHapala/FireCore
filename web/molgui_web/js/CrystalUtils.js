import { Vec3 } from "../../common_js/Vec3.js";
import { Mat3 } from "../../common_js/Mat3.js";
import { EditableMolecule } from "./EditableMolecule.js";

export function parseFraction(x) {
    if (typeof x === 'number') return x;
    if (typeof x !== 'string') throw new Error(`parseFraction: expected number|string, got ${typeof x}`);
    const s = x.trim();
    if (!s) throw new Error('parseFraction: empty string');
    const i = s.indexOf('/');
    if (i >= 0) {
        const a = parseFloat(s.slice(0, i));
        const b = parseFloat(s.slice(i + 1));
        if (!(isFinite(a) && isFinite(b)) || Math.abs(b) < 1e-16) throw new Error(`parseFraction: bad fraction '${x}'`);
        return a / b;
    }
    const v = parseFloat(s);
    if (!isFinite(v)) throw new Error(`parseFraction: bad number '${x}'`);
    return v;
}

export function wrapFrac(x) {
    let y = x % 1.0;
    if (y < 0) y += 1.0;
    return y;
}

function _gridKey(ix, iy, iz) { return `${ix}|${iy}|${iz}`; }

function _dedupInsertOrGet(pos, Z, q, buckets, atoms, tol2, h, fnAddAtom = null) {
    const ix = Math.floor(pos.x / h) | 0;
    const iy = Math.floor(pos.y / h) | 0;
    const iz = Math.floor(pos.z / h) | 0;
    for (let dz = -1; dz <= 1; dz++) {
        for (let dy = -1; dy <= 1; dy++) {
            for (let dx = -1; dx <= 1; dx++) {
                const key = _gridKey(ix + dx, iy + dy, iz + dz);
                const list = buckets.get(key);
                if (!list) continue;
                for (let k = 0; k < list.length; k++) {
                    const ia = list[k] | 0;
                    const a = atoms[ia];
                    if ((a.Z | 0) !== (Z | 0)) {
                        const d0x = a.pos.x - pos.x;
                        const d0y = a.pos.y - pos.y;
                        const d0z = a.pos.z - pos.z;
                        const r02 = d0x * d0x + d0y * d0y + d0z * d0z;
                        if (r02 <= tol2) throw new Error(`dedup: collapsing atoms of different elements Z=${Z} and Z=${a.Z} at r=${Math.sqrt(r02)}`);
                        continue;
                    }
                    const dx0 = a.pos.x - pos.x;
                    const dy0 = a.pos.y - pos.y;
                    const dz0 = a.pos.z - pos.z;
                    const r2 = dx0 * dx0 + dy0 * dy0 + dz0 * dz0;
                    if (r2 <= tol2) {
                        if (q !== undefined && q !== null) {
                            const aq = (a.charge !== undefined) ? +a.charge : 0.0;
                            if (Math.abs(aq - (+q)) > 1e-6) throw new Error(`dedup: collapsing atoms with different charges q=${q} vs qKeep=${aq} (Z=${Z})`);
                        }
                        return { keepIndex: ia, keepId: a.id };
                    }
                }
            }
        }
    }
    const key0 = _gridKey(ix, iy, iz);
    let list0 = buckets.get(key0);
    if (!list0) { list0 = []; buckets.set(key0, list0); }
    const idNew = fnAddAtom ? fnAddAtom(Z) : null;
    if (idNew !== null) {
        const iaNew = atoms.length - 1;
        list0.push(iaNew);
        return { keepIndex: iaNew, keepId: idNew };
    }
    return { keepIndex: -1, keepId: -1 };
}

export function dedupFracSitesByTolA(sites, lvec, tolA = 0.1) {
    const tol = +tolA;
    if (!(tol > 0)) throw new Error('dedupFracSitesByTolA: tolA must be >0');
    if (!Array.isArray(sites)) throw new Error('dedupFracSitesByTolA: sites must be array');
    if (!lvec || lvec.length !== 3) throw new Error('dedupFracSitesByTolA: lvec must be Vec3[3]');
    const tol2 = tol * tol;
    const h = tol;
    const buckets = new Map();
    const out = [];
    const tmp = new Vec3();
    const pos = new Vec3();
    for (let i = 0; i < sites.length; i++) {
        const s = sites[i];
        const Z = EditableMolecule.asZ(s.element);
        fracToCart([+s.x, +s.y, +s.z], lvec, pos);
        const ix = Math.floor(pos.x / h) | 0;
        const iy = Math.floor(pos.y / h) | 0;
        const iz = Math.floor(pos.z / h) | 0;
        let found = -1;
        for (let dz = -1; dz <= 1 && found < 0; dz++) {
            for (let dy = -1; dy <= 1 && found < 0; dy++) {
                for (let dx = -1; dx <= 1 && found < 0; dx++) {
                    const key = _gridKey(ix + dx, iy + dy, iz + dz);
                    const list = buckets.get(key);
                    if (!list) continue;
                    for (let k = 0; k < list.length; k++) {
                        const j = list[k] | 0;
                        const t = out[j];
                        const Zt = EditableMolecule.asZ(t.element);
                        fracToCart([+t.x, +t.y, +t.z], lvec, tmp);
                        const dx0 = tmp.x - pos.x;
                        const dy0 = tmp.y - pos.y;
                        const dz0 = tmp.z - pos.z;
                        const r2 = dx0 * dx0 + dy0 * dy0 + dz0 * dz0;
                        if (r2 > tol2) continue;
                        if (Zt !== Z) throw new Error(`dedupFracSitesByTolA: collapsing different elements '${s.element}' vs '${t.element}' (r=${Math.sqrt(r2)})`);
                        const qs = (s.q !== undefined && s.q !== null) ? +s.q : 0.0;
                        const qt = (t.q !== undefined && t.q !== null) ? +t.q : 0.0;
                        if (Math.abs(qs - qt) > 1e-6) throw new Error(`dedupFracSitesByTolA: collapsing different charges q=${qs} vs qKeep=${qt} (el=${s.element})`);
                        found = j;
                        break;
                    }
                }
            }
        }
        if (found >= 0) continue;
        const key0 = _gridKey(ix, iy, iz);
        let list0 = buckets.get(key0);
        if (!list0) { list0 = []; buckets.set(key0, list0); }
        list0.push(out.length);
        out.push({ element: s.element, x: +s.x, y: +s.y, z: +s.z, q: (s.q !== undefined && s.q !== null) ? +s.q : 0.0 });
    }
    return out;
}

export function dedupMolAtomsByTolA(mol, tolA = 0.1, opts = {}) {
    const tol = +tolA;
    if (!(tol > 0)) throw new Error('dedupMolAtomsByTolA: tolA must be >0');
    if (!mol || !mol.atoms || !mol.removeAtomById) throw new Error('dedupMolAtomsByTolA: mol must be EditableMolecule-like');
    const bPrint = !!opts.bPrint;
    const bError = !!opts.bError;
    const tol2 = tol * tol;
    const h = tol;
    const buckets = new Map();
    const removeIds = [];
    const dupToKeep = new Map();

    const atoms = mol.atoms;
    for (let i = 0; i < atoms.length; i++) {
        const a = atoms[i];
        if (!a) continue;
        if (dupToKeep.has(a.id)) continue;

        const p = a.pos;
        const ix = Math.floor(p.x / h) | 0;
        const iy = Math.floor(p.y / h) | 0;
        const iz = Math.floor(p.z / h) | 0;
        let keep = -1;
        let bestR2 = 1e300;

        for (let dz = -1; dz <= 1; dz++) {
            for (let dy = -1; dy <= 1; dy++) {
                for (let dx = -1; dx <= 1; dx++) {
                    const key = _gridKey(ix + dx, iy + dy, iz + dz);
                    const list = buckets.get(key);
                    if (!list) continue;
                    for (let k = 0; k < list.length; k++) {
                        const j = list[k] | 0;
                        const b = atoms[j];
                        if (!b) continue;
                        if ((b.Z | 0) !== (a.Z | 0)) continue;
                        const dx0 = b.pos.x - p.x;
                        const dy0 = b.pos.y - p.y;
                        const dz0 = b.pos.z - p.z;
                        const r2 = dx0 * dx0 + dy0 * dy0 + dz0 * dz0;
                        if (r2 > tol2) continue;
                        const aq = (a.charge !== undefined) ? +a.charge : 0.0;
                        const bq = (b.charge !== undefined) ? +b.charge : 0.0;
                        if (Math.abs(aq - bq) > 1e-6) throw new Error(`dedupMolAtomsByTolA: collapsing atoms with different charges q=${aq} vs qKeep=${bq} (Z=${a.Z})`);
                        if (r2 < bestR2) { bestR2 = r2; keep = j; }
                    }
                }
            }
        }

        if (keep >= 0) {
            const keepAtom = atoms[keep];
            if (!keepAtom) throw new Error('dedupMolAtomsByTolA: internal error (keep atom missing)');
            dupToKeep.set(a.id, keepAtom.id);
            removeIds.push(a.id);
            if (bPrint) console.log(`dedupMolAtomsByTolA dup id=${a.id} keepId=${keepAtom.id} i=${i} keepI=${keep} Z=${a.Z} r=${Math.sqrt(bestR2)}`);
            if (bError) throw new Error(`dedupMolAtomsByTolA: found duplicate atoms id=${a.id} keepId=${keepAtom.id} r=${Math.sqrt(bestR2)}`);
            continue;
        }

        const key0 = _gridKey(ix, iy, iz);
        let list0 = buckets.get(key0);
        if (!list0) { list0 = []; buckets.set(key0, list0); }
        list0.push(i);
    }

    for (let k = 0; k < removeIds.length; k++) mol.removeAtomById(removeIds[k]);
    return { nRemoved: removeIds.length, dupToKeep };
}

export function latticeVectorsFromParams(params) {
    if (!params) throw new Error('latticeVectorsFromParams: params required');
    const a = +params.a;
    const b = (params.b !== undefined) ? +params.b : a;
    const c = (params.c !== undefined) ? +params.c : a;
    const alphaDeg = (params.alpha !== undefined) ? +params.alpha : 90;
    const betaDeg = (params.beta !== undefined) ? +params.beta : 90;
    const gammaDeg = (params.gamma !== undefined) ? +params.gamma : 90;
    if (!(a > 0 && b > 0 && c > 0)) throw new Error(`latticeVectorsFromParams: invalid a,b,c = ${a},${b},${c}`);

    const toRad = Math.PI / 180.0;
    const alpha = alphaDeg * toRad;
    const beta = betaDeg * toRad;
    const gamma = gammaDeg * toRad;
    const ca = Math.cos(alpha), cb = Math.cos(beta), cg = Math.cos(gamma);
    const sg = Math.sin(gamma);
    if (Math.abs(sg) < 1e-12) throw new Error(`latticeVectorsFromParams: sin(gamma) too small gamma=${gammaDeg}`);

    const va = new Vec3(a, 0, 0);
    const vb = new Vec3(b * cg, b * sg, 0);
    const cx = c * cb;
    const cy = c * (ca - cb * cg) / sg;
    const cz2 = c * c - cx * cx - cy * cy;
    if (cz2 < -1e-10) throw new Error(`latticeVectorsFromParams: invalid angles -> cz^2=${cz2}`);
    const vc = new Vec3(cx, cy, Math.sqrt(Math.max(0, cz2)));
    return [va, vb, vc];
}

export function fracToCart(p, lvec, out = new Vec3()) {
    if (!lvec || lvec.length !== 3) throw new Error('fracToCart: lvec must be Vec3[3]');
    out.set(0, 0, 0);
    out.addMul(lvec[0], p[0]);
    out.addMul(lvec[1], p[1]);
    out.addMul(lvec[2], p[2]);
    return out;
}

export function cartToFrac(p, lvec, out = [0, 0, 0]) {
    if (!lvec || lvec.length !== 3) throw new Error('cartToFrac: lvec must be Vec3[3]');
    const a = lvec[0], b = lvec[1], c = lvec[2];
    const bxc = new Vec3().setCross(b, c);
    const cxa = new Vec3().setCross(c, a);
    const axb = new Vec3().setCross(a, b);
    const V = a.dot(bxc);
    if (Math.abs(V) < 1e-12) throw new Error('cartToFrac: singular cell');
    const invV = 1.0 / V;
    out[0] = (p.x * bxc.x + p.y * bxc.y + p.z * bxc.z) * invV;
    out[1] = (p.x * cxa.x + p.y * cxa.y + p.z * cxa.z) * invV;
    out[2] = (p.x * axb.x + p.y * axb.y + p.z * axb.z) * invV;
    return out;
}

export function parseLatticeText(txt) {
    const lines = String(txt || '').replace(/\r/g, '').split('\n').map(s => s.trim()).filter(s => s && !s.startsWith('#'));
    if (lines.length < 3) throw new Error('parseLatticeText: expected at least 3 lines (a,b,c vectors)');
    const v = [];
    for (let i = 0; i < 3; i++) {
        const toks = lines[i].split(/\s+/).filter(t => t.length);
        if (toks.length < 3) throw new Error(`parseLatticeText: line ${i + 1} has <3 numbers`);
        const x = parseFloat(toks[0]), y = parseFloat(toks[1]), z = parseFloat(toks[2]);
        if (!(isFinite(x) && isFinite(y) && isFinite(z))) throw new Error(`parseLatticeText: bad numbers on line ${i + 1}`);
        v.push(new Vec3(x, y, z));
    }
    return v;
}

export function formatLatticeText(lvec) {
    if (!lvec || lvec.length !== 3) throw new Error('formatLatticeText: lvec must be Vec3[3]');
    const f = (x) => {
        if (!Number.isFinite(x)) return String(x);
        let s = x.toFixed(6);
        if (!s.includes('.')) s += '.000000';
        return s;
    };
    const pad = (s, w = 12) => (s.length >= w) ? s : (' '.repeat(w - s.length) + s);
    const row = (v) => `${pad(f(v.x))} ${pad(f(v.y))} ${pad(f(v.z))}`;
    return `${row(lvec[0])}\n${row(lvec[1])}\n${row(lvec[2])}\n`;
}

export function parseSitesTextXYZ(txt, coordMode = 'frac') {
    const lines = String(txt || '').replace(/\r/g, '').split('\n');
    const out = [];
    for (let il = 0; il < lines.length; il++) {
        const l0 = lines[il].trim();
        if (!l0 || l0.startsWith('#')) continue;
        const toks = l0.split(/\s+/).filter(t => t.length);
        if (toks.length < 4) throw new Error(`parseSitesTextXYZ: line ${il + 1} expected: El x y z [q]`);
        const element = toks[0];
        const x = parseFloat(toks[1]), y = parseFloat(toks[2]), z = parseFloat(toks[3]);
        if (!(isFinite(x) && isFinite(y) && isFinite(z))) throw new Error(`parseSitesTextXYZ: bad numbers on line ${il + 1}`);
        let q = 0.0;
        if (toks.length >= 5) {
            q = parseFloat(toks[4]);
            if (!isFinite(q)) throw new Error(`parseSitesTextXYZ: bad charge on line ${il + 1}`);
        }
        out.push({ element, x, y, z, q, coordMode });
    }
    return out;
}

export function formatSitesTextXYZFrac(sites) {
    if (!Array.isArray(sites)) throw new Error('formatSitesTextXYZFrac: sites must be array');
    const f = (x) => {
        if (!Number.isFinite(x)) return String(x);
        let s = x.toFixed(6);
        if (!s.includes('.')) s += '.000000';
        return s;
    };
    const pad = (s, w = 11) => (s.length >= w) ? s : (' '.repeat(w - s.length) + s);
    const hasQ = sites.some(a => (a.q !== undefined && a.q !== null && Math.abs(+a.q) > 1e-12));
    const padQ = (s, w = 10) => (s.length >= w) ? s : (' '.repeat(w - s.length) + s);
    let s = '';
    for (const a of sites) {
        const el = String(a.element || '').padEnd(3, ' ');
        const line = `${el} ${pad(f(+a.x))} ${pad(f(+a.y))} ${pad(f(+a.z))}`;
        if (hasQ) {
            const q = (a.q !== undefined && a.q !== null) ? +a.q : 0.0;
            s += `${line} ${padQ(f(q))}\n`;
        } else {
            s += `${line}\n`;
        }
    }
    return s;
}

export function formatSitesTextXYZCart(sitesFrac, lvec) {
    if (!Array.isArray(sitesFrac)) throw new Error('formatSitesTextXYZCart: sitesFrac must be array');
    if (!lvec || lvec.length !== 3) throw new Error('formatSitesTextXYZCart: lvec must be Vec3[3]');
    const f = (x) => {
        if (!Number.isFinite(x)) return String(x);
        let s = x.toFixed(6);
        if (!s.includes('.')) s += '.000000';
        return s;
    };
    const pad = (s, w = 11) => (s.length >= w) ? s : (' '.repeat(w - s.length) + s);
    const p = new Vec3();
    const hasQ = sitesFrac.some(a => (a.q !== undefined && a.q !== null && Math.abs(+a.q) > 1e-12));
    const padQ = (s, w = 10) => (s.length >= w) ? s : (' '.repeat(w - s.length) + s);
    let s = '';
    for (const a of sitesFrac) {
        fracToCart([+a.x, +a.y, +a.z], lvec, p);
        const el = String(a.element || '').padEnd(3, ' ');
        const line = `${el} ${pad(f(p.x))} ${pad(f(p.y))} ${pad(f(p.z))}`;
        if (hasQ) {
            const q = (a.q !== undefined && a.q !== null) ? +a.q : 0.0;
            s += `${line} ${padQ(f(q))}\n`;
        } else {
            s += `${line}\n`;
        }
    }
    return s;
}

export function parseSymOpsText(txt) {
    const lines = String(txt || '').replace(/\r/g, '').split('\n');
    const out = [];
    for (let il = 0; il < lines.length; il++) {
        const l0 = lines[il].trim();
        if (!l0 || l0.startsWith('#')) continue;
        out.push(parseSymOpXYZ(l0));
    }
    return out;
}

export function formatSymOpsText(symOpsStrs) {
    if (!Array.isArray(symOpsStrs)) return '';
    let s = '';
    for (const op of symOpsStrs) s += `${_stripQuotes(op)}\n`;
    return s;
}

function _stripQuotes(s) {
    const t = (s !== undefined && s !== null) ? String(s).trim() : '';
    if (t.length >= 2) {
        const a = t.charCodeAt(0);
        const b = t.charCodeAt(t.length - 1);
        if ((a === 39 && b === 39) || (a === 34 && b === 34)) return t.slice(1, -1);
    }
    return t;
}

function _tokenizeCIFLine(line) {
    const re = /'[^']*'|"[^"]*"|\S+/g;
    const out = [];
    let m;
    while ((m = re.exec(line)) !== null) out.push(m[0]);
    return out;
}

export function parseCIF(cifText) {
    const lines = String(cifText).replace(/\r/g, '').split('\n');
    const data = { tags: {}, loops: [] };
    let i = 0;
    while (i < lines.length) {
        let line = lines[i];
        if (line === undefined) break;
        line = line.trim();
        if (!line || line.startsWith('#')) { i++; continue; }
        if (line.toLowerCase().startsWith('data_')) { i++; continue; }
        if (line.toLowerCase() === 'loop_') {
            i++;
            const headers = [];
            while (i < lines.length) {
                const l = (lines[i] || '').trim();
                if (!l || l.startsWith('#')) { i++; continue; }
                if (!l.startsWith('_')) break;
                headers.push(l.split(/\s+/)[0]);
                i++;
            }
            if (headers.length === 0) throw new Error('parseCIF: loop_ without headers');
            const rows = [];
            while (i < lines.length) {
                const raw = lines[i] || '';
                const l = raw.trim();
                if (!l || l.startsWith('#')) { i++; if (!l) break; continue; }
                const ll = l.toLowerCase();
                if (l.startsWith('_') || ll === 'loop_' || ll.startsWith('data_')) break;
                const toks = _tokenizeCIFLine(raw).map(_stripQuotes);
                if (toks.length === 0) { i++; continue; }
                if (toks.length < headers.length) throw new Error(`parseCIF: loop row has ${toks.length} tokens, expected >= ${headers.length}`);
                rows.push(toks.slice(0, headers.length));
                i++;
            }
            data.loops.push({ headers, rows });
            continue;
        }
        if (line.startsWith('_')) {
            const toks = _tokenizeCIFLine(lines[i]);
            const key = toks[0];
            const val = toks.slice(1).join(' ');
            data.tags[key] = _stripQuotes(val);
            i++;
            continue;
        }
        i++;
    }
    return data;
}

function _parseSymExpr(expr) {
    const s = _stripQuotes(expr).replace(/\s+/g, '');
    let cx = 0, cy = 0, cz = 0, c0 = 0;
    const parts = s.match(/[+-]?[^+-]+/g) || [];
    for (const p of parts) {
        if (!p) continue;
        if (p === 'x' || p === '+x') cx += 1;
        else if (p === '-x') cx -= 1;
        else if (p === 'y' || p === '+y') cy += 1;
        else if (p === '-y') cy -= 1;
        else if (p === 'z' || p === '+z') cz += 1;
        else if (p === '-z') cz -= 1;
        else c0 += parseFraction(p);
    }
    return { c: [cx, cy, cz], t: c0 };
}

export function parseSymOpXYZ(opStr) {
    const s = _stripQuotes(opStr);
    const parts = s.split(',');
    if (parts.length !== 3) throw new Error(`parseSymOpXYZ: expected 3 comma-separated parts, got '${opStr}'`);
    const ex = _parseSymExpr(parts[0]);
    const ey = _parseSymExpr(parts[1]);
    const ez = _parseSymExpr(parts[2]);
    const m = new Int8Array(9);
    m[0] = ex.c[0]; m[1] = ex.c[1]; m[2] = ex.c[2];
    m[3] = ey.c[0]; m[4] = ey.c[1]; m[5] = ey.c[2];
    m[6] = ez.c[0]; m[7] = ez.c[1]; m[8] = ez.c[2];
    const t = [wrapFrac(ex.t), wrapFrac(ey.t), wrapFrac(ez.t)];
    return { m, t };
}

export function applySymmetryOpsFracSites(sites, symOps, opts = {}) {
    const tol = (opts.tol !== undefined) ? +opts.tol : 1e-6;
    if (!(tol > 0)) throw new Error('applySymmetryOpsFracSites: tol must be >0');
    if (!Array.isArray(sites)) throw new Error('applySymmetryOpsFracSites: sites must be array');
    if (!Array.isArray(symOps) || symOps.length === 0) throw new Error('applySymmetryOpsFracSites: symOps must be non-empty array');
    const inv = 1.0 / tol;
    const out = [];
    const seen = new Set();
    for (const a of sites) {
        const el = a.element;
        const x = +a.x, y = +a.y, z = +a.z;
        const q = (a.q !== undefined && a.q !== null) ? +a.q : 0.0;
        for (const op of symOps) {
            const m = op.m;
            const t = op.t;
            const fx = wrapFrac(m[0] * x + m[1] * y + m[2] * z + t[0]);
            const fy = wrapFrac(m[3] * x + m[4] * y + m[5] * z + t[1]);
            const fz = wrapFrac(m[6] * x + m[7] * y + m[8] * z + t[2]);
            const ix = Math.round(fx * inv) | 0;
            const iy = Math.round(fy * inv) | 0;
            const iz = Math.round(fz * inv) | 0;
            const key = `${el}|${ix}|${iy}|${iz}`;
            if (seen.has(key)) continue;
            seen.add(key);
            out.push({ element: el, x: fx, y: fy, z: fz, q });
        }
    }
    return out;
}

function _normalizeElementSymbol(s) {
    const t = _stripQuotes(s);
    if (!t) return t;
    const m = t.match(/[A-Za-z]+/);
    if (!m) return t;
    const r = m[0];
    if (r.length === 1) return r.toUpperCase();
    return r[0].toUpperCase() + r.slice(1).toLowerCase();
}

export function cifToCrystalData(cifText) {
    const cif = parseCIF(cifText);
    const tags = cif.tags;
    const lat = {
        a: parseFloat(tags._cell_length_a),
        b: parseFloat(tags._cell_length_b),
        c: parseFloat(tags._cell_length_c),
        alpha: parseFloat(tags._cell_angle_alpha),
        beta: parseFloat(tags._cell_angle_beta),
        gamma: parseFloat(tags._cell_angle_gamma)
    };
    if (!(lat.a > 0 && lat.b > 0 && lat.c > 0)) throw new Error('cifToCrystalData: missing/invalid cell lengths');
    const spg = (tags._symmetry_Int_Tables_number !== undefined) ? (parseInt(tags._symmetry_Int_Tables_number) | 0) : 0;

    let symOpStrs = null;
    let siteHeaders = null;
    let siteRows = null;
    for (const loop of cif.loops) {
        const hs = loop.headers;
        if (hs.includes('_symmetry_equiv_pos_as_xyz')) symOpStrs = loop.rows.map(r => r[hs.indexOf('_symmetry_equiv_pos_as_xyz')]);
        if (hs.includes('_space_group_symop_operation_xyz')) symOpStrs = loop.rows.map(r => r[hs.indexOf('_space_group_symop_operation_xyz')]);
        if (hs.includes('_atom_site_fract_x') && hs.includes('_atom_site_fract_y') && hs.includes('_atom_site_fract_z')) {
            siteHeaders = hs;
            siteRows = loop.rows;
        }
    }
    if (!siteRows) throw new Error('cifToCrystalData: missing atom site loop (_atom_site_fract_x/y/z)');
    const ix = siteHeaders.indexOf('_atom_site_fract_x');
    const iy = siteHeaders.indexOf('_atom_site_fract_y');
    const iz = siteHeaders.indexOf('_atom_site_fract_z');
    const itype = siteHeaders.indexOf('_atom_site_type_symbol');
    if (ix < 0 || iy < 0 || iz < 0) throw new Error('cifToCrystalData: atom site loop missing fract columns');
    const sites = [];
    for (const r of siteRows) {
        const el = _normalizeElementSymbol((itype >= 0) ? r[itype] : '');
        if (!el) continue;
        sites.push({ element: el, x: parseFraction(r[ix]), y: parseFraction(r[iy]), z: parseFraction(r[iz]) });
    }
    const symOps = (symOpStrs && symOpStrs.length) ? symOpStrs.map(parseSymOpXYZ) : [];
    return { lattice: lat, spg, symOps, symOpStrs: symOpStrs || [], sites };
}

export function cellDataFromFracSites(lvec, sites) {
    if (!lvec || lvec.length !== 3) throw new Error('cellDataFromFracSites: lvec must be Vec3[3]');
    if (!Array.isArray(sites) || sites.length === 0) throw new Error('cellDataFromFracSites: sites required');
    const n = sites.length | 0;
    const basisPos = new Float32Array(n * 3);
    const basisTypes = new Uint8Array(n);
    const basisCharges = new Float32Array(n);
    const tmp = new Vec3();
    for (let i = 0; i < n; i++) {
        const a = sites[i];
        const Z = EditableMolecule.asZ(a.element);
        basisTypes[i] = Z;
        basisCharges[i] = (a.q !== undefined && a.q !== null) ? +a.q : 0.0;
        fracToCart([+a.x, +a.y, +a.z], lvec, tmp);
        const i3 = i * 3;
        basisPos[i3] = tmp.x;
        basisPos[i3 + 1] = tmp.y;
        basisPos[i3 + 2] = tmp.z;
    }
    return { lvec, basisPos, basisTypes, basisCharges };
}

export function crystalSitesToCellData(crystal, opts = {}) {
    if (!crystal) throw new Error('crystalSitesToCellData: crystal required');
    const lvec = latticeVectorsFromParams(crystal.lattice);
    const sites = crystal.sites;
    if (!Array.isArray(sites) || sites.length === 0) throw new Error('crystalSitesToCellData: no sites');
    const n = sites.length | 0;
    const basisPos = new Float32Array(n * 3);
    const basisTypes = new Uint8Array(n);
    const tmp = new Vec3();
    for (let i = 0; i < n; i++) {
        const a = sites[i];
        const Z = EditableMolecule.asZ(a.element);
        basisTypes[i] = Z;
        fracToCart([a.x, a.y, a.z], lvec, tmp);
        const i3 = i * 3;
        basisPos[i3] = tmp.x;
        basisPos[i3 + 1] = tmp.y;
        basisPos[i3 + 2] = tmp.z;
    }
    return { lvec, basisPos, basisTypes };
}

export function genCrystalFromCIF(cifText, params = {}) {
    const nRep = (params.nRep !== undefined) ? params.nRep : [1, 1, 1];
    const origin = (params.origin !== undefined) ? params.origin : new Vec3(0, 0, 0);
    const applySym = (params.applySymmetry !== undefined) ? !!params.applySymmetry : false;
    const dedupSym = (params.dedupSymmetry !== undefined) ? !!params.dedupSymmetry : false;
    const dedupTol = (params.dedupTol !== undefined) ? +params.dedupTol : 0.1;
    const crystal = cifToCrystalData(cifText);
    if (applySym) {
        if (!crystal.symOps || crystal.symOps.length === 0) throw new Error('genCrystalFromCIF: no symmetry operations in CIF');
        crystal.sites = applySymmetryOpsFracSites(crystal.sites, crystal.symOps, { tol: (params.tol !== undefined) ? params.tol : undefined });
    }
    let cell = crystalSitesToCellData(crystal);
    if (dedupSym) {
        const sitesDedup = dedupFracSitesByTolA(crystal.sites, cell.lvec, dedupTol);
        crystal.sites = sitesDedup;
        cell = crystalSitesToCellData(crystal);
    }
    let mol = null;
    if (params.slab) {
        if (params.buildBonds) throw new Error('genCrystalFromCIF: buildBonds not supported with slab cut');
        const hkl = params.slab.hkl;
        if (!hkl || hkl.length < 3) throw new Error('genCrystalFromCIF: slab.hkl must be [h,k,l]');
        const cmin = +params.slab.cmin;
        const cmax = +params.slab.cmax;
        if (!(cmax > cmin)) throw new Error('genCrystalFromCIF: slab requires cmax>cmin');
        const b = reciprocalLattice(cell.lvec);
        const n = new Vec3().setLincomb3(hkl[0] | 0, b[0], hkl[1] | 0, b[1], hkl[2] | 0, b[2]);
        const ln = n.normalize();
        if (!(ln > 0)) throw new Error('genCrystalFromCIF: slab normal is zero');
        mol = genReplicatedCellSlab({ lvec: cell.lvec, basisPos: cell.basisPos, basisTypes: cell.basisTypes, nRep, origin, nHat: n, cmin, cmax, dedup: !!params.dedupReplicate, dedupTol });
    } else {
        mol = genReplicatedCell({ lvec: cell.lvec, basisPos: cell.basisPos, basisTypes: cell.basisTypes, nRep, origin, buildBonds: !!params.buildBonds, mmParams: params.mmParams, bondTol: params.bondTol, dedup: !!params.dedupReplicate, dedupTol });
    }
    mol.lvec = [cell.lvec[0].clone().mulScalar(nRep[0] | 0), cell.lvec[1].clone().mulScalar(nRep[1] | 0), cell.lvec[2].clone().mulScalar(nRep[2] | 0)];
    return { mol, lvec: mol.lvec, spg: crystal.spg, lattice: crystal.lattice };
}

export function mpJsonToCellData(mpJson, opts = {}) { throw new Error('mpJsonToCellData: MP JSON support was removed from MolGUI Web; use CIF/unit-cell editor instead'); }

export function genCrystalFromMPJson(mpJson, params = {}) { throw new Error('genCrystalFromMPJson: MP JSON support was removed from MolGUI Web; use CIF/unit-cell editor instead'); }

export function genReplicatedCellSlab(params = {}) {
    const { lvec, basisPos, basisTypes, basisCharges = null, nRep = [1, 1, 1], origin = new Vec3(0, 0, 0), nHat = null, cmin = 0.0, cmax = 0.0, dedup = false, dedupTol = 0.1 } = params;
    if (!lvec || !basisPos || !basisTypes) throw new Error('genReplicatedCellSlab: lvec, basisPos, basisTypes required');
    const na = nRep[0] | 0;
    const nb = nRep[1] | 0;
    const nc = nRep[2] | 0;
    if (na <= 0 || nb <= 0 || nc <= 0) throw new Error(`genReplicatedCellSlab: invalid nRep=${nRep}`);
    if (!(cmax > cmin)) throw new Error('genReplicatedCellSlab: requires cmax>cmin');
    if (!(lvec[0] instanceof Vec3) || !(lvec[1] instanceof Vec3) || !(lvec[2] instanceof Vec3)) throw new Error('genReplicatedCellSlab: lvec must be Vec3[3]');
    if (!(origin instanceof Vec3)) throw new Error('genReplicatedCellSlab: origin must be Vec3');
    if (!(nHat instanceof Vec3)) throw new Error('genReplicatedCellSlab: nHat must be Vec3');

    const n2 = nHat.x * nHat.x + nHat.y * nHat.y + nHat.z * nHat.z;
    if (Math.abs(n2 - 1.0) > 1e-6) throw new Error('genReplicatedCellSlab: nHat must be normalized');

    const nBasis = basisTypes.length | 0;
    const basisIsFlat = (basisPos.length | 0) === (nBasis * 3);
    const basisIsVec3 = (basisPos.length | 0) === nBasis && (basisPos[0] instanceof Vec3);
    if (!basisIsFlat && !basisIsVec3) throw new Error('genReplicatedCellSlab: basisPos must be flat array length 3*nBasis or Vec3[] length nBasis');
    if (basisCharges && ((basisCharges.length | 0) !== nBasis)) throw new Error('genReplicatedCellSlab: basisCharges length mismatch');

    const a = lvec[0], b = lvec[1], c = lvec[2];
    const da = a.dot(nHat);
    const db = b.dot(nHat);
    const dc = c.dot(nHat);
    const dMin = (da < 0 ? da : 0) + (db < 0 ? db : 0) + (dc < 0 ? dc : 0);
    const dMax = (da > 0 ? da : 0) + (db > 0 ? db : 0) + (dc > 0 ? dc : 0);

    const mol = new EditableMolecule();
    const tol = +dedupTol;
    const doDedup = !!dedup;
    const tol2 = tol * tol;
    const h = tol;
    const buckets = doDedup ? new Map() : null;
    const s = new Vec3();
    for (let iz = 0; iz < nc; iz++) {
        for (let iy = 0; iy < nb; iy++) {
            for (let ix = 0; ix < na; ix++) {
                s.setV(origin);
                s.addMul(a, ix);
                s.addMul(b, iy);
                s.addMul(c, iz);
                const c0 = s.dot(nHat);
                const cellCmin = c0 + dMin;
                const cellCmax = c0 + dMax;
                if (cellCmax < cmin || cellCmin > cmax) continue;
                for (let ib = 0; ib < nBasis; ib++) {
                    const Z = basisTypes[ib] | 0;
                    const q = basisCharges ? +basisCharges[ib] : 0.0;
                    let x = 0, y = 0, z = 0;
                    if (basisIsVec3) {
                        const p = basisPos[ib];
                        x = p.x; y = p.y; z = p.z;
                    } else {
                        const b3 = ib * 3;
                        x = basisPos[b3]; y = basisPos[b3 + 1]; z = basisPos[b3 + 2];
                    }
                    const cx = s.x + x;
                    const cy = s.y + y;
                    const cz = s.z + z;
                    const cc = cx * nHat.x + cy * nHat.y + cz * nHat.z;
                    if (cc < cmin || cc > cmax) continue;
                    if (doDedup) {
                        const pos = new Vec3(cx, cy, cz);
                        const fnAdd = (Z_) => {
                            const id = mol.addAtom(pos.x, pos.y, pos.z, Z_);
                            const ia = mol.getAtomIndex(id);
                            if (ia >= 0) mol.atoms[ia].charge = q;
                            return id;
                        };
                        const r = _dedupInsertOrGet(pos, Z, q, buckets, mol.atoms, tol2, h, fnAdd);
                        if (r.keepId < 0) throw new Error('genReplicatedCellSlab: dedup internal error');
                    } else {
                        const id = mol.addAtom(cx, cy, cz, Z);
                        const ia = mol.getAtomIndex(id);
                        if (ia >= 0) mol.atoms[ia].charge = q;
                    }
                }
            }
        }
    }
    mol.lvec = [a.clone().mulScalar(na), b.clone().mulScalar(nb), c.clone().mulScalar(nc)];
    return mol;
}

export function genReplicatedCellCutPlanes(params = {}) {
    const { lvec, basisPos, basisTypes, basisCharges = null, nRep = [1, 1, 1], origin = new Vec3(0, 0, 0), planes = null, planeMode = 'ang', centered = false, dedup = false, dedupTol = 0.1 } = params;
    if (!lvec || !basisPos || !basisTypes) throw new Error('genReplicatedCellCutPlanes: lvec, basisPos, basisTypes required');
    if (!Array.isArray(planes) || planes.length === 0) throw new Error('genReplicatedCellCutPlanes: planes must be non-empty array');
    const na = nRep[0] | 0;
    const nb = nRep[1] | 0;
    const nc = nRep[2] | 0;
    if (na <= 0 || nb <= 0 || nc <= 0) throw new Error(`genReplicatedCellCutPlanes: invalid nRep=${nRep}`);
    if (!(lvec[0] instanceof Vec3) || !(lvec[1] instanceof Vec3) || !(lvec[2] instanceof Vec3)) throw new Error('genReplicatedCellCutPlanes: lvec must be Vec3[3]');
    if (!(origin instanceof Vec3)) throw new Error('genReplicatedCellCutPlanes: origin must be Vec3');
    if (planeMode !== 'ang' && planeMode !== 'frac') throw new Error("genReplicatedCellCutPlanes: planeMode must be 'ang' or 'frac'");

    const nBasis = basisTypes.length | 0;
    const basisIsFlat = (basisPos.length | 0) === (nBasis * 3);
    const basisIsVec3 = (basisPos.length | 0) === nBasis && (basisPos[0] instanceof Vec3);
    if (!basisIsFlat && !basisIsVec3) throw new Error('genReplicatedCellCutPlanes: basisPos must be flat array length 3*nBasis or Vec3[] length nBasis');
    if (basisCharges && ((basisCharges.length | 0) !== nBasis)) throw new Error('genReplicatedCellCutPlanes: basisCharges length mismatch');

    const planeDefs = [];
    for (const p of planes) {
        if (!p || !(p.n instanceof Vec3)) throw new Error('genReplicatedCellCutPlanes: each plane requires Vec3 n');
        const cmin = +p.cmin;
        const cmax = +p.cmax;
        if (!(cmax > cmin)) throw new Error('genReplicatedCellCutPlanes: requires cmax>cmin for each plane');
        const n = p.n.clone();
        if (planeMode === 'ang') {
            const ln = n.normalize();
            if (!(ln > 0)) throw new Error('genReplicatedCellCutPlanes: plane normal is zero');
        } else {
            const n2 = n.x * n.x + n.y * n.y + n.z * n.z;
            if (!(n2 > 0)) throw new Error('genReplicatedCellCutPlanes: plane normal is zero');
        }
        planeDefs.push({ n, cmin, cmax });
    }

    const a = lvec[0], b = lvec[1], c = lvec[2];
    const mol = new EditableMolecule();
    const tol = +dedupTol;
    const doDedup = !!dedup;
    const tol2 = tol * tol;
    const h = tol;
    const buckets = doDedup ? new Map() : null;
    const s = new Vec3();

    const ix0 = centered ? -na : 0;
    const ix1 = centered ? na : (na - 1);
    const iy0 = centered ? -nb : 0;
    const iy1 = centered ? nb : (nb - 1);
    const iz0 = centered ? -nc : 0;
    const iz1 = centered ? nc : (nc - 1);
    const na_ = centered ? (2 * na + 1) : na;
    const nb_ = centered ? (2 * nb + 1) : nb;
    const nc_ = centered ? (2 * nc + 1) : nc;

    for (let iz = iz0; iz <= iz1; iz++) {
        for (let iy = iy0; iy <= iy1; iy++) {
            for (let ix = ix0; ix <= ix1; ix++) {
                const icell = ((iz - iz0) * (nb_ * na_) + (iy - iy0) * na_ + (ix - ix0)) | 0; // dense cell index
                s.setV(origin);
                s.addMul(a, ix);
                s.addMul(b, iy);
                s.addMul(c, iz);
                for (let ib = 0; ib < nBasis; ib++) {
                    const Z = basisTypes[ib] | 0;
                    const q = basisCharges ? +basisCharges[ib] : 0.0;
                    let x = 0, y = 0, z = 0;
                    if (basisIsVec3) {
                        const p = basisPos[ib];
                        x = p.x; y = p.y; z = p.z;
                    } else {
                        const b3 = ib * 3;
                        x = basisPos[b3]; y = basisPos[b3 + 1]; z = basisPos[b3 + 2];
                    }
                    const cx = s.x + x;
                    const cy = s.y + y;
                    const cz = s.z + z;
                    let ok = true;
                    for (const pl of planeDefs) {
                        const cc = cx * pl.n.x + cy * pl.n.y + cz * pl.n.z;
                        if (cc < pl.cmin || cc > pl.cmax) { ok = false; break; }
                    }
                    if (!ok) continue;
                    if (doDedup) {
                        const pos = new Vec3(cx, cy, cz);
                        const fnAdd = (Z_) => {
                            const id = mol.addAtom(pos.x, pos.y, pos.z, Z_);
                            const ia = mol.getAtomIndex(id);
                            if (ia >= 0) {
                                mol.atoms[ia].charge = q;
                                mol.atoms[ia].cellIndex = icell;
                            }
                            return id;
                        };
                        const r = _dedupInsertOrGet(pos, Z, q, buckets, mol.atoms, tol2, h, fnAdd);
                        if (r.keepId < 0) throw new Error('genReplicatedCellCutPlanes: dedup internal error');
                    } else {
                        const id = mol.addAtom(cx, cy, cz, Z);
                        const ia = mol.getAtomIndex(id);
                        if (ia >= 0) {
                            mol.atoms[ia].charge = q;
                            mol.atoms[ia].cellIndex = icell;
                        }
                    }
                }
            }
        }
    }
    const sx = centered ? (2 * na + 1) : na;
    const sy = centered ? (2 * nb + 1) : nb;
    const sz = centered ? (2 * nc + 1) : nc;
    mol.lvec = [a.clone().mulScalar(sx), b.clone().mulScalar(sy), c.clone().mulScalar(sz)];
    return mol;
}

function _basisToVec3Array(basisPos, nBasis) {
    const basisIsVec3 = (basisPos.length | 0) === nBasis && (basisPos[0] instanceof Vec3);
    if (basisIsVec3) return basisPos;
    const out = new Array(nBasis);
    for (let i = 0; i < nBasis; i++) {
        const i3 = i * 3;
        out[i] = new Vec3(basisPos[i3], basisPos[i3 + 1], basisPos[i3 + 2]);
    }
    return out;
}

function _computeBasisBonds(lvec, basisPosVec3, basisTypes, mmParams, bondTol = 0.2) {
    if (!mmParams || !mmParams.getBondL0) throw new Error('_computeBasisBonds: mmParams with getBondL0 required');
    const tol = (bondTol !== undefined) ? +bondTol : 0.2;
    if (!(tol >= 0)) throw new Error('_computeBasisBonds: bondTol must be >=0');
    const nBasis = basisTypes.length | 0;
    const a = lvec[0], b = lvec[1], c = lvec[2];
    const off = [];
    for (let iz = -1; iz <= 1; iz++) {
        for (let iy = -1; iy <= 1; iy++) {
            for (let ix = -1; ix <= 1; ix++) {
                if (ix === 0 && iy === 0 && iz === 0) continue;
                if (iz < 0) continue;
                if (iz === 0 && iy < 0) continue;
                if (iz === 0 && iy === 0 && ix < 0) continue;
                off.push([ix, iy, iz]);
            }
        }
    }
    off.push([0, 0, 0]);

    const bonds = [];
    const d = new Vec3();
    const t = new Vec3();
    for (const [ox, oy, oz] of off) {
        t.set(0, 0, 0);
        if (ox) t.addMul(a, ox);
        if (oy) t.addMul(b, oy);
        if (oz) t.addMul(c, oz);
        for (let i = 0; i < nBasis; i++) {
            const zi = basisTypes[i] | 0;
            const pi = basisPosVec3[i];
            const j0 = (ox === 0 && oy === 0 && oz === 0) ? (i + 1) : 0;
            for (let j = j0; j < nBasis; j++) {
                const zj = basisTypes[j] | 0;
                const bt = mmParams.getBondL0(zi, zj);
                if (!bt) continue;
                const l0 = bt.l0;
                const rMax = l0 * (1.0 + tol);
                const rMax2 = rMax * rMax;
                const pj = basisPosVec3[j];
                d.setV(pj);
                d.add(t);
                d.sub(pi);
                const r2 = d.x * d.x + d.y * d.y + d.z * d.z;
                if (r2 <= rMax2) bonds.push({ i, j, ox, oy, oz });
            }
        }
    }
    return bonds;
}

export function genReplicatedCell(params = {}) {
    const { lvec, basisPos, basisTypes, basisCharges = null, nRep = [1, 1, 1], origin = new Vec3(0, 0, 0), buildBonds = false, mmParams = null, bondTol = 0.2, dedup = false, dedupTol = 0.1 } = params;
    if (!lvec || !basisPos || !basisTypes) throw new Error('genReplicatedCell: lvec, basisPos, basisTypes required');
    const na = nRep[0] | 0;
    const nb = nRep[1] | 0;
    const nc = nRep[2] | 0;
    if (na <= 0 || nb <= 0 || nc <= 0) throw new Error(`genReplicatedCell: invalid nRep=${nRep}`);

    if (!(lvec[0] instanceof Vec3) || !(lvec[1] instanceof Vec3) || !(lvec[2] instanceof Vec3)) throw new Error('genReplicatedCell: lvec must be Vec3[3]');
    if (!(origin instanceof Vec3)) throw new Error('genReplicatedCell: origin must be Vec3');

    const nBasis = basisTypes.length | 0;
    const basisIsFlat = (basisPos.length | 0) === (nBasis * 3);
    const basisIsVec3 = (basisPos.length | 0) === nBasis && (basisPos[0] instanceof Vec3);
    if (!basisIsFlat && !basisIsVec3) throw new Error('genReplicatedCell: basisPos must be flat array length 3*nBasis or Vec3[] length nBasis');
    if (basisCharges && ((basisCharges.length | 0) !== nBasis)) throw new Error('genReplicatedCell: basisCharges length mismatch');

    const a = lvec[0], b = lvec[1], c = lvec[2];

    const mol = new EditableMolecule();
    const ids = (buildBonds ? new Array((na * nb * nc) * nBasis) : null);
    const tol = +dedupTol;
    const doDedup = !!dedup;
    const tol2 = tol * tol;
    const h = tol;
    const buckets = doDedup ? new Map() : null;
    const s = new Vec3();
    for (let iz = 0; iz < nc; iz++) {
        for (let iy = 0; iy < nb; iy++) {
            for (let ix = 0; ix < na; ix++) {
                const icell = ((iz * nb + iy) * na + ix) | 0;
                s.setV(origin);
                s.addMul(a, ix);
                s.addMul(b, iy);
                s.addMul(c, iz);
                for (let ib = 0; ib < nBasis; ib++) {
                    const Z = basisTypes[ib] | 0;
                    const q = basisCharges ? +basisCharges[ib] : 0.0;
                    let x = 0, y = 0, z = 0;
                    if (basisIsVec3) {
                        const p = basisPos[ib];
                        x = p.x; y = p.y; z = p.z;
                    } else {
                        const b3 = ib * 3;
                        x = basisPos[b3]; y = basisPos[b3 + 1]; z = basisPos[b3 + 2];
                    }
                    let id = -1;
                    if (doDedup) {
                        const pos = new Vec3(s.x + x, s.y + y, s.z + z);
                        const fnAdd = (Z_) => {
                            const id_ = mol.addAtom(pos.x, pos.y, pos.z, Z_);
                            const ia = mol.getAtomIndex(id_);
                            if (ia >= 0) { mol.atoms[ia].charge = q; mol.atoms[ia].cellIndex = icell; }
                            return id_;
                        };
                        const r = _dedupInsertOrGet(pos, Z, q, buckets, mol.atoms, tol2, h, fnAdd);
                        id = r.keepId;
                        if (id < 0) throw new Error('genReplicatedCell: dedup internal error');
                    } else {
                        id = mol.addAtom(s.x + x, s.y + y, s.z + z, Z);
                        const ia = mol.getAtomIndex(id);
                        if (ia >= 0) { mol.atoms[ia].charge = q; mol.atoms[ia].cellIndex = icell; }
                    }
                    if (ids) {
                        const ic = ((iz * nb + iy) * na + ix) * nBasis + ib;
                        ids[ic] = id;
                    }
                }
            }
        }
    }

    if (buildBonds) {
        if (!mmParams) throw new Error('genReplicatedCell: buildBonds requires mmParams');
        const basisPosVec3 = _basisToVec3Array(basisPos, nBasis);
        const bBonds = _computeBasisBonds(lvec, basisPosVec3, basisTypes, mmParams, bondTol);
        for (let iz = 0; iz < nc; iz++) {
            for (let iy = 0; iy < nb; iy++) {
                for (let ix = 0; ix < na; ix++) {
                    const ic0 = ((iz * nb + iy) * na + ix);
                    for (const bb of bBonds) {
                        const jx = ix + (bb.ox | 0);
                        const jy = iy + (bb.oy | 0);
                        const jz = iz + (bb.oz | 0);
                        if (jx < 0 || jx >= na) continue;
                        if (jy < 0 || jy >= nb) continue;
                        if (jz < 0 || jz >= nc) continue;
                        const ic1 = ((jz * nb + jy) * na + jx);
                        const aId = ids[ic0 * nBasis + (bb.i | 0)];
                        const bId = ids[ic1 * nBasis + (bb.j | 0)];
                        mol.addBond(aId, bId);
                    }
                }
            }
        }
    }
    mol.lvec = [a.clone().mulScalar(na), b.clone().mulScalar(nb), c.clone().mulScalar(nc)];
    return mol;
}

export function genNaClStep(params = {}) {
    const a = (params.a !== undefined) ? +params.a : (5.6413 / 2);
    const nx = (params.nx !== undefined) ? (params.nx | 0) : 13;
    const ny = (params.ny !== undefined) ? (params.ny | 0) : 12;
    const nz = (params.nz !== undefined) ? (params.nz | 0) : 3;
    if (nx <= 0 || ny <= 0 || nz <= 0) throw new Error(`genNaClStep: invalid nx,ny,nz = ${nx},${ny},${nz}`);

    const Lx = a * nx;
    const LxStep = Lx * 0.5;
    const ax = -1.0 / nx;

    const mol = new EditableMolecule();
    for (let iz = 0; iz < nz; iz++) {
        for (let ix = 0; ix < nx; ix++) {
            for (let iy = 0; iy < ny; iy++) {
                let x = ix * a;
                const y = iy * a;
                let z = iz * a;
                let i = ix + iy + iz;
                const z_ = z + ax * x;
                x -= ax * z;
                z = z_;
                if (x > LxStep) {
                    z += a;
                    i += 1;
                }
                const Z = ((i & 1) === 0) ? 11 : 17;
                mol.addAtom(x, y, z, Z);
            }
        }
    }
    mol.lvec = [new Vec3(nx * a, 0, 0), new Vec3(0, ny * a, 0), new Vec3(0, 0, nz * a)];
    return mol;
}

export function reciprocalLattice(lvec) {
    if (!lvec || lvec.length !== 3) throw new Error('reciprocalLattice: lvec must be length 3');
    const a = lvec[0], b = lvec[1], c = lvec[2];
    if (!(a instanceof Vec3) || !(b instanceof Vec3) || !(c instanceof Vec3)) throw new Error('reciprocalLattice: lvec must be Vec3[3]');

    const axb = new Vec3().setCross(a, b);
    const bxc = new Vec3().setCross(b, c);
    const cxa = new Vec3().setCross(c, a);
    const V = a.dot(bxc);
    if (Math.abs(V) < 1e-12) throw new Error('reciprocalLattice: singular cell');
    const s = (2 * Math.PI) / V;
    bxc.mulScalar(s);
    cxa.mulScalar(s);
    axb.mulScalar(s);
    return [bxc, cxa, axb];
}

export function rotationAlignVectorToZ(v) {
    if (!(v instanceof Vec3)) throw new Error('rotationAlignVectorToZ: v must be Vec3');
    return Mat3.alignVectorToZ(v);
}

export function rotateLvec(lvec, R) {
    if (!lvec || lvec.length !== 3) throw new Error('rotateLvec: lvec must be Vec3[3]');
    if (!(R instanceof Mat3)) throw new Error('rotateLvec: R must be Mat3');
    return [R.mulVec(lvec[0]), R.mulVec(lvec[1]), R.mulVec(lvec[2])];
}

export function rotateMoleculeInPlace(mol, R) {
    if (!(R instanceof Mat3)) throw new Error('rotateMoleculeInPlace: R must be Mat3');
    const n = mol.atoms.length;
    const tmp = new Vec3();
    for (let i = 0; i < n; i++) {
        const p = mol.atoms[i].pos;
        tmp.setV(p);
        R.mulVec(tmp, p);
    }
    mol.dirtyGeom = true;
    mol.dirtyExport = true;
}
