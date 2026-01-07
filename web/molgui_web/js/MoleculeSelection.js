import { EditableMolecule } from "./EditableMolecule.js";

/// Parse a count-set string like "{1,2}" into a Set of ints.
export function parseCountSet(s) {
    const t = String(s || '').trim();
    if (!t.startsWith('{') || !t.endsWith('}')) throw new Error(`selectQuery: expected count set {..}, got '${t}'`);
    const inner = t.slice(1, -1).trim();
    if (!inner) return new Set();
    const parts = inner.split(',');
    const out = new Set();
    for (let i = 0; i < parts.length; i++) {
        const x = parseInt(parts[i].trim(), 10);
        if (!isFinite(x)) throw new Error(`selectQuery: invalid count '${parts[i]}'`);
        out.add(x | 0);
    }
    return out;
}

/// Compile a token-set matcher for elements/atom-types (or '*' wildcard).
export function compileTokenSetToMatcher(mmParams, tokSetStr, asZ) {
    const s = String(tokSetStr || '').trim();
    if (!s) throw new Error('selectQuery: empty token set');
    if (s === '*') return { bAny: true, zSet: null, tSet: null, fnMatch: (_a) => true };
    const toks = s.split('|').map(x => x.trim()).filter(x => x.length > 0);
    if (toks.length === 0) throw new Error('selectQuery: empty token set');
    const zSet = new Set();
    const tSet = new Set();
    for (let i = 0; i < toks.length; i++) {
        const tok = toks[i];
        if (tok === '*') return { bAny: true, zSet: null, tSet: null, fnMatch: (_a) => true };
        if (tok.indexOf('_') >= 0) {
            const it = mmParams.atomTypeNameToIndex(tok);
            if (it < 0) throw new Error(`selectQuery: unknown atom type '${tok}'`);
            tSet.add(it | 0);
        } else {
            const z = asZ(tok);
            if (!(z > 0)) throw new Error(`selectQuery: invalid element '${tok}'`);
            zSet.add(z | 0);
        }
    }
    const fnMatch = (a) => {
        if (!a) return false;
        const iz = a.Z | 0;
        if (zSet.size && zSet.has(iz)) return true;
        const it = (a.atype !== undefined) ? (a.atype | 0) : -1;
        if (tSet.size && it >= 0 && tSet.has(it)) return true;
        return false;
    };
    return { bAny: false, zSet, tSet, fnMatch };
}

/// Compile a full selection query (atom token + neighbor count constraints).
export function compileSelectQuerySpec(mmParams, q, asZ) {
    const s = String(q || '').trim();
    if (!s) throw new Error('selectQuery: empty query');
    const parts = s.split(/\s+/).filter(x => x.length > 0);
    if (parts.length === 0) throw new Error('selectQuery: empty query');
    const atomTok = parts[0];
    const atomMatcher = compileTokenSetToMatcher(mmParams, atomTok, asZ);
    const constraints = [];
    for (let i = 1; i < parts.length; i++) {
        const p = parts[i];
        if (!p) continue;
        if (p.startsWith('n{')) {
            const j = p.indexOf('}');
            if (j < 0) throw new Error(`selectQuery: missing '}' in '${p}'`);
            const neiSetStr = p.slice(2, j + 1);
            const rest = p.slice(j + 1);
            if (!rest.startsWith('={')) throw new Error(`selectQuery: expected '={' after n{..}, got '${p}'`);
            const cntSetStr = rest.slice(1);
            const cntSet = parseCountSet(cntSetStr);
            const neiInner = neiSetStr.slice(2, -1).trim();
            const neiMatcher = compileTokenSetToMatcher(mmParams, neiInner || '*', asZ);
            constraints.push({ neiMatcher, cntSet, src: p });
            continue;
        }
        if (p.startsWith('deg')) {
            const rest = p.slice(3);
            if (!rest.startsWith('={')) throw new Error(`selectQuery: expected deg={..}, got '${p}'`);
            const cntSet = parseCountSet(rest.slice(1));
            const neiMatcher = compileTokenSetToMatcher(mmParams, '*', asZ);
            constraints.push({ neiMatcher, cntSet, src: p });
            continue;
        }
        throw new Error(`selectQuery: cannot parse token '${p}'`);
    }
    return { q: s, atomMatcher, constraints };
}

/// Apply compiled selection query to a molecule-like object (atoms/bonds/selection).
export function applySelectQuery(mol, compiled, opts = {}) {
    const mode = (opts.mode !== undefined) ? String(opts.mode) : 'replace';
    const bPrint = !!opts.bPrint;
    if (!compiled || !compiled.atomMatcher || !compiled.constraints) throw new Error('applySelectQuery: invalid compiled query');
    const atomMatch = compiled.atomMatcher.fnMatch;
    const constraints = compiled.constraints;
    const sel = mol.selection;
    if (mode === 'replace') sel.clear();
    const nAtoms = mol.atoms.length | 0;
    const counts = new Int32Array(constraints.length | 0);
    let nHit = 0;
    for (let ia = 0; ia < nAtoms; ia++) {
        const a = mol.atoms[ia];
        if (!a) continue;
        if (!atomMatch(a)) continue;
        for (let k = 0; k < counts.length; k++) counts[k] = 0;
        const bs = a.bonds;
        for (let ib0 = 0; ib0 < bs.length; ib0++) {
            const ib = bs[ib0] | 0;
            const bnd = mol.bonds[ib];
            if (!bnd) continue;
            bnd.ensureIndices(mol);
            const jb = bnd.other(ia);
            if (jb < 0) continue;
            const nb = mol.atoms[jb];
            if (!nb) continue;
            for (let ic = 0; ic < constraints.length; ic++) {
                if (constraints[ic].neiMatcher.fnMatch(nb)) counts[ic]++;
            }
        }
        let ok = true;
        for (let ic = 0; ic < constraints.length; ic++) {
            const cs = constraints[ic].cntSet;
            if (cs.size && !cs.has(counts[ic] | 0)) { ok = false; break; }
        }
        if (!ok) continue;
        const id = a.id;
        if (mode === 'subtract') sel.delete(id);
        else sel.add(id);
        nHit++;
        if (bPrint) console.log(`selectQuery hit id=${id} i=${ia} Z=${a.Z}`);
    }
    mol.dirtyExport = true;
    return { nHit };
}

/// Install selection helpers onto EditableMolecule (uses mol.mmParams and cls.asZ).
export function installMoleculeSelectionMethods(cls = EditableMolecule) {
    console.log('[installMoleculeSelectionMethods] installing on', cls && cls.name);
    if (!cls || !cls.prototype) throw new Error('installMoleculeSelectionMethods: invalid class');
    cls.prototype.compileSelectQuerySpec = function (q, opts = {}) {
        const mmParams = opts.mmParams || this.mmParams;
        if (!mmParams) throw new Error('compileSelectQuerySpec: mmParams missing (set mol.mmParams or pass opts.mmParams)');
        const asZ = opts.asZ || cls.asZ;
        if (!asZ) throw new Error('compileSelectQuerySpec: asZ not available on class');
        return compileSelectQuerySpec(mmParams, q, asZ);
    };
    cls.prototype.applySelectQuery = function (compiled, opts = {}) {
        return applySelectQuery(this, compiled, opts);
    };
    // Override static stubs so static calls also work after install
    cls.compileSelectQuery = function (q, mmParams, asZFn) {
        const asZ = asZFn || cls.asZ;
        return cls.prototype.compileSelectQuerySpec.call({ mmParams }, q, { mmParams, asZ });
    };
    cls.applySelectQuery = function (compiled, opts = {}) {
        // Use the instance method; 'this' should be an EditableMolecule instance when called as instance
        if (this && this !== cls) {
            return cls.prototype.applySelectQuery.call(this, compiled, opts);
        }
        // Fallback: require an instance via opts.mol
        if (opts && opts.mol) return cls.prototype.applySelectQuery.call(opts.mol, compiled, opts);
        throw new Error('applySelectQuery: instance required');
    };
}
