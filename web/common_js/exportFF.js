/// Build topology NPZ from mol2 + relaxed xyz: MMFFL topology + collision workgroups + exclusions.
import fs from 'node:fs';
import path from 'node:path';
import { buildMMFFLTopology, packLinearTopologyForGPU, enumerateGraphDistances, exportTopologySpringList } from '../molgui_webgpu/MMFFLTopology.js';
import { writeTopologyNpzFile } from '../molgui_webgpu/LinearizedTopologyNpz.js';
import { exportStickViewerHTML, exportStickViewerNpzLoaderHTML } from '../molgui_webgpu/LinearizedTopologyViewer.js';
import { buildCollisionWorkgroups, buildExclIcol_1_2_3 } from './CollisionWorkgroups.js';
import { readXyzPositions } from './npzIO.js';
import { loadMMParamsFromDir, loadMolFromMol2, loadMol, applyPositions } from './MolIO.js';
import { Z_TO_SYMBOL } from '../molgui_webgpu/MoleculeIO.js';

/// Build topology NPZ from in-memory mol: MMFFL packing + optional collision workgroups / surface AABBs.
export function buildTopologyNpzFromMol({ mol, mmParams, outNpzPath = null, maxNeighbors = 48, addAngle = true, addDihedral = true, K12 = 0, K13 = 0, K14 = 0, groupCap = 32, exportSurfaceGroups = true, source = '', extra = null }) {
    if (!mol || !mol.atoms || mol.atoms.length === 0) throw new Error('buildTopologyNpzFromMol: empty mol');
    if (!mmParams) throw new Error('buildTopologyNpzFromMol: mmParams required');
    const topo = buildMMFFLTopology(mol, mmParams, {
        type_source: 'table',
        add_angle: addAngle,
        add_dihedral: addDihedral,
        add_pi: false,
        add_epair: false,
        K12, K13, K14,
    });
    const packing = packLinearTopologyForGPU(topo, { maxNeighbors });
    const nAtoms = topo.n_real | 0;
    const meta = {
        source: source || '',
        n_bond: topo.bonds_primary.length | 0,
        n_angle: topo.bonds_angle.length | 0,
        n_dihedral: (topo.bonds_dihedral ? topo.bonds_dihedral.length : 0) | 0,
    };
    let collisionExtra = null;
    if (exportSurfaceGroups) {
        const posFlat = new Float64Array(nAtoms * 3);
        for (let i = 0; i < nAtoms; i++) {
            const q = topo.apos[i];
            posFlat[i * 3 + 0] = +q[0];
            posFlat[i * 3 + 1] = +q[1];
            posFlat[i * 3 + 2] = +q[2];
        }
        const radius = new Float64Array(nAtoms);
        for (let i = 0; i < nAtoms; i++) {
            const a = mol.atoms[i];
            const at = mmParams.getAtomTypeForAtom(a);
            const r = (at && at.RvdW > 0) ? at.RvdW : ((at && at.element && at.element.RvdW > 0) ? at.element.RvdW : 1.5);
            if (!(r > 0)) throw new Error(`buildTopologyNpzFromMol: invalid RvdW=${r} at atom ${i}`);
            radius[i] = +r;
        }
        radius._shape = [nAtoms];
        const wg = buildCollisionWorkgroups({ pos: posFlat, mol, radius, groupCap, fillFactor: 0.8 });
        const ex2 = buildExclIcol_1_2_3({ mol, icol: wg.icol, EXCL_MAX: 16, ipbc: 0 });
        wg.group_atoms._shape = [wg.nGroups, wg.groupCap];
        wg.group_nAtoms._shape = [wg.nGroups];
        wg.icol._shape = [nAtoms];
        wg.icolGroup._shape = [nAtoms];
        wg.bbox_min._shape = [wg.nGroups, 3];
        wg.bbox_max._shape = [wg.nGroups, 3];
        ex2.excl._shape = [nAtoms, ex2.EXCL_MAX];
        ex2.excl_count._shape = [nAtoms];
        collisionExtra = {
            radius,
            icol: wg.icol,
            icolGroup: wg.icolGroup,
            group_atoms: wg.group_atoms,
            group_nAtoms: wg.group_nAtoms,
            group_bbox_min: wg.bbox_min,
            group_bbox_max: wg.bbox_max,
            excl_icol: ex2.excl,
            excl_count: ex2.excl_count,
            excl_max: new Int32Array([ex2.EXCL_MAX | 0]),
            group_cap: new Int32Array([wg.groupCap | 0]),
            n_groups: new Int32Array([wg.nGroups | 0]),
        };
    }
    const mergedExtra = collisionExtra ? { ...collisionExtra, ...(extra || {}) } : extra;
    if (outNpzPath) writeTopologyNpzFile(fs, outNpzPath, topo, packing, mol, meta, mergedExtra);
    return { natoms: nAtoms, meta, packing, topo, mol, extra: mergedExtra };
}

/// Full topology NPZ export: angles + dihedrals + surface group AABBs (matches build_linearized_topology.mjs defaults).
export function buildTopologyNpzFull({ mol, mmParams, outNpzPath, maxNeighbors = 48, addAngle = true, addDihedral = true, K12 = 0, K13 = 0, K14 = 0, groupCap = 32, exportSurfaceGroups = true, source = '', extra = null }) {
    return buildTopologyNpzFromMol({ mol, mmParams, outNpzPath, maxNeighbors, addAngle, addDihedral, K12, K13, K14, groupCap, exportSurfaceGroups, source, extra });
}

export function buildTopologyNpz({ mol2Path, relaxedXyzPath, outNpzPath, mm = null, maxNeighbors = 48, addAngle = false, addDihedral = false, groupCap = 32 }) {
    const mmParams = mm || loadMMParamsFromDir();
    const mol = loadMolFromMol2(mol2Path, mmParams);
    if (relaxedXyzPath) {
        const { pos } = readXyzPositions(fs, relaxedXyzPath);
        applyPositions(mol, pos);
    }
    return buildTopologyNpzFromMol({ mol, mmParams, outNpzPath, maxNeighbors, addAngle, addDihedral, K12: 0, K13: 0, K14: 0, groupCap, exportSurfaceGroups: true, source: mol2Path });
}

/// Build primary adjacency list from mol bonds with bond params (l0, K).
export function buildPrimaryAdj(mol, mm, typeNames) {
    const n = mol.atoms.length;
    const adj = new Array(n);
    for (let i = 0; i < n; i++) adj[i] = [];
    for (let i = 0; i < mol.bonds.length; i++) {
        const b = mol.bonds[i];
        b.ensureIndices(mol);
        const a = b.a | 0, c = b.b | 0;
        const bp = mm.getBondParams(typeNames[a], typeNames[c]);
        const l0 = +bp.l0, K = +bp.k;
        adj[a].push([c, l0, K]);
        adj[c].push([a, l0, K]);
    }
    return adj;
}

/// Validate that spring graph distances match expected (1-2=1, 1-3=2, 1-4=3).
export function validateSpringGraphDistances(topo, bondsAdj1) {
    const gdist = enumerateGraphDistances(bondsAdj1, 6);
    const pairDist = (a, b) => {
        const n = bondsAdj1.length;
        const dist = new Int32Array(n);
        dist.fill(-1);
        const q = [a | 0];
        dist[a] = 0;
        let qi = 0;
        while (qi < q.length) {
            const u = q[qi++];
            const du = dist[u] | 0;
            for (const nb of bondsAdj1[u] || []) {
                const v = nb[0] | 0;
                if (dist[v] >= 0) continue;
                dist[v] = du + 1;
                q.push(v);
            }
        }
        return dist[b | 0] | 0;
    };
    const errors = [];
    for (let k = 0; k < topo.bonds_primary.length; k++) {
        const e = topo.bonds_primary[k];
        const d = pairDist(e[0], e[1]);
        if (d !== 1) errors.push(`bond (${e[0]},${e[1]}) graph_dist=${d} expected 1`);
    }
    for (const e of topo.bonds_angle) {
        const d = pairDist(e[0], e[1]);
        if (d !== 2) errors.push(`angle spring (${e[0]},${e[1]}) graph_dist=${d} expected 2`);
    }
    for (const e of topo.bonds_dihedral) {
        const d = pairDist(e[0], e[1]);
        if (d !== 3) errors.push(`dihedral spring (${e[0]},${e[1]}) graph_dist=${d} expected 3 (MMFF 1-4)`);
    }
    return { gdist, errors };
}

/// Build full topology with validation, viewer export, and NPZ output (unified from build_linearized_topology.mjs).
export async function buildTopologyFull({ inputPath, outDir, mm = null, addAngle = true, addDihedral = true, K12 = 0, K13 = 0, K14 = 0, maxNeighbors = 48, graphDist = true, exportNpz = true, exportJson = false, htmlViewer = true, npzViewer = true }) {
    const mmParams = mm || loadMMParamsFromDir();
    const mol = loadMol(inputPath, mmParams);
    const topoOpts = { type_source: 'table', add_angle: addAngle, add_dihedral: addDihedral, add_pi: false, add_epair: false, K12, K13, K14 };
    const topo = buildMMFFLTopology(mol, mmParams, topoOpts);
    const nBond = topo.bonds_primary.length | 0;
    const nAng = topo.bonds_angle.length | 0;
    const nDih = (topo.bonds_dihedral ? topo.bonds_dihedral.length : 0) | 0;
    const typeNames = [];
    for (let i = 0; i < topo.n_real; i++) typeNames[i] = mmParams.resolveTypeNameTable(Z_TO_SYMBOL[mol.atoms[i].Z] || 'X');
    const bondsAdj1 = buildPrimaryAdj(mol, mmParams, typeNames);
    let graphReport = null;
    if (graphDist) {
        const { gdist, errors } = validateSpringGraphDistances(topo, bondsAdj1);
        graphReport = { hist_undirected_pairs: gdist.hist, per_atom_max_dist: Array.from(gdist.perAtomMax), validation_errors: errors };
        if (errors.length > 0) throw new Error(`graph distance validation failed (${errors.length} errors)`);
    }
    const packing = packLinearTopologyForGPU(topo, { maxNeighbors });
    const base = path.basename(inputPath, path.extname(inputPath));
    const meta = { source: inputPath, n_bond: nBond, n_angle: nAng, n_dihedral: nDih };
    const results = { topo, packing, mol, meta, graphReport, base, nBond, nAng, nDih };
    if (exportNpz) {
        const npzPath = path.join(outDir, `${base}_topology.npz`);
        writeTopologyNpzFile(fs, npzPath, topo, packing, mol, meta);
        results.npzPath = npzPath;
    }
    if (exportJson) {
        const springs = exportTopologySpringList(topo);
        const topoJson = { meta: { ...meta, natoms: topo.n_real | 0, K12_override: K12, K13_override: K13, K14_override: K14, max_neighbors_required: packing.maxDegree, max_neighbors_budget: maxNeighbors, graph_dist: graphReport }, springs };
        const jsonPath = path.join(outDir, `${base}_topology.json`);
        fs.writeFileSync(jsonPath, JSON.stringify(topoJson, null, 2));
        results.jsonPath = jsonPath;
    }
    if (htmlViewer) {
        const html = exportStickViewerHTML(topo, `${base} MMFFL`);
        const htmlPath = path.join(outDir, `${base}_debug_viewer.html`);
        fs.writeFileSync(htmlPath, html);
        results.htmlPath = htmlPath;
    }
    if (npzViewer) {
        const npzHtml = exportStickViewerNpzLoaderHTML(`${base} MMFFL`, base);
        const npzHtmlPath = path.join(outDir, `${base}_npz_viewer.html`);
        fs.writeFileSync(npzHtmlPath, npzHtml);
        results.npzHtmlPath = npzHtmlPath;
    }
    return results;
}
