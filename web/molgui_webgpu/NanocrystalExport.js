/// Unified nanocrystal NPZ export: crystal geometry + full topology + surface group AABBs.
import fs from 'node:fs';
import path from 'node:path';
import { molToCrystalArrays, writeCrystalNpz } from '../common_js/npzIO.js';
import { buildTopologyNpzFull } from '../common_js/exportFF.js';

/// exportNanocrystalBundle(mol, mmParams, opts) → { crystalNpz, topologyNpz, meta, metaPath }
export function exportNanocrystalBundle(mol, mmParams, opts = {}) {
    const {
        outDir,
        id = 'crystal',
        genParams = {},
        exportTopology = true,
        exportSurfaceGroups = true,
        addAngle = true,
        addDihedral = true,
        K12 = 0, K13 = 0, K14 = 0,
        maxNeighbors = 48,
        groupCap = 32,
        defectsMeta = {},
        writeMeta = true,
        crystalName = '01_crystal.npz',
        topologyName = '03_topology.npz',
    } = opts;
    if (!outDir) throw new Error('exportNanocrystalBundle: outDir required');
    if (!mol || !mol.atoms || mol.atoms.length === 0) throw new Error('exportNanocrystalBundle: empty mol');
    if (!mmParams) throw new Error('exportNanocrystalBundle: mmParams required');
    fs.mkdirSync(outDir, { recursive: true });

    const genJson = JSON.stringify(genParams);
    const arrays = molToCrystalArrays(mol);
    const crystalNpz = path.join(outDir, crystalName);
    writeCrystalNpz(fs, crystalNpz, { ...arrays, gen_params: genJson });

    let topologyNpz = null;
    let topoMeta = null;
    if (exportTopology) {
        topologyNpz = path.join(outDir, topologyName);
        const defectsJson = {
            nCollapsed: defectsMeta.nCollapsed ?? 0,
            nInserted: defectsMeta.nInserted ?? 0,
            nFused: defectsMeta.nFused ?? 0,
            insertProb: defectsMeta.insertProb ?? genParams.insertProb ?? 0,
            collapseProb: defectsMeta.collapseProb ?? genParams.collapseProb ?? 0,
            fuseProb: defectsMeta.fuseProb ?? genParams.fuseProb ?? 0,
        };
        const defectsBytes = new Uint8Array(new TextEncoder().encode(JSON.stringify(defectsJson)));
        defectsBytes._shape = [defectsBytes.length];
        const schema_version = new Int32Array([1]);
        const { meta } = buildTopologyNpzFull({
            mol, mmParams, outNpzPath: topologyNpz,
            addAngle, addDihedral, K12, K13, K14,
            maxNeighbors, groupCap, exportSurfaceGroups,
            source: id,
            extra: { schema_version, defects_json: defectsBytes },
        });
        topoMeta = meta;
    }

    const meta = {
        id,
        natoms: arrays.natoms,
        nbonds: arrays.nbonds,
        gen_params: genParams,
        defects: {
            nCollapsed: defectsMeta.nCollapsed ?? 0,
            nInserted: defectsMeta.nInserted ?? 0,
            nFused: defectsMeta.nFused ?? 0,
            insertProb: defectsMeta.insertProb ?? genParams.insertProb ?? 0,
            collapseProb: defectsMeta.collapseProb ?? genParams.collapseProb ?? 0,
            fuseProb: defectsMeta.fuseProb ?? genParams.fuseProb ?? 0,
        },
        crystal_npz: crystalName,
        topology_npz: topologyNpz ? topologyName : null,
        topology: topoMeta,
    };
    const metaPath = path.join(outDir, 'meta.json');
    if (writeMeta) fs.writeFileSync(metaPath, JSON.stringify(meta, null, 2));
    return { crystalNpz, topologyNpz, meta, metaPath };
}
