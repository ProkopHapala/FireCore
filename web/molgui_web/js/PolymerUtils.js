import { EditableMolecule } from "./EditableMolecule.js";
import { Vec3 } from "../../common_js/Vec3.js";


export function assemblePolymerFromTokens(tokens, monomers, opts = {}) {
    const _0 = (opts._0 !== undefined) ? (opts._0 | 0) : 1;
    if (!tokens || tokens.length === 0) throw new Error('assemblePolymerFromTokens: tokens is empty');

    const mol = new EditableMolecule();

    const pos = new Vec3(0, 0, 0);
    let prev = null;
    let prevTailId = -1;

    for (let i = 0; i < tokens.length; i++) {
        const key = tokens[i];
        if (!key) continue;
        const rec = monomers[key];
        if (!rec || !rec.parsed) throw new Error(`assemblePolymerFromTokens: missing monomer '${key}' or its parsed data`);
        if (!rec.anchors || rec.anchors.length < 2) throw new Error(`assemblePolymerFromTokens: monomer '${key}' missing anchors [head,tail] (1-based)`);

        if (!prev) {
            const ids = mol.appendParsedSystem(rec.parsed, { pos });
            const iTail = ((rec.anchors[1] | 0) - _0);
            if (iTail < 0 || iTail >= ids.length) throw new Error(`assemblePolymerFromTokens: tail anchor out of range for '${key}' iTail=${iTail}`);
            prevTailId = ids[iTail];
            prev = rec;
        } else {
            const lv = (rec.parsed.lvec && rec.parsed.lvec[1]) ? rec.parsed.lvec[1] : null;
            if (lv) {
                if (!(lv instanceof Vec3)) throw new Error('assemblePolymerFromTokens: parsed.lvec must be Vec3[3]');
                pos.add(lv);
            }
            const ids = mol.appendParsedSystem(rec.parsed, { pos });
            const iHead = ((rec.anchors[0] | 0) - _0);
            const iTail = ((rec.anchors[1] | 0) - _0);
            if (iHead < 0 || iHead >= ids.length) throw new Error(`assemblePolymerFromTokens: head anchor out of range for '${key}' iHead=${iHead}`);
            if (iTail < 0 || iTail >= ids.length) throw new Error(`assemblePolymerFromTokens: tail anchor out of range for '${key}' iTail=${iTail}`);
            const headId = ids[iHead];
            mol.addBond(prevTailId, headId);
            prevTailId = ids[iTail];
            prev = rec;
        }
    }

    return mol;
}


