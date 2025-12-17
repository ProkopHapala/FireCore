import { MoleculeSystem } from '../js/MoleculeSystem.js';

function getArg(name, defVal = null) {
    const argv = process.argv;
    const i = argv.indexOf(name);
    if (i < 0) return defVal;
    const v = argv[i + 1];
    if (v === undefined) return defVal;
    return v;
}

function hasFlag(name) { return process.argv.indexOf(name) >= 0; }

function parseMap(str) {
    // format: D=path:head,tail;T=path:head,tail
    const out = {};
    if (!str) return out;
    const items = str.split(';').map(s => s.trim()).filter(Boolean);
    for (const it of items) {
        const [key, rest] = it.split('=');
        if (!key || !rest) throw new Error(`Bad map item '${it}'`);
        const [pathPart, anchPart] = rest.split(':');
        if (!pathPart || !anchPart) throw new Error(`Bad map item '${it}' (need path:head,tail)`);
        const [h, t] = anchPart.split(',').map(x => parseInt(x));
        if (!(h > 0 && t > 0)) throw new Error(`Bad anchors in '${it}'`);
        out[key] = { path: pathPart, anchors: [h, t] };
    }
    return out;
}

function usage() {
    console.log('Usage: node polymerize_sequence.mjs --seq DDDD_DDDD --map "D=../../tests/tAttach/backbones/DANA_CG_2inv.mol2:1,4;P=../../tests/tAttach/PNA_CG.mol2:2,19" --out out.xyz');
    console.log('Notes: anchors are 1-based atom indices [head,tail] like in tests/tAttach/polymerize.py');
}

if (hasFlag('--help') || hasFlag('-h')) { usage(); process.exit(0); }

const seq = getArg('--seq', null);
const mapStr = getArg('--map', null);
const outPath = getArg('--out', 'sequence.xyz');
const keepLvs = hasFlag('--lvs');

if (!seq || !mapStr) {
    usage();
    throw new Error('Missing --seq or --map');
}

const fs = await import('node:fs');
const mp = parseMap(mapStr);

const monomers = {};
let lvec0 = null;
for (const key of Object.keys(mp)) {
    const rec = mp[key];
    const txt = fs.readFileSync(rec.path, 'utf8');
    const parsed = MoleculeSystem.parseMol2(txt);
    monomers[key] = { parsed, anchors: rec.anchors };
    if (!lvec0) lvec0 = parsed.lvec;
}

const sys = MoleculeSystem.assemblePolymerSequence(seq, monomers, { _0: 1, capacity: 200000 });
const xyz = sys.toXYZString({ lvec: keepLvs ? lvec0 : null });
fs.writeFileSync(outPath, xyz, 'utf8');
console.log(`Wrote ${outPath} atoms=${sys.nAtoms} bonds=${sys.bonds.length}`);
