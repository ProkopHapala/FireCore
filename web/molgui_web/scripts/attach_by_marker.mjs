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

function usage() {
    console.log('Usage: node attach_by_marker.mjs --backbone <backbone.mol2> --group <group.mol2> --markerX Se --markerY Cl --out out.xyz');
    console.log('Example: node attach_by_marker.mjs --backbone ../../tests/tAttach/porphironoids/PorhQuad_4SeCl_30deg.mol2 --group ../../tests/tAttach/endgroups/guanine-SeCl.mol2 --markerX Se --markerY Cl --out PorhQuad_G_30deg_from_js.xyz');
}

if (hasFlag('--help') || hasFlag('-h')) { usage(); process.exit(0); }

const backbonePath = getArg('--backbone', null);
const groupPath = getArg('--group', null);
const markerX = getArg('--markerX', 'Xe');
const markerY = getArg('--markerY', 'He');
const outPath = getArg('--out', 'attached.xyz');
const keepMol2Cell = hasFlag('--lvs');

if (!backbonePath || !groupPath) {
    usage();
    throw new Error('Missing --backbone or --group');
}

const fs = await import('node:fs');

const backboneText = fs.readFileSync(backbonePath, 'utf8');
const groupText = fs.readFileSync(groupPath, 'utf8');

const backboneParsed = MoleculeSystem.parseMol2(backboneText);
const groupParsed = MoleculeSystem.parseMol2(groupText);

const sys = new MoleculeSystem(Math.max(1024, backboneParsed.types.length * 2));
sys.appendParsedSystem(backboneParsed);

sys.attachGroupByMarker(groupParsed, markerX, markerY);

const xyz = sys.toXYZString({ lvec: keepMol2Cell ? backboneParsed.lvec : null });
fs.writeFileSync(outPath, xyz, 'utf8');
console.log(`Wrote ${outPath}  atoms=${sys.nAtoms}  bonds=${sys.bonds.length}`);
