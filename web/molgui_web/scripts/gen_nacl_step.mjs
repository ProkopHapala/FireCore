import { MoleculeSystem } from '../js/MoleculeSystem.js';

// use like this: node gen_nacl_step.mjs --out NaCl_cubic.xyz --nx 13 --ny 12 --nz 3 --a 2.82065 --Q0 0.7

function getArg(name, defVal = null) {
    const argv = process.argv;
    const i = argv.indexOf(name);
    if (i < 0) return defVal;
    const v = argv[i + 1];
    if (v === undefined) return defVal;
    return v;
}

function hasFlag(name) {
    return process.argv.indexOf(name) >= 0;
}

function fmt(x) {
    return (Number.isFinite(x) ? x : 0).toFixed(5);
}

function genXYZ({ pos, types, qs, nTot, lvec }) {
    const a = [];
    a.push(String(nTot));
    a.push(`lvs ${lvec[0][0]} ${lvec[0][1]} ${lvec[0][2]}   ${lvec[1][0]} ${lvec[1][1]} ${lvec[1][2]}   ${lvec[2][0]} ${lvec[2][1]} ${lvec[2][2]}`);
    for (let i = 0; i < nTot; i++) {
        const i3 = i * 3;
        const t = types[i];
        const S = (t === 11) ? 'Na' : (t === 17) ? 'Cl' : 'X';
        const q = qs ? qs[i] : 0.0;
        a.push(`${S} ${fmt(pos[i3])} ${fmt(pos[i3 + 1])} ${fmt(pos[i3 + 2])} ${q}`);
    }
    return a.join('\n') + '\n';
}

const nx = parseInt(getArg('--nx', '13'));
const ny = parseInt(getArg('--ny', '12'));
const nz = parseInt(getArg('--nz', '3'));
const a0 = parseFloat(getArg('--a', String(5.6413 / 2)));
const Q0 = parseFloat(getArg('--Q0', '0.7'));
const out = getArg('--out', 'NaCl_cubic.xyz');
const toStdout = hasFlag('--stdout');

const res = MoleculeSystem.genNaClStep({ a: a0, nx, ny, nz, Q0 });
const xyz = genXYZ(res);

if (toStdout) {
    process.stdout.write(xyz);
} else {
    const fs = await import('node:fs');
    fs.writeFileSync(out, xyz, 'utf8');
    console.log(`Wrote ${out} (${res.nTot} atoms)`);
}
