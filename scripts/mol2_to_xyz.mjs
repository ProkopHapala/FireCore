#!/usr/bin/env node
/** Convert .mol2 → .xyz using EditableMolecule (headless). */
import fs from 'node:fs';

import { EditableMolecule } from '../web/molgui_webgpu/EditableMolecule.js';
import { toXYZString, installMoleculeIOMethods } from '../web/molgui_webgpu/MoleculeIO.js';

installMoleculeIOMethods(EditableMolecule);

const inp = process.argv[2];
const outp = process.argv[3];
if (!inp || !outp) {
    console.error('usage: node scripts/mol2_to_xyz.mjs <in.mol2> <out.xyz>');
    process.exit(1);
}
const parsed = EditableMolecule.parseMol2(fs.readFileSync(inp, 'utf8'));
const mol = new EditableMolecule();
mol.appendParsedSystem(parsed);
fs.writeFileSync(outp, toXYZString(mol, {}));
console.log(`[mol2_to_xyz] ${inp} -> ${outp} atoms=${mol.atoms.length}`);
