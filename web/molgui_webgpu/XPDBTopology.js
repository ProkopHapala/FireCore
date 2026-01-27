import { MMParams } from "./MMParams.js";

/// DEPRECATED: GUI-only lightweight topology builder. Do not use for parity.
/// Use buildMMFFLTopology/buildXPDBInputsFromMol (MMFFLTopology.js) or buildXPDBInputsFromXYZArgs (dump_xpdb_topology.mjs) instead.
/// This file is kept temporarily as a fallback while the GUI migrates to the shared path.
/// Helper to build bondsAdj format for XPDB_WebGPU from EditableMolecule + MMParams
/// bondsAdj[i] = [[neighbor_idx, rest_length, stiffness], ...]
export function buildXPDBTopology(mol, mmParams, opts = {}) {
    if (!mol || !mmParams) throw new Error('buildXPDBTopology: mol and mmParams required');

    const nAtoms = mol.nAtoms || mol.atoms?.length || 0;
    if (nAtoms === 0) throw new Error('buildXPDBTopology: no atoms in molecule');

    const includeAngleConstraints = (opts.includeAngleConstraints !== undefined) ? !!opts.includeAngleConstraints : true;
    const maxBonds = opts.maxBonds || 16;
    const defaultL = opts.defaultL || 1.5;
    const defaultK = opts.defaultK || 100.0;

    console.log(`[XPDBTopology] Building topology: atoms=${nAtoms} maxBonds=${maxBonds} includeAngles=${includeAngleConstraints}`);

    // Initialize adjacency lists: first-order (real) and second-order (angles), plus merged output
    const bondsAdj1 = new Array(nAtoms); // real bonds only
    const bondsAdj2 = new Array(nAtoms); // angle-derived virtual bonds
    for (let i = 0; i < nAtoms; i++) { bondsAdj1[i] = []; bondsAdj2[i] = []; }

    // Helper to get atom type name
    const getAtomTypeName = (atom) => {
        if (atom.typeName) return atom.typeName;
        const at = mmParams.getAtomTypeForAtom(atom);
        return at ? at.name : '*';
    };

    // Helper to add bond entry into a target adjacency list (mutates target)
    const addBondEntry = (adj, i, j, l0, k, label='') => {
        if (i < 0 || i >= nAtoms || j < 0 || j >= nAtoms) {
            throw new Error(`addBondEntry: out of range i=${i} j=${j} nAtoms=${nAtoms}`);
        }
        if (adj[i].length === (maxBonds - 1)) {
            console.log(`[XPDBTopology][DEBUG] addBondEntry about to reach maxBonds: atom=${i} nextNeighbor=${j} maxBonds=${maxBonds} ${label}`);
        }
        if (adj[i].length >= maxBonds) {
            console.warn(`[XPDBTopology][WARN] addBondEntry: atom ${i} exceeds max bonds (${maxBonds}), skipping neighbor ${j} ${label}`);
            return;
        }
        adj[i].push([j, l0, k]);
    };

    // Helper to find existing bond entry in provided adjacency
    const findBondEntry = (adj, i, j) => {
        const arr = adj[i];
        for (let k = 0; k < arr.length; k++) {
            if (arr[k][0] === j) return k;
        }
        return -1;
    };

    // 1. Process real bonds
    const bonds = mol.bonds || [];
    let realConstraintsCount = 0;

    console.log(`[XPDBTopology] Processing ${bonds.length} real bonds...`);
    for (const bond of bonds) {
        if (typeof bond.ensureIndices === 'function') {
            bond.ensureIndices(mol);
        }
        const i = (Number.isInteger(bond.a)) ? bond.a : bond.i;
        const j = (Number.isInteger(bond.b)) ? bond.b : bond.j;
        if (!Number.isInteger(i) || !Number.isInteger(j)) {
            throw new Error(`buildXPDBTopology: bond ${bond.id ?? 'unknown'} missing atom indices`);
        }
        if (i < 0 || i >= nAtoms || j < 0 || j >= nAtoms) {
            throw new Error(`buildXPDBTopology: bond ${bond.id ?? 'unknown'} out of range i=${i} j=${j}`);
        }

        const atomA = mol.atoms[i];
        const atomB = mol.atoms[j];
        if (!atomA || !atomB) continue;

        const zA = atomA.Z || 0;
        const zB = atomB.Z || 0;
        const bondData = mmParams.getBondL0(zA, zB);
        if (!bondData) {
            const enA = mmParams.getElementName ? mmParams.getElementName(zA) : zA;
            const enB = mmParams.getElementName ? mmParams.getElementName(zB) : zB;
            throw new Error(`buildXPDBTopology: missing bond params for pair ${i}(${enA})-${j}(${enB}) z=(${zA},${zB})`);
        }

        let l0 = bondData ? bondData.l0 : defaultL;
        let k = bondData ? bondData.k : defaultK;

        addBondEntry(bondsAdj1, i, j, l0, k, '[real]');
        addBondEntry(bondsAdj1, j, i, l0, k, '[real]');
        realConstraintsCount++;
    }

    console.log(`[XPDBTopology] Real bonds processed: ${realConstraintsCount} constraints`);

    // 2. Process angles to create virtual bonds (kept separate from real bonds)
    let virtualConstraintsCount = 0;
    let angleTriplesProcessed = 0;
    let angleParamsFound = 0;
    let angleParamsMissing = 0;
    let angleBondsAdded = 0;
    let angleBondsAccumulated = 0;

    if (includeAngleConstraints) {
        console.log(`[XPDBTopology] Processing angle-derived constraints...`);

        for (let b = 0; b < nAtoms; b++) {
            const neighs = bondsAdj1[b];
            if (!neighs || neighs.length < 2) continue;

            const atomB = mol.atoms[b];
            const nameB = getAtomTypeName(atomB);

            // Enumerate all pairs of neighbors to form angles A-B-C
            for (let ni = 0; ni < neighs.length; ni++) {
                for (let nj = ni + 1; nj < neighs.length; nj++) {
                    let a = neighs[ni][0];
                    let c = neighs[nj][0];
                    if (a < 0 || a >= nAtoms || c < 0 || c >= nAtoms) continue;

                    // Canonical ordering to avoid duplicate A-C virtual bonds from different triples.
                    // (We still add both directions a->c and c->a below.)
                    if (c < a) { const t = a; a = c; c = t; }

                    const atomA = mol.atoms[a];
                    const atomC = mol.atoms[c];
                    if (!atomA || !atomC) continue;

                    const nameA = getAtomTypeName(atomA);
                    const nameC = getAtomTypeName(atomC);

                    angleTriplesProcessed++;

                    const angleParams = mmParams.getAngleParams(nameA, nameB, nameC);
                    if (!angleParams) {
                        const enA = mmParams.getElementName ? mmParams.getElementName(atomA.Z || 0) : atomA.Z;
                        const enB = mmParams.getElementName ? mmParams.getElementName(atomB.Z || 0) : atomB.Z;
                        const enC = mmParams.getElementName ? mmParams.getElementName(atomC.Z || 0) : atomC.Z;
                        throw new Error(`buildXPDBTopology: missing angle params for ${enA}-${enB}-${enC} names=(${nameA},${nameB},${nameC}) at center ${b}`);
                    }

                    angleParamsFound++;

                    // Get bond lengths for AB and BC
                    let lab = defaultL, lbc = defaultL;
                    const idxAB = findBondEntry(bondsAdj1, a, b);
                    const idxBC = findBondEntry(bondsAdj1, b, c);
                    if (idxAB >= 0) lab = bondsAdj1[a][idxAB][1];
                    if (idxBC >= 0) lbc = bondsAdj1[b][idxBC][1];

                    // Convert angle to distance constraint
                    const virtBond = mmParams.convertAngleToDistance(
                        lab, lbc, angleParams.ang0, angleParams.k
                    );

                    // Add virtual bond A-C (accumulate if exists)
                    const idxAC = findBondEntry(bondsAdj2, a, c);
                    if (idxAC >= 0) {
                        bondsAdj2[a][idxAC][2] += virtBond.stiffness;
                        angleBondsAccumulated++;
                        const idxCA = findBondEntry(bondsAdj2, c, a);
                        if (idxCA >= 0) bondsAdj2[c][idxCA][2] += virtBond.stiffness;
                    } else {
                        addBondEntry(bondsAdj2, a, c, virtBond.restLength, virtBond.stiffness, '[angle]');
                        addBondEntry(bondsAdj2, c, a, virtBond.restLength, virtBond.stiffness, '[angle]');
                        angleBondsAdded++;
                        virtualConstraintsCount++;
                    }
                }
            }
        }

        console.log(`[XPDBTopology] Angle processing: triples=${angleTriplesProcessed} found=${angleParamsFound} missing=${angleParamsMissing} added=${angleBondsAdded} accumulated=${angleBondsAccumulated} virtual=${virtualConstraintsCount}`);
    } else {
        console.log(`[XPDBTopology] Angle-derived constraints SKIPPED`);
    }

    // 3. Merge first- and second-order neighbors respecting maxBonds, keeping real bonds first
    const bondsAdjMerged = bondsAdj1.map((list, i) => list.slice());
    let mergedSkipped = 0;
    for (let i = 0; i < nAtoms; i++) {
        const targ = bondsAdjMerged[i];
        const extras = bondsAdj2[i];
        for (const b of extras) {
            if (targ.length >= maxBonds) { mergedSkipped++; continue; }
            targ.push(b);
        }
    }

    const totalConstraints = realConstraintsCount + virtualConstraintsCount;
    console.log(`[XPDBTopology] Summary: atoms=${nAtoms} real=${realConstraintsCount} virtual=${virtualConstraintsCount} total=${totalConstraints} mergedSkipped=${mergedSkipped}`);

    return {
        bondsAdj: bondsAdjMerged,
        stats: {
            atoms: nAtoms,
            realConstraints: realConstraintsCount,
            virtualConstraints: virtualConstraintsCount,
            mergedSkipped,
            totalConstraints
        }
    };
}

/// Helper to compute atom radii from MMParams
export function buildAtomRadii(mol, mmParams) {
    const nAtoms = mol.nAtoms || mol.atoms?.length || 0;
    const radii = new Float32Array(nAtoms);

    for (let i = 0; i < nAtoms; i++) {
        const atom = mol.atoms[i];
        if (!atom) continue;

        const at = mmParams.getAtomTypeForAtom(atom);
        radii[i] = (at && at.RvdW > 0) ? at.RvdW : 1.0;
    }

    return radii;
}

/// Helper to compute atom masses (approximate from atomic number)
export function buildAtomMasses(mol) {
    const nAtoms = mol.nAtoms || mol.atoms?.length || 0;
    const masses = new Float32Array(nAtoms);

    // Simple approximation: mass ≈ atomic number (not accurate but sufficient for PD)
    for (let i = 0; i < nAtoms; i++) {
        const atom = mol.atoms[i];
        if (!atom) continue;
        masses[i] = (atom.Z || 1) * 1.0;
    }

    return masses;
}

/// Helper to compute max radius from molecule
export function getMaxRadius(mol, mmParams, fallback = 1.0) {
    let rMax = fallback;
    for (let i = 0; i < mol.atoms.length; i++) {
        const atom = mol.atoms[i];
        if (!atom) continue;
        const at = mmParams.getAtomTypeForAtom(atom);
        if (at && at.RvdW > rMax) rMax = at.RvdW;
    }
    return rMax;
}
