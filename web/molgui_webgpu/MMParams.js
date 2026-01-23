/// Simple element record with covalent/van der Waals radii and QEq props.
class ElementType {
    constructor(name) {
        this.name = name;
        this.iZ = 0;
        this.neval = 0;
        this.valence = 0;
        this.piMax = 0;
        this.color = 0xFFFFFF;
        this.Rcov = 0.0;
        this.RvdW = 0.0;
        this.EvdW = 0.0;
        this.Quff = 0.0;
        this.Uuff = 0.0;
        this.Vuff = 0.0;
        // QEq
        this.bQEq = false;
        this.Eaff = 0.0;
        this.Ehard = 0.0;
        this.Ra = 0.0;
        this.eta = 0.0;
    }
}


/// Atom-type record with valence, radii, MMFF params, and links to element/parent/epair.
class AtomType {
    constructor(name) {
        this.name = name;
        this.parent_name = "*";
        this.element_name = "";
        this.epair_name = "";
        this.valence = 0;
        this.nepair = 0;
        this.npi = 0;
        this.sym = 0;
        this.Ruff = 0.0;
        this.RvdW = 0.0;
        this.EvdW = 0.0;
        this.Qbase = 0.0;
        this.Hb = 0.0;
        // MMFF
        this.bMMFF = false;
        this.Ass = 0.0;
        this.Asp = 0.0;
        this.Kss = 0.0;
        this.Ksp = 0.0;
        this.Kep = 0.0;
        this.Kpp = 0.0;

        // References
        this.element = null; // ElementType object
        this.parent = null;  // AtomType object
        this.epair = null;   // AtomType object

        // Cached from element for convenience
        this.iZ = 0;
        this.color = 0xFFFFFF;
    }
}

/// Container for element/atom/bond/angle type tables with parsing helpers.
export class MMParams {
    constructor() {
        this.elementTypes = {}; // name -> ElementType
        this.atomTypes = {};    // name -> AtomType
        this.byAtomicNumber = {}; // iZ -> ElementType
        this.bondTypes = {};    // key -> { l0, k, order, a, b }
        this.angleTypes = [];   // array of { a1, a2, a3, ang0, k } with wildcard support
    }

    /// Build cached name↔index maps on atomTypes if missing.
    ensureTypeIndex() {
        if (!this.atomTypes) throw new Error('ensureTypeIndex: atomTypes required');
        if (this._atomTypeNameToIndex && this._atomTypeIndexToName) return;
        const names = Object.keys(this.atomTypes).sort();
        const n2i = {};
        const i2n = names.slice();
        for (let i = 0; i < i2n.length; i++) n2i[i2n[i]] = i;
        this._atomTypeNameToIndex = n2i;
        this._atomTypeIndexToName = i2n;
    }

    /// Resolve atom-type name to index (or -1 if missing).
    atomTypeNameToIndex(name) {
        this.ensureTypeIndex();
        const i = this._atomTypeNameToIndex[name];
        return (i === undefined) ? -1 : (i | 0);
    }

    /// Resolve atom-type index to name (or null).
    atomTypeIndexToName(i) {
        this.ensureTypeIndex();
        const n = this._atomTypeIndexToName[i | 0];
        return n ? n : null;
    }

    /// Normalize type/element token to atom-type record with name, index, and iZ.
    resolveTypeOrElementToAtomType(nameOrEl) {
        const s = String(nameOrEl || '').trim();
        if (!s) throw new Error('resolveTypeOrElementToAtomType: empty type');
        if (this.atomTypes && this.atomTypes[s]) {
            const at = this.atomTypes[s];
            return { name: s, atype: this.atomTypeNameToIndex(s), iZ: at.iZ | 0 };
        }
        if (this.elementTypes && this.elementTypes[s]) {
            const el = this.elementTypes[s];
            const iZ = el.iZ | 0;
            if (this.atomTypes && this.atomTypes[s]) {
                const at = this.atomTypes[s];
                return { name: s, atype: this.atomTypeNameToIndex(s), iZ: at.iZ | 0 };
            }
            if (this.atomTypes) {
                for (const k in this.atomTypes) {
                    const at = this.atomTypes[k];
                    if (at && at.element_name === s) return { name: k, atype: this.atomTypeNameToIndex(k), iZ };
                }
            }
            return { name: s, atype: -1, iZ };
        }
        throw new Error(`Unknown cap type/element '${s}'`);
    }

    /// Resolve atom-type object for a given atom (atype or element fallback).
    getAtomTypeForAtom(atom) {
        if (!atom) throw new Error('getAtomTypeForAtom: atom required');
        if (atom.atype !== undefined && (atom.atype | 0) >= 0) {
            const name = this.atomTypeIndexToName(atom.atype | 0);
            if (name && this.atomTypes && this.atomTypes[name]) return this.atomTypes[name];
        }
        if (this.byAtomicNumber) {
            const et = this.byAtomicNumber[atom.Z | 0];
            if (et && this.atomTypes && this.atomTypes[et.name]) return this.atomTypes[et.name];
        }
        if (this.atomTypes && this.atomTypes['*']) return this.atomTypes['*'];
        throw new Error(`getAtomTypeForAtom: cannot resolve atom type for Z=${atom.Z}`);
    }

    /// Compute squared bond cutoff from covalent radii (fallback to default) for a Z pair.
    bondCutoff2(z1, z2, defaultRcut2, bondFactor, stats) {
        if (stats) stats.nBondCut2 = (stats.nBondCut2 | 0) + 1;
        let r1 = 0.7;
        let r2 = 0.7;
        if (this.byAtomicNumber) {
            const e1 = this.byAtomicNumber[z1 | 0];
            const e2 = this.byAtomicNumber[z2 | 0];
            if (e1 && (e1.Rcov > 0)) r1 = e1.Rcov;
            if (e2 && (e2.Rcov > 0)) r2 = e2.Rcov;
            if (e1 && e2) {
                const rSum = (r1 + r2) * bondFactor;
                if (stats) stats.nRcovCutoff = (stats.nRcovCutoff | 0) + 1;
                return rSum * rSum;
            }
        }
        if (stats) stats.nDefaultCutoff = (stats.nDefaultCutoff | 0) + 1;
        return defaultRcut2;
    }

    /// Find maximum covalent radius present in atom list (fallback when unknown).
    maxRcovFromAtoms(atoms, fallback = 0.7) {
        let r = +fallback;
        if (this.byAtomicNumber) {
            for (let i = 0; i < atoms.length; i++) {
                const z = atoms[i].Z | 0;
                const et = this.byAtomicNumber[z];
                if (et && (et.Rcov > r)) r = et.Rcov;
            }
        }
        return r;
    }

    /// Estimate bond length from covalent radii (fallback constant).
    bondLengthEstimate(zA, zB, bondFactor = 1.1, fallback = 1.0) {
        if (this.byAtomicNumber) {
            const e1 = this.byAtomicNumber[zA | 0];
            const e2 = this.byAtomicNumber[zB | 0];
            const r1 = (e1 && (e1.Rcov > 0)) ? e1.Rcov : 0.0;
            const r2 = (e2 && (e2.Rcov > 0)) ? e2.Rcov : 0.0;
            if (r1 > 0 && r2 > 0) return (r1 + r2) * bondFactor;
        }
        return +fallback;
    }

    async loadResources(elementPath, atomPath, bondPath = null, anglePath = null) {
        try {
            const eRes = await fetch(elementPath);
            const eText = await eRes.text();
            this.parseElementTypes(eText);

            const aRes = await fetch(atomPath);
            const aText = await aRes.text();
            this.parseAtomTypes(aText);

            if (bondPath) {
                const bRes = await fetch(bondPath);
                const bText = await bRes.text();
                this.parseBondTypes(bText);
            }

            if (anglePath) {
                const angRes = await fetch(anglePath);
                const angText = await angRes.text();
                this.parseAngleTypes(angText);
            }

            window.logger.info("MMParams loaded successfully.");
        } catch (e) {
            window.logger.error("Failed to load MMParams: " + e);
        }
    }

    parseElementTypes(content) {
        const lines = content.split('\n');
        this.elementTypes = {};
        this.byAtomicNumber = {};

        for (let line of lines) {
            line = line.trim();
            if (!line || line.startsWith('#')) continue;

            const parts = line.split(/\s+/);
            if (parts.length < 12) continue;

            const name = parts[0];
            const et = new ElementType(name);

            try {
                et.iZ = parseInt(parts[1]);
                et.neval = parseInt(parts[2]);
                et.valence = parseInt(parts[3]);
                et.piMax = parseInt(parts[4]);

                let colStr = parts[5];
                if (colStr.startsWith('0x')) {
                    et.color = parseInt(colStr, 16);
                } else {
                    et.color = parseInt(colStr);
                }

                et.Rcov = parseFloat(parts[6]);
                et.RvdW = parseFloat(parts[7]);
                et.EvdW = parseFloat(parts[8]);
                et.Quff = parseFloat(parts[9]);
                et.Uuff = parseFloat(parts[10]);
                et.Vuff = parseFloat(parts[11]);

                if (parts.length >= 16) {
                    et.bQEq = true;
                    et.Eaff = parseFloat(parts[12]);
                    et.Ehard = parseFloat(parts[13]);
                    et.Ra = parseFloat(parts[14]);
                    et.eta = parseFloat(parts[15]);
                }

                this.elementTypes[name] = et;
                this.byAtomicNumber[et.iZ] = et;

            } catch (e) {
                console.warn("Error parsing element line:", line, e);
            }
        }
        window.logger.info(`Parsed ${Object.keys(this.elementTypes).length} element types.`);
        return this.elementTypes;
    }

    parseAtomTypes(content) {
        const lines = content.split('\n');
        this.atomTypes = {};

        // Pass 1: Create objects
        for (let line of lines) {
            line = line.trim();
            if (!line || line.startsWith('#')) continue;
            if (line.startsWith('*')) continue; // Skip header line starting with *

            const parts = line.split(/\s+/);
            if (parts.length < 5) continue;

            const name = parts[0];
            const at = new AtomType(name);
            at.parent_name = parts[1];
            at.element_name = parts[2];
            at.epair_name = parts[3];

            try {
                at.valence = parseInt(parts[4]);
                at.nepair = parseInt(parts[5]);
                at.npi = parseInt(parts[6]);
                at.sym = parseInt(parts[7]);
                at.Ruff = parseFloat(parts[8]);
                at.RvdW = parseFloat(parts[9]);
                at.EvdW = parseFloat(parts[10]);
                at.Qbase = parseFloat(parts[11]);
                at.Hb = parseFloat(parts[12]);

                if (parts.length >= 19) {
                    at.bMMFF = true;
                    at.Ass = parseFloat(parts[13]);
                    at.Asp = parseFloat(parts[14]);
                    at.Kss = parseFloat(parts[15]);
                    at.Ksp = parseFloat(parts[16]);
                    at.Kep = parseFloat(parts[17]);
                    at.Kpp = parseFloat(parts[18]);
                }

                this.atomTypes[name] = at;

            } catch (e) {
                console.warn("Error parsing atom type line:", line, e);
            }
        }

        // Pass 2: Resolve references
        for (const key in this.atomTypes) {
            const at = this.atomTypes[key];

            // Element
            if (this.elementTypes[at.element_name]) {
                at.element = this.elementTypes[at.element_name];
                at.iZ = at.element.iZ;
                at.color = at.element.color;
            }

            // Parent
            if (at.parent_name !== '*' && this.atomTypes[at.parent_name]) {
                at.parent = this.atomTypes[at.parent_name];
            }

            // Epair
            if (at.epair_name !== '*' && this.atomTypes[at.epair_name]) {
                at.epair = this.atomTypes[at.epair_name];
            }
        }

        window.logger.info(`Parsed ${Object.keys(this.atomTypes).length} atom types.`);
        return this.atomTypes;
    }

    parseBondTypes(content) {
        const lines = String(content).replace(/\r/g, '').split('\n');
        this.bondTypes = {};
        let nOk = 0;
        for (let line of lines) {
            line = line.trim();
            if (!line || line.startsWith('#')) continue;
            const parts = line.split(/\s+/);
            if (parts.length < 5) continue;
            const a = parts[0];
            const b = parts[1];
            const order = parseInt(parts[2]) | 0;
            const l0 = parseFloat(parts[3]);
            const k = parseFloat(parts[4]);
            if (!(order > 0) || !(l0 > 0)) continue;

            const za = this.atomTypes[a] ? (this.atomTypes[a].iZ | 0) : 0;
            const zb = this.atomTypes[b] ? (this.atomTypes[b].iZ | 0) : 0;
            if (!(za > 0) || !(zb > 0)) continue;

            const z1 = Math.min(za, zb) | 0;
            const z2 = Math.max(za, zb) | 0;
            const key = `${z1}|${z2}|${order}`;
            const prev = this.bondTypes[key];
            if (!prev || l0 < prev.l0) {
                this.bondTypes[key] = { l0, k, order, a, b, z1, z2 };
            }
            nOk++;
        }
        window.logger.info(`Parsed ${Object.keys(this.bondTypes).length} bond types (from ${nOk} lines).`);
        return this.bondTypes;
    }

    getBondL0(zA, zB) {
        const z1 = Math.min(zA | 0, zB | 0) | 0;
        const z2 = Math.max(zA | 0, zB | 0) | 0;
        if (!(z1 > 0) || !(z2 > 0)) return null;
        let best = null;
        for (const key in this.bondTypes) {
            const bt = this.bondTypes[key];
            if ((bt.z1 | 0) !== z1 || (bt.z2 | 0) !== z2) continue;
            if (!best) best = bt;
            else if ((bt.order | 0) === 1 && (best.order | 0) !== 1) best = bt;
            else if ((bt.order | 0) === (best.order | 0) && bt.l0 < best.l0) best = bt;
            else if ((best.order | 0) !== 1 && bt.l0 < best.l0) best = bt;
        }
        return best ? best : null;
    }

    parseAngleTypes(content) {
        const lines = String(content).replace(/\r/g, '').split('\n');
        this.angleTypes = [];
        for (const line of lines) {
            const trimmed = line.trim();
            if (!trimmed || trimmed.startsWith('#')) continue;
            const parts = trimmed.split(/\s+/);
            if (parts.length < 5) continue;
            const a1 = parts[0];
            const a2 = parts[1];
            const a3 = parts[2];
            const ang0 = parseFloat(parts[3]);
            const k = parseFloat(parts[4]);
            if (!(ang0 > 0) || !(k > 0)) continue;
            this.angleTypes.push({ a1, a2, a3, ang0, k });
        }
        window.logger.info(`Parsed ${this.angleTypes.length} angle types.`);
        return this.angleTypes;
    }

    getAngleParams(nameA, nameB, nameC) {
        const na = nameA || '*';
        const nb = nameB || '*';
        const nc = nameC || '*';
        for (const at of this.angleTypes) {
            const matchA = (at.a1 === '*') || (at.a1 === na);
            const matchB = (at.a2 === '*') || (at.a2 === nb);
            const matchC = (at.a3 === '*') || (at.a3 === nc);
            if (matchA && matchB && matchC) {
                return { ang0: at.ang0, k: at.k };
            }
        }
        return null;
    }

    convertAngleToDistance(lab, lbc, angleDeg, kAng) {
        const theta = angleDeg * (Math.PI / 180.0);
        const lacSq = (lab * lab) + (lbc * lbc) - (2 * lab * lbc * Math.cos(theta));
        const lac = Math.sqrt(lacSq);
        const sinTheta = Math.sin(theta);
        let kLin = 0;
        if (Math.abs(sinTheta) > 1e-5) {
            const factor = lac / (lab * lbc * sinTheta);
            kLin = kAng * (factor * factor);
        } else {
            if (window.logger && window.logger.warn) {
                window.logger.warn(`Angle ${angleDeg} deg is too close to 0 or 180 degrees. Stiffness set to 0 to avoid explosion.`);
            }
            kLin = 0;
        }
        return { restLength: lac, stiffness: kLin };
    }

    // --- Helpers for Renderer ---

    getColor(iZ) {
        const et = this.byAtomicNumber[iZ];
        if (et) {
            // Convert integer color to [r, g, b] normalized
            const r = ((et.color >> 16) & 255) / 255.0;
            const g = ((et.color >> 8) & 255) / 255.0;
            const b = (et.color & 255) / 255.0;
            return [r, g, b];
        }
        return [1.0, 0.0, 1.0]; // Default Magenta
    }

    getRadius(iZ) {
        const et = this.byAtomicNumber[iZ];
        if (et) return et.RvdW;
        return 1.5; // Default
    }

    getElementName(iZ) {
        const et = this.byAtomicNumber[iZ];
        return et ? et.name : "X";
    }
}
