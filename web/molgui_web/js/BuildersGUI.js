import * as CrystalUtils from './CrystalUtils.js';
import * as PolymerUtils from './PolymerUtils.js';
import { Vec3 } from '../../common_js/Vec3.js';
import { EditableMolecule } from './EditableMolecule.js';
import { GUIutils } from '../../common_js/GUIutils.js';

export class BuildersGUI {
    constructor(gui) {
        this.gui = gui;
    }

    addSubstrateSection(sidebar) {
        const gui = this.gui;
        gui.createSection(sidebar, 'Builder: Substrate', (container) => {
            const rowPreset = GUIutils.row(container);
            GUIutils.span(rowPreset, 'Preset: ');
            const selPreset = GUIutils.selectList(rowPreset, ['NaCl(step)', 'NaCl(rocksalt)', 'KBr(rocksalt)', 'MgO(rocksalt)', 'CaF2(fluorite)', 'CaCO3(todo)'], null, null, { flexGrow: '1' });

            const rowA = GUIutils.row(container);
            GUIutils.span(rowA, 'a(Å): ');
            const inpA = GUIutils.num(rowA, 2.82, { step: '0.01' }, { width: '70px', flexGrow: '0' });

            const MP_PREPARED_CRYSTALS = {
                'C(diamond) mp-66':      { path: '../../cpp/common_resources/crystals/C_diamant_mp-66.json', lattice: { a: 3.567, b: 3.567, c: 3.567, alpha: 90, beta: 90, gamma: 90 } },
                'Si(diamond) mp-149':    { path: '../../cpp/common_resources/crystals/Si_diamond_mp-149.json', lattice: { a: 5.431, b: 5.431, c: 5.431, alpha: 90, beta: 90, gamma: 90 } },
                'CaF2(fluorite) mp-2741': { path: '../../cpp/common_resources/crystals/CaF2_mp-2741.json', lattice: { a: 5.4623, b: 5.4623, c: 5.4623, alpha: 90, beta: 90, gamma: 90 } },
                'CaCO3(calcite) mp-3953 (todo)': { path: '../../cpp/common_resources/crystals/CaCO3_mp-3953.json', lattice: { a: 4.989, b: 4.989, c: 17.062, alpha: 90, beta: 90, gamma: 120 } },
                'CaF2(Pnma) mp-10464 (todo)': { path: '../../cpp/common_resources/crystals/CaF2_mp-10464.json', lattice: null },
            };

            const PRESETS = {
                'NaCl(step)':       { a0: 5.6413 / 2 },
                'NaCl(rocksalt)':   { a0: 5.6413 },
                'KBr(rocksalt)':    { a0: 6.60 },
                'MgO(rocksalt)':    { a0: 4.212 },
                'CaF2(fluorite)':   { a0: 5.4623 },
                'CaCO3(todo)':      { a0: 4.99 }
            };

            const updatePresetDefaults = () => {
                const p = PRESETS[selPreset.value];
                if (p && (p.a0 !== undefined)) inpA.value = String(p.a0);
            };
            selPreset.onchange = () => { updatePresetDefaults(); };
            updatePresetDefaults();

            const rowN = GUIutils.row(container);
            const mkInt = (row, label, val) => (GUIutils.span(row, label, { marginRight: '2px' }), GUIutils.num(row, val, { step: '1' }, { width: '44px', flexGrow: '0' }));
            const inpNx = mkInt(rowN, 'nx', 10);
            const inpNy = mkInt(rowN, 'ny', 10);
            const inpNz = mkInt(rowN, 'nz', 3);

            GUIutils.el(container, 'hr', null, { borderColor: '#444', margin: '8px 0' });

            const rowMP = GUIutils.row(container);
            GUIutils.span(rowMP, 'MP JSON: ');
            const mpKeys = Object.keys(MP_PREPARED_CRYSTALS);
            const selMP = GUIutils.selectList(rowMP, ['(none)', ...mpKeys], null, null, { flexGrow: '1' });
            const btnLoadMP = GUIutils.btn(rowMP, 'Load', null, { marginLeft: '4px', flexGrow: '0' });

            const btnFileMP = GUIutils.btn(container, 'Load crystal JSON file...', null, { marginTop: '4px' });
            const inpFileMP = GUIutils.input(container, { type: 'file', attrs: { accept: '.json,application/json' } }, { display: 'none' });

            GUIutils.el(container, 'hr', null, { borderColor: '#444', margin: '8px 0' });

            const CIF_PREPARED = {
                'C(diamond) sym':    { path: '../../cpp/common_resources/crystals/C_diamond_sym.cif', defaultApplySymmetry: true },
                'C(diamond) nosym':  { path: '../../cpp/common_resources/crystals/C_diamond_nosym.cif', defaultApplySymmetry: false },
                'CaF2 sym':          { path: '../../cpp/common_resources/crystals/CaF2-sym.cif', defaultApplySymmetry: true },
                'CaF2 nosym':        { path: '../../cpp/common_resources/crystals/CaF2-nosym.cif', defaultApplySymmetry: false },
                'CaCO3 sym (Pnma)':  { path: '../../cpp/common_resources/crystals/CaCO3-sym.cif', defaultApplySymmetry: true },
                'CaCO3 nosym':       { path: '../../cpp/common_resources/crystals/CaCO3-nosym.cif', defaultApplySymmetry: false },
                'Si sym':            { path: '../../cpp/common_resources/crystals/Si-sym.cif', defaultApplySymmetry: true },
                'Si nosym':          { path: '../../cpp/common_resources/crystals/Si-nosym.cif', defaultApplySymmetry: false },
            };

            const rowCIF = GUIutils.row(container);
            GUIutils.span(rowCIF, 'CIF: ');
            const cifKeys = Object.keys(CIF_PREPARED);
            const selCIF = GUIutils.selectList(rowCIF, ['(none)', ...cifKeys], null, null, { flexGrow: '1' });
            const btnLoadCIF = GUIutils.btn(rowCIF, 'Load', null, { marginLeft: '4px', flexGrow: '0' });

            const chkSym = GUIutils.labelCheck(container, 'Apply symmetry ops', false, null, { marginTop: '4px' }).input;
            const chkBonds = GUIutils.labelCheck(container, 'Build bonds (BondTypes)', false, null, { marginTop: '4px' }).input;
            const btnFileCIF = GUIutils.btn(container, 'Load CIF file...', null, { marginTop: '4px' });
            const inpFileCIF = GUIutils.input(container, { type: 'file', attrs: { accept: '.cif,.CIF,text/plain' } }, { display: 'none' });

            const chkMiller = GUIutils.labelCheck(container, 'Orient by Miller (h k l) -> z', false, null, { marginTop: '6px' }).input;

            const rowHKL = GUIutils.row(container, { marginTop: '4px' });
            const inpH = mkInt(rowHKL, 'h', 1);
            const inpK = mkInt(rowHKL, 'k', 0);
            const inpL = mkInt(rowHKL, 'l', 0);

            const chkSlab = GUIutils.labelCheck(container, 'Slab cut (cmin/cmax along n)', false, null, { marginTop: '6px' }).input;
            const rowSlab = GUIutils.row(container, { marginTop: '4px' });
            GUIutils.span(rowSlab, 'cmin', { marginRight: '2px' });
            const inpCmin = GUIutils.num(rowSlab, 0.0, { step: '0.1' }, { width: '70px', flexGrow: '0' });
            GUIutils.span(rowSlab, 'cmax', { marginLeft: '6px', marginRight: '2px' });
            const inpCmax = GUIutils.num(rowSlab, 10.0, { step: '0.1' }, { width: '70px', flexGrow: '0' });
            GUIutils.span(rowSlab, 'Å', { marginLeft: '4px', color: '#aaa' });

            const btnRow = GUIutils.row(container, { marginTop: '6px' });
            const btnReplace = GUIutils.btn(btnRow, 'Generate (Replace)');
            const btnAppend = GUIutils.btn(btnRow, 'Generate (Append)', null, { marginLeft: '4px' });

            const getHKL = () => [parseInt(inpH.value), parseInt(inpK.value), parseInt(inpL.value)];

            const applyMiller = (data) => {
                if (!chkMiller.checked) return data;
                const [h, k, l] = getHKL();
                if ((h | 0) === 0 && (k | 0) === 0 && (l | 0) === 0) throw new Error('HKL cannot be (0,0,0)');
                if (!data.lvec) throw new Error('applyMiller: missing lvec');
                const b = CrystalUtils.reciprocalLattice(data.lvec);
                const n = new Vec3().setLincomb3(h, b[0], k, b[1], l, b[2]);
                const R = CrystalUtils.rotationAlignVectorToZ(n);
                return { ...data, Rmiller: R, lvec: CrystalUtils.rotateLvec(data.lvec, R) };
            };

            const getSlab = () => {
                if (!chkSlab.checked) return null;
                const cmin = +inpCmin.value;
                const cmax = +inpCmax.value;
                if (!(cmax > cmin)) throw new Error('Slab cut requires cmax>cmin');
                return { hkl: getHKL(), cmin, cmax };
            };

            const buildData = () => {
                const preset = selPreset.value;
                const a = +inpA.value;
                const nx = parseInt(inpNx.value) | 0;
                const ny = parseInt(inpNy.value) | 0;
                const nz = parseInt(inpNz.value) | 0;
                if (!(a > 0)) throw new Error('a must be >0');
                if (nx <= 0 || ny <= 0 || nz <= 0) throw new Error('nx,ny,nz must be >0');
                const slab = getSlab();
                if (slab && chkBonds.checked) throw new Error('Build bonds is not supported with slab cut');
                let data = null;
                if (preset === 'NaCl(step)') {
                    const mol = CrystalUtils.genNaClStep({ a, nx, ny, nz });
                    data = { mol, lvec: mol.lvec || null };
                } else if (preset.endsWith('(rocksalt)')) {
                    const zA = (preset.startsWith('NaCl')) ? 11 : (preset.startsWith('KBr') ? 19 : 12);
                    const zB = (preset.startsWith('NaCl')) ? 17 : (preset.startsWith('KBr') ? 35 : 8);
                    const lvec = [new Vec3(a, 0, 0), new Vec3(0, a, 0), new Vec3(0, 0, a)];
                    const basis = [
                        [0, 0, 0, zA], [0.5 * a, 0.5 * a, 0, zA], [0.5 * a, 0, 0.5 * a, zA], [0, 0.5 * a, 0.5 * a, zA],
                        [0.5 * a, 0, 0, zB], [0, 0.5 * a, 0, zB], [0, 0, 0.5 * a, zB], [0.5 * a, 0.5 * a, 0.5 * a, zB]
                    ];
                    const basisPos = new Float32Array(basis.length * 3);
                    const basisTypes = new Uint8Array(basis.length);
                    for (let i = 0; i < basis.length; i++) {
                        basisPos[i * 3] = basis[i][0];
                        basisPos[i * 3 + 1] = basis[i][1];
                        basisPos[i * 3 + 2] = basis[i][2];
                        basisTypes[i] = basis[i][3];
                    }
                    if (slab) {
                        const b = CrystalUtils.reciprocalLattice(lvec);
                        const n = new Vec3().setLincomb3(slab.hkl[0] | 0, b[0], slab.hkl[1] | 0, b[1], slab.hkl[2] | 0, b[2]);
                        const ln = n.normalize();
                        if (!(ln > 0)) throw new Error('Slab normal is zero');
                        data = { mol: CrystalUtils.genReplicatedCellSlab({ lvec, basisPos, basisTypes, nRep: [nx, ny, nz], origin: new Vec3(0, 0, 0), nHat: n, cmin: slab.cmin, cmax: slab.cmax }), lvec };
                    } else {
                        data = { mol: CrystalUtils.genReplicatedCell({ lvec, basisPos, basisTypes, nRep: [nx, ny, nz], origin: new Vec3(0, 0, 0) }), lvec };
                    }
                } else if (preset === 'CaF2(fluorite)') {
                    const zCa = 20;
                    const zF = 9;
                    const lvec = [new Vec3(a, 0, 0), new Vec3(0, a, 0), new Vec3(0, 0, a)];
                    const basis = [
                        [0, 0, 0, zCa], [0, 0.5 * a, 0.5 * a, zCa], [0.5 * a, 0, 0.5 * a, zCa], [0.5 * a, 0.5 * a, 0, zCa],
                        [0.25 * a, 0.25 * a, 0.25 * a, zF], [0.25 * a, 0.75 * a, 0.75 * a, zF], [0.75 * a, 0.25 * a, 0.75 * a, zF], [0.75 * a, 0.75 * a, 0.25 * a, zF],
                        [0.75 * a, 0.75 * a, 0.75 * a, zF], [0.75 * a, 0.25 * a, 0.25 * a, zF], [0.25 * a, 0.75 * a, 0.25 * a, zF], [0.25 * a, 0.25 * a, 0.75 * a, zF]
                    ];
                    const basisPos = new Float32Array(basis.length * 3);
                    const basisTypes = new Uint8Array(basis.length);
                    for (let i = 0; i < basis.length; i++) {
                        basisPos[i * 3] = basis[i][0];
                        basisPos[i * 3 + 1] = basis[i][1];
                        basisPos[i * 3 + 2] = basis[i][2];
                        basisTypes[i] = basis[i][3];
                    }
                    if (slab) {
                        const b = CrystalUtils.reciprocalLattice(lvec);
                        const n = new Vec3().setLincomb3(slab.hkl[0] | 0, b[0], slab.hkl[1] | 0, b[1], slab.hkl[2] | 0, b[2]);
                        const ln = n.normalize();
                        if (!(ln > 0)) throw new Error('Slab normal is zero');
                        data = { mol: CrystalUtils.genReplicatedCellSlab({ lvec, basisPos, basisTypes, nRep: [nx, ny, nz], origin: new Vec3(0, 0, 0), nHat: n, cmin: slab.cmin, cmax: slab.cmax }), lvec };
                    } else {
                        data = { mol: CrystalUtils.genReplicatedCell({ lvec, basisPos, basisTypes, nRep: [nx, ny, nz], origin: new Vec3(0, 0, 0) }), lvec };
                    }
                } else if (preset === 'CaCO3(todo)') {
                    throw new Error('CaCO3 preset not implemented yet (non-trivial basis). Use a custom unit cell (e.g. extended XYZ with Lattice=...) once supported.');
                } else {
                    throw new Error(`Preset not implemented: ${preset}`);
                }
                return applyMiller(data);
            };

            const applyToScene = (data, mode) => {
                if (mode === 'replace') gui.system.clear();
                const mol = data.mol;
                if (!mol) throw new Error('applyToScene: missing mol');
                if (data.Rmiller) CrystalUtils.rotateMoleculeInPlace(mol, data.Rmiller);
                const idsNew = new Map();
                for (const a of mol.atoms) {
                    const id = gui.system.addAtom(a.pos.x, a.pos.y, a.pos.z, a.Z);
                    idsNew.set(a.id, id);
                }
                if (mol.bonds && mol.bonds.length) {
                    for (const b of mol.bonds) {
                        const aId = idsNew.get(b.aId);
                        const cId = idsNew.get(b.bId);
                        if (aId === undefined || cId === undefined) throw new Error('applyToScene: bond endpoint missing after copy');
                        gui.system.addBond(aId, cId);
                    }
                }
                gui.renderer.update();
                gui.requestRender();
                window.logger.info(`Substrate generated: atoms=${gui.system.nAtoms}`);
            };

            btnReplace.onclick = () => { try { applyToScene(buildData(), 'replace'); } catch (e) { window.logger.error(String(e)); throw e; } };
            btnAppend.onclick = () => { try { applyToScene(buildData(), 'append'); } catch (e) { window.logger.error(String(e)); throw e; } };

            const buildFromMPJson = (mpJson) => {
                const nx = parseInt(inpNx.value) | 0;
                const ny = parseInt(inpNy.value) | 0;
                const nz = parseInt(inpNz.value) | 0;
                if (nx <= 0 || ny <= 0 || nz <= 0) throw new Error('nx,ny,nz must be >0');
                const slab = getSlab();
                const mm = (window.app && window.app.mmParams) ? window.app.mmParams : null;
                if (slab && chkBonds.checked) throw new Error('Build bonds is not supported with slab cut');
                if (chkBonds.checked && !mm) throw new Error('Build bonds requested but window.app.mmParams is not available');
                const sel = selMP.value;
                const preset = (sel && sel !== '(none)') ? MP_PREPARED_CRYSTALS[sel] : null;
                const lat = preset ? preset.lattice : null;
                if (!lat) throw new Error('MP JSON: missing lattice params for this structure; add to MP_PREPARED_CRYSTALS (BuildersGUI) or supply lattice UI');
                const data0 = CrystalUtils.genCrystalFromMPJson(mpJson, { lattice: lat, nRep: [nx, ny, nz], origin: new Vec3(0, 0, 0), buildBonds: chkBonds.checked, mmParams: mm, slab });
                return applyMiller(data0);
            };

            btnLoadMP.onclick = async () => {
                try {
                    const key = selMP.value;
                    if (!key || key === '(none)') throw new Error('Select a prepared structure');
                    const preset = MP_PREPARED_CRYSTALS[key];
                    if (!preset || !preset.path) throw new Error('Bad preset');
                    const res = await fetch(preset.path);
                    if (!res.ok) throw new Error(`Failed to fetch ${preset.path} (HTTP ${res.status})`);
                    const mpJson = await res.json();
                    applyToScene(buildFromMPJson(mpJson), 'replace');
                } catch (e) { window.logger.error(String(e)); throw e; }
            };

            btnFileMP.onclick = () => inpFileMP.click();
            inpFileMP.onchange = async (e) => {
                try {
                    if (e.target.files.length <= 0) return;
                    const file = e.target.files[0];
                    const txt = await gui.readFileText(file);
                    const mpJson = JSON.parse(txt);
                    applyToScene(buildFromMPJson(mpJson), 'replace');
                    inpFileMP.value = '';
                } catch (err) { window.logger.error(String(err)); throw err; }
            };

            const buildFromCIFText = (cifText) => {
                const nx = parseInt(inpNx.value) | 0;
                const ny = parseInt(inpNy.value) | 0;
                const nz = parseInt(inpNz.value) | 0;
                if (nx <= 0 || ny <= 0 || nz <= 0) throw new Error('nx,ny,nz must be >0');
                const slab = getSlab();
                const mm = (window.app && window.app.mmParams) ? window.app.mmParams : null;
                if (slab && chkBonds.checked) throw new Error('Build bonds is not supported with slab cut');
                if (chkBonds.checked && !mm) throw new Error('Build bonds requested but window.app.mmParams is not available');
                const data0 = CrystalUtils.genCrystalFromCIF(cifText, { applySymmetry: chkSym.checked, nRep: [nx, ny, nz], origin: new Vec3(0, 0, 0), buildBonds: chkBonds.checked, mmParams: mm, slab });
                return applyMiller(data0);
            };

            const updateCIFDefaults = () => {
                const key = selCIF.value;
                const p = (key && key !== '(none)') ? CIF_PREPARED[key] : null;
                if (p) chkSym.checked = !!p.defaultApplySymmetry;
            };
            selCIF.onchange = () => updateCIFDefaults();
            updateCIFDefaults();

            btnLoadCIF.onclick = async () => {
                try {
                    const key = selCIF.value;
                    if (!key || key === '(none)') throw new Error('Select a prepared CIF');
                    const preset = CIF_PREPARED[key];
                    if (!preset || !preset.path) throw new Error('Bad CIF preset');
                    const res = await fetch(preset.path);
                    if (!res.ok) throw new Error(`Failed to fetch ${preset.path} (HTTP ${res.status})`);
                    const txt = await res.text();
                    applyToScene(buildFromCIFText(txt), 'replace');
                } catch (e) { window.logger.error(String(e)); throw e; }
            };

            btnFileCIF.onclick = () => inpFileCIF.click();
            inpFileCIF.onchange = async (e) => {
                try {
                    if (e.target.files.length <= 0) return;
                    const file = e.target.files[0];
                    const txt = await gui.readFileText(file);
                    applyToScene(buildFromCIFText(txt), 'replace');
                    inpFileCIF.value = '';
                } catch (err) { window.logger.error(String(err)); throw err; }
            };
        }, { collapsible: true, open: false });
    }

    addPolymersSection(sidebar) {
        const gui = this.gui;
        gui.createSection(sidebar, 'Builder: Polymers', (container) => {
            const mkLbl = (t) => {
                const d = GUIutils.div(container, null, { marginTop: '4px', fontSize: '0.9em', color: '#ccc' });
                d.textContent = t;
                return d;
            };

            mkLbl('Monomer library (JSON + mol2 files)');

            const btnLoadLib = GUIutils.btn(container, 'Load Library JSON...');
            const inpLib = GUIutils.input(container, { type: 'file', attrs: { accept: '.json,application/json' } }, { display: 'none' });
            const lblLib = GUIutils.div(container, null, { fontSize: '0.85em', color: '#aaa', marginTop: '2px' });
            lblLib.textContent = 'No library loaded';

            btnLoadLib.onclick = () => inpLib.click();
            inpLib.onchange = async (e) => {
                if (e.target.files.length <= 0) return;
                const file = e.target.files[0];
                const txt = await gui.readFileText(file);
                gui.monomerLib = JSON.parse(txt);
                if (!gui.monomerLib || !gui.monomerLib.monomers) throw new Error('Library JSON missing monomers[]');
                lblLib.textContent = `Library: ${gui.monomerLib.name || file.name} (monomers=${gui.monomerLib.monomers.length})`;
                window.logger.info(`Loaded library JSON: ${file.name}`);
                inpLib.value = '';
            };

            const btnLoadMol2 = GUIutils.btn(container, 'Load mol2 geometries...', null, { marginTop: '4px' });
            const inpMol2 = GUIutils.input(container, { type: 'file', attrs: { accept: '.mol2' } }, { display: 'none' });
            inpMol2.multiple = true;
            const lblMol2 = GUIutils.div(container, null, { fontSize: '0.85em', color: '#aaa', marginTop: '2px' });
            lblMol2.textContent = 'mol2 loaded: 0';

            btnLoadMol2.onclick = () => inpMol2.click();
            inpMol2.onchange = async (e) => {
                if (e.target.files.length <= 0) return;
                for (const file of e.target.files) {
                    const txt = await gui.readFileText(file);
                    const parsed = EditableMolecule.parseMol2(txt);
                    gui.monomerGeoms.set(file.name, parsed);
                }
                lblMol2.textContent = `mol2 loaded: ${gui.monomerGeoms.size}`;
                window.logger.info(`Loaded mol2 geometries: ${e.target.files.length}`);
                inpMol2.value = '';
            };

            mkLbl('Sequence');

            const inpSeq = GUIutils.textInput(container, '', { placeholder: 'e.g. D3Gly6A' });

            const rowSeqBtn = GUIutils.row(container, { marginTop: '4px' });
            const btnSeqReplace = GUIutils.btn(rowSeqBtn, 'Build (Replace)');
            const btnSeqAppend = GUIutils.btn(rowSeqBtn, 'Build (Append)', null, { marginLeft: '4px' });

            const buildMonomerMap = () => {
                if (!gui.monomerLib) throw new Error('No monomer library loaded');
                const map = {};
                for (const m of gui.monomerLib.monomers) {
                    if (!m.id) throw new Error('Monomer missing id');
                    if (!m.file) throw new Error(`Monomer '${m.id}' missing file`);
                    const parsed = gui.monomerGeoms.get(m.file);
                    if (!parsed) throw new Error(`Monomer '${m.id}' missing loaded geometry file '${m.file}'`);
                    if (!m.anchors || !m.anchors.head || !m.anchors.tail) throw new Error(`Monomer '${m.id}' missing anchors.head/tail`);
                    if (m.anchors.head.type !== 'index' || m.anchors.tail.type !== 'index') throw new Error(`Monomer '${m.id}': only anchors type='index' supported (v0)`);
                    map[m.id] = { parsed, anchors: [m.anchors.head.value, m.anchors.tail.value] };
                    if (m.aliases) for (const a of m.aliases) map[a] = map[m.id];
                }
                return map;
            };

            const applyPolymerToScene = (poly, mode) => {
                if (!poly || !poly.atoms || !poly.bonds) throw new Error('applyPolymerToScene: poly must be EditableMolecule');
                if (mode === 'replace') gui.system.clear();

                const idsNew = new Map();
                for (const a of poly.atoms) {
                    const id = gui.system.addAtom(a.pos.x, a.pos.y, a.pos.z, a.Z);
                    idsNew.set(a.id, id);
                }
                for (const b of poly.bonds) {
                    const aId = idsNew.get(b.aId);
                    const cId = idsNew.get(b.bId);
                    if (aId === undefined || cId === undefined) throw new Error('applyPolymerToScene: bond endpoint missing after copy');
                    gui.system.addBond(aId, cId);
                }

                gui.renderer.update();
                gui.requestRender();
                window.logger.info(`Polymer built: atoms=${poly.atoms.length} bonds=${poly.bonds.length}`);
            };

            const doBuildSeq = (mode) => {
                const seq = inpSeq.value.trim();
                if (!seq) throw new Error('Empty sequence');
                const tokens = gui.parseSequenceTokens(seq);
                const monomers = buildMonomerMap();
                const poly = PolymerUtils.assemblePolymerFromTokens(tokens, monomers, { _0: 1 });
                applyPolymerToScene(poly, mode);
            };

            btnSeqReplace.onclick = () => { try { doBuildSeq('replace'); } catch (e) { window.logger.error(String(e)); throw e; } };
            btnSeqAppend.onclick = () => { try { doBuildSeq('append'); } catch (e) { window.logger.error(String(e)); throw e; } };

            GUIutils.el(container, 'hr', null, { borderColor: '#444', margin: '8px 0' });

            mkLbl('End-group attachment');

            const btnLoadEnd = GUIutils.btn(container, 'Load Endgroup mol2...');

            const inpEnd = GUIutils.input(container, { type: 'file', attrs: { accept: '.mol2' } }, { display: 'none' });

            const lblEnd = GUIutils.div(container, null, { fontSize: '0.85em', color: '#aaa', marginTop: '2px' });
            lblEnd.textContent = 'No endgroup loaded';

            btnLoadEnd.onclick = () => inpEnd.click();
            inpEnd.onchange = async (e) => {
                if (e.target.files.length <= 0) return;
                const file = e.target.files[0];
                const txt = await gui.readFileText(file);
                gui.endGroupParsed = EditableMolecule.parseMol2(txt);
                lblEnd.textContent = `Endgroup: ${file.name} (atoms=${gui.endGroupParsed.types.length})`;
                window.logger.info(`Loaded endgroup mol2: ${file.name}`);
                inpEnd.value = '';
            };

            const rowMarker = GUIutils.row(container, { marginTop: '4px' });
            const mkSmall = (lab, val) => (GUIutils.span(rowMarker, lab), GUIutils.textInput(rowMarker, val, null, { width: '48px', flexGrow: '0', marginLeft: '2px' }));
            const inpMX = mkSmall('X', 'Se');
            const inpMY = mkSmall('Y', 'Cl');

            const btnAttachMarker = GUIutils.btn(container, 'Attach by markers (all sites)', () => {
                if (!gui.endGroupParsed) throw new Error('No endgroup loaded');
                gui.system.attachGroupByMarker(gui.endGroupParsed, inpMX.value.trim(), inpMY.value.trim());
                gui.renderer.update();
                gui.requestRender();
            }, { marginTop: '4px' });

            GUIutils.el(container, 'hr', null, { borderColor: '#444', margin: '8px 0' });

            mkLbl('Attach by picked direction (cap atom)');

            const rowIds = GUIutils.row(container);
            GUIutils.span(rowIds, 'cap ');
            const inpCap = GUIutils.num(rowIds, '', { placeholder: 'cap id', step: '1' }, { width: '56px', flexGrow: '0' });
            GUIutils.span(rowIds, ' back ');
            const inpBack = GUIutils.num(rowIds, '', { placeholder: 'back id', step: '1' }, { width: '56px', flexGrow: '0' });

            GUIutils.btn(container, 'Use selection -> cap/back', () => {
                const ids = Array.from(gui.system.selection).sort((a, b) => a - b);
                if (ids.length <= 0) throw new Error('Nothing selected');
                inpCap.value = String(ids[0]);
                inpBack.value = (ids.length > 1) ? String(ids[1]) : '';
            }, { marginTop: '4px' });

            const rowGeom = GUIutils.row(container, { marginTop: '4px' });
            GUIutils.span(rowGeom, 'bond ');
            const inpBond = GUIutils.num(rowGeom, 1.50, { step: '0.01' }, { width: '62px', flexGrow: '0' });
            GUIutils.span(rowGeom, ' twist° ');
            const inpTw = GUIutils.num(rowGeom, 0, { step: '1' }, { width: '62px', flexGrow: '0' });

            const rowUp = GUIutils.row(container, { marginTop: '4px' });
            GUIutils.span(rowUp, 'up ');
            const mkF = (val) => GUIutils.num(rowUp, val, { step: '0.1' }, { width: '44px', flexGrow: '0' });
            const upx = mkF(0);
            const upy = mkF(0);
            const upz = mkF(1);

            const rowRefs = GUIutils.row(container, { marginTop: '4px' });
            const mkIntRef = (lbl, val) => (GUIutils.span(rowRefs, lbl + ' '), GUIutils.num(rowRefs, val, { step: '1' }, { width: '44px', flexGrow: '0' }));
            const gA = mkIntRef('gA', 1);
            const gF = mkIntRef('gF', 2);
            const gU = mkIntRef('gU', 0);

            GUIutils.btn(container, 'Attach (direction)', () => {
                if (!gui.endGroupParsed) throw new Error('No endgroup loaded');
                const cap = parseInt(inpCap.value);
                if (!(cap >= 0)) throw new Error('cap id missing');
                const back = inpBack.value.trim() ? parseInt(inpBack.value) : undefined;
                const up = [parseFloat(upx.value), parseFloat(upy.value), parseFloat(upz.value)];
                const params = {
                    backAtom: (back !== undefined) ? back : undefined,
                    bondLen: parseFloat(inpBond.value),
                    up,
                    twistDeg: parseFloat(inpTw.value),
                    groupAnchor: parseInt(gA.value),
                    groupForwardRef: parseInt(gF.value),
                    groupUpRef: parseInt(gU.value)
                };
                gui.system.attachParsedByDirection(cap, gui.endGroupParsed, params);
                gui.renderer.update();
                gui.requestRender();
            }, { marginTop: '4px' });
        }, { collapsible: true, open: false });
    }
}
