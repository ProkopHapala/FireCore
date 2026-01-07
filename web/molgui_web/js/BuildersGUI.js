import * as CrystalUtils from './CrystalUtils.js';
import * as PolymerUtils from './MoleculeUtils.js';
import { Vec3 } from '../../common_js/Vec3.js';
import { EditableMolecule } from './EditableMolecule.js';
import { MoleculeRenderer, PackedMolecule } from './MoleculeRenderer.js';
import { buildCrystalCellBucketsFromMol } from '../../common_js/Buckets.js';
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
            rowPreset.style.display = 'none';

            const rowA = GUIutils.row(container);
            GUIutils.span(rowA, 'a(Å): ');
            const inpA = GUIutils.num(rowA, 2.82, { step: '0.01' }, { width: '70px', flexGrow: '0' });

            const rowRep = GUIutils.row(container, { marginTop: '4px' });
            const chkRep = GUIutils.labelCheck(rowRep, 'Visual Replicas', false, (e) => {
                if (window.app) {
                    const lat = window.app.getLattice('substrate');
                    lat.show = e.target.checked;
                    window.app.updateReplicas('substrate');
                }
            }).input;

            const rowNRep = GUIutils.row(container, { marginTop: '4px' });
            GUIutils.span(rowNRep, 'nRep: ');
            const inpNxRep = GUIutils.num(rowNRep, 1, { min: 0, step: 1 }, { width: '40px' });
            const inpNyRep = GUIutils.num(rowNRep, 1, { min: 0, step: 1 }, { width: '40px', marginLeft: '5px' });
            const inpNzRep = GUIutils.num(rowNRep, 1, { min: 0, step: 1 }, { width: '40px', marginLeft: '5px' });

            const updateSubstrateRep = () => {
                if (window.app) {
                    const lat = window.app.getLattice('substrate');
                    lat.nrep.x = parseInt(inpNxRep.value) || 0;
                    lat.nrep.y = parseInt(inpNyRep.value) || 0;
                    lat.nrep.z = parseInt(inpNzRep.value) || 0;
                    const a = parseFloat(inpA.value) || 2.82;
                    const nx = parseInt(inpNx.value) || 10;
                    const ny = parseInt(inpNy.value) || 10;
                    const nz = parseInt(inpNz.value) || 3;
                    lat.lvec[0].set(a * nx, 0, 0);
                    lat.lvec[1].set(0, a * ny, 0);
                    lat.lvec[2].set(0, 0, a * nz);
                    window.app.updateReplicas('substrate');
                }
            };
            inpNxRep.onchange = updateSubstrateRep;
            inpNyRep.onchange = updateSubstrateRep;
            inpNzRep.onchange = updateSubstrateRep;
            inpA.onchange = updateSubstrateRep;


            // const PRESETS = {
            //     'NaCl(step)':       { a0: 5.6413 / 2 },
            //     'NaCl(rocksalt)':   { a0: 5.6413 },
            //     'KBr(rocksalt)':    { a0: 6.60 },
            //     'MgO(rocksalt)':    { a0: 4.212 },
            //     'CaF2(fluorite)':   { a0: 5.4623 },
            //     'CaCO3(todo)':      { a0: 4.99 }
            // };

            // const updatePresetDefaults = () => {
            //     const p = PRESETS[selPreset.value];
            //     if (p && (p.a0 !== undefined)) inpA.value = String(p.a0);
            // };
            // selPreset.onchange = () => { updatePresetDefaults(); };
            // updatePresetDefaults();

            const rowN = GUIutils.row(container);
            const mkInt = (row, label, val) => (GUIutils.span(row, label, { marginRight: '2px' }), GUIutils.num(row, val, { step: '1' }, { width: '44px', flexGrow: '0' }));
            const inpNx = mkInt(rowN, 'nx', 10);
            const inpNy = mkInt(rowN, 'ny', 10);
            const inpNz = mkInt(rowN, 'nz', 3);

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
            const btnFileCIF = GUIutils.btn(rowCIF, 'Load CIF', null, { marginLeft: '4px', flexGrow: '0' });
            const inpFileCIF = GUIutils.input(container, { type: 'file', attrs: { accept: '.cif,.CIF,text/plain' } }, { display: 'none' });

            const chkSym = GUIutils.labelCheck(container, 'sym', false, null, { marginTop: '4px' }).input;
            const chkDedupSym = GUIutils.labelCheck(container, 'dedup sym', true, null, { marginTop: '4px' }).input;
            const chkDedupRep = GUIutils.labelCheck(container, 'dedup rep', true, null, { marginTop: '4px' }).input;
            const rowDedupTol = GUIutils.row(container, { marginTop: '4px' });
            GUIutils.span(rowDedupTol, 'dedup tol', { marginRight: '4px' });
            const inpDedupTol = GUIutils.num(rowDedupTol, 0.1, { step: '0.01' }, { width: '70px', flexGrow: '0' });
            GUIutils.span(rowDedupTol, 'Å', { marginLeft: '4px', color: '#aaa' });
            const chkBonds = GUIutils.labelCheck(container, 'bonds', false, null, { marginTop: '4px' }).input;

            GUIutils.el(container, 'hr', null, { borderColor: '#444', margin: '8px 0' });

            const rowCellMode = GUIutils.row(container);
            GUIutils.span(rowCellMode, 'Cell coords: ', { marginRight: '4px' });
            const selCellMode = GUIutils.selectList(rowCellMode, ['fractional', 'cartesian'], 'fractional', null, { flexGrow: '1' });
            let _cellModePrev = 'fractional';
            const chkPreview = GUIutils.labelCheck(container, 'Preview', true, null, { marginTop: '4px' }).input;
            const chkAutoSync = GUIutils.labelCheck(container, 'Auto-sync', true, null, { marginTop: '4px' }).input;

            const taLat = GUIutils.textArea(container, '', { display: 'block', height: '72px', marginTop: '4px', placeholder: 'Lattice vectors (3 lines):\nax ay az\nbx by bz\ncx cy cz' });
            const taAtoms = GUIutils.textArea(container, '', { display: 'block', height: '140px', marginTop: '4px', placeholder: 'Atoms in unit cell (XYZ-like):\nEl x y z\n...' });
            const taSym = GUIutils.textArea(container, '', { display: 'block', height: '90px', marginTop: '4px', placeholder: 'Symmetry operations (one per line):\nx,y,z\n-x,-y,-z\n...' });

            const rowCellBtns = GUIutils.row(container, { marginTop: '6px' });
            const btnCellPreview = GUIutils.btn(rowCellBtns, 'Preview cell', null, { flexGrow: '1' });
            const btnCellGenerateReplace = GUIutils.btn(rowCellBtns, 'Generate from cell (Replace)', null, { marginLeft: '4px', flexGrow: '1' });
            const btnCellGenerateAppend = GUIutils.btn(rowCellBtns, '... Append', null, { marginLeft: '4px', flexGrow: '0' });

            const lblCellStatus = GUIutils.div(container, null, { marginTop: '4px', fontSize: '0.85em', color: '#aaa' });
            lblCellStatus.textContent = '';

            // --- Unit-cell preview renderer (non-editable) ---
            let previewMol = null;
            let previewPacked = null;
            let previewRenderer = null;
            let previewPlanes = null;
            let previewCellBox = null;
            const ensurePreview = () => {
                if (previewRenderer) return;
                if (!window.app || !window.app.scene || !window.app.shaders || !window.app.mmParams) throw new Error('Unit-cell preview requires window.app.scene/shaders/mmParams');
                previewMol = new EditableMolecule();
                previewPacked = new PackedMolecule(1024);
                previewRenderer = new MoleculeRenderer(window.app.scene, previewPacked, window.app.shaders, window.app.mmParams, previewMol);
                if (previewRenderer.atomMesh) previewRenderer.atomMesh.renderOrder = 10;
                if (previewRenderer.bondLines) previewRenderer.bondLines.renderOrder = 9;
            };

            const ensurePreviewPlanes = () => {
                if (previewPlanes) return;
                if (!window.app || !window.app.scene) throw new Error('Preview planes require window.app.scene');
                const THREE = window.THREE;
                if (!THREE) throw new Error('Preview planes require window.THREE');
                const geom = new THREE.BufferGeometry();
                geom.setAttribute('position', new THREE.BufferAttribute(new Float32Array(0), 3));
                const mat = new THREE.LineBasicMaterial({ color: 0x44ff88, transparent: true, opacity: 0.9 });
                previewPlanes = new THREE.LineSegments(geom, mat);
                previewPlanes.renderOrder = 11;
                window.app.scene.add(previewPlanes);
            };

            const ensurePreviewCellBox = () => {
                if (previewCellBox) return;
                if (!window.app || !window.app.scene) throw new Error('Preview cell box requires window.app.scene');
                const THREE = window.THREE;
                if (!THREE) throw new Error('Preview cell box requires window.THREE');
                const geom = new THREE.BufferGeometry();
                geom.setAttribute('position', new THREE.BufferAttribute(new Float32Array(0), 3));
                const mat = new THREE.LineBasicMaterial({ color: 0x4488ff, transparent: true, opacity: 0.7 });
                previewCellBox = new THREE.LineSegments(geom, mat);
                previewCellBox.renderOrder = 8;
                window.app.scene.add(previewCellBox);
            };

            const setPreviewVisible = (b) => {
                const v = !!b;
                if (previewRenderer) {
                    if (previewRenderer.atomMesh) previewRenderer.atomMesh.visible = v;
                    if (previewRenderer.bondLines) previewRenderer.bondLines.visible = v;
                }
                if (previewPlanes) previewPlanes.visible = v;
                if (previewCellBox) previewCellBox.visible = v;
                gui.requestRender();
            };

            const fillCellTextAreasFromCIF = (cifText) => {
                const crystal = CrystalUtils.cifToCrystalData(cifText);
                const lvec = CrystalUtils.latticeVectorsFromParams(crystal.lattice);
                inpA.value = String(+crystal.lattice.a);
                const a0 = +inpA.value;
                if (!(a0 > 0)) throw new Error('CIF: invalid a');
                const lvecN = [lvec[0].clone().mulScalar(1.0 / a0), lvec[1].clone().mulScalar(1.0 / a0), lvec[2].clone().mulScalar(1.0 / a0)];
                taLat.value = CrystalUtils.formatLatticeText(lvecN);
                taSym.value = CrystalUtils.formatSymOpsText(crystal.symOpStrs || []);
                if (selCellMode.value === 'cartesian') {
                    taAtoms.value = CrystalUtils.formatSitesTextXYZCart(crystal.sites, lvec);
                } else {
                    taAtoms.value = CrystalUtils.formatSitesTextXYZFrac(crystal.sites);
                }
                lblCellStatus.textContent = `Loaded cell: sites=${crystal.sites.length} symOps=${(crystal.symOpStrs || []).length}`;
                if (chkPreview.checked) updatePreviewFromEditor();
            };

            const parseUnitCellEditor = () => {
                const dedupTol = +inpDedupTol.value;
                if (!(dedupTol > 0)) throw new Error('Dedup tol must be >0');
                const a0 = +inpA.value;
                if (!(a0 > 0)) throw new Error('a must be >0');
                let lvecN = null;
                try {
                    lvecN = CrystalUtils.parseLatticeText(taLat.value);
                } catch (e) {
                    const s = (taLat.value || '').trim();
                    if (s.length === 0) {
                        taLat.value = '1 0 0\n0 1 0\n0 0 1\n';
                        lvecN = CrystalUtils.parseLatticeText(taLat.value);
                    } else {
                        throw e;
                    }
                }
                const lvec = [lvecN[0].clone().mulScalar(a0), lvecN[1].clone().mulScalar(a0), lvecN[2].clone().mulScalar(a0)];
                const mode = selCellMode.value;
                const atomsTxt = CrystalUtils.parseSitesTextXYZ(taAtoms.value, mode);
                const sitesFrac = [];
                const p = new Vec3();
                const f = [0, 0, 0];
                if (mode === 'cartesian') {
                    for (const a of atomsTxt) {
                        p.set(+a.x, +a.y, +a.z);
                        CrystalUtils.cartToFrac(p, lvec, f);
                        sitesFrac.push({ element: a.element, x: f[0], y: f[1], z: f[2], q: +a.q });
                    }
                } else {
                    for (const a of atomsTxt) sitesFrac.push({ element: a.element, x: +a.x, y: +a.y, z: +a.z, q: +a.q });
                }

                let sites = sitesFrac;
                const symOps = CrystalUtils.parseSymOpsText(taSym.value);
                if (chkSym.checked) {
                    if (!symOps || symOps.length === 0) throw new Error('Apply symmetry ops enabled but no sym ops provided');
                    sites = CrystalUtils.applySymmetryOpsFracSites(sites, symOps, { tol: 1e-6 });
                }
                if (chkDedupSym.checked) {
                    sites = CrystalUtils.dedupFracSitesByTolA(sites, lvec, dedupTol);
                }

                const cell = CrystalUtils.cellDataFromFracSites(lvec, sites);
                return { lvec, sitesFrac: sites, cell };
            };

            const updatePreviewPlanesFromEditor = () => {
                if (!chkPreview.checked) return;
                ensurePreviewPlanes();
                const THREE = window.THREE;
                const { cell } = parseUnitCellEditor();
                const planes = getPlanes();
                const slab = planes ? null : getSlab();

                const pls = [];
                const b = CrystalUtils.reciprocalLattice(cell.lvec);
                if (planes) {
                    for (const p of planes.planes) {
                        const n = new Vec3().setLincomb3(p.h, b[0], p.k, b[1], p.l, b[2]);
                        pls.push({ n, cmin: p.cmin, cmax: p.cmax, mode: planes.planeMode });
                    }
                } else if (slab) {
                    const n = new Vec3().setLincomb3(slab.hkl[0] | 0, b[0], slab.hkl[1] | 0, b[1], slab.hkl[2] | 0, b[2]);
                    pls.push({ n, cmin: slab.cmin, cmax: slab.cmax, mode: 'ang' });
                }

                if (pls.length === 0) {
                    previewPlanes.geometry.setAttribute('position', new THREE.BufferAttribute(new Float32Array(0), 3));
                    previewPlanes.geometry.computeBoundingSphere();
                    return;
                }

                const L = Math.max(cell.lvec[0].norm(), cell.lvec[1].norm(), cell.lvec[2].norm());
                const R = 0.6 * L;
                const verts = [];
                const tmp = new Vec3();
                for (const pl of pls) {
                    const n = pl.n.clone();
                    if (pl.mode === 'ang') n.normalize();
                    const aN = Math.abs(n.x), bN = Math.abs(n.y), cN = Math.abs(n.z);
                    const ref = (aN < 0.9) ? new Vec3(1, 0, 0) : new Vec3(0, 1, 0);
                    const u = new Vec3().setCross(n, ref);
                    const lu = u.normalize();
                    if (!(lu > 0)) continue;
                    const v = new Vec3().setCross(n, u);
                    v.normalize();
                    for (const d0 of [pl.cmin, pl.cmax]) {
                        const c0 = tmp.setV(n).mulScalar(d0);
                        const p00 = new Vec3().setV(c0).addMul(u, +R).addMul(v, +R);
                        const p10 = new Vec3().setV(c0).addMul(u, -R).addMul(v, +R);
                        const p11 = new Vec3().setV(c0).addMul(u, -R).addMul(v, -R);
                        const p01 = new Vec3().setV(c0).addMul(u, +R).addMul(v, -R);
                        const ps = [p00, p10, p11, p01];
                        for (let i = 0; i < 4; i++) {
                            const a = ps[i];
                            const b = ps[(i + 1) & 3];
                            verts.push(a.x, a.y, a.z, b.x, b.y, b.z);
                        }
                    }
                }
                previewPlanes.geometry.setAttribute('position', new THREE.BufferAttribute(new Float32Array(verts), 3));
                previewPlanes.geometry.computeBoundingSphere();
            };

            const updatePreviewCellBoxFromEditor = () => {
                if (!chkPreview.checked) return;
                ensurePreviewCellBox();
                const THREE = window.THREE;
                const { cell } = parseUnitCellEditor();
                const a = cell.lvec[0], b = cell.lvec[1], c = cell.lvec[2];
                const o = new Vec3(0, 0, 0);
                const p100 = new Vec3().setV(a);
                const p010 = new Vec3().setV(b);
                const p001 = new Vec3().setV(c);
                const p110 = new Vec3().setV(a).add(b);
                const p101 = new Vec3().setV(a).add(c);
                const p011 = new Vec3().setV(b).add(c);
                const p111 = new Vec3().setV(a).add(b).add(c);
                const edges = [
                    o, p100, o, p010, o, p001,
                    p100, p110, p100, p101,
                    p010, p110, p010, p011,
                    p001, p101, p001, p011,
                    p110, p111, p101, p111, p011, p111,
                ];
                const verts = new Float32Array(edges.length * 3);
                for (let i = 0; i < edges.length; i++) {
                    const v = edges[i];
                    const i3 = i * 3;
                    verts[i3] = v.x; verts[i3 + 1] = v.y; verts[i3 + 2] = v.z;
                }
                previewCellBox.geometry.setAttribute('position', new THREE.BufferAttribute(verts, 3));
                previewCellBox.geometry.computeBoundingSphere();
            };

            const updatePreviewFromEditor = () => {
                if (!chkPreview.checked) return;
                ensurePreview();
                const { cell, sitesFrac } = parseUnitCellEditor();
                previewMol.clear();
                const n = cell.basisTypes.length | 0;
                for (let i = 0; i < n; i++) {
                    const i3 = i * 3;
                    const id = previewMol.addAtom(cell.basisPos[i3], cell.basisPos[i3 + 1], cell.basisPos[i3 + 2], cell.basisTypes[i]);
                    const ia = previewMol.getAtomIndex(id);
                    if (ia >= 0 && cell.basisCharges) previewMol.atoms[ia].charge = +cell.basisCharges[i];
                }
                previewMol.lvec = cell.lvec;
                previewRenderer.update();
                gui.renderer.update();
                gui.requestRender();
                lblCellStatus.textContent = `Preview: atoms=${n} (sites=${sitesFrac.length})`;
                updatePreviewPlanesFromEditor();
                updatePreviewCellBoxFromEditor();
            };

            let _syncTimer = null;
            const scheduleAutoSync = () => {
                if (!chkPreview.checked) return;
                if (!chkAutoSync.checked) return;
                if (_syncTimer) clearTimeout(_syncTimer);
                _syncTimer = setTimeout(() => {
                    _syncTimer = null;
                    try { updatePreviewFromEditor(); } catch (e) { lblCellStatus.textContent = String(e); throw e; }
                }, 250);
            };

            chkPreview.onchange = () => {
                try {
                    if (chkPreview.checked) {
                        ensurePreview();
                        ensurePreviewPlanes();
                        ensurePreviewCellBox();
                        setPreviewVisible(true);
                        updatePreviewFromEditor();
                    } else {
                        setPreviewVisible(false);
                    }
                } catch (e) { lblCellStatus.textContent = String(e); throw e; }
            };

            taLat.oninput = () => scheduleAutoSync();
            taAtoms.oninput = () => scheduleAutoSync();
            taSym.oninput = () => scheduleAutoSync();
            chkSym.onchange = () => scheduleAutoSync();
            chkDedupSym.onchange = () => scheduleAutoSync();
            selCellMode.onchange = () => {
                try {
                    const lvec = CrystalUtils.parseLatticeText(taLat.value);
                    const modeNew = selCellMode.value;
                    const modeOld = _cellModePrev;
                    const sites = CrystalUtils.parseSitesTextXYZ(taAtoms.value, modeOld);
                    if (selCellMode.value === 'cartesian') {
                        // switch frac->cart for display
                        const sitesFrac = sites.map(a => ({ element: a.element, x: +a.x, y: +a.y, z: +a.z }));
                        taAtoms.value = CrystalUtils.formatSitesTextXYZCart(sitesFrac, lvec);
                    } else {
                        // switch cart->frac for display
                        const sitesFrac = [];
                        const p = new Vec3();
                        const f = [0, 0, 0];
                        for (const a of sites) {
                            p.set(+a.x, +a.y, +a.z);
                            CrystalUtils.cartToFrac(p, lvec, f);
                            sitesFrac.push({ element: a.element, x: f[0], y: f[1], z: f[2] });
                        }
                        taAtoms.value = CrystalUtils.formatSitesTextXYZFrac(sitesFrac);
                    }
                    _cellModePrev = modeNew;
                } catch (e) { lblCellStatus.textContent = String(e); throw e; }
                scheduleAutoSync();
            };

            btnCellPreview.onclick = () => { try { chkPreview.checked = true; updatePreviewFromEditor(); } catch (e) { window.logger.error(String(e)); throw e; } };

            const buildFromUnitCellEditor = (mode) => {
                const nx = parseInt(inpNx.value) | 0;
                const ny = parseInt(inpNy.value) | 0;
                const nz = parseInt(inpNz.value) | 0;
                const dedupTol = +inpDedupTol.value;
                if (!(dedupTol > 0)) throw new Error('Dedup tol must be >0');
                if (nx <= 0 || ny <= 0 || nz <= 0) throw new Error('nx,ny,nz must be >0');
                const planes = getPlanes();
                const slab = planes ? null : getSlab();
                if ((slab || planes) && chkBonds.checked) throw new Error('Build bonds is not supported with cut planes');
                if (planes) lblCellStatus.textContent = `Cut planes: n=${planes.planes.length} mode=${planes.planeMode} centered=${planes.centered}`;
                else lblCellStatus.textContent = slab ? `Slab: HKL=(${slab.hkl[0]},${slab.hkl[1]},${slab.hkl[2]}) c=[${slab.cmin},${slab.cmax}]` : 'Bulk replication';

                const mm = (window.app && window.app.mmParams) ? window.app.mmParams : null;
                if (chkBonds.checked && !mm) throw new Error('Build bonds requested but window.app.mmParams is not available');

                const { cell } = parseUnitCellEditor();
                let mol = null;
                if (planes) {
                    const b = CrystalUtils.reciprocalLattice(cell.lvec);
                    const pls = [];
                    for (const p of planes.planes) {
                        const n = new Vec3().setLincomb3(p.h, b[0], p.k, b[1], p.l, b[2]);
                        pls.push({ n, cmin: p.cmin, cmax: p.cmax });
                    }
                    mol = CrystalUtils.genReplicatedCellCutPlanes({ lvec: cell.lvec, basisPos: cell.basisPos, basisTypes: cell.basisTypes, basisCharges: cell.basisCharges, nRep: [nx, ny, nz], origin: new Vec3(0, 0, 0), planes: pls, planeMode: planes.planeMode, centered: planes.centered, dedup: chkDedupRep.checked, dedupTol });
                } else if (slab) {
                    const b = CrystalUtils.reciprocalLattice(cell.lvec);
                    const n = new Vec3().setLincomb3(slab.hkl[0] | 0, b[0], slab.hkl[1] | 0, b[1], slab.hkl[2] | 0, b[2]);
                    const ln = n.normalize();
                    if (!(ln > 0)) throw new Error('Slab normal is zero');
                    mol = CrystalUtils.genReplicatedCellSlab({ lvec: cell.lvec, basisPos: cell.basisPos, basisTypes: cell.basisTypes, basisCharges: cell.basisCharges, nRep: [nx, ny, nz], origin: new Vec3(0, 0, 0), nHat: n, cmin: slab.cmin, cmax: slab.cmax, dedup: chkDedupRep.checked, dedupTol });
                } else {
                    mol = CrystalUtils.genReplicatedCell({ lvec: cell.lvec, basisPos: cell.basisPos, basisTypes: cell.basisTypes, basisCharges: cell.basisCharges, nRep: [nx, ny, nz], origin: new Vec3(0, 0, 0), buildBonds: chkBonds.checked, mmParams: mm, dedup: chkDedupRep.checked, dedupTol });
                }
                // Charge neutrality report
                let qTot = 0.0;
                for (const a of mol.atoms) qTot += (a.charge !== undefined && a.charge !== null) ? +a.charge : 0.0;
                lblCellStatus.textContent = `${lblCellStatus.textContent} | atoms=${mol.atoms.length} | Q=${qTot}`;
                if (Math.abs(qTot) > 1e-3) window.logger.error(`Total charge not neutral: Q=${qTot}`);
                else window.logger.info(`Total charge Q=${qTot}`);
                const data0 = { mol, lvec: mol.lvec || null, bucketCells: { na: nx, nb: ny, nc: nz, lvecCell: cell.lvec } };
                applyToScene(applyMiller(data0), mode);
            };

            btnCellGenerateReplace.onclick = () => { try { buildFromUnitCellEditor('replace'); } catch (e) { window.logger.error(String(e)); throw e; } };
            btnCellGenerateAppend.onclick = () => { try { buildFromUnitCellEditor('append'); } catch (e) { window.logger.error(String(e)); throw e; } };

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

            const chkPlanes = GUIutils.labelCheck(container, 'Cut by planes', false, null, { marginTop: '6px' }).input;
            const rowPlanesOpts = GUIutils.row(container, { marginTop: '4px' });
            const chkPlaneFrac = GUIutils.labelCheck(rowPlanesOpts, 'c fractional', false, null, { marginRight: '8px' }).input;
            const chkCentered = GUIutils.labelCheck(rowPlanesOpts, 'centered rep', true, null).input;
            const rowPlaneBtns = GUIutils.row(container, { marginTop: '4px' });
            const btnAdd100 = GUIutils.btn(rowPlaneBtns, '+{100}', null, { flexGrow: '1' });
            const btnAdd110 = GUIutils.btn(rowPlaneBtns, '+{110}', null, { marginLeft: '4px', flexGrow: '1' });
            const btnAdd111 = GUIutils.btn(rowPlaneBtns, '+{111}', null, { marginLeft: '4px', flexGrow: '1' });
            const taPlanes = GUIutils.textArea(container, '', { display: 'block', height: '90px', marginTop: '4px', placeholder: 'Planes (per line):\nh k l cmin cmax\n(or) h k l c  => symmetric [-c,+c]' });

            taPlanes.oninput = () => scheduleAutoSync();
            chkSlab.onchange = () => scheduleAutoSync();
            chkPlanes.onchange = () => scheduleAutoSync();
            chkPlaneFrac.onchange = () => scheduleAutoSync();
            chkCentered.onchange = () => scheduleAutoSync();

            const btnRow = GUIutils.row(container, { marginTop: '6px' });
            const btnReplace = GUIutils.btn(btnRow, 'Generate (Replace)');
            const btnAppend = GUIutils.btn(btnRow, 'Generate (Append)', null, { marginLeft: '4px' });
            // UI cleanup: hide old generation buttons (generation should go through unit-cell editor)
            btnRow.style.display = 'none';

            const getHKL = () => [parseInt(inpH.value), parseInt(inpK.value), parseInt(inpL.value)];

            const applyMiller = (data) => {
                if (!chkMiller.checked) return data;
                const [h, k, l] = getHKL();
                if ((h | 0) === 0 && (k | 0) === 0 && (l | 0) === 0) throw new Error('HKL cannot be (0,0,0)');
                if (!data.lvec) throw new Error('applyMiller: missing lvec');
                const b = CrystalUtils.reciprocalLattice(data.lvec);
                const n = new Vec3().setLincomb3(h, b[0], k, b[1], l, b[2]);
                const R = CrystalUtils.rotationAlignVectorToZ(n);
                const out = { ...data, Rmiller: R, lvec: CrystalUtils.rotateLvec(data.lvec, R) };
                if (data.bucketCells && data.bucketCells.lvecCell) out.bucketCells = { ...data.bucketCells, lvecCell: CrystalUtils.rotateLvec(data.bucketCells.lvecCell, R) };
                return out;
            };

            const getSlab = () => {
                if (!chkSlab.checked) return null;
                const cmin = +inpCmin.value;
                const cmax = +inpCmax.value;
                if (!(cmax > cmin)) throw new Error('Slab cut requires cmax>cmin');
                const hkl = getHKL();
                if (!Number.isFinite(hkl[0]) || !Number.isFinite(hkl[1]) || !Number.isFinite(hkl[2])) throw new Error('Slab cut requires integer HKL');
                if (((hkl[0] | 0) === 0) && ((hkl[1] | 0) === 0) && ((hkl[2] | 0) === 0)) throw new Error('Slab cut requires non-zero HKL');
                return { hkl, cmin, cmax };
            };

            const getPlanes = () => {
                if (!chkPlanes.checked) return null;
                const lines = (taPlanes.value || '').split(/\r?\n/);
                const out = [];
                for (const ln of lines) {
                    const s = ln.trim();
                    if (!s) continue;
                    const ws = s.split(/\s+/);
                    if (ws.length !== 4 && ws.length !== 5) throw new Error(`Planes: expected 4 or 5 columns, got '${ln}'`);
                    const h = +ws[0], k = +ws[1], l = +ws[2];
                    if (!Number.isFinite(h) || !Number.isFinite(k) || !Number.isFinite(l)) throw new Error(`Planes: invalid hkl in '${ln}'`);
                    let cmin = 0.0, cmax = 0.0;
                    if (ws.length === 4) {
                        const c = +ws[3];
                        if (!Number.isFinite(c)) throw new Error(`Planes: invalid c in '${ln}'`);
                        cmin = -c; cmax = +c;
                    } else {
                        cmin = +ws[3]; cmax = +ws[4];
                        if (!Number.isFinite(cmin) || !Number.isFinite(cmax)) throw new Error(`Planes: invalid cmin/cmax in '${ln}'`);
                    }
                    if (!(cmax > cmin)) throw new Error(`Planes: requires cmax>cmin in '${ln}'`);
                    out.push({ h, k, l, cmin, cmax });
                }
                if (out.length === 0) throw new Error('Planes: no planes specified');
                return { planes: out, planeMode: chkPlaneFrac.checked ? 'frac' : 'ang', centered: !!chkCentered.checked };
            };

            const _appendPlaneFamily = (hkls) => {
                chkPlanes.checked = true;
                const c = +inpCmax.value;
                if (!(c > 0)) throw new Error('Planes: set cmax>0 (used as symmetric c)');
                const existing = new Set();
                const lines0 = (taPlanes.value || '').split(/\r?\n/);
                for (const ln of lines0) {
                    const s = ln.trim();
                    if (!s) continue;
                    const ws = s.split(/\s+/);
                    if (ws.length < 3) continue;
                    const h = +ws[0], k = +ws[1], l = +ws[2];
                    if (!Number.isFinite(h) || !Number.isFinite(k) || !Number.isFinite(l)) continue;
                    let hs = h, ks = k, ls = l;
                    if (hs < 0 || (hs === 0 && ks < 0) || (hs === 0 && ks === 0 && ls < 0)) { hs = -hs; ks = -ks; ls = -ls; }
                    existing.add(`${hs},${ks},${ls}`);
                }
                const out = [];
                for (const v of hkls) {
                    let h = v[0], k = v[1], l = v[2];
                    if (h < 0 || (h === 0 && k < 0) || (h === 0 && k === 0 && l < 0)) { h = -h; k = -k; l = -l; }
                    const key = `${h},${k},${l}`;
                    if (existing.has(key)) continue;
                    existing.add(key);
                    out.push(`${h} ${k} ${l} ${c}`);
                }
                if (out.length === 0) return;
                const sep = (taPlanes.value && taPlanes.value.trim().length > 0) ? '\n' : '';
                taPlanes.value = `${taPlanes.value}${sep}${out.join('\n')}\n`;
            };

            btnAdd100.onclick = () => { try { _appendPlaneFamily([[1, 0, 0], [0, 1, 0], [0, 0, 1]]); } catch (e) { window.logger.error(String(e)); throw e; } };
            btnAdd110.onclick = () => { try { _appendPlaneFamily([[1, 1, 0], [1, 0, 1], [0, 1, 1], [1, -1, 0], [1, 0, -1], [0, 1, -1]]); } catch (e) { window.logger.error(String(e)); throw e; } };
            btnAdd111.onclick = () => { try { _appendPlaneFamily([[1, 1, 1], [1, 1, -1], [1, -1, 1], [-1, 1, 1]]); } catch (e) { window.logger.error(String(e)); throw e; } };

            const buildData = () => {
                const preset = selPreset.value;
                const a = +inpA.value;
                const nx = parseInt(inpNx.value) | 0;
                const ny = parseInt(inpNy.value) | 0;
                const nz = parseInt(inpNz.value) | 0;
                const dedupTol = +inpDedupTol.value;
                if (!(dedupTol > 0)) throw new Error('Dedup tol must be >0');
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
                        data = { mol: CrystalUtils.genReplicatedCellSlab({ lvec, basisPos, basisTypes, nRep: [nx, ny, nz], origin: new Vec3(0, 0, 0), nHat: n, cmin: slab.cmin, cmax: slab.cmax, dedup: chkDedupRep.checked, dedupTol }), lvec };
                    } else {
                        data = { mol: CrystalUtils.genReplicatedCell({ lvec, basisPos, basisTypes, nRep: [nx, ny, nz], origin: new Vec3(0, 0, 0), dedup: chkDedupRep.checked, dedupTol }), lvec };
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
                        data = { mol: CrystalUtils.genReplicatedCellSlab({ lvec, basisPos, basisTypes, nRep: [nx, ny, nz], origin: new Vec3(0, 0, 0), nHat: n, cmin: slab.cmin, cmax: slab.cmax, dedup: chkDedupRep.checked, dedupTol }), lvec };
                    } else {
                        data = { mol: CrystalUtils.genReplicatedCell({ lvec, basisPos, basisTypes, nRep: [nx, ny, nz], origin: new Vec3(0, 0, 0), dedup: chkDedupRep.checked, dedupTol }), lvec };
                    }
                } else if (preset === 'CaCO3(todo)') {
                    throw new Error('CaCO3 preset not implemented yet (non-trivial basis). Use a custom unit cell (e.g. extended XYZ with Lattice=...) once supported.');
                } else {
                    throw new Error(`Preset not implemented: ${preset}`);
                }
                return applyMiller(data);
            };

            const applyToScene = (data, mode) => {
                const app = window.app;
                if (!app) return;
                const sys = app.systems.substrate;
                const rend = app.renderers.substrate;

                if (mode === 'replace') {
                    sys.clear();
                    sys.id2atom.clear();
                    sys.id2bond.clear();
                    sys.lastAtomId = 0;
                    sys.lastBondId = 0;
                    sys.atoms.length = 0;
                    sys.bonds.length = 0;
                }
                const mol = data.mol;
                if (!mol) throw new Error('applyToScene: missing mol');
                if (data.Rmiller) CrystalUtils.rotateMoleculeInPlace(mol, data.Rmiller);
                const idsNew = new Map();
                for (const a of mol.atoms) {
                    const id = sys.addAtom(a.pos.x, a.pos.y, a.pos.z, a.Z);
                    idsNew.set(a.id, id);
                    const ia = sys.getAtomIndex(id);
                    if (ia >= 0 && (a.cellIndex !== undefined && a.cellIndex !== null)) sys.atoms[ia].cellIndex = a.cellIndex | 0;
                }
                if (mol.bonds && mol.bonds.length) {
                    for (const b of mol.bonds) {
                        const aId = idsNew.get(b.aId);
                        const cId = idsNew.get(b.bId);
                        if (aId === undefined || cId === undefined) throw new Error('applyToScene: bond endpoint missing after copy');
                        sys.addBond(aId, cId);
                    }
                }

                if (data.lvec) {
                    sys.lvec = [data.lvec[0].clone(), data.lvec[1].clone(), data.lvec[2].clone()];
                }

                if (data.bucketCells && data.bucketCells.lvecCell) {
                    try {
                        const bg = buildCrystalCellBucketsFromMol(sys, data.bucketCells.na | 0, data.bucketCells.nb | 0, data.bucketCells.nc | 0, data.bucketCells.lvecCell, new Vec3(0, 0, 0));
                        app.lastBucketGraph = bg;
                        if (typeof app.updateBucketOverlay === 'function') app.updateBucketOverlay();
                    } catch (e) {
                        window.logger.error(`Failed to build bucket graph: ${String(e)}`);
                        throw e;
                    }
                }

                rend.update();
                app.requestRender();
                if (app.gui) app.gui.updateSelectionUI();
                window.logger.info(`Substrate generated: atoms=${sys.atoms.length}`);
            };

            btnReplace.onclick = () => { try { applyToScene(buildData(), 'replace'); } catch (e) { window.logger.error(String(e)); throw e; } };
            btnAppend.onclick = () => { try { applyToScene(buildData(), 'append'); } catch (e) { window.logger.error(String(e)); throw e; } };

            const buildFromCIFText = (cifText) => {
                const nx = parseInt(inpNx.value) | 0;
                const ny = parseInt(inpNy.value) | 0;
                const nz = parseInt(inpNz.value) | 0;
                const dedupTol = +inpDedupTol.value;
                if (!(dedupTol > 0)) throw new Error('Dedup tol must be >0');
                if (nx <= 0 || ny <= 0 || nz <= 0) throw new Error('nx,ny,nz must be >0');
                const slab = getSlab();
                const mm = (window.app && window.app.mmParams) ? window.app.mmParams : null;
                if (slab && chkBonds.checked) throw new Error('Build bonds is not supported with slab cut');
                if (chkBonds.checked && !mm) throw new Error('Build bonds requested but window.app.mmParams is not available');
                const data0 = CrystalUtils.genCrystalFromCIF(cifText, { applySymmetry: chkSym.checked, dedupSymmetry: chkDedupSym.checked, dedupReplicate: chkDedupRep.checked, dedupTol, nRep: [nx, ny, nz], origin: new Vec3(0, 0, 0), buildBonds: chkBonds.checked, mmParams: mm, slab });
                return applyMiller(data0);
            };

            const updateCIFDefaults = () => {
                const key = selCIF.value;
                const p = (key && key !== '(none)') ? CIF_PREPARED[key] : null;
                if (p) chkSym.checked = !!p.defaultApplySymmetry;
            };
            updateCIFDefaults();

            selCIF.onchange = async () => {
                try {
                    updateCIFDefaults();
                    const key = selCIF.value;
                    if (!key || key === '(none)') return;
                    const preset = CIF_PREPARED[key];
                    if (!preset || !preset.path) throw new Error('Bad CIF preset');
                    const res = await fetch(preset.path);
                    if (!res.ok) throw new Error(`Failed to fetch ${preset.path} (HTTP ${res.status})`);
                    const txt = await res.text();
                    fillCellTextAreasFromCIF(txt);
                    if (chkAutoSync.checked) updatePreviewFromEditor();
                } catch (e) { window.logger.error(String(e)); throw e; }
            };

            btnFileCIF.onclick = () => inpFileCIF.click();
            inpFileCIF.onchange = async (e) => {
                try {
                    if (e.target.files.length <= 0) return;
                    const file = e.target.files[0];
                    const txt = await gui.readFileText(file);
                    fillCellTextAreasFromCIF(txt);
                    if (chkAutoSync.checked) updatePreviewFromEditor();
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

            mkLbl('Examples (from tests/tAttach)');

            const EXAMPLES = {
                backbones: {
                    'PorhQuad_4SeCl': '../../tests/tAttach/porphironoids/PorhQuad_4SeCl.mol2',
                    'PorhQuad_4SeCl_20deg': '../../tests/tAttach/porphironoids/PorhQuad_4SeCl_20deg.mol2',
                    'PorhQuad_4SeCl_30deg': '../../tests/tAttach/porphironoids/PorhQuad_4SeCl_30deg.mol2',
                    'diphenatroline_CO_subs': '../../tests/tAttach/porphironoids/diphenatroline_CO_subs.mol2',
                    'diphenatroline_NH_subs': '../../tests/tAttach/porphironoids/diphenatroline_NH_subs.mol2',
                    'porphirin_3sites_subs': '../../tests/tAttach/porphironoids/porphirin_3sites_subs.mol2',
                    'porphirin_subs2': '../../tests/tAttach/porphironoids/porphirin_subs2.mol2',
                    'tisite_subs': '../../tests/tAttach/porphironoids/tisite_subs.mol2',
                },
                endgroups: {
                    'guanine-SeCl': { path: '../../tests/tAttach/endgroups/guanine-SeCl.mol2', markers: { X: 'Se', Y: 'Cl' } },
                    'guanine-AlCl': { path: '../../tests/tAttach/endgroups/guanine-AlCl.mol2', markers: { X: 'Al', Y: 'Cl' } },
                    'guanine-SF': { path: '../../tests/tAttach/endgroups/guanine-SF.mol2', markers: { X: 'S', Y: 'F' } },
                },
                monomers: {
                    D: { file: '../../tests/tAttach/backbones/DANA_CG_2inv.mol2', anchors: [1, 4] },
                    T: { file: '../../tests/tAttach/TNA_CG-.mol2', anchors: [6, 5] },
                    P: { file: '../../tests/tAttach/PNA_CG.mol2', anchors: [2, 19] },
                },
                sequences: {
                    'DDDD_DDDD': 'DDDD_DDDD',
                    'PPPPPPP': 'PPPPPPP',
                    'TPTTDD': 'TPTTDD',
                },
                // Default marker pair used only if a preset does not specify its own.
                markers: { X: 'Se', Y: 'Cl' },
            };

            const rowEx0 = GUIutils.row(container);
            GUIutils.span(rowEx0, 'Backbone: ');
            const selBackbone = GUIutils.selectList(rowEx0, Object.keys(EXAMPLES.backbones), null, null, { flexGrow: '1' });
            const btnLoadBackbone = GUIutils.btn(rowEx0, 'Load', null, { marginLeft: '4px', flexGrow: '0' });

            const rowEx1 = GUIutils.row(container, { marginTop: '4px' });
            GUIutils.span(rowEx1, 'Endgroup: ');
            const selEnd = GUIutils.selectList(rowEx1, Object.keys(EXAMPLES.endgroups), null, null, { flexGrow: '1' });
            const btnLoadEnd = GUIutils.btn(rowEx1, 'Load', null, { marginLeft: '4px', flexGrow: '0' });

            const rowEx2 = GUIutils.row(container, { marginTop: '4px' });
            const btnRunMarkerExample = GUIutils.btn(rowEx2, 'Run marker-attach example', null, { flexGrow: '1' });

            const rowEx3 = GUIutils.row(container, { marginTop: '6px' });
            GUIutils.span(rowEx3, 'Seq: ');
            const selSeq = GUIutils.selectList(rowEx3, Object.keys(EXAMPLES.sequences), null, null, { flexGrow: '1' });
            const btnLoadMonoSet = GUIutils.btn(rowEx3, 'Load monomers', null, { marginLeft: '4px', flexGrow: '0' });
            const btnBuildSeq = GUIutils.btn(rowEx3, 'Build', null, { marginLeft: '4px', flexGrow: '0' });

            const rowDbg = GUIutils.row(container, { marginTop: '4px' });
            const btnSaveMol2 = GUIutils.btn(rowDbg, 'Save MOL2', null, { flexGrow: '0' });
            const btnLogBonds = GUIutils.btn(rowDbg, 'Log longest bonds', null, { marginLeft: '4px', flexGrow: '0' });

            const fetchMol2 = async (path) => {
                const res = await fetch(path);
                if (!res.ok) throw new Error(`Failed to fetch ${path} (HTTP ${res.status})`);
                const txt = await res.text();
                return EditableMolecule.parseMol2(txt);
            };

            const applyParsedReplace = (parsed) => {
                gui.system.clear();
                gui.system.appendParsedSystem(parsed);
                gui.renderer.update();
                gui.requestRender();
            };

            btnSaveMol2.onclick = () => {
                try {
                    gui.saveMol2File('molgui_export.mol2');
                } catch (e) { window.logger.error(String(e)); throw e; }
            };

            btnLogBonds.onclick = () => {
                try {
                    const mol = gui.system;
                    const out = [];
                    for (let ib = 0; ib < mol.bonds.length; ib++) {
                        const b = mol.bonds[ib];
                        b.ensureIndices(mol);
                        const a = mol.atoms[b.a];
                        const c = mol.atoms[b.b];
                        const dx = a.pos.x - c.pos.x;
                        const dy = a.pos.y - c.pos.y;
                        const dz = a.pos.z - c.pos.z;
                        const r = Math.sqrt(dx * dx + dy * dy + dz * dz);
                        out.push({ r, ib, aId: a.id, bId: c.id, aZ: a.Z, bZ: c.Z, ax: a.pos.x, ay: a.pos.y, az: a.pos.z, bx: c.pos.x, by: c.pos.y, bz: c.pos.z });
                    }
                    out.sort((u, v) => v.r - u.r);
                    const n = Math.min(20, out.length);
                    window.logger.info(`Longest bonds (top ${n} of ${out.length}):`);
                    for (let i = 0; i < n; i++) {
                        const e = out[i];
                        window.logger.info(`#${i} ib=${e.ib} r=${e.r.toFixed(3)} aId=${e.aId}(Z=${e.aZ}) bId=${e.bId}(Z=${e.bZ}) A=(${e.ax.toFixed(3)},${e.ay.toFixed(3)},${e.az.toFixed(3)}) B=(${e.bx.toFixed(3)},${e.by.toFixed(3)},${e.bz.toFixed(3)})`);
                    }
                } catch (e) { window.logger.error(String(e)); throw e; }
            };

            btnLoadBackbone.onclick = async () => {
                try {
                    const key = selBackbone.value;
                    const path = EXAMPLES.backbones[key];
                    if (!path) throw new Error('Backbone preset missing path');
                    const parsed = await fetchMol2(path);
                    applyParsedReplace(parsed);
                    window.logger.info(`Loaded backbone preset '${key}' atoms=${parsed.types.length}`);
                } catch (e) { window.logger.error(String(e)); throw e; }
            };

            btnLoadEnd.onclick = async () => {
                try {
                    const key = selEnd.value;
                    const rec = EXAMPLES.endgroups[key];
                    const path = (typeof rec === 'string') ? rec : (rec ? rec.path : null);
                    if (!path) throw new Error('Endgroup preset missing path');
                    gui.endGroupParsed = await fetchMol2(path);
                    window.logger.info(`Loaded endgroup preset '${key}' atoms=${gui.endGroupParsed.types.length}`);
                } catch (e) { window.logger.error(String(e)); throw e; }
            };

            btnRunMarkerExample.onclick = async () => {
                try {
                    const backKey = selBackbone.value;
                    const endKey = selEnd.value;
                    const backPath = EXAMPLES.backbones[backKey];
                    const endRec = EXAMPLES.endgroups[endKey];
                    const endPath = (typeof endRec === 'string') ? endRec : (endRec ? endRec.path : null);
                    if (!backPath || !endPath) throw new Error('Example selection missing path');
                    const backbone = await fetchMol2(backPath);
                    const endg = await fetchMol2(endPath);
                    applyParsedReplace(backbone);
                    gui.endGroupParsed = endg;
                    const mkBack = EXAMPLES.markers;
                    const mkGroup = (endRec && typeof endRec === 'object' && endRec.markers) ? endRec.markers : EXAMPLES.markers;
                    if (!mkBack || !mkBack.X || !mkBack.Y) throw new Error('Example preset missing backbone markers');
                    if (!mkGroup || !mkGroup.X || !mkGroup.Y) throw new Error('Example preset missing endgroup markers');
                    try {
                        const same = (mkBack.X === mkGroup.X) && (mkBack.Y === mkGroup.Y);
                        if (same) gui.system.attachGroupByMarker(gui.endGroupParsed, mkBack.X, mkBack.Y);
                        else gui.system.attachGroupByMarker(gui.endGroupParsed, mkBack.X, mkBack.Y, { groupMarkerX: mkGroup.X, groupMarkerY: mkGroup.Y });
                    } catch (err) {
                        const msg = String(err);
                        if (msg.includes('group must have exactly one marker pair, got 0')) {
                            throw new Error(`Marker attach failed: endgroup '${endKey}' does not contain marker pair X='${mkGroup.X}' Y='${mkGroup.Y}'`);
                        }
                        throw err;
                    }
                    gui.renderer.update();
                    gui.requestRender();
                    window.logger.info(`Marker-attach example done: backbone='${backKey}' endgroup='${endKey}'`);
                } catch (e) { window.logger.error(String(e)); throw e; }
            };

            btnLoadMonoSet.onclick = async () => {
                try {
                    gui.monomerLib = { name: 'tAttach polymerize.py', monomers: [] };
                    gui.monomerGeoms.clear();
                    for (const [id, rec] of Object.entries(EXAMPLES.monomers)) {
                        const parsed = await fetchMol2(rec.file);
                        const fname = rec.file.split('/').slice(-1)[0];
                        gui.monomerGeoms.set(fname, parsed);
                        gui.monomerLib.monomers.push({
                            id,
                            file: fname,
                            anchors: { head: { type: 'index', value: rec.anchors[0] }, tail: { type: 'index', value: rec.anchors[1] } }
                        });
                    }
                    window.logger.info(`Loaded monomer preset set: monomers=${gui.monomerLib.monomers.length} mol2=${gui.monomerGeoms.size}`);
                } catch (e) { window.logger.error(String(e)); throw e; }
            };

            btnBuildSeq.onclick = async () => {
                try {
                    const seq = EXAMPLES.sequences[selSeq.value];
                    if (!seq) throw new Error('Sequence preset missing');

                    if (!gui.monomerLib) await btnLoadMonoSet.onclick();

                    const tokens = gui.parseSequenceTokens(seq);
                    const monomers = (() => {
                        if (!gui.monomerLib) throw new Error('No monomer library loaded');
                        const map = {};
                        for (const m of gui.monomerLib.monomers) {
                            const parsed = gui.monomerGeoms.get(m.file);
                            if (!parsed) throw new Error(`Missing geometry for '${m.id}' file='${m.file}'`);
                            map[m.id] = { parsed, anchors: [m.anchors.head.value, m.anchors.tail.value] };
                        }
                        return map;
                    })();

                    const poly = PolymerUtils.assemblePolymerFromTokens(tokens, monomers, { _0: 1 });
                    gui.system.clear();
                    const idsNew = new Map();
                    for (const a of poly.atoms) {
                        const id = gui.system.addAtom(a.pos.x, a.pos.y, a.pos.z, a.Z);
                        idsNew.set(a.id, id);
                    }
                    for (const b of poly.bonds) {
                        const aId = idsNew.get(b.aId);
                        const cId = idsNew.get(b.bId);
                        if (aId === undefined || cId === undefined) throw new Error('btnBuildSeq: bond endpoint missing after copy');
                        gui.system.addBond(aId, cId);
                    }
                    gui.renderer.update();
                    gui.requestRender();
                    window.logger.info(`Built polymer preset sequence '${seq}' atoms=${gui.system.nAtoms}`);
                } catch (e) { window.logger.error(String(e)); throw e; }
            };

            GUIutils.el(container, 'hr', null, { borderColor: '#444', margin: '8px 0' });

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

            const btnLoadEndFile = GUIutils.btn(container, 'Load Endgroup mol2...');

            const inpEnd = GUIutils.input(container, { type: 'file', attrs: { accept: '.mol2' } }, { display: 'none' });

            const lblEnd = GUIutils.div(container, null, { fontSize: '0.85em', color: '#aaa', marginTop: '2px' });
            lblEnd.textContent = 'No endgroup loaded';

            btnLoadEndFile.onclick = () => inpEnd.click();
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
