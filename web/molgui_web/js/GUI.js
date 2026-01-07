import * as CrystalUtils from './CrystalUtils.js';
import * as PolymerUtils from './MoleculeUtils.js';
import { Vec3 } from '../../common_js/Vec3.js';
import { GUIutils } from '../../common_js/GUIutils.js';
import { EditableMolecule } from './EditableMolecule.js';
import { installMoleculeIOMethods } from './MoleculeIO.js';
import { installMoleculeUtilsMethods } from './MoleculeUtils.js';
import { installMoleculeSelectionMethods } from './MoleculeSelection.js';
import { BuildersGUI } from './BuildersGUI.js';

// Ensure IO and utils methods are installed on EditableMolecule before any static use.
installMoleculeIOMethods(EditableMolecule);
installMoleculeUtilsMethods(EditableMolecule);
installMoleculeSelectionMethods(EditableMolecule);

export class GUI {
    constructor(system, renderer) {
        this.system = system;
        this.renderer = renderer;
        this.monomerLib = null;
        this.monomerGeoms = new Map();
        this.endGroupParsed = null;
        this.buildersGUI = new BuildersGUI(this);
        this.init();
    }

    async readFileText(file) {
        if (!file) throw new Error('readFileText: file is null');
        return await new Promise((resolve, reject) => {
            const reader = new FileReader();
            reader.onload = (ev) => resolve(ev.target.result);
            reader.onerror = (e) => reject(e);
            reader.readAsText(file);
        });
    }

    requestRender() {
        if (window.app && window.app.requestRender) window.app.requestRender();
    }

    parseSequenceTokens(seq) {
        const s = (seq !== undefined && seq !== null) ? String(seq) : '';
        const tokens = [];
        let i = 0;
        while (i < s.length) {
            const c = s[i];
            if (c === '_' || c === '-' || c === ' ' || c === '\t' || c === '\n' || c === '\r') { i++; continue; }
            if (!(c >= 'A' && c <= 'Z')) throw new Error(`parseSequenceTokens: expected token starting with [A-Z] at i=${i}, got '${c}'`);
            let j = i + 1;
            while (j < s.length) {
                const cj = s[j];
                const isAlpha = ((cj >= 'A' && cj <= 'Z') || (cj >= 'a' && cj <= 'z'));
                if (!isAlpha) break;
                j++;
            }
            const tok0 = s.slice(i, j);
            i = j;
            let k = i;
            while (k < s.length) {
                const ck = s[k];
                if (!(ck >= '0' && ck <= '9')) break;
                k++;
            }
            let n = 1;
            if (k > i) {
                n = parseInt(s.slice(i, k), 10) | 0;
                if (!(n > 0)) throw new Error(`parseSequenceTokens: invalid repeat count for token '${tok0}'`);
                i = k;
            }

            // If token is ALL-UPPERCASE and has no explicit count, interpret it as
            // a sequence of single-letter tokens (so DDDD_DDDD works).
            // If it HAS a count, keep it intact (so PNA10 works).
            let bAllUpper = true;
            for (let ii = 0; ii < tok0.length; ii++) {
                const ch = tok0[ii];
                if (!(ch >= 'A' && ch <= 'Z')) { bAllUpper = false; break; }
            }
            if (bAllUpper && n === 1 && tok0.length > 1) {
                for (let ii = 0; ii < tok0.length; ii++) tokens.push(tok0[ii]);
            } else {
                for (let rep = 0; rep < n; rep++) tokens.push(tok0);
            }
        }
        return tokens;
    }

    saveXYZFile(filename = 'molecule.xyz') {
        const content = this.system.toXYZString();
        const blob = new Blob([content], { type: 'text/plain' });
        const url = URL.createObjectURL(blob);
        const a = document.createElement('a');
        a.href = url;
        a.download = filename;
        a.click();
        URL.revokeObjectURL(url);
        window.logger.info(`Saved ${filename}`);
    }

    saveMol2File(filename = 'molecule.mol2') {
        const content = this.system.toMol2String({ name: filename.replace(/\.mol2$/i, '') });
        const blob = new Blob([content], { type: 'text/plain' });
        const url = URL.createObjectURL(blob);
        const a = document.createElement('a');
        a.href = url;
        a.download = filename;
        a.click();
        URL.revokeObjectURL(url);
        window.logger.info(`Saved ${filename}`);
    }

    loadXYZString(content) {
        const parsed = EditableMolecule.parseXYZ(content);
        this.system.clear();
        this.system.appendParsedSystem(parsed);
        const mm = (window.app && window.app.mmParams) ? window.app.mmParams : null;
        this.system.recalculateBonds(mm);
        this.renderer.update();
        this.requestRender();
        window.logger.info(`Loaded XYZ: atoms=${parsed.types.length}`);
    }

    init() {
        // Create Sidebar
        const sidebar = document.createElement('div');
        sidebar.id = 'gui-sidebar';
        document.body.appendChild(sidebar);

        // Create Resizer
        const resizer = document.createElement('div');
        resizer.id = 'gui-resizer';
        document.body.appendChild(resizer);

        // Resizer Logic
        let isResizing = false;
        resizer.addEventListener('mousedown', (e) => {
            isResizing = true;
            document.body.style.cursor = 'col-resize';
            e.preventDefault(); // Prevent selection
        });

        window.addEventListener('mousemove', (e) => {
            if (!isResizing) return;
            const newWidth = e.clientX;
            if (newWidth < 150 || newWidth > 600) return; // Limits

            sidebar.style.width = `${newWidth}px`;
            resizer.style.left = `${newWidth}px`;

            const canvasContainer = document.getElementById('canvas-container');
            if (canvasContainer) {
                const resizerWidth = 5;
                canvasContainer.style.left = `${newWidth + resizerWidth}px`;
                canvasContainer.style.width = `calc(100vw - ${newWidth + resizerWidth}px)`;
            }

            // Trigger Resize for Renderer
            if (window.app) {
                window.app.onWindowResize();
            }
        });

        window.addEventListener('mouseup', () => {
            if (isResizing) {
                isResizing = false;
                document.body.style.cursor = 'default';
            }
        });

        // --- Section: Selection ---
        this.createSection(sidebar, 'Selection', (container) => {
            const row = GUIutils.row(container);
            GUIutils.span(row, 'Count: ', { fontSize: '0.9em', marginRight: '5px' });
            this.lblCount = GUIutils.span(row, '0', { fontWeight: 'bold' });
            this.inpSelection = GUIutils.textInput(container, '', { placeholder: 'IDs (e.g. 1,5)', onchange: (e) => this.onSelectionInputChange(e.target.value) });

            const rowSelAll = GUIutils.row(container, { marginTop: '4px' });
            GUIutils.btn(rowSelAll, 'Select All', () => {
                this.system.selectAll();
                this.updateSelectionUI();
                this.renderer.update();
                if (window.app && window.app.editor) window.app.editor.updateGizmo();
                this.requestRender();
            }, { flexGrow: '1' });

            GUIutils.el(container, 'hr', null, { borderColor: '#444', margin: '10px 0' });

            GUIutils.div(container, null, { fontSize: '0.9em', marginBottom: '5px' }).textContent = 'Selection Query:';
            this.inpSelQuery = GUIutils.textInput(container, '', { placeholder: 'e.g. N|C n{F|Br|Cl}={1,2}' });

            const rowQEx = GUIutils.row(container, { marginTop: '5px' });
            GUIutils.span(rowQEx, 'Examples: ', { fontSize: '0.9em', marginRight: '5px' });
            const selQEx = GUIutils.el(rowQEx, 'select', { className: 'gui-select' }, { flexGrow: '1' });
            const qExamples = [
                { value: '', text: '-- pick example --' },
                { value: 'C', text: 'All carbon (C)' },
                { value: 'O', text: 'All oxygen (O)' },
                { value: 'H', text: 'All hydrogen (H)' },
                { value: 'C|N|O', text: 'Backbone (C,N,O)' },
                { value: 'C deg={4}', text: 'sp3-like carbon: C with degree 4' },
                { value: 'N deg={3}', text: 'Nitrogen degree 3' },
                { value: 'O deg={1}', text: 'Terminal oxygen: O degree 1' },
                { value: 'C n{H}={3}', text: 'Methyl carbon: C with 3 H neighbors' },
                { value: 'N n{H}={2}', text: 'Amine-like N: N with 2 H neighbors' },
                { value: 'N|O n{C}={1}', text: 'Hetero attached to one carbon' },
                { value: 'C n{F|Cl|Br|I}={1,2,3}', text: 'Halogenated carbon (1-3 halogens)' }
            ];
            GUIutils.setSelectOptions(selQEx, qExamples, { selectedValue: '', selectFirst: true });
            selQEx.onchange = (e) => {
                const v = (e && e.target) ? String(e.target.value) : '';
                if (!v) return;
                if (this.inpSelQuery) this.inpSelQuery.value = v;
            };

            const rowQBtns = GUIutils.row(container, { marginTop: '5px' });
            const applyQ = (mode) => {
                if (!(window.app && window.app.mmParams)) throw new Error('Selection Query: window.app.mmParams is missing');
                const q = (this.inpSelQuery && this.inpSelQuery.value) ? String(this.inpSelQuery.value).trim() : '';
                if (!q) throw new Error('Selection Query: query is empty');
                const cq = EditableMolecule.compileSelectQuery(q, window.app.mmParams);
                this.system.applySelectQuery(cq, { mode });
                this.updateSelectionUI();
                if (window.app && window.app.editor) window.app.editor.updateGizmo();
                if (window.app && window.app.molRenderer) window.app.molRenderer.updateSelection();
                this.renderer.update();
                this.requestRender();
            };
            GUIutils.btn(rowQBtns, 'Replace', () => applyQ('replace'), { marginRight: '4px' });
            GUIutils.btn(rowQBtns, 'Add', () => applyQ('add'), { marginRight: '4px' });
            GUIutils.btn(rowQBtns, 'Subtract', () => applyQ('subtract') );
        });

        // --- Section: View ---
        this.createSection(sidebar, 'View', (container) => {
            // Zoom Slider
            const zoomRow = GUIutils.row(container, { flexDirection: 'column', alignItems: 'flex-start' });
            GUIutils.el(zoomRow, 'label', { text: 'Zoom Level (Log10)' }, { fontSize: '0.9em' });
            GUIutils.range(zoomRow, 1.0, -2.0, 3.0, 0.1, (e) => { const val = parseFloat(e.target.value); this.setZoom(Math.pow(10, val)); }, { width: '100%' });

            // View Buttons
            const viewRow = GUIutils.row(container, { marginTop: '10px' });

            const views = [
                { name: 'XY', pos: [0, 0, 20], up: [0, 1, 0] },
                { name: 'XZ', pos: [0, 20, 0], up: [0, 0, -1] },
                { name: 'YZ', pos: [20, 0, 0], up: [0, 1, 0] }
            ];

            views.forEach(v => {
                GUIutils.btn(viewRow, v.name, () => this.setView(v.pos, v.up), { marginRight: '2px' });
            });

            // Axis Toggle
            const axisRow = GUIutils.row(container, { marginTop: '10px' });
            GUIutils.labelCheck(axisRow, 'Show Axes', true, (e) => { if (window.app && window.app.molRenderer) window.app.molRenderer.toggleAxes(e.target.checked); this.requestRender(); });

            const renderRow = GUIutils.row(container, { marginTop: '6px' });
            GUIutils.labelCheck(renderRow, 'Continuous Render (Animate)', false, (e) => {
                if (window.app && window.app.setContinuousRender) window.app.setContinuousRender(e.target.checked);
                this.requestRender();
            });

            // Make on-demand the explicit default (checkbox unchecked)
            if (window.app && window.app.setContinuousRender) window.app.setContinuousRender(false);

            // Label Mode Dropdown
            const labelRow = GUIutils.row(container, { marginTop: '10px' });
            GUIutils.span(labelRow, 'Labels: ');
            const labelSel = GUIutils.el(labelRow, 'select', { className: 'gui-select' }, { flexGrow: '1' });

            const modes = [
                { value: 'none',    text: 'None'      },
                { value: 'id',      text: 'Atom ID'   },
                { value: 'element', text: 'Element'   },
                { value: 'type',    text: 'Atom Type' }
            ];

            GUIutils.setSelectOptions(labelSel, modes, { selectedValue: 'none', selectFirst: true });

            labelSel.onchange = (e) => {
                if (window.app && window.app.molRenderer) {
                    window.app.molRenderer.setLabelMode(e.target.value);
                }
                this.requestRender();
            };

            // Label Color & Size
            const styleRow = GUIutils.row(container, { marginTop: '5px' });
            GUIutils.span(styleRow, 'Color: ');
            const colorInput = GUIutils.input(styleRow, { type: 'color', value: '#ffffff' }, { flexGrow: '0', width: '40px' });
            GUIutils.span(styleRow, ' Size: ', { marginLeft: '10px' });
            const sizeInput = GUIutils.num(styleRow, 0.5, { step: '0.1', min: '0.1' }, { width: '50px' });

            const updateLabelStyle = () => {
                if (window.app && window.app.molRenderer) {
                    window.app.molRenderer.setLabelStyle(colorInput.value, sizeInput.value);
                }
                this.requestRender();
            };

            colorInput.oninput = updateLabelStyle;
            sizeInput.oninput = updateLabelStyle;

        });

        // --- Section: Gizmo ---
        this.createSection(sidebar, 'Gizmo', (container) => {
            // Enable Checkbox
            GUIutils.labelCheck(container, 'Enable Gizmo', true, (e) => {
                if (window.app && window.app.editor) window.app.editor.toggleGizmo(e.target.checked);
            });

            // Lock Selection Checkbox
            GUIutils.labelCheck(container, 'Lock Selection', false, (e) => {
                if (window.app && window.app.editor) {
                    window.app.editor.selectionLocked = e.target.checked;
                    window.logger.info(`Selection ${e.target.checked ? 'Locked' : 'Unlocked'}`);
                }
            }, { marginTop: '5px' });

            // Modes
            const modes = ['translate', 'rotate', 'scale'];
            const modeContainer = GUIutils.div(container, null, { marginTop: '10px' });

            modes.forEach(mode => {
                const mLabel = GUIutils.el(modeContainer, 'label', { className: 'gui-checkbox-label' }, { marginBottom: '5px' });
                const radio = GUIutils.el(mLabel, 'input', { type: 'radio', attrs: { name: 'gizmo-mode' }, value: mode });
                radio.checked = (mode === 'translate');
                radio.onchange = () => { if (window.app && window.app.editor) window.app.editor.setGizmoMode(mode); };
                GUIutils.span(mLabel, mode.charAt(0).toUpperCase() + mode.slice(1));
            });
        });

        // --- Section: Structure ---
        this.createSection(sidebar, 'Structure', (container) => {
            // Recalculate Bonds
            GUIutils.btn(container, 'Recalculate Bonds', () => { if (window.app && window.app.editor) window.app.editor.recalculateBonds(); });

            const rowBondMode = GUIutils.row(container, { marginTop: '6px' });
            GUIutils.span(rowBondMode, 'Bond mode: ', { marginRight: '4px', fontSize: '0.9em' });
            const selBondMode = GUIutils.selectList(rowBondMode, ['brute', 'bucketNeighbors', 'bucketAllPairsAABB'], null, null, { flexGrow: '1' });
            selBondMode.value = (window.app && window.app.bondRecalcMode) ? String(window.app.bondRecalcMode) : 'brute';
            selBondMode.onchange = () => {
                if (window.app) window.app.bondRecalcMode = String(selBondMode.value);
            };

            const rowBuckets = GUIutils.row(container, { marginTop: '4px' });
            const chkBuckets = GUIutils.labelCheck(rowBuckets, 'Show buckets', false, null).input;
            chkBuckets.checked = !!(window.app && window.app.showBucketBoxes);
            chkBuckets.onchange = () => {
                if (window.app) {
                    window.app.showBucketBoxes = !!chkBuckets.checked;
                    if (typeof window.app.updateBucketOverlay === 'function') window.app.updateBucketOverlay();
                }
            };

            const rowBucketAuto = GUIutils.row(container, { marginTop: '4px' });
            const chkBucketAuto = GUIutils.labelCheck(rowBucketAuto, 'Auto-update buckets', true, null).input;
            chkBucketAuto.checked = (window.app && window.app.autoUpdateBuckets !== undefined) ? !!window.app.autoUpdateBuckets : true;
            chkBucketAuto.onchange = () => {
                if (window.app) window.app.autoUpdateBuckets = !!chkBucketAuto.checked;
            };

            const rowBucketLines = GUIutils.row(container, { marginTop: '4px' });
            const chkBucketLines = GUIutils.labelCheck(rowBucketLines, 'Show atom→bucket lines', false, null).input;
            chkBucketLines.checked = !!(window.app && window.app.showBucketAtomLines);
            chkBucketLines.onchange = () => {
                if (window.app) {
                    window.app.showBucketAtomLines = !!chkBucketLines.checked;
                    if (typeof window.app.refreshBucketDebug === 'function') window.app.refreshBucketDebug();
                    else if (typeof window.app.updateBucketOverlay === 'function') window.app.updateBucketOverlay();
                }
            };

            GUIutils.el(container, 'hr', null, { borderColor: '#444', margin: '10px 0' });

            // Passivation
            GUIutils.div(container, null, { fontSize: '0.9em', marginBottom: '5px' }).textContent = 'Passivation:';
            const rowCapType = GUIutils.row(container);
            GUIutils.span(rowCapType, 'Cap type: ', { marginRight: '4px', fontSize: '0.9em' });
            this.inpCapType = GUIutils.textInput(rowCapType, 'H', { placeholder: 'H or atom type (e.g. C_3)' });
            this.inpCapType.style.flexGrow = '1';

            const rowPassivBtns = GUIutils.row(container, { marginTop: '5px' });
            GUIutils.btn(rowPassivBtns, 'Add Caps', () => {
                if (!(window.app && window.app.mmParams)) throw new Error('Add Caps: window.app.mmParams is missing');
                const onlySelection = (this.system && this.system.selection && this.system.selection.size > 0);
                const cap = (this.inpCapType && this.inpCapType.value) ? String(this.inpCapType.value).trim() : 'H';
                if (!cap) throw new Error('Add Caps: cap type is empty');
                try {
                    this.system.addCappingAtoms(window.app.mmParams, cap, { onlySelection, bBond: true });
                } catch (e) {
                    if (window.logger && window.logger.error) window.logger.error(String(e && e.stack ? e.stack : e));
                    throw e;
                }
                this.renderer.update();
                this.requestRender();
            }, { marginRight: '4px' });

            GUIutils.btn(rowPassivBtns, 'Add EPairs', () => {
                if (!(window.app && window.app.mmParams)) throw new Error('Add EPairs: window.app.mmParams is missing');
                const onlySelection = (this.system && this.system.selection && this.system.selection.size > 0);
                try {
                    this.system.addExplicitEPairs(window.app.mmParams, { onlySelection, bBond: true });
                } catch (e) {
                    if (window.logger && window.logger.error) window.logger.error(String(e && e.stack ? e.stack : e));
                    throw e;
                }
                this.renderer.update();
                this.requestRender();
            });

            GUIutils.el(container, 'hr', null, { borderColor: '#444', margin: '10px 0' });

            // Add Atom Controls
            GUIutils.div(container, null, { fontSize: '0.9em', marginBottom: '5px' }).textContent = 'Add Atom Settings:';

            // Element Dropdown
            const rowEl = GUIutils.row(container);
            GUIutils.span(rowEl, 'Element: ');
            this.selElement = GUIutils.el(rowEl, 'select', { className: 'gui-select', onchange: (e) => this.onElementChange(e.target.value) }, { flexGrow: '1' });

            // Atom Type Dropdown
            const rowType = GUIutils.row(container, { marginTop: '5px' });
            GUIutils.span(rowType, 'Type: ');
            this.selAtomType = GUIutils.el(rowType, 'select', { className: 'gui-select', onchange: (e) => this.onAtomTypeChange(e.target.value) }, { flexGrow: '1' });

            // Populate initially (delayed to ensure MMParams loaded)
            setTimeout(() => this.populateStructureControls(), 1000);
        });

        // --- Section: Geometry ---
        this.createSection(sidebar, 'Geometry', (container) => {
            // Load
            const btnLoad = GUIutils.btn(container, 'Load XYZ File...');
            const fileInput = GUIutils.input(document.body, { type: 'file', attrs: { accept: '.xyz' } }, { display: 'none' });
            btnLoad.onclick = () => fileInput.click();
            fileInput.onchange = (e) => {
                if (e.target.files.length > 0) {
                    const file = e.target.files[0];
                    const reader = new FileReader();
                    reader.onload = (ev) => this.loadXYZString(ev.target.result);
                    reader.readAsText(file);
                    fileInput.value = '';
                }
            };

            // Save
            GUIutils.btn(container, 'Save XYZ', () => this.saveXYZFile(), { marginTop: '5px' });

            // Clear
            GUIutils.btn(container, 'Clear Scene', () => {
                this.system.clear();
                this.renderer.update();
                this.requestRender();
                window.logger.info("Scene cleared.");
            }, { marginTop: '5px' });

            GUIutils.el(container, 'hr', null, { borderColor: '#444', margin: '10px 0' });

            // Examples (mol/mol2 preferred, xyz fallback)
            const examples = [
                { name: 'H2O',             path: '../../cpp/common_resources/mol/H2O.mol2', fmt: 'mol2' },
                { name: 'NH3',             path: '../../cpp/common_resources/mol/NH3.mol2', fmt: 'mol2' },
                { name: 'CH4',             path: '../../cpp/common_resources/xyz/CH4.xyz',  fmt: 'xyz' },
                { name: 'Formaldehyde',    path: '../../cpp/common_resources/mol/formaldehyde.mol2', fmt: 'mol2' },
                { name: 'Formic Acid',     path: '../../cpp/common_resources/mol/formic_acid.mol2', fmt: 'mol2' },
                { name: 'Methanol',        path: '../../cpp/common_resources/mol/methanol.mol2',   fmt: 'mol2' },
                { name: 'HCN',             path: '../../cpp/common_resources/mol/HCN.mol2',        fmt: 'mol2' },
                { name: 'Benzene',         path: '../../cpp/common_resources/xyz/benzene.mol2',    fmt: 'mol2' },
                { name: 'Pyridine',        path: '../../cpp/common_resources/xyz/pyridine.xyz',    fmt: 'xyz' },
                { name: 'Pentacene',       path: '../../cpp/common_resources/xyz/pentacene.xyz',   fmt: 'xyz' },
                { name: 'PTCDA',           path: '../../cpp/common_resources/xyz/PTCDA.xyz',       fmt: 'xyz' },
                { name: 'Porphirin',       path: '../../cpp/common_resources/mol/porphirin.mol2',  fmt: 'mol2' },
                { name: 'Thymine',         path: '../../cpp/common_resources/mol/thymine.mol2',    fmt: 'mol2' },
                { name: 'Adenine',         path: '../../cpp/common_resources/xyz/adenine.xyz',     fmt: 'xyz' },
                { name: 'Cytosine',        path: '../../cpp/common_resources/xyz/citosine.xyz',    fmt: 'xyz' },
                { name: 'Uracil',          path: '../../cpp/common_resources/xyz/uracil.xyz',      fmt: 'xyz' },
                { name: 'Adamantane',      path: '../../cpp/common_resources/mol/adamantane.mol2', fmt: 'mol2' },
                { name: 'Si10',            path: '../../cpp/common_resources/mol/Si10.mol2',       fmt: 'mol2' },
                { name: 'Xylitol',         path: '../../cpp/common_resources/mol/xylitol.mol2',    fmt: 'mol2' },
                { name: 'NaCl 1x1 L3',     path: '../../cpp/common_resources/xyz/NaCl_1x1_L3.xyz', fmt: 'xyz' },

                // agregates
                { name: 'Adenine-Thymine pair',  path: '../../cpp/common_resources/xyz/adenine-thymine.xyz', fmt: 'xyz' },
                { name: 'Guanine-Cytosine pair', path: '../../cpp/common_resources/xyz/guanine-cytosine.xyz', fmt: 'xyz' },
            ];

            GUIutils.div(container, null, { fontSize: '0.9em', marginBottom: '5px' }).textContent = 'Examples:';
            const rowEx = GUIutils.row(container);
            const selEx = GUIutils.el(rowEx, 'select', { className: 'gui-select' }, { flexGrow: '1' });
            GUIutils.setSelectOptions(selEx, [{ value: '', text: '-- pick example --', selected: true }, ...examples.map((e, i) => ({ value: String(i), text: e.name }))], { selectFirst: true });
            GUIutils.btn(rowEx, 'Add', async () => {
                const idx = parseInt(selEx.value);
                if (isNaN(idx) || idx < 0 || idx >= examples.length) return;
                const ex = examples[idx];
                try {
                    const resp = await fetch(ex.path);
                    if (!resp.ok) throw new Error(`Fetch failed: ${resp.status} ${resp.statusText}`);
                    const txt = await resp.text();
                    const parsed = (ex.fmt === 'mol' || ex.fmt === 'mol2') ? EditableMolecule.parseMol2(txt) : EditableMolecule.parseXYZ(txt);
                    const n0 = this.system.atoms.length;
                    this.system.appendParsedSystem(parsed);
                    // Auto-bond XYZ imports (mol/mol2 already have bonds)
                    if (ex.fmt === 'xyz') {
                        const mm = (window.app && window.app.mmParams) ? window.app.mmParams : null;
                        this.system.recalculateBonds(mm);
                    }
                    const n1 = this.system.atoms.length;
                    this.system.clearSelection();
                    for (let i = n0; i < n1; i++) this.system.select(this.system.atoms[i].id, 'add');
                    this.renderer.update();
                    if (window.app && window.app.editor) {
                        window.app.editor.updateGizmo();
                        window.app.editor.toggleGizmo(true);
                    }
                    this.requestRender();
                    window.logger.info(`Added example '${ex.name}' atoms=${n1 - n0}`);
                } catch (err) {
                    window.logger.error(String(err));
                }
            }, { marginLeft: '5px' });

            // Separator
            GUIutils.el(container, 'hr', null, { borderColor: '#444', margin: '10px 0' });
        });

        this.createSection(sidebar, 'XYZ (Edit)', (container) => {
            const btnToggle = GUIutils.btn(container, 'Show/Hide XYZ');
            const textArea = GUIutils.textArea(container, '', { height: '150px', placeholder: 'Paste XYZ content here...' });

            btnToggle.onclick = () => {
                const isHidden = textArea.style.display === 'none';
                textArea.style.display = isHidden ? 'block' : 'none';
                if (isHidden) {
                    // Populate with current system state if showing
                    textArea.value = this.system.toXYZString();
                }
            };

            GUIutils.btn(container, 'Apply XYZ', () => { if (textArea.value.trim()) this.loadXYZString(textArea.value); }, { marginTop: '5px' });
        });

        this.buildersGUI.addSubstrateSection(sidebar);
        this.buildersGUI.addPolymersSection(sidebar);

        // --- Section: User Script ---
        this.createSection(sidebar, 'User Script', (container) => {
            GUIutils.div(container, null, { fontSize: '0.9em', marginBottom: '4px' }).textContent = 'JavaScript (async supported):';
            const defaultScript = `
// PTCDA on NaCl step
molecule.clear();
substrate.clear();
substrate.build_substrate('NaCl', { size:[13,12,3], step_edge:[0,0,1] });
mol.load("../../cpp/common_resources/xyz/PTCDA.xyz");
mol.move([0, 0, 10]);
mol.move([12, 0, 0]);
mol.rotate([1, 0, 0], 60.0);
mol.rotate([0, 0, 1], 30.0);
molecule.addLvec({lvec: [ [18.0, 9.0, 0.0], [0.0, 16.0, 0.0], [0.0, 0.0, 18.0] ]});
molecule.setViewReplicas({ show:true, showBox:true, nrep:[1,1,0] });
`.trim();

            const sampleScripts = [
                { id: 'blank',         name: '-- custom --', code: '' },
                { id: 'ptcda_nacl',    name: 'PTCDA on NaCl step', code: defaultScript },
                { id: 'benzene_clean', name: 'Benzene on clean NaCl', code: `
// Benzene on clean NaCl
molecule.clear();
substrate.clear();
substrate.build_substrate('NaCl', { size:[10,10,2] });
mol.load("../../cpp/common_resources/mol/benzene.mol2");
mol.move([0, 0, 8]);
mol.rotate([0, 0, 1], 30.0);
`.trim() }
                ,
                { id: 'adamantane_collapse', name: 'Collapse bridge in adamantane', code: `
// Collapse one CH2 bridge in adamantane using modular helpers
molecule.clear();
molecule.load("../../cpp/common_resources/mol/adamantane.mol2");

// Select bridge carbons (heavy-heavy with optional H2). Prefer CH2, else degree=2 heavy.
let nsel = select_bridge_candidates({ requireH2: true });
if (nsel === 0) nsel = select_bridge_candidates({ requireH2: false });
if (nsel === 0) throw new Error("No bridge atoms found for collapse");

// Pick one from selection (random by default) and collapse it
const picked = pick_random_selection();
collapse_bridge_at({ id: picked });
`.trim() }
                ,
                { id: 'adamantane_collapse_all', name: 'Collapse ALL bridges in adamantane', code: `
// Collapse all CH2-like bridges in adamantane
molecule.clear();
molecule.load("../../cpp/common_resources/mol/adamantane.mol2");
const nCollapsed = collapse_all_bridges();
logger.info("Collapsed bridges: " + nCollapsed);
`.trim() }
                ,
                { id: 'adamantane_insert_bridge', name: 'Insert bridge between CH carbons', code: `
// Insert a CH2 bridge into adamantane between two CH carbons
molecule.clear();
molecule.load("../../cpp/common_resources/mol/tetrahedrane.mol2");

// Pick a random bond between carbons with >=3 heavy neighbors and >=1 hydrogen (approx CH)
const newC = insert_bridge_random({ minHeavy: 3, minHyd: 1 });
logger.info("Inserted bridge carbon id=" + newC);
`.trim() }
            ];

            const rowSel = GUIutils.row(container, { marginBottom: '5px' });
            GUIutils.span(rowSel, 'Sample: ', { marginRight: '4px' });
            const selSamples = GUIutils.el(rowSel, 'select', { className: 'gui-select' }, { flexGrow: '1' });
            GUIutils.setSelectOptions(selSamples, sampleScripts.map(s => ({ value: s.id, text: s.name })), { selectedValue: 'ptcda_nacl', selectFirst: true });

            let ta = null;
            selSamples.onchange = (e) => {
                const id = e && e.target ? e.target.value : '';
                const sample = sampleScripts.find(s => s.id === id);
                if (sample && ta) ta.value = sample.code;
            };

            ta = GUIutils.textArea(container, defaultScript, { 
                width: '100%', height: '180px', fontSize: '11px', fontFamily: 'monospace', backgroundColor: '#222', color: '#eee', display: 'block'
            });
            
            const rowBtns = GUIutils.row(container, { marginTop: '5px' });
            GUIutils.btn(rowBtns, 'Run Script', async () => {
                if (window.app && window.app.scriptRunner) {
                    await window.app.scriptRunner.run(ta.value);
                    this.updateSelectionUI();
                }
            }, { flexGrow: '1' });
        }, { collapsible: true });

        // --- Section: Replicas / Lattice ---
        // Store refs to keep refresh callable from main.js
        const latticeUI = {};

        this.createSection(sidebar, 'Replicas / Lattice', (container) => {
            const rowSystem = GUIutils.row(container);
            GUIutils.span(rowSystem, 'System: ');
            const rowLvec = GUIutils.row(container, { marginTop: '5px' });
            GUIutils.span(rowLvec, 'Lattice: ', { verticalAlign: 'top' });
            const taLvec = GUIutils.textArea(rowLvec, '', { height: '60px', flexGrow: '1', fontSize: '0.8em', fontFamily: 'monospace' });
            latticeUI.taLvec = taLvec;

            const updateLvecFromTA = () => {
                const name = selLattice.value;
                const lat = window.app.getLattice(name);
                const lines = taLvec.value.trim().split('\n');
                if (lines.length >= 3) {
                    const lvec = lines.slice(0, 3).map(line => {
                        const parts = line.trim().split(/\s+/).map(Number);
                        return new Vec3(parts[0] || 0, parts[1] || 0, parts[2] || 0);
                    });
                    lat.lvec = lvec;
                    const sys = window.app.systems[name];
                    if (sys) sys.lvec = [lvec[0].clone(), lvec[1].clone(), lvec[2].clone()];
                    window.app.updateReplicas(name);
                }
            };
            taLvec.onchange = updateLvecFromTA;

            const updateLvecTA = () => {
                const name = selLattice.value;
                const lat = window.app.getLattice(name);
                if (lat && lat.lvec) {
                    taLvec.value = lat.lvec.map(v => `${v.x.toFixed(4)} ${v.y.toFixed(4)} ${v.z.toFixed(4)}`).join('\n');
                } else {
                    taLvec.value = '';
                }
            };

            const selLattice = GUIutils.selectList(rowSystem, ['molecule', 'substrate'], 'molecule', (val) => {
                if (!window.app || !window.app.getLattice) return;
                const lat = window.app.getLattice(val);
                if (!lat) return;
                chkShow.checked = !!lat.show;
                chkShowBox.checked = !!lat.showBox;
                inpNx.value = lat.nrep.x;
                inpNy.value = lat.nrep.y;
                inpNz.value = lat.nrep.z;
                updateLvecTA();
            }, { flexGrow: '1' });
            latticeUI.selLattice = selLattice;

            const rowShow = GUIutils.row(container, { marginTop: '5px' });
            const chkShow = GUIutils.labelCheck(rowShow, 'Show Replicas', false, (e) => {
                if (window.app) {
                    const name = selLattice.value;
                    const lat = window.app.getLattice ? window.app.getLattice(name) : null;
                    if (!lat) return;
                    lat.show = e.target.checked;
                    window.app.updateReplicas(name);
                }
            }).input;
            latticeUI.chkShow = chkShow;

            const rowShowBox = GUIutils.row(container, { marginTop: '2px' });
            const chkShowBox = GUIutils.labelCheck(rowShowBox, 'Show Lattice Box', false, (e) => {
                if (window.app) {
                    const name = selLattice.value;
                    const lat = window.app.getLattice ? window.app.getLattice(name) : null;
                    if (!lat) return;
                    lat.showBox = e.target.checked;
                    window.app.updateLatticeBox(name);
                }
            }).input;
            latticeUI.chkShowBox = chkShowBox;

            const rowN = GUIutils.row(container, { marginTop: '5px' });
            GUIutils.span(rowN, 'nRep: ');
            const mkInt = (val, cb) => {
                const inp = GUIutils.num(rowN, val, { min: 0, step: 1 }, { width: '40px', marginLeft: '5px' });
                inp.onchange = cb;
                return inp;
            };
            const updateN = () => {
                if (window.app) {
                    const name = selLattice.value;
                    const lat = window.app.getLattice ? window.app.getLattice(name) : null;
                    if (!lat) return;
                    lat.nrep.x = parseInt(inpNx.value) || 0;
                    lat.nrep.y = parseInt(inpNy.value) || 0;
                    lat.nrep.z = parseInt(inpNz.value) || 0;
                    window.app.updateReplicas(name);
                }
            };
            const inpNx = mkInt(1, updateN);
            const inpNy = mkInt(1, updateN);
            const inpNz = mkInt(1, updateN);
            latticeUI.inpNx = inpNx;
            latticeUI.inpNy = inpNy;
            latticeUI.inpNz = inpNz;

            const rowBtns = GUIutils.row(container, { marginTop: '5px' });
            GUIutils.btn(rowBtns, 'Bake Replicas', () => {
                if (window.app) window.app.bakeReplicas(selLattice.value);
            }, { flexGrow: '1' });
        });

        // Expose refresh helper for external callers (main.js) to keep UI in sync
        this.refreshLatticeControls = (name = 'molecule') => {
            if (!window.app || !window.app.getLattice) return;
            const lat = window.app.getLattice(name);
            if (!lat || !latticeUI.selLattice) return;
            latticeUI.selLattice.value = name;
            if (latticeUI.chkShow) latticeUI.chkShow.checked = !!lat.show;
            if (latticeUI.chkShowBox) latticeUI.chkShowBox.checked = !!lat.showBox;
            if (latticeUI.inpNx) latticeUI.inpNx.value = lat.nrep.x;
            if (latticeUI.inpNy) latticeUI.inpNy.value = lat.nrep.y;
            if (latticeUI.inpNz) latticeUI.inpNz.value = lat.nrep.z;
            if (latticeUI.taLvec && lat.lvec) {
                latticeUI.taLvec.value = lat.lvec.map(v => `${v.x.toFixed(4)} ${v.y.toFixed(4)} ${v.z.toFixed(4)}`).join('\n');
            }
        };

        // --- Section: Parameters ---
        this.createSection(sidebar, 'Parameters', (container) => {
            // Element Types
            this.createParamControl(container, 'Element Types', 'common_resources/ElementTypes.dat',
                (content) => {
                    if (window.app && window.app.mmParams) {
                        window.app.mmParams.parseElementTypes(content);
                        window.app.molRenderer.updateStructure();
                    }
                }
            );

            // Atom Types
            this.createParamControl(container, 'Atom Types', 'common_resources/AtomTypes.dat',
                (content) => {
                    if (window.app && window.app.mmParams) {
                        window.app.mmParams.parseAtomTypes(content);
                    }
                }
            );
        });

        // --- Section: Help ---
        this.createSection(sidebar, 'Help', (container) => {
            const toggleHelp = () => {
                const help = document.getElementById('help-overlay');
                if (!help) return;
                help.style.display = chk.input.checked ? 'block' : 'none';
            };

            const chk = GUIutils.labelCheck(container, 'Show Help', true, toggleHelp);

            // Initialize on first load
            toggleHelp();
        });

        // --- Section: System Log ---
        this.createSection(sidebar, 'System Log', (container) => {
            // Verbosity
            const row = GUIutils.row(container);
            GUIutils.span(row, 'Verbosity: ', { fontSize: '0.9em' });
            const input = GUIutils.num(row, window.VERBOSITY_LEVEL, { min: 0, max: 4, step: '1', onchange: (e) => {
                window.VERBOSITY_LEVEL = parseInt(e.target.value);
                window.logger.info(`Verbosity set to ${window.VERBOSITY_LEVEL}`);
            } }, { width: '50px', flexGrow: '0' });

            GUIutils.btn(row, 'Clear', () => { if (window.logger) window.logger.clear(); }, { marginLeft: '5px', flexGrow: '0' });

            // Log Output
            const logOut = GUIutils.div(container, 'gui-log-output');

            if (window.logger) window.logger.setContainer(logOut);
        });

        // Final pass to apply current defaults and ensure scene/render is in sync after async init
        const finalizeInit = () => {
            const app = window.app;
            if (!app) return;
            if (app.molRenderer) {
                app.molRenderer.updateSelection();
                app.molRenderer.update();
            }
            if (app.requestRender) app.requestRender();
        };
        // Defer so window.app is set by main.js before we touch it
        setTimeout(finalizeInit, 0);
    }

    createSection(parent, title, contentFn, opts = null) {
        const section = GUIutils.div(null, 'gui-section');

        const titleEl = GUIutils.div(section, 'gui-section-title');
        titleEl.textContent = title;

        const content = GUIutils.div(section);

        const collapsible = opts && opts.collapsible;
        const openDefault = !(opts && (opts.open === false));
        if (collapsible) {
            titleEl.style.cursor = 'pointer';
            content.style.display = openDefault ? 'block' : 'none';
            titleEl.onclick = () => {
                const isOpen = content.style.display !== 'none';
                content.style.display = isOpen ? 'none' : 'block';
            };
        }

        parent.appendChild(section);
        contentFn(content);
        return { section, content, titleEl };
    }

    updateSelectionCount() {
        if (!this.lblCount) return;
        const n = (this.system && this.system.selection) ? this.system.selection.size : 0;
        this.lblCount.textContent = String(n);
    }

    updateSelectionUI() {
        this.updateSelectionCount();
        if (this.inpSelection && this.system && this.system.selection) {
            const ids = Array.from(this.system.selection);
            ids.sort((a, b) => a - b);
            this.inpSelection.value = ids.join(',');
        }
    }

    onSelectionInputChange(value) {
        const parts = value.split(',');
        this.system.clearSelection();
        let count = 0;
        for (const part of parts) {
            const id = parseInt(part.trim());
            if (!isNaN(id) && id >= 0 && id < this.system.nAtoms) {
                this.system.select(id, 'add');
                count++;
            }
        }
        this.renderer.update();
        this.requestRender();
        if (this.onSelectionChanged) this.onSelectionChanged();
        this.updateSelectionUI();
        // Also update editor gizmo if possible
        if (window.app && window.app.editor) {
            window.app.editor.updateGizmo();
        }
        window.logger.info(`Selection set from input: ${count} atoms.`);
    }

    setZoom(zoom) {
        if (!(window.app && window.app.camera && window.app.controls)) return;
        const cam = window.app.camera;
        cam.zoom = zoom;
        cam.updateProjectionMatrix();
        window.app.controls.update();
    }

    createParamControl(container, label, defaultPath, onApply) {
        const wrapper = GUIutils.div(null, null, { marginBottom: '10px', borderBottom: '1px solid #444', paddingBottom: '5px' });

        const btnToggle = GUIutils.btn(wrapper, `Show/Hide ${label}`);

        const textArea = GUIutils.textArea(wrapper, '', { placeholder: `Paste ${label} content here...` });

        // Load Default (Fetch)
        btnToggle.onclick = async () => {
            const isHidden = textArea.style.display === 'none';
            textArea.style.display = isHidden ? 'block' : 'none';
            if (isHidden && !textArea.value) {
                try {
                    const res = await fetch(defaultPath);
                    if (res.ok) {
                        textArea.value = await res.text();
                    }
                } catch (e) {
                    console.warn("Failed to fetch default params:", e);
                }
            }
        };

        // Load File Button
        const btnLoad = GUIutils.btn(wrapper, 'Load File...', null, { marginTop: '5px' });

        const fileInput = GUIutils.input(wrapper, { type: 'file' }, { display: 'none' });

        btnLoad.onclick = () => fileInput.click();
        fileInput.onchange = (e) => {
            if (e.target.files.length > 0) {
                const reader = new FileReader();
                reader.onload = (ev) => {
                    textArea.value = ev.target.result;
                    textArea.style.display = 'block'; // Show it
                };
                reader.readAsText(e.target.files[0]);
                fileInput.value = '';
            }
        };

        // Apply Button
        GUIutils.btn(wrapper, 'Apply', () => {
            if (textArea.value.trim()) {
                onApply(textArea.value);
                window.logger.info(`${label} updated.`);
            }
        }, { marginTop: '5px' });

        container.appendChild(wrapper);
    }

    setView(pos, up) {
        if (window.app && window.app.camera && window.app.controls) {
            const cam = window.app.camera;
            cam.position.set(pos[0], pos[1], pos[2]);
            cam.up.set(up[0], up[1], up[2]);
            cam.lookAt(0, 0, 0);
            cam.updateProjectionMatrix();
            window.app.controls.update();
            this.requestRender();
        }
    }

    populateStructureControls() {
        if (!window.app || !window.app.mmParams) return;
        const params = window.app.mmParams;

        // Populate Elements
        const elements = Object.values(params.elementTypes).sort((a, b) => a.iZ - b.iZ);

        GUIutils.setSelectOptions(this.selElement, elements.map(el => ({ value: el.iZ, text: `${el.name} (${el.iZ})`, selected: (el.name === 'C') })), { selectFirst: true });

        // Trigger update for Atom Types
        this.onElementChange(this.selElement.value);
    }

    onElementChange(iZ) {
        iZ = parseInt(iZ);
        if (!window.app || !window.app.mmParams) return;
        const params = window.app.mmParams;

        // Update Editor
        if (window.app.editor) {
            window.app.editor.selectedElement = iZ;
        }

        // Filter Atom Types by Element Name
        const el = params.byAtomicNumber[iZ];
        if (!el) return;
        const elName = el.name;

        this.selAtomType.innerHTML = '';
        const types = Object.values(params.atomTypes).filter(t => t.element_name === elName);

        if (types.length === 0) {
            GUIutils.setSelectOptions(this.selAtomType, [{ value: elName, text: elName, selected: true }], { selectFirst: true });
        } else {
            GUIutils.setSelectOptions(this.selAtomType, types.map(t => ({ value: t.name, text: t.name })), { selectFirst: true });
        }

        // Trigger Atom Type Update
        if (this.selAtomType.options.length > 0) {
            this.onAtomTypeChange(this.selAtomType.options[0].value);
        }
    }

    onAtomTypeChange(typeName) {
        if (window.app && window.app.editor) {
            window.app.editor.selectedAtomType = typeName;
            window.logger.info(`Selected Atom Type: ${typeName}`);
        }
    }

}
