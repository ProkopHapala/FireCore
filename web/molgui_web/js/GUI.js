import { MoleculeSystem } from './MoleculeSystem.js';

export class GUI {
    constructor(io) {
        this.io = io;
        this.monomerLib = null;
        this.monomerGeoms = new Map();
        this.endGroupParsed = null;
        this.init();
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
            const row = document.createElement('div');
            row.className = 'gui-row';

            const lbl = document.createElement('span');
            lbl.textContent       = 'Count: ';
            lbl.style.fontSize    = '0.9em';
            lbl.style.marginRight = '5px';
            row.appendChild(lbl);

            this.lblCount = document.createElement('span');
            this.lblCount.textContent      = '0';
            this.lblCount.style.fontWeight = 'bold';
            row.appendChild(this.lblCount);

            container.appendChild(row);

            this.inpSelection = document.createElement('input');
            this.inpSelection.type        = 'text';
            this.inpSelection.className   = 'gui-input';
            this.inpSelection.placeholder = 'IDs (e.g. 1,5)';
            this.inpSelection.onchange    = (e) => this.onSelectionInputChange(e.target.value);
            container.appendChild(this.inpSelection);
        });

        // --- Section: View ---
        this.createSection(sidebar, 'View', (container) => {
            // Zoom Slider
            const zoomRow = document.createElement('div');
            zoomRow.className           = 'gui-row';
            zoomRow.style.flexDirection = 'column';
            zoomRow.style.alignItems    = 'flex-start';

            const zoomLabel = document.createElement('label');
            zoomLabel.textContent   = 'Zoom Level (Log10)';
            zoomLabel.style.fontSize  = '0.9em';
            zoomRow.appendChild(zoomLabel);

            const zoomSlider = document.createElement('input');
            zoomSlider.type  = 'range';
            zoomSlider.min   = '-2.0';
            zoomSlider.max   = '3.0';
            zoomSlider.step  = '0.1';
            zoomSlider.value = '1.0'; // Default
            zoomSlider.style.width = '100%';
            zoomSlider.oninput = (e) => {
                const val = parseFloat(e.target.value);
                this.setZoom(Math.pow(10, val));
            };
            zoomRow.appendChild(zoomSlider);
            container.appendChild(zoomRow);

            // View Buttons
            const viewRow           = document.createElement('div');
            viewRow.className       = 'gui-row';
            viewRow.style.marginTop = '10px';

            const views = [
                { name: 'XY', pos: [0, 0, 20], up: [0, 1, 0] },
                { name: 'XZ', pos: [0, 20, 0], up: [0, 0, -1] },
                { name: 'YZ', pos: [20, 0, 0], up: [0, 1, 0] }
            ];

            views.forEach(v => {
                const btn = document.createElement('button');
                btn.textContent       = v.name;
                btn.className         = 'gui-btn';
                btn.style.marginRight = '2px';
                btn.onclick           = () => this.setView(v.pos, v.up);
                viewRow.appendChild(btn);
            });
            container.appendChild(viewRow);

            // Axis Toggle
            const axisRow = document.createElement('div');
            axisRow.className = 'gui-row';
            axisRow.style.marginTop = '10px';

            const axisLabel = document.createElement('label');
            axisLabel.className = 'gui-checkbox-label';
            const axisChk       = document.createElement('input');
            axisChk.type        = 'checkbox';
            axisChk.checked     = false;
            axisChk.onchange    = (e) => {
                if (window.app && window.app.molRenderer) {
                    window.app.molRenderer.toggleAxes(e.target.checked);
                }
            };
            axisLabel.appendChild(axisChk);
            axisLabel.appendChild(document.createTextNode('Show Axes'));
            axisRow.appendChild(axisLabel);
            container.appendChild(axisRow);

            // Label Mode Dropdown
            const labelRow = document.createElement('div');
            labelRow.className = 'gui-row';
            labelRow.style.marginTop = '10px';

            const labelLbl = document.createElement('span');
            labelLbl.textContent = 'Labels: ';
            labelRow.appendChild(labelLbl);

            const labelSel = document.createElement('select');
            labelSel.className = 'gui-select';
            labelSel.style.flexGrow = '1';

            const modes = [
                { value: 'none',    text: 'None'      },
                { value: 'id',      text: 'Atom ID'   },
                { value: 'element', text: 'Element'   },
                { value: 'type',    text: 'Atom Type' }
            ];

            modes.forEach(m => {
                const opt = document.createElement('option');
                opt.value = m.value;
                opt.textContent = m.text;
                labelSel.appendChild(opt);
            });

            labelSel.onchange = (e) => {
                if (window.app && window.app.molRenderer) {
                    window.app.molRenderer.setLabelMode(e.target.value);
                }
            };
            labelRow.appendChild(labelSel);
            container.appendChild(labelRow);

            // Label Color & Size
            const styleRow = document.createElement('div');
            styleRow.className = 'gui-row';
            styleRow.style.marginTop = '5px';

            const colorLbl = document.createElement('span');
            colorLbl.textContent = 'Color: ';
            styleRow.appendChild(colorLbl);

            const colorInput = document.createElement('input');
            colorInput.type = 'color';
            colorInput.value = '#ffffff';
            colorInput.style.flexGrow = '0';
            colorInput.style.width = '40px';
            styleRow.appendChild(colorInput);

            const sizeLbl = document.createElement('span');
            sizeLbl.textContent = ' Size: ';
            sizeLbl.style.marginLeft = '10px';
            styleRow.appendChild(sizeLbl);

            const sizeInput = document.createElement('input');
            sizeInput.type = 'number';
            sizeInput.value = '0.5';
            sizeInput.step = '0.1';
            sizeInput.min = '0.1';
            sizeInput.style.width = '50px';
            styleRow.appendChild(sizeInput);

            const updateLabelStyle = () => {
                if (window.app && window.app.molRenderer) {
                    window.app.molRenderer.setLabelStyle(colorInput.value, sizeInput.value);
                }
            };

            colorInput.oninput = updateLabelStyle;
            sizeInput.oninput = updateLabelStyle;

            container.appendChild(styleRow);
        });

        // --- Section: Gizmo ---
        this.createSection(sidebar, 'Gizmo', (container) => {
            // Enable Checkbox
            const label = document.createElement('label');
            label.className = 'gui-checkbox-label';
            const chk = document.createElement('input');
            chk.type = 'checkbox';
            chk.checked = true; // Default on
            chk.onchange = (e) => {
                if (window.app && window.app.editor) {
                    window.app.editor.toggleGizmo(e.target.checked);
                }
            };
            label.appendChild(chk);
            label.appendChild(document.createTextNode('Enable Gizmo'));
            container.appendChild(label);

            // Lock Selection Checkbox
            const labelLock = document.createElement('label');
            labelLock.className = 'gui-checkbox-label';
            labelLock.style.marginTop = '5px';
            const chkLock = document.createElement('input');
            chkLock.type = 'checkbox';
            chkLock.checked = false;
            chkLock.onchange = (e) => {
                if (window.app && window.app.editor) {
                    window.app.editor.selectionLocked = e.target.checked;
                    window.logger.info(`Selection ${e.target.checked ? 'Locked' : 'Unlocked'}`);
                }
            };
            labelLock.appendChild(chkLock);
            labelLock.appendChild(document.createTextNode('Lock Selection'));
            container.appendChild(labelLock);

            // Modes
            const modes = ['translate', 'rotate', 'scale'];
            const modeContainer = document.createElement('div');
            modeContainer.style.marginTop = '10px';

            modes.forEach(mode => {
                const mLabel = document.createElement('label');
                mLabel.className = 'gui-checkbox-label';
                mLabel.style.marginBottom = '5px';

                const radio = document.createElement('input');
                radio.type = 'radio';
                radio.name = 'gizmo-mode';
                radio.value = mode;
                radio.checked = (mode === 'translate');
                radio.onchange = () => {
                    if (window.app && window.app.editor) {
                        window.app.editor.setGizmoMode(mode);
                    }
                };

                mLabel.appendChild(radio);
                mLabel.appendChild(document.createTextNode(mode.charAt(0).toUpperCase() + mode.slice(1)));
                modeContainer.appendChild(mLabel);
            });
            container.appendChild(modeContainer);
        });

        // --- Section: Structure ---
        this.createSection(sidebar, 'Structure', (container) => {
            // Recalculate Bonds
            const btnRecalc = document.createElement('button');
            btnRecalc.textContent = 'Recalculate Bonds';
            btnRecalc.className = 'gui-btn';
            btnRecalc.onclick = () => {
                if (window.app && window.app.editor) {
                    window.app.editor.recalculateBonds();
                }
            };
            container.appendChild(btnRecalc);

            const hr = document.createElement('hr');
            hr.style.borderColor = '#444';
            hr.style.margin = '10px 0';
            container.appendChild(hr);

            // Add Atom Controls
            const lblAdd = document.createElement('div');
            lblAdd.textContent = 'Add Atom Settings:';
            lblAdd.style.fontSize = '0.9em';
            lblAdd.style.marginBottom = '5px';
            container.appendChild(lblAdd);

            // Element Dropdown
            const rowEl = document.createElement('div');
            rowEl.className = 'gui-row';
            const lblEl = document.createElement('span');
            lblEl.textContent = 'Element: ';
            rowEl.appendChild(lblEl);

            this.selElement = document.createElement('select');
            this.selElement.className = 'gui-select';
            this.selElement.style.flexGrow = '1';
            this.selElement.onchange = (e) => this.onElementChange(e.target.value);
            rowEl.appendChild(this.selElement);
            container.appendChild(rowEl);

            // Atom Type Dropdown
            const rowType = document.createElement('div');
            rowType.className = 'gui-row';
            rowType.style.marginTop = '5px';
            const lblType = document.createElement('span');
            lblType.textContent = 'Type: ';
            rowType.appendChild(lblType);

            this.selAtomType = document.createElement('select');
            this.selAtomType.className = 'gui-select';
            this.selAtomType.style.flexGrow = '1';
            this.selAtomType.onchange = (e) => this.onAtomTypeChange(e.target.value);
            rowType.appendChild(this.selAtomType);
            container.appendChild(rowType);

            // Populate initially (delayed to ensure MMParams loaded)
            setTimeout(() => this.populateStructureControls(), 1000);
        });

        // --- Section: Geometry ---
        this.createSection(sidebar, 'Geometry', (container) => {
            // Load
            const btnLoad = document.createElement('button');
            btnLoad.textContent = 'Load XYZ File...';
            btnLoad.className = 'gui-btn';
            container.appendChild(btnLoad);

            const fileInput = document.createElement('input');
            fileInput.type = 'file';
            fileInput.accept = '.xyz';
            fileInput.style.display = 'none';
            document.body.appendChild(fileInput);

            btnLoad.onclick = () => fileInput.click();
            fileInput.onchange = (e) => {
                if (e.target.files.length > 0) {
                    this.io.loadXYZ(e.target.files[0]);
                    fileInput.value = '';
                }
            };

            // Save
            const btnSave = document.createElement('button');
            btnSave.textContent = 'Save XYZ';
            btnSave.className = 'gui-btn';
            btnSave.style.marginTop = '5px';
            btnSave.onclick = () => this.io.saveFile();
            container.appendChild(btnSave);

            // Clear
            const btnClear = document.createElement('button');
            btnClear.textContent = 'Clear Scene';
            btnClear.className = 'gui-btn';
            btnClear.style.marginTop = '5px';
            btnClear.onclick = () => {
                this.io.system.clear();
                this.io.renderer.update();
                window.logger.info("Scene cleared.");
            };
            container.appendChild(btnClear);

            // Separator
            const hr = document.createElement('hr');
            hr.style.borderColor = '#444';
            hr.style.margin = '10px 0';
            container.appendChild(hr);

            // Manual Edit
            const btnToggle = document.createElement('button');
            btnToggle.textContent = 'Edit XYZ Manually';
            btnToggle.className = 'gui-btn';
            container.appendChild(btnToggle);

            const textArea = document.createElement('textarea');
            textArea.className = 'gui-textarea';
            textArea.style.display = 'none';
            textArea.style.height = '150px';
            textArea.style.width = '100%';
            textArea.style.marginTop = '5px';
            textArea.style.fontSize = '0.8em';
            textArea.style.fontFamily = 'monospace';
            textArea.placeholder = 'Paste XYZ content here...';
            container.appendChild(textArea);

            btnToggle.onclick = () => {
                const isHidden = textArea.style.display === 'none';
                textArea.style.display = isHidden ? 'block' : 'none';
                if (isHidden) {
                    // Populate with current system state if showing
                    textArea.value = this.io.generateXYZ();
                }
            };

            const btnApply = document.createElement('button');
            btnApply.textContent = 'Apply XYZ';
            btnApply.className = 'gui-btn';
            btnApply.style.marginTop = '5px';
            btnApply.onclick = () => {
                if (textArea.value.trim()) {
                    this.io.loadXYZString(textArea.value);
                }
            };
            container.appendChild(btnApply);
        });

        // --- Section: Builder - Crystal/Substrate ---
        this.createSection(sidebar, 'Builder: Substrate', (container) => {
            const rowPreset = document.createElement('div');
            rowPreset.className = 'gui-row';
            const lblPreset = document.createElement('span');
            lblPreset.textContent = 'Preset: ';
            rowPreset.appendChild(lblPreset);
            const selPreset = document.createElement('select');
            selPreset.className = 'gui-select';
            selPreset.style.flexGrow = '1';
            ['NaCl(step)', 'NaCl(rocksalt)', 'KBr(rocksalt)', 'MgO(rocksalt)', 'CaF2(fluorite)', 'CaCO3(todo)'].forEach(name => {
                const opt = document.createElement('option');
                opt.value = name;
                opt.textContent = name;
                selPreset.appendChild(opt);
            });
            rowPreset.appendChild(selPreset);
            container.appendChild(rowPreset);

            const rowA = document.createElement('div');
            rowA.className = 'gui-row';
            const lblA = document.createElement('span');
            lblA.textContent = 'a(Å): ';
            rowA.appendChild(lblA);
            const inpA = document.createElement('input');
            inpA.type = 'number';
            inpA.step = '0.01';
            inpA.value = '2.82';
            inpA.className = 'gui-input';
            inpA.style.width = '70px';
            inpA.style.flexGrow = '0';
            rowA.appendChild(inpA);
            container.appendChild(rowA);

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

            const rowN = document.createElement('div');
            rowN.className = 'gui-row';
            const mkInt = (label, val) => {
                const sp = document.createElement('span');
                sp.textContent = label;
                sp.style.marginRight = '2px';
                const inp = document.createElement('input');
                inp.type = 'number';
                inp.step = '1';
                inp.value = String(val);
                inp.className = 'gui-input';
                inp.style.width = '44px';
                inp.style.flexGrow = '0';
                rowN.appendChild(sp);
                rowN.appendChild(inp);
                return inp;
            };
            const inpNx = mkInt('nx', 10);
            const inpNy = mkInt('ny', 10);
            const inpNz = mkInt('nz', 3);
            container.appendChild(rowN);

            const chkMillerLbl = document.createElement('label');
            chkMillerLbl.className = 'gui-checkbox-label';
            chkMillerLbl.style.marginTop = '6px';
            const chkMiller = document.createElement('input');
            chkMiller.type = 'checkbox';
            chkMiller.checked = false;
            chkMillerLbl.appendChild(chkMiller);
            chkMillerLbl.appendChild(document.createTextNode('Orient by Miller (h k l) -> z'));
            container.appendChild(chkMillerLbl);

            const rowHKL = document.createElement('div');
            rowHKL.className = 'gui-row';
            rowHKL.style.marginTop = '4px';
            const mkIntHKL = (label, val) => {
                const sp = document.createElement('span');
                sp.textContent = label;
                sp.style.marginRight = '2px';
                const inp = document.createElement('input');
                inp.type = 'number';
                inp.step = '1';
                inp.value = String(val);
                inp.className = 'gui-input';
                inp.style.width = '44px';
                inp.style.flexGrow = '0';
                rowHKL.appendChild(sp);
                rowHKL.appendChild(inp);
                return inp;
            };
            const inpH = mkIntHKL('h', 1);
            const inpK = mkIntHKL('k', 0);
            const inpL = mkIntHKL('l', 0);
            container.appendChild(rowHKL);

            const btnRow = document.createElement('div');
            btnRow.className = 'gui-row';
            btnRow.style.marginTop = '6px';
            const btnReplace = document.createElement('button');
            btnReplace.textContent = 'Generate (Replace)';
            btnReplace.className = 'gui-btn';
            const btnAppend = document.createElement('button');
            btnAppend.textContent = 'Generate (Append)';
            btnAppend.className = 'gui-btn';
            btnAppend.style.marginLeft = '4px';
            btnRow.appendChild(btnReplace);
            btnRow.appendChild(btnAppend);
            container.appendChild(btnRow);

            const getHKL = () => [parseInt(inpH.value), parseInt(inpK.value), parseInt(inpL.value)];

            const applyMiller = (data) => {
                if (!chkMiller.checked) return data;
                const [h, k, l] = getHKL();
                if ((h | 0) === 0 && (k | 0) === 0 && (l | 0) === 0) throw new Error('HKL cannot be (0,0,0)');
                if (!data.lvec) throw new Error('applyMiller: missing lvec');
                const b = MoleculeSystem.reciprocalLattice(data.lvec);
                const n = [h * b[0][0] + k * b[1][0] + l * b[2][0], h * b[0][1] + k * b[1][1] + l * b[2][1], h * b[0][2] + k * b[1][2] + l * b[2][2]];
                const R = MoleculeSystem.rotationAlignVectorToZ(n);
                const pos = MoleculeSystem.rotatePosArray(data.pos, R);
                const lvec = MoleculeSystem.rotateLvec(data.lvec, R);
                return { ...data, pos, lvec };
            };

            const buildData = () => {
                const preset = selPreset.value;
                const a = +inpA.value;
                const nx = parseInt(inpNx.value) | 0;
                const ny = parseInt(inpNy.value) | 0;
                const nz = parseInt(inpNz.value) | 0;
                if (!(a > 0)) throw new Error('a must be >0');
                if (nx <= 0 || ny <= 0 || nz <= 0) throw new Error('nx,ny,nz must be >0');
                let data = null;
                if (preset === 'NaCl(step)') {
                    data = MoleculeSystem.genNaClStep({ a, nx, ny, nz, Q0: 0.7 });
                } else if (preset.endsWith('(rocksalt)')) {
                    const zA = (preset.startsWith('NaCl')) ? 11 : (preset.startsWith('KBr') ? 19 : 12);
                    const zB = (preset.startsWith('NaCl')) ? 17 : (preset.startsWith('KBr') ? 35 : 8);
                    const lvec = [[a, 0, 0], [0, a, 0], [0, 0, a]];
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
                    data = MoleculeSystem.genReplicatedCell({ lvec, basisPos, basisTypes, nRep: [nx, ny, nz] });
                } else if (preset === 'CaF2(fluorite)') {
                    const zCa = 20;
                    const zF = 9;
                    const lvec = [[a, 0, 0], [0, a, 0], [0, 0, a]];
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
                    data = MoleculeSystem.genReplicatedCell({ lvec, basisPos, basisTypes, nRep: [nx, ny, nz] });
                } else if (preset === 'CaCO3(todo)') {
                    throw new Error('CaCO3 preset not implemented yet (non-trivial basis). Use a custom unit cell (e.g. extended XYZ with Lattice=...) once supported.');
                } else {
                    throw new Error(`Preset not implemented: ${preset}`);
                }
                return applyMiller(data);
            };

            const applyToScene = (data, mode) => {
                if (mode === 'replace') this.io.system.clear();
                const offset = this.io.system.nAtoms;
                this.io.system.addAtomsFromArrays(data.pos, data.types);
                if (data.bonds && data.bonds.length) {
                    for (const [a0, b0] of data.bonds) this.io.system.bonds.push([a0 + offset, b0 + offset]);
                    this.io.system.updateNeighborList();
                }
                this.io.renderer.update();
                window.logger.info(`Substrate generated: atoms=${this.io.system.nAtoms}`);
            };

            btnReplace.onclick = () => {
                try { applyToScene(buildData(), 'replace'); } catch (e) { window.logger.error(String(e)); throw e; }
            };
            btnAppend.onclick = () => {
                try { applyToScene(buildData(), 'append'); } catch (e) { window.logger.error(String(e)); throw e; }
            };
        }, { collapsible: true, open: false });

        // --- Section: Builder - Polymers / Attachment ---
        this.createSection(sidebar, 'Builder: Polymers', (container) => {
            const mkLbl = (t) => { const d = document.createElement('div'); d.textContent = t; d.style.marginTop = '4px'; d.style.fontSize = '0.9em'; d.style.color = '#ccc'; return d; };

            container.appendChild(mkLbl('Monomer library (JSON + mol2 files)'));

            const btnLoadLib = document.createElement('button');
            btnLoadLib.textContent = 'Load Library JSON...';
            btnLoadLib.className = 'gui-btn';
            container.appendChild(btnLoadLib);

            const inpLib = document.createElement('input');
            inpLib.type = 'file';
            inpLib.accept = '.json,application/json';
            inpLib.style.display = 'none';
            container.appendChild(inpLib);

            const lblLib = document.createElement('div');
            lblLib.textContent = 'No library loaded';
            lblLib.style.fontSize = '0.85em';
            lblLib.style.color = '#aaa';
            lblLib.style.marginTop = '2px';
            container.appendChild(lblLib);

            btnLoadLib.onclick = () => inpLib.click();
            inpLib.onchange = async (e) => {
                if (e.target.files.length <= 0) return;
                const file = e.target.files[0];
                const txt = await this.readFileText(file);
                this.monomerLib = JSON.parse(txt);
                if (!this.monomerLib || !this.monomerLib.monomers) throw new Error('Library JSON missing monomers[]');
                lblLib.textContent = `Library: ${this.monomerLib.name || file.name} (monomers=${this.monomerLib.monomers.length})`;
                window.logger.info(`Loaded library JSON: ${file.name}`);
                inpLib.value = '';
            };

            const btnLoadMol2 = document.createElement('button');
            btnLoadMol2.textContent = 'Load mol2 geometries...';
            btnLoadMol2.className = 'gui-btn';
            btnLoadMol2.style.marginTop = '4px';
            container.appendChild(btnLoadMol2);

            const inpMol2 = document.createElement('input');
            inpMol2.type = 'file';
            inpMol2.accept = '.mol2';
            inpMol2.multiple = true;
            inpMol2.style.display = 'none';
            container.appendChild(inpMol2);

            const lblMol2 = document.createElement('div');
            lblMol2.textContent = 'mol2 loaded: 0';
            lblMol2.style.fontSize = '0.85em';
            lblMol2.style.color = '#aaa';
            lblMol2.style.marginTop = '2px';
            container.appendChild(lblMol2);

            btnLoadMol2.onclick = () => inpMol2.click();
            inpMol2.onchange = async (e) => {
                if (e.target.files.length <= 0) return;
                for (const file of e.target.files) {
                    const txt = await this.readFileText(file);
                    const parsed = MoleculeSystem.parseMol2(txt);
                    this.monomerGeoms.set(file.name, parsed);
                }
                lblMol2.textContent = `mol2 loaded: ${this.monomerGeoms.size}`;
                window.logger.info(`Loaded mol2 geometries: ${e.target.files.length}`);
                inpMol2.value = '';
            };

            container.appendChild(mkLbl('Sequence'));

            const inpSeq = document.createElement('input');
            inpSeq.type = 'text';
            inpSeq.className = 'gui-input';
            inpSeq.placeholder = 'e.g. D3Gly6A';
            container.appendChild(inpSeq);

            const rowSeqBtn = document.createElement('div');
            rowSeqBtn.className = 'gui-row';
            rowSeqBtn.style.marginTop = '4px';
            const btnSeqReplace = document.createElement('button');
            btnSeqReplace.textContent = 'Build (Replace)';
            btnSeqReplace.className = 'gui-btn';
            const btnSeqAppend = document.createElement('button');
            btnSeqAppend.textContent = 'Build (Append)';
            btnSeqAppend.className = 'gui-btn';
            btnSeqAppend.style.marginLeft = '4px';
            rowSeqBtn.appendChild(btnSeqReplace);
            rowSeqBtn.appendChild(btnSeqAppend);
            container.appendChild(rowSeqBtn);

            const buildMonomerMap = () => {
                if (!this.monomerLib) throw new Error('No monomer library loaded');
                const map = {};
                for (const m of this.monomerLib.monomers) {
                    if (!m.id) throw new Error('Monomer missing id');
                    if (!m.file) throw new Error(`Monomer '${m.id}' missing file`);
                    const parsed = this.monomerGeoms.get(m.file);
                    if (!parsed) throw new Error(`Monomer '${m.id}' missing loaded geometry file '${m.file}'`);
                    if (!m.anchors || !m.anchors.head || !m.anchors.tail) throw new Error(`Monomer '${m.id}' missing anchors.head/tail`);
                    if (m.anchors.head.type !== 'index' || m.anchors.tail.type !== 'index') throw new Error(`Monomer '${m.id}': only anchors type='index' supported (v0)`);
                    map[m.id] = { parsed, anchors: [m.anchors.head.value, m.anchors.tail.value] };
                    if (m.aliases) for (const a of m.aliases) map[a] = map[m.id];
                }
                return map;
            };

            const applyPolymerToScene = (poly, mode) => {
                const n = poly.nAtoms;
                const pos = poly.pos.subarray(0, n * 3);
                const types = poly.types.subarray(0, n);
                const parsed = { pos, types, bonds: poly.bonds };
                if (mode === 'replace') this.io.system.clear();
                const off = this.io.system.appendParsedSystem(parsed);
                if (off !== 0) this.io.system.updateNeighborList();
                this.io.renderer.update();
                window.logger.info(`Polymer built: atoms=${n} bonds=${poly.bonds.length}`);
            };

            const doBuildSeq = (mode) => {
                const seq = inpSeq.value.trim();
                if (!seq) throw new Error('Empty sequence');
                const tokens = this.parseSequenceTokens(seq);
                const monomers = buildMonomerMap();
                const poly = MoleculeSystem.assemblePolymerFromTokens(tokens, monomers, { _0: 1, capacity: 200000 });
                applyPolymerToScene(poly, mode);
            };

            btnSeqReplace.onclick = () => { try { doBuildSeq('replace'); } catch (e) { window.logger.error(String(e)); throw e; } };
            btnSeqAppend.onclick = () => { try { doBuildSeq('append'); } catch (e) { window.logger.error(String(e)); throw e; } };

            const hr = document.createElement('hr');
            hr.style.borderColor = '#444';
            hr.style.margin = '8px 0';
            container.appendChild(hr);

            container.appendChild(mkLbl('End-group attachment'));

            const btnLoadEnd = document.createElement('button');
            btnLoadEnd.textContent = 'Load Endgroup mol2...';
            btnLoadEnd.className = 'gui-btn';
            container.appendChild(btnLoadEnd);

            const inpEnd = document.createElement('input');
            inpEnd.type = 'file';
            inpEnd.accept = '.mol2';
            inpEnd.style.display = 'none';
            container.appendChild(inpEnd);

            const lblEnd = document.createElement('div');
            lblEnd.textContent = 'No endgroup loaded';
            lblEnd.style.fontSize = '0.85em';
            lblEnd.style.color = '#aaa';
            lblEnd.style.marginTop = '2px';
            container.appendChild(lblEnd);

            btnLoadEnd.onclick = () => inpEnd.click();
            inpEnd.onchange = async (e) => {
                if (e.target.files.length <= 0) return;
                const file = e.target.files[0];
                const txt = await this.readFileText(file);
                this.endGroupParsed = MoleculeSystem.parseMol2(txt);
                lblEnd.textContent = `Endgroup: ${file.name} (atoms=${this.endGroupParsed.types.length})`;
                window.logger.info(`Loaded endgroup mol2: ${file.name}`);
                inpEnd.value = '';
            };

            const rowMarker = document.createElement('div');
            rowMarker.className = 'gui-row';
            rowMarker.style.marginTop = '4px';
            const mkSmall = (lab, val) => {
                const sp = document.createElement('span');
                sp.textContent = lab;
                const inp = document.createElement('input');
                inp.type = 'text';
                inp.value = val;
                inp.className = 'gui-input';
                inp.style.width = '48px';
                inp.style.flexGrow = '0';
                inp.style.marginLeft = '2px';
                rowMarker.appendChild(sp);
                rowMarker.appendChild(inp);
                return inp;
            };
            const inpMX = mkSmall('X', 'Se');
            const inpMY = mkSmall('Y', 'Cl');
            container.appendChild(rowMarker);

            const btnAttachMarker = document.createElement('button');
            btnAttachMarker.textContent = 'Attach by markers (all sites)';
            btnAttachMarker.className = 'gui-btn';
            btnAttachMarker.style.marginTop = '4px';
            btnAttachMarker.onclick = () => {
                if (!this.endGroupParsed) throw new Error('No endgroup loaded');
                this.io.system.attachGroupByMarker(this.endGroupParsed, inpMX.value.trim(), inpMY.value.trim());
                this.io.renderer.update();
            };
            container.appendChild(btnAttachMarker);

            const hr2 = document.createElement('hr');
            hr2.style.borderColor = '#444';
            hr2.style.margin = '8px 0';
            container.appendChild(hr2);

            container.appendChild(mkLbl('Attach by picked direction (cap atom)'));

            const rowIds = document.createElement('div');
            rowIds.className = 'gui-row';
            const inpCap = document.createElement('input');
            inpCap.type = 'number';
            inpCap.className = 'gui-input';
            inpCap.placeholder = 'cap id';
            inpCap.style.width = '56px';
            inpCap.style.flexGrow = '0';
            const inpBack = document.createElement('input');
            inpBack.type = 'number';
            inpBack.className = 'gui-input';
            inpBack.placeholder = 'back id';
            inpBack.style.width = '56px';
            inpBack.style.flexGrow = '0';
            rowIds.appendChild(document.createTextNode('cap '));
            rowIds.appendChild(inpCap);
            rowIds.appendChild(document.createTextNode(' back '));
            rowIds.appendChild(inpBack);
            container.appendChild(rowIds);

            const btnUseSel = document.createElement('button');
            btnUseSel.textContent = 'Use selection -> cap/back';
            btnUseSel.className = 'gui-btn';
            btnUseSel.style.marginTop = '4px';
            btnUseSel.onclick = () => {
                const ids = Array.from(this.io.system.selection).sort((a, b) => a - b);
                if (ids.length <= 0) throw new Error('Nothing selected');
                inpCap.value = String(ids[0]);
                inpBack.value = (ids.length > 1) ? String(ids[1]) : '';
            };
            container.appendChild(btnUseSel);

            const rowGeom = document.createElement('div');
            rowGeom.className = 'gui-row';
            rowGeom.style.marginTop = '4px';
            const inpBond = document.createElement('input');
            inpBond.type = 'number';
            inpBond.step = '0.01';
            inpBond.value = '1.50';
            inpBond.className = 'gui-input';
            inpBond.style.width = '62px';
            inpBond.style.flexGrow = '0';
            const inpTw = document.createElement('input');
            inpTw.type = 'number';
            inpTw.step = '1';
            inpTw.value = '0';
            inpTw.className = 'gui-input';
            inpTw.style.width = '62px';
            inpTw.style.flexGrow = '0';
            rowGeom.appendChild(document.createTextNode('bond '));
            rowGeom.appendChild(inpBond);
            rowGeom.appendChild(document.createTextNode(' twist° '));
            rowGeom.appendChild(inpTw);
            container.appendChild(rowGeom);

            const rowUp = document.createElement('div');
            rowUp.className = 'gui-row';
            rowUp.style.marginTop = '4px';
            const mkF = (val) => { const i = document.createElement('input'); i.type = 'number'; i.step = '0.1'; i.value = String(val); i.className = 'gui-input'; i.style.width = '44px'; i.style.flexGrow = '0'; return i; };
            const upx = mkF(0); const upy = mkF(0); const upz = mkF(1);
            rowUp.appendChild(document.createTextNode('up '));
            rowUp.appendChild(upx); rowUp.appendChild(upy); rowUp.appendChild(upz);
            container.appendChild(rowUp);

            const rowRefs = document.createElement('div');
            rowRefs.className = 'gui-row';
            rowRefs.style.marginTop = '4px';
            const mkInt = (lbl, val) => { rowRefs.appendChild(document.createTextNode(lbl + ' ')); const i = document.createElement('input'); i.type = 'number'; i.step = '1'; i.value = String(val); i.className = 'gui-input'; i.style.width = '44px'; i.style.flexGrow = '0'; rowRefs.appendChild(i); return i; };
            const gA = mkInt('gA', 1);
            const gF = mkInt('gF', 2);
            const gU = mkInt('gU', 0);
            container.appendChild(rowRefs);

            const btnAttachDir = document.createElement('button');
            btnAttachDir.textContent = 'Attach (direction)';
            btnAttachDir.className = 'gui-btn';
            btnAttachDir.style.marginTop = '4px';
            btnAttachDir.onclick = () => {
                if (!this.endGroupParsed) throw new Error('No endgroup loaded');
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
                this.io.system.attachParsedByDirection(cap, this.endGroupParsed, params);
                this.io.renderer.update();
            };
            container.appendChild(btnAttachDir);
        }, { collapsible: true, open: false });

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
                        // Atom types might not affect rendering directly yet (unless we use them for something)
                        // But good to update.
                    }
                }
            );
        });

        // --- Section: Help ---
        this.createSection(sidebar, 'Help', (container) => {
            const label = document.createElement('label');
            label.className = 'gui-checkbox-label';

            const chk = document.createElement('input');
            chk.type = 'checkbox';
            chk.checked = true; // Help visible by default

            const toggleHelp = () => {
                const help = document.getElementById('help-overlay');
                if (!help) return;
                help.style.display = chk.checked ? 'block' : 'none';
            };

            chk.onchange = toggleHelp;
            label.appendChild(chk);
            label.appendChild(document.createTextNode('Show Help'));
            container.appendChild(label);

            // Initialize on first load
            toggleHelp();
        });

        // --- Section: System Log ---
        this.createSection(sidebar, 'System Log', (container) => {
            // Verbosity
            const row = document.createElement('div');
            row.className = 'gui-row';

            const lbl = document.createElement('span');
            lbl.textContent = 'Verbosity: ';
            lbl.style.fontSize = '0.9em';
            row.appendChild(lbl);

            const input = document.createElement('input');
            input.type = 'number';
            input.className = 'gui-input';
            input.style.width = '50px';
            input.style.flexGrow = '0';
            input.value = window.VERBOSITY_LEVEL;
            input.min = 0;
            input.max = 4;
            input.onchange = (e) => {
                window.VERBOSITY_LEVEL = parseInt(e.target.value);
                window.logger.info(`Verbosity set to ${window.VERBOSITY_LEVEL}`);
            };
            row.appendChild(input);

            const btnClear = document.createElement('button');
            btnClear.textContent = 'Clear';
            btnClear.className = 'gui-btn';
            btnClear.style.marginLeft = '5px';
            btnClear.style.flexGrow = '0';
            btnClear.onclick = () => {
                if (window.logger) window.logger.clear();
            };
            row.appendChild(btnClear);
            container.appendChild(row);

            // Log Output
            const logOut = document.createElement('div');
            logOut.className = 'gui-log-output';
            container.appendChild(logOut);

            if (window.logger) window.logger.setContainer(logOut);
        });
    }

    createSection(parent, title, contentFn, opts = null) {
        const section = document.createElement('div');
        section.className = 'gui-section';

        const titleEl = document.createElement('div');
        titleEl.className = 'gui-section-title';
        titleEl.textContent = title;
        section.appendChild(titleEl);

        const content = document.createElement('div');
        section.appendChild(content);

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
        const n = (this.io && this.io.system && this.io.system.selection) ? this.io.system.selection.size : 0;
        this.lblCount.textContent = String(n);
    }

    onSelectionInputChange(value) {
        const parts = value.split(',');
        this.io.system.clearSelection();
        let count = 0;
        for (const part of parts) {
            const id = parseInt(part.trim());
            if (!isNaN(id) && id >= 0 && id < this.io.system.nAtoms) {
                this.io.system.select(id, 'add');
                count++;
            }
        }
        this.io.renderer.update();
        if (this.onSelectionChanged) this.onSelectionChanged();
        // Also update editor gizmo if possible
        if (window.app && window.app.editor) {
            window.app.editor.updateGizmo();
        }
        window.logger.info(`Selection set from input: ${count} atoms.`);
    }

    setZoom(zoomValue) {
        if (window.app && window.app.camera) {
            const cam = window.app.camera;
            // For OrthographicCamera, zoom is controlled by 'zoom' property
            cam.zoom = zoomValue;
            cam.updateProjectionMatrix();
        }
    }

    createParamControl(container, label, defaultPath, onApply) {
        const wrapper = document.createElement('div');
        wrapper.style.marginBottom = '10px';
        wrapper.style.borderBottom = '1px solid #444';
        wrapper.style.paddingBottom = '5px';

        const btnToggle = document.createElement('button');
        btnToggle.textContent = `Show/Hide ${label}`;
        btnToggle.className = 'gui-btn';
        wrapper.appendChild(btnToggle);

        const textArea = document.createElement('textarea');
        textArea.className = 'gui-textarea';
        textArea.style.display = 'none';
        textArea.style.height = '100px';
        textArea.style.width = '100%';
        textArea.style.marginTop = '5px';
        textArea.style.fontSize = '0.8em';
        textArea.style.fontFamily = 'monospace';
        textArea.placeholder = `Paste ${label} content here...`;
        wrapper.appendChild(textArea);

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
        const btnLoad = document.createElement('button');
        btnLoad.textContent = 'Load File...';
        btnLoad.className = 'gui-btn';
        btnLoad.style.marginTop = '5px';
        wrapper.appendChild(btnLoad);

        const fileInput = document.createElement('input');
        fileInput.type = 'file';
        fileInput.style.display = 'none';
        wrapper.appendChild(fileInput);

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
        const btnApply = document.createElement('button');
        btnApply.textContent = 'Apply';
        btnApply.className = 'gui-btn';
        btnApply.style.marginTop = '5px';
        btnApply.onclick = () => {
            if (textArea.value.trim()) {
                onApply(textArea.value);
                window.logger.info(`${label} updated.`);
            }
        };
        wrapper.appendChild(btnApply);

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
        }
    }

    populateStructureControls() {
        if (!window.app || !window.app.mmParams) return;
        const params = window.app.mmParams;

        // Populate Elements
        this.selElement.innerHTML = '';
        const elements = Object.values(params.elementTypes).sort((a, b) => a.iZ - b.iZ);

        // Add common ones first if needed, or just all
        for (const el of elements) {
            const opt = document.createElement('option');
            opt.value = el.iZ;
            opt.textContent = `${el.name} (${el.iZ})`;
            if (el.name === 'C') opt.selected = true;
            this.selElement.appendChild(opt);
        }

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
            // Fallback if no specific types
            const opt = document.createElement('option');
            opt.value = elName;
            opt.textContent = elName;
            this.selAtomType.appendChild(opt);
        } else {
            for (const t of types) {
                const opt = document.createElement('option');
                opt.value = t.name;
                opt.textContent = t.name;
                this.selAtomType.appendChild(opt);
            }
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
