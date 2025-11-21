class GUI {
    constructor(io) {
        this.io = io;
        this.init();
    }

    init() {
        // Create Sidebar
        const sidebar = document.createElement('div');
        sidebar.id = 'gui-sidebar';
        document.body.appendChild(sidebar);

        // --- Section: File ---
        this.createSection(sidebar, 'File', (container) => {
            // Load
            const btnLoad = document.createElement('button');
            btnLoad.textContent = 'Load XYZ';
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
            btnSave.onclick = () => this.io.saveFile();
            container.appendChild(btnSave);

            // Clear
            const btnClear = document.createElement('button');
            btnClear.textContent = 'Clear Scene';
            btnClear.className = 'gui-btn';
            btnClear.onclick = () => {
                this.io.system.clear();
                this.io.renderer.update();
                // We need to clear editor selection too, but we don't have direct ref to editor here easily.
                // But system.clear() clears data. Editor updateGizmo checks selection size.
                // Ideally we should trigger a global "SceneCleared" event.
                // For now, let's rely on the fact that system is empty.
                window.logger.info("Scene cleared.");
            };
            container.appendChild(btnClear);
        });

        // --- Section: Selection ---
        this.createSection(sidebar, 'Selection', (container) => {
            const row = document.createElement('div');
            row.className = 'gui-row';

            const lbl = document.createElement('span');
            lbl.textContent = 'Count: ';
            lbl.style.fontSize = '0.9em';
            lbl.style.marginRight = '5px';
            row.appendChild(lbl);

            this.lblCount = document.createElement('span');
            this.lblCount.textContent = '0';
            this.lblCount.style.fontWeight = 'bold';
            row.appendChild(this.lblCount);

            container.appendChild(row);

            this.inpSelection = document.createElement('input');
            this.inpSelection.type = 'text';
            this.inpSelection.className = 'gui-input';
            this.inpSelection.placeholder = 'IDs (e.g. 1,5)';
            this.inpSelection.onchange = (e) => this.onSelectionInputChange(e.target.value);
            container.appendChild(this.inpSelection);
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

        // --- Section: View ---
        this.createSection(sidebar, 'View', (container) => {
            const row = document.createElement('div');
            row.className = 'gui-row';

            const btnZoomIn = document.createElement('button');
            btnZoomIn.textContent = 'Zoom In (+)';
            btnZoomIn.className = 'gui-btn';
            btnZoomIn.onclick = () => this.zoomCamera(-1);
            row.appendChild(btnZoomIn);

            const btnZoomOut = document.createElement('button');
            btnZoomOut.textContent = 'Zoom Out (-)';
            btnZoomOut.className = 'gui-btn';
            btnZoomOut.onclick = () => this.zoomCamera(1);
            row.appendChild(btnZoomOut);

            container.appendChild(row);
        });

        // --- Section: System Log ---
        this.createSection(sidebar, 'System Log', (container) => {
            // Verbosity
            const row = document.createElement('div');
            row.className = 'gui-row';

            const lbl = document.createElement('span');
            lbl.textContent = 'Level: ';
            lbl.style.fontSize = '0.9em';
            row.appendChild(lbl);

            const sel = document.createElement('select');
            sel.className = 'gui-input';
            sel.style.width = 'auto';
            sel.style.flexGrow = '1';
            ['DEBUG', 'INFO', 'WARN', 'ERROR'].forEach(l => {
                const opt = document.createElement('option');
                opt.value = l;
                opt.textContent = l;
                if (l === 'INFO') opt.selected = true;
                sel.appendChild(opt);
            });
            sel.onchange = (e) => {
                if (window.logger) window.logger.setLevel(e.target.value);
            };
            row.appendChild(sel);

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

        // --- Section: Help ---
        this.createSection(sidebar, 'Help', (container) => {
            const btnHelp = document.createElement('button');
            btnHelp.textContent = 'Show Controls (?)';
            btnHelp.className = 'gui-btn';
            btnHelp.onclick = () => {
                const help = document.getElementById('help-overlay');
                help.style.display = help.style.display === 'none' ? 'block' : 'none';
            };
            container.appendChild(btnHelp);
        });
    }

    createSection(parent, title, contentFn) {
        const section = document.createElement('div');
        section.className = 'gui-section';

        const titleEl = document.createElement('div');
        titleEl.className = 'gui-section-title';
        titleEl.textContent = title;
        section.appendChild(titleEl);

        contentFn(section);
        parent.appendChild(section);
    }

    updateSelectionCount() {
        const count = this.io.system.selection.size;
        const ids = Array.from(this.io.system.selection).sort((a, b) => a - b).join(', ');
        if (this.lblCount) this.lblCount.textContent = count;
        if (this.inpSelection) this.inpSelection.value = ids;
        window.logger.info(`Selection Updated: ${count} atoms`);
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

    zoomCamera(delta) {
        if (window.app && window.app.camera) {
            const cam = window.app.camera;
            // Move camera along its view vector
            const dir = new THREE.Vector3();
            cam.getWorldDirection(dir);
            cam.position.addScaledVector(dir, -delta * 2); // Reverse delta because moving forward is negative? No, addScaledVector(dir, dist). 
            // Actually, let's just use OrbitControls dolly if available, or simple translate.
            // Simple translate:
            // cam.translateZ(delta * 2); 
            // But OrbitControls might fight this.

            // Better: Update OrbitControls target distance?
            // Or just let OrbitControls handle it via dollyIn/dollyOut if accessible.
            // But they are internal.

            // Let's just move the camera and let OrbitControls re-sync on next update.
            cam.translateZ(delta * 2);
        }
    }
}
