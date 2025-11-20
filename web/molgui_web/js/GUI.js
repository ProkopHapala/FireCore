class GUI {
    constructor(io) {
        this.io = io;
        this.init();
    }

    init() {
        // Create Menu Bar
        const menu = document.createElement('div');
        menu.id = 'gui-menu';
        document.body.appendChild(menu);

        // Load Button
        const btnLoad = document.createElement('button');
        btnLoad.textContent = 'Load XYZ';
        btnLoad.className = 'gui-btn';
        menu.appendChild(btnLoad);

        // Hidden File Input
        const fileInput = document.createElement('input');
        fileInput.type = 'file';
        fileInput.accept = '.xyz';
        fileInput.style.display = 'none';
        document.body.appendChild(fileInput);

        btnLoad.onclick = () => fileInput.click();
        fileInput.onchange = (e) => {
            if (e.target.files.length > 0) {
                this.io.loadXYZ(e.target.files[0]);
                fileInput.value = ''; // Reset
            }
        };

        // Save Button
        const btnSave = document.createElement('button');
        btnSave.textContent = 'Save XYZ';
        btnSave.className = 'gui-btn';
        btnSave.onclick = () => this.io.saveFile();
        menu.appendChild(btnSave);

        // Clear Button
        const btnClear = document.createElement('button');
        btnClear.textContent = 'Clear';
        btnClear.className = 'gui-btn';
        btnClear.onclick = () => {
            this.io.system.clear();
            this.io.renderer.update();
            window.logger.info("Scene cleared.");
        };
        menu.appendChild(btnClear);

        // Selection Count / Input
        const lblSel = document.createElement('span');
        lblSel.textContent = ' Sel: ';
        lblSel.style.marginLeft = '10px';
        lblSel.style.fontSize = '0.9em';
        menu.appendChild(lblSel);

        this.inpSelection = document.createElement('input');
        this.inpSelection.type = 'text';
        this.inpSelection.className = 'gui-input';
        this.inpSelection.placeholder = 'IDs (e.g. 1,5)';
        this.inpSelection.style.width = '100px';
        this.inpSelection.onchange = (e) => this.onSelectionInputChange(e.target.value);
        menu.appendChild(this.inpSelection);

        // Help Button
        const btnHelp = document.createElement('button');
        btnHelp.textContent = '?';
        btnHelp.className = 'gui-btn';
        btnHelp.style.marginLeft = 'auto'; // Push to right
        btnHelp.onclick = () => {
            const help = document.getElementById('help-overlay');
            help.style.display = help.style.display === 'none' ? 'block' : 'none';
        };
        menu.appendChild(btnHelp);
    }

    updateSelectionCount() {
        const count = this.io.system.selection.size;
        const ids = Array.from(this.io.system.selection).sort((a, b) => a - b).join(', ');
        this.inpSelection.value = ids;
        window.logger.info(`Selection Updated: ${count} atoms (${ids})`);
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
        // Trigger external update if needed (renderer update usually handles visual, but we need to trigger selection renderer)
        // We can access main app via global or pass callback? 
        // For now, let's assume the system change triggers what's needed if we hook it up right.
        // Actually, we need to trigger selectionRenderer update.
        // The GUI doesn't have access to selectionRenderer directly, but IO has system.
        // We need a way to notify main.
        if (this.onSelectionChanged) this.onSelectionChanged();

        window.logger.info(`Selection set from input: ${count} atoms.`);
    }
}
