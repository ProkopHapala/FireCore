class MoleculeSystem {
    constructor(capacity = 1000) {
        this.capacity = capacity;
        this.nAtoms = 0;

        // Structure of Arrays (SoA)
        this.pos = new Float32Array(this.capacity * 3); // x, y, z
        this.types = new Uint8Array(this.capacity);     // element type (atomic number)

        // Bonds (Simple list for now, can be optimized later)
        this.bonds = []; // Array of [id1, id2]

        // Dirty flags
        this.isDirty = true;

        // Selection
        this.selection = new Set();
    }

    select(id, mode = 'replace') {
        if (mode === 'replace') {
            this.selection.clear();
            this.selection.add(id);
        } else if (mode === 'add') {
            this.selection.add(id);
        } else if (mode === 'subtract') {
            this.selection.delete(id);
        }
        this.isDirty = true; // Trigger redraw to update colors
    }

    clearSelection() {
        this.selection.clear();
        this.isDirty = true;
    }

    addAtom(x, y, z, type) {
        if (this.nAtoms >= this.capacity) {
            this.resize(this.capacity * 2);
        }

        const i = this.nAtoms;
        this.pos[i * 3] = x;
        this.pos[i * 3 + 1] = y;
        this.pos[i * 3 + 2] = z;
        this.types[i] = type;

        this.nAtoms++;
        this.isDirty = true;
        return i;
    }

    addBond(id1, id2) {
        this.bonds.push([id1, id2]);
        this.isDirty = true;
    }

    resize(newCapacity) {
        window.logger.info(`Resizing MoleculeSystem to ${newCapacity}`);
        const newPos = new Float32Array(newCapacity * 3);
        const newTypes = new Uint8Array(newCapacity);

        newPos.set(this.pos);
        newTypes.set(this.types);

        this.pos = newPos;
        this.types = newTypes;
        this.capacity = newCapacity;
    }

    clear() {
        this.nAtoms = 0;
        this.bonds = [];
        this.isDirty = true;
    }
}
