import { MeshRenderer } from '../../common_js/MeshRenderer.js';
import { Draw3D } from '../../common_js/Draw3D.js';
import { Logger } from '../../common_js/Logger.js';

export class PackedMolecule {
    constructor(capacity = 1024) {
        this.capacity = capacity | 0;
        this.nAtoms = 0;
        this.pos = new Float32Array(this.capacity * 3);
        this.types = new Uint8Array(this.capacity);
        this.bonds = [];
        this.selection = new Set();
        this.atomIds = new Int32Array(this.capacity);
        this.isDirty = true;
    }

    resize(newCapacity) {
        newCapacity = newCapacity | 0;
        if (newCapacity <= (this.capacity | 0)) return;
        const newPos = new Float32Array(newCapacity * 3);
        const newTypes = new Uint8Array(newCapacity);
        const newAtomIds = new Int32Array(newCapacity);
        newPos.set(this.pos);
        newTypes.set(this.types);
        newAtomIds.set(this.atomIds);
        this.pos = newPos;
        this.types = newTypes;
        this.atomIds = newAtomIds;
        this.capacity = newCapacity;
        this.isDirty = true;
    }

    setBonds(bonds) {
        this.bonds = bonds || [];
        this.isDirty = true;
    }
}

export class MoleculeRenderer extends MeshRenderer {
    constructor(scene, system, shaders, mmParams, source = null) {
        super(scene, shaders, system.capacity);
        this.system = system;
        this.source = source;
        this.mmParams = mmParams;

        this.axesHelper = null;

        this.labelMode = 'none';

        // selection is handled by MeshRenderer
    }

    update() {
        // Main update loop
        if (this.source && (this.source.dirtyExport || this.source.dirtyTopo || this.source.dirtyGeom)) {
            this.source.exportToMoleculeSystem(this.system);
        }
        if (this.system && (this.system.capacity | 0) > (this.capacity | 0)) {
            this.ensureCapacity(this.system.capacity);
        }
        if (this.system.isDirty) {
            this.updateStructure();
            this.updatePositions(); // Structure change implies position update too
            this.system.isDirty = false;
        }

        // Update Labels (visibility check handled in parent)
        // But we need to update content if mode changes or structure changes.
        // Here we just ensure visibility logic if needed.
    }

    updatePositions() {
        super.updatePositions(this.system.pos, this.system.nAtoms);
    }

    updateStructure() {
        const nAtoms = this.system.nAtoms;

        // Update Particles
        this.updateParticles(nAtoms, (i) => {
            const type = this.system.types[i];
            // Use MMParams
            let col = [1, 0, 1]; // Default magenta
            if (this.mmParams) {
                col = this.mmParams.getColor(type);
            }
            return col;
        }, (i) => {
            const type = this.system.types[i];
            let radius = 1.0;
            if (this.mmParams) {
                radius = this.mmParams.getRadius(type) * 0.4;
            }
            return radius;
        });

        // Update Bonds
        this.updateBonds(this.system.bonds);

        this.updateSelection();
        this.updateLabelsContent(); // Update labels when structure changes

        if (window.VERBOSITY_LEVEL >= Logger.INFO) {
            // window.logger.info(`Renderer structure updated: ${ nAtoms } atoms.`);
        }
    }

    updateSelection() {
        if (this.source && this.source.dirtyExport) {
            this.source.exportToMoleculeSystem(this.system);
        }
        const selectedIDs = Array.from(this.system.selection);
        super.updateSelection(selectedIDs);
    }

    updateLabelsContent() {
        if (!this.labelMesh || this.labelMode === 'none') return;

        const nAtoms = this.system.nAtoms;

        this.updateLabels((i) => {
            if (this.labelMode === 'id') {
                return i.toString();
            } else if (this.labelMode === 'element') {
                const type = this.system.types[i];
                if (this.mmParams && this.mmParams.byAtomicNumber[type]) {
                    return this.mmParams.byAtomicNumber[type].name;
                }
                return type.toString();
            } else if (this.labelMode === 'type') {
                return "?"; // Placeholder as per previous logic
            }
            return "";
        }, nAtoms);
    }

    toggleAxes(visible) {
        if (!this.axesHelper) {
            this.axesHelper = new THREE.AxesHelper(5); // 5 units length
            this.scene.add(this.axesHelper);
        }
        this.axesHelper.visible = visible;
    }

    setLabelMode(mode) {
        if (this.labelMode !== mode) {
            this.labelMode = mode;
            if (this.labelMesh) {
                this.labelMesh.visible = (mode !== 'none');
                if (this.labelMesh.visible) {
                    this.updateLabelsContent();
                }
            }
        }
    }

    setLabelStyle(colorHex, scale) {
        if (this.labelMesh) {
            if (colorHex) {
                this.labelMesh.material.uniforms.uColor.value.set(colorHex);
            }
            if (scale !== undefined && scale !== null) {
                this.labelMesh.material.uniforms.uScale.value = parseFloat(scale);
            }
        }
    }
}

