
class MoleculeRenderer {
    constructor(scene, system, shaders) {
        this.scene = scene;
        this.system = system;
        this.shaders = shaders; // Expected to contain { atom: {vertex, fragment}, bond: {vertex, fragment}, selection: {vertex, fragment} }

        this.atomMesh = null;
        this.bondLines = null;
        this.selectionMesh = null;

        this.posTexture = null;
        this.posData = null; // Float32Array for texture
        this.texSize = 1024; // Fixed width for simplicity
        this.texHeight = 0;

        // CPK Colors (Simple map)
        this.colors = {
            1: [1.0, 1.0, 1.0], // H
            6: [0.2, 0.2, 0.2], // C
            7: [0.1, 0.3, 1.0], // N
            8: [1.0, 0.1, 0.1], // O
            16: [1.0, 1.0, 0.3], // S
            // Default
            0: [1.0, 0.0, 1.0]
        };

        this.init();
    }

    init() {
        // --- 1. Position Texture ---
        // Allocate texture large enough for max capacity
        // We use a square texture or just Nx1?
        // Max capacity is usually e.g. 10000. 1024x1024 is 1M atoms. Plenty.
        // Let's use 1024 width.
        this.texHeight = Math.ceil(this.system.capacity / this.texSize);
        const size = this.texSize * this.texHeight * 4; // RGBA float
        this.posData = new Float32Array(size);

        this.posTexture = new THREE.DataTexture(
            this.posData,
            this.texSize,
            this.texHeight,
            THREE.RGBAFormat,
            THREE.FloatType
        );
        this.posTexture.minFilter = THREE.NearestFilter;
        this.posTexture.magFilter = THREE.NearestFilter;
        this.posTexture.needsUpdate = true;

        const commonUniforms = {
            uPosTex: { value: this.posTexture },
            uTexSize: { value: new THREE.Vector2(this.texSize, this.texHeight) }
        };

        // --- 2. Atoms (InstancedMesh) ---
        // Geometry: A simple plane (quad)
        const atomGeo = new THREE.PlaneBufferGeometry(1, 1);
        this.atomMesh = Draw3D.createTextureBasedInstancedMesh(
            this.system.capacity,
            atomGeo,
            this.shaders.atom,
            commonUniforms
        );

        // Add extra attributes for color and scale
        const colors = new Float32Array(this.system.capacity * 3);
        const scales = new Float32Array(this.system.capacity);
        this.atomMesh.geometry.setAttribute('instanceColor', new THREE.InstancedBufferAttribute(colors, 3));
        this.atomMesh.geometry.setAttribute('instanceScale', new THREE.InstancedBufferAttribute(scales, 1));

        this.scene.add(this.atomMesh);

        // --- 3. Bonds (LineSegments) ---
        const maxBonds = this.system.capacity * 4;
        this.bondLines = Draw3D.createTextureBasedLineSegments(
            this.system.capacity,
            maxBonds,
            this.shaders.bond,
            {
                ...commonUniforms,
                uColor: { value: new THREE.Vector4(1.0, 1.0, 1.0, 1.0) } // White
            }
        );
        this.scene.add(this.bondLines);

        // --- 4. Selection (InstancedMesh - Rings) ---
        this.selectionBaseGeo = Draw3D.createOctSphereGeometry(16, 1.0);
        this.selectionMesh = Draw3D.createTextureBasedSelectionLines(
            this.system.capacity,
            this.selectionBaseGeo,
            this.shaders.selection,
            {
                ...commonUniforms,
                uColor: { value: new THREE.Vector4(1.0, 1.0, 0.0, 0.8) } // Yellow
            }
        );
        this.scene.add(this.selectionMesh);
    }

    update() {
        // Main update loop
        // If system is dirty, we assume topology or colors changed -> Full Update
        // If only positions moved (handled by Gizmo), we call updatePositions() explicitly
        // But for safety, if system.isDirty is true, we do everything.

        if (this.system.isDirty) {
            this.updateStructure();
            this.updatePositions(); // Structure change implies position update too
            this.system.isDirty = false;
        }
    }

    updatePositions() {
        const nAtoms = this.system.nAtoms;
        // Update Texture
        for (let i = 0; i < nAtoms; i++) {
            this.posData[i * 4] = this.system.pos[i * 3];
            this.posData[i * 4 + 1] = this.system.pos[i * 3 + 1];
            this.posData[i * 4 + 2] = this.system.pos[i * 3 + 2];
            this.posData[i * 4 + 3] = 1.0; // w
        }
        this.posTexture.needsUpdate = true;

        // Selection might need update if it depends on positions? 
        // No, selection depends on IDs which index into texture. 
        // So if texture updates, selection updates visually.
    }

    updateStructure() {
        const nAtoms = this.system.nAtoms;
        this.atomMesh.count = nAtoms;

        const colorAttr = this.atomMesh.geometry.getAttribute('instanceColor');
        const scaleAttr = this.atomMesh.geometry.getAttribute('instanceScale');

        for (let i = 0; i < nAtoms; i++) {
            // Color & Scale
            const type = this.system.types[i];
            let col = this.colors[type] || this.colors[0];
            colorAttr.setXYZ(i, col[0], col[1], col[2]);

            let scale = 1.0;
            if (type === 1) scale = 0.6;
            else if (type === 6) scale = 1.4;
            else if (type === 8) scale = 1.3;

            scaleAttr.setX(i, scale);
        }
        colorAttr.needsUpdate = true;
        scaleAttr.needsUpdate = true;

        // Update Bonds
        const bonds = this.system.bonds;
        const bondIDAttr = this.bondLines.geometry.getAttribute('aAtomID');
        let ptr = 0;
        for (let i = 0; i < bonds.length; i++) {
            bondIDAttr.setX(ptr++, bonds[i][0]);
            bondIDAttr.setX(ptr++, bonds[i][1]);
        }
        bondIDAttr.needsUpdate = true;
        this.bondLines.geometry.setDrawRange(0, bonds.length * 2);

        this.updateSelection();

        if (window.VERBOSITY_LEVEL >= Logger.INFO) {
            // window.logger.info(`Renderer structure updated: ${ nAtoms } atoms.`);
        }
    }

    updateSelection() {
        const selectedIDs = Array.from(this.system.selection);
        const count = selectedIDs.length;

        // Update instance count for InstancedBufferGeometry
        this.selectionMesh.geometry.instanceCount = count;

        const idAttr = this.selectionMesh.geometry.getAttribute('aAtomID');
        for (let i = 0; i < count; i++) {
            idAttr.setX(i, selectedIDs[i]);
        }
        idAttr.needsUpdate = true;
    }
}

