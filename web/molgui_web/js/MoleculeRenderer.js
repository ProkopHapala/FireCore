class MoleculeRenderer {
    constructor(scene, system) {
        this.scene = scene;
        this.system = system;

        this.atomMesh = null;
        this.bondLines = null;

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
        // --- Atoms (InstancedMesh) ---
        // Geometry: A simple plane (quad)
        const geometry = new THREE.PlaneBufferGeometry(1, 1);

        // Material: Custom Shader
        const material = new THREE.ShaderMaterial({
            vertexShader: AtomVertexShader,
            fragmentShader: AtomFragmentShader,
            uniforms: {},
            side: THREE.DoubleSide
        });

        this.atomMesh = new THREE.InstancedMesh(geometry, material, this.system.capacity);
        this.atomMesh.instanceMatrix.setUsage(THREE.DynamicDrawUsage); // We will update this often

        // Add extra attributes for color and scale
        // We need to allocate buffers for the max capacity
        const colors = new Float32Array(this.system.capacity * 3);
        const scales = new Float32Array(this.system.capacity);

        this.atomMesh.geometry.setAttribute('instanceColor', new THREE.InstancedBufferAttribute(colors, 3));
        this.atomMesh.geometry.setAttribute('instanceScale', new THREE.InstancedBufferAttribute(scales, 1));

        this.scene.add(this.atomMesh);

        // --- Bonds (LineSegments) ---
        // Initial empty geometry
        const lineGeo = new THREE.BufferGeometry();
        const lineMat = new THREE.LineBasicMaterial({ color: 0xffffff });
        this.bondLines = new THREE.LineSegments(lineGeo, lineMat);
        this.scene.add(this.bondLines);
    }

    update() {
        if (!this.system.isDirty) return;

        const nAtoms = this.system.nAtoms;
        this.atomMesh.count = nAtoms;

        const dummy = new THREE.Object3D();
        const colorAttr = this.atomMesh.geometry.getAttribute('instanceColor');
        const scaleAttr = this.atomMesh.geometry.getAttribute('instanceScale');

        for (let i = 0; i < nAtoms; i++) {
            // 1. Position
            dummy.position.set(
                this.system.pos[i * 3],
                this.system.pos[i * 3 + 1],
                this.system.pos[i * 3 + 2]
            );
            dummy.updateMatrix();
            this.atomMesh.setMatrixAt(i, dummy.matrix);

            // 2. Color & Scale
            const type = this.system.types[i];
            let col = this.colors[type] || this.colors[0];

            colorAttr.setXYZ(i, col[0], col[1], col[2]);

            // Simple radius mapping
            let scale = 1.0; // Default diameter
            if (type === 1) scale = 0.6; // H
            else if (type === 6) scale = 1.4; // C
            else if (type === 8) scale = 1.3; // O

            scaleAttr.setX(i, scale);
        }

        this.atomMesh.instanceMatrix.needsUpdate = true;
        colorAttr.needsUpdate = true;
        scaleAttr.needsUpdate = true;

        // --- Update Bonds ---
        // Rebuild vertex buffer for lines
        const bonds = this.system.bonds;
        const positions = [];

        for (let i = 0; i < bonds.length; i++) {
            const id1 = bonds[i][0];
            const id2 = bonds[i][1];

            positions.push(
                this.system.pos[id1 * 3], this.system.pos[id1 * 3 + 1], this.system.pos[id1 * 3 + 2],
                this.system.pos[id2 * 3], this.system.pos[id2 * 3 + 1], this.system.pos[id2 * 3 + 2]
            );
        }

        this.bondLines.geometry.setAttribute('position', new THREE.Float32BufferAttribute(positions, 3));

        this.system.isDirty = false;
        window.logger.info(`Renderer updated: ${nAtoms} atoms, ${bonds.length} bonds.`);
    }
}
