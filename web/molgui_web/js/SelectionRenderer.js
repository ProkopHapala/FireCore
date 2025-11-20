class SelectionRenderer {
    constructor(scene, system) {
        this.scene = scene;
        this.system = system;
        this.mesh = null;

        // Pre-calculate circle vertices (unit radius)
        this.circleVertices = this.generateCircleVertices();

        this.init();
    }

    generateCircleVertices() {
        const vertices = [];
        const segments = 32;

        // XY Circle
        for (let i = 0; i < segments; i++) {
            const theta1 = (i / segments) * Math.PI * 2;
            const theta2 = ((i + 1) / segments) * Math.PI * 2;
            vertices.push(Math.cos(theta1), Math.sin(theta1), 0);
            vertices.push(Math.cos(theta2), Math.sin(theta2), 0);
        }

        // XZ Circle
        for (let i = 0; i < segments; i++) {
            const theta1 = (i / segments) * Math.PI * 2;
            const theta2 = ((i + 1) / segments) * Math.PI * 2;
            vertices.push(Math.cos(theta1), 0, Math.sin(theta1));
            vertices.push(Math.cos(theta2), 0, Math.sin(theta2));
        }

        // YZ Circle
        for (let i = 0; i < segments; i++) {
            const theta1 = (i / segments) * Math.PI * 2;
            const theta2 = ((i + 1) / segments) * Math.PI * 2;
            vertices.push(0, Math.cos(theta1), Math.sin(theta1));
            vertices.push(0, Math.cos(theta2), Math.sin(theta2));
        }

        return vertices;
    }

    init() {
        const geometry = new THREE.BufferGeometry();
        const material = new THREE.LineBasicMaterial({
            color: 0xffff00, // Yellow/Gold
            depthTest: false, // Always visible on top? Maybe true is better for 3D feel
            depthWrite: false,
            transparent: true,
            opacity: 0.8
        });

        this.mesh = new THREE.LineSegments(geometry, material);
        // Ensure it renders on top if we want, or just normally
        // this.mesh.renderOrder = 999; 

        this.scene.add(this.mesh);
    }

    update() {
        const selectedIDs = Array.from(this.system.selection);
        if (selectedIDs.length === 0) {
            this.mesh.visible = false;
            return;
        }

        this.mesh.visible = true;

        // Rebuild geometry
        // Total vertices = selected * circleVertices.length
        const totalVertices = selectedIDs.length * this.circleVertices.length;
        const positions = new Float32Array(totalVertices);

        let ptr = 0;
        const radius = 1.0; // Increased radius to be clearly visible outside atoms

        for (const id of selectedIDs) {
            const x = this.system.pos[id * 3];
            const y = this.system.pos[id * 3 + 1];
            const z = this.system.pos[id * 3 + 2];

            // Get atom type to adjust radius? 
            // For now fixed radius is fine, or maybe scale by atom size?
            // Let's use fixed slightly larger radius for consistency.

            for (let i = 0; i < this.circleVertices.length; i += 3) {
                positions[ptr++] = x + this.circleVertices[i] * radius;
                positions[ptr++] = y + this.circleVertices[i + 1] * radius;
                positions[ptr++] = z + this.circleVertices[i + 2] * radius;
            }
        }

        this.mesh.geometry.setAttribute('position', new THREE.BufferAttribute(positions, 3));
        this.mesh.geometry.attributes.position.needsUpdate = true;

        // If count changed significantly, we might need to recreate buffer if it grew too much
        // But BufferAttribute handles resizing if we pass a new array usually?
        // Actually setAttribute with new buffer is safest.
    }
}
