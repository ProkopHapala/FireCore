class MolGUIApp {
    constructor() {
        this.container = document.getElementById('canvas-container');
        this.scene = null;
        this.camera = null;
        this.renderer = null;
        this.controls = null;
    }

    init() {
        window.logger.info("Initializing MolGUI...");

        // 1. Scene
        this.scene = new THREE.Scene();
        this.scene.background = new THREE.Color(0x222222);

        // 2. Camera
        this.camera = new THREE.PerspectiveCamera(60, window.innerWidth / window.innerHeight, 0.1, 1000);
        this.camera.position.set(5, 5, 5);

        // 3. Renderer
        this.renderer = new THREE.WebGLRenderer({
            antialias: true,
            powerPreference: "high-performance"
        });
        this.renderer.setSize(window.innerWidth, window.innerHeight);
        this.renderer.setPixelRatio(window.devicePixelRatio);
        this.container.appendChild(this.renderer.domElement);

        // 4. Controls
        if (THREE.OrbitControls) {
            this.controls = new THREE.OrbitControls(this.camera, this.renderer.domElement);
            this.controls.enableDamping = false; // Disable inertia for raw performance feel
            this.controls.dampingFactor = 0.05;

            // Unity Style Controls:
            // LMB: Selection (Handled by Editor) -> Set to NULL here
            // MMB: Zoom/Dolly (Standard)
            // RMB: Rotate (Standard)
            // Shift + RMB: Pan (Custom)

            this.controls.mouseButtons = {
                LEFT: null, // Disable Left Click for Camera
                MIDDLE: THREE.MOUSE.DOLLY,
                RIGHT: THREE.MOUSE.ROTATE
            };

            // Listen for Shift Key to toggle Pan Mode on RMB
            window.addEventListener('keydown', (e) => {
                if (e.key === 'Shift') {
                    this.controls.mouseButtons.RIGHT = THREE.MOUSE.PAN;
                }
            });

            window.addEventListener('keyup', (e) => {
                if (e.key === 'Shift') {
                    this.controls.mouseButtons.RIGHT = THREE.MOUSE.ROTATE;
                }
            });
        } else {
            window.logger.error("OrbitControls not loaded!");
        }

        // 5. Molecule System
        this.system = new MoleculeSystem();
        this.molRenderer = new MoleculeRenderer(this.scene, this.system);

        // 6. I/O and GUI
        this.io = new IO(this.system, this.molRenderer);
        this.gui = new GUI(this.io);

        // 8. Selection Renderer
        this.selectionRenderer = new SelectionRenderer(this.scene, this.system);

        // 7. Editor (Selection, Gizmo)
        this.editor = new Editor(this.scene, this.camera, this.renderer, this.system, this.molRenderer, this.selectionRenderer);
        this.editor.onSelectionChange = () => {
            this.gui.updateSelectionCount();
        };

        // 9. Shortcut Manager
        this.shortcuts = new ShortcutManager(this.editor);

        // Hook into Editor or System updates?
        // Ideally System should dispatch events, but we can just patch the update method or call it manually.
        // For now, let's monkey-patch the renderer update or just call it in animate loop?
        // Better: Pass it to Editor? Or make Editor call a global update?
        // Let's make Editor call it.
        this.editor.onSelectionChange = () => {
            this.selectionRenderer.update();
            this.gui.updateSelectionCount();
        };

        // Also hook GUI input back to renderer
        this.gui.onSelectionChanged = () => {
            this.selectionRenderer.update();
        };

        // Create Test Molecule (Water: H-O-H)
        // O at 0,0,0
        const o = this.system.addAtom(0, 0, 0, 8);
        // H at 0.8, 0.6, 0
        const h1 = this.system.addAtom(0.8, 0.6, 0, 1);
        // H at -0.8, 0.6, 0
        const h2 = this.system.addAtom(-0.8, 0.6, 0, 1);

        this.system.addBond(o, h1);
        this.system.addBond(o, h2);

        // Add some more atoms to test performance/visuals
        // Methane-like structure nearby
        const c = this.system.addAtom(3, 0, 0, 6);
        const h3 = this.system.addAtom(3, 1, 0, 1);
        const h4 = this.system.addAtom(3, -0.5, 0.8, 1);
        const h5 = this.system.addAtom(3, -0.5, -0.8, 1);

        this.system.addBond(c, h3);
        this.system.addBond(c, h4);
        this.system.addBond(c, h5);

        this.molRenderer.update();

        // Lights
        const ambientLight = new THREE.AmbientLight(0xffffff, 0.5);
        this.scene.add(ambientLight);
        const dirLight = new THREE.DirectionalLight(0xffffff, 1);
        dirLight.position.set(5, 10, 7);
        this.scene.add(dirLight);

        // Events
        window.addEventListener('resize', this.onWindowResize.bind(this));

        // Start Loop
        this.animate();

        window.logger.info("Initialization Complete.");
    }

    onWindowResize() {
        this.camera.aspect = window.innerWidth / window.innerHeight;
        this.camera.updateProjectionMatrix();
        this.renderer.setSize(window.innerWidth, window.innerHeight);
    }

    animate() {
        requestAnimationFrame(this.animate.bind(this));

        if (this.controls) this.controls.update();

        this.renderer.render(this.scene, this.camera);
    }
}

// Start
window.onload = () => {
    const app = new MolGUIApp();
    window.app = app; // Expose for debugging
    app.init();
};
