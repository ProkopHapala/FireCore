class Editor {
    constructor(scene, camera, renderer, system, molRenderer) {
        this.scene = scene;
        this.camera = camera;
        this.renderer = renderer;
        this.system = system;
        this.molRenderer = molRenderer;

        this.raycaster = new THREE.Raycaster();
        this.mouse = new THREE.Vector2();

        this.isDragging = false;
        this.startPos = { x: 0, y: 0 };
        this.selectionBox = document.createElement('div');
        this.selectionBox.className = 'selection-box';
        document.body.appendChild(this.selectionBox);

        window.editorInitialized = true;
        this.initEvents();
    }

    initEvents() {
        const canvas = this.renderer.domElement;
        console.log(`[Editor] InitEvents. Window Size: ${window.innerWidth}x${window.innerHeight}`);

        // Use capture: true for ALL events to ensure we intercept them before OrbitControls/Canvas
        // We attach to window for move/up to handle dragging outside canvas
        // Also attach mousedown to window to be absolutely sure
        // Use POINTER events to be more robust
        document.addEventListener('pointerdown', this.onMouseDown.bind(this), { capture: true });
        window.addEventListener('pointermove', this.onMouseMove.bind(this), { capture: true });
        window.addEventListener('pointerup', this.onMouseUp.bind(this), { capture: true });

        window.editorEventsAttached = true;
    }

    onMouseDown(e) {
        console.log(`[Editor] MouseDown: Button=${e.button} Client=(${e.clientX}, ${e.clientY})`);
        window.lastMouseDownEvent = {
            button: e.button,
            clientX: e.clientX,
            clientY: e.clientY,
            time: Date.now(),
            target: e.target.tagName
        };
        // if (e.button !== 0) return; // ALLOW ALL BUTTONS FOR DEBUGGING

        this.isDragging = true;
        this.startPos.x = e.clientX;
        this.startPos.y = e.clientY;

        // Reset box
        this.selectionBox.style.display = 'none';
        this.selectionBox.style.width = '0';
        this.selectionBox.style.height = '0';
    }

    onMouseMove(e) {
        if (!this.isDragging) return;

        const dx = e.clientX - this.startPos.x;
        const dy = e.clientY - this.startPos.y;

        // If moved enough, show box (Marquee)
        if (Math.abs(dx) > 5 || Math.abs(dy) > 5) {
            this.selectionBox.style.display = 'block';

            const x = dx < 0 ? e.clientX : this.startPos.x;
            const y = dy < 0 ? e.clientY : this.startPos.y;

            this.selectionBox.style.left = x + 'px';
            this.selectionBox.style.top = y + 'px';
            this.selectionBox.style.width = Math.abs(dx) + 'px';
            this.selectionBox.style.height = Math.abs(dy) + 'px';
        }
    }

    onMouseUp(e) {
        console.log(`[Editor] MouseUp: Button=${e.button} Dragging=${this.isDragging}`);
        window.lastMouseUpEvent = {
            button: e.button,
            clientX: e.clientX,
            clientY: e.clientY,
            time: Date.now(),
            wasDragging: this.isDragging
        };
        if (!this.isDragging) return;
        this.isDragging = false;
        this.selectionBox.style.display = 'none';

        const dx = e.clientX - this.startPos.x;
        const dy = e.clientY - this.startPos.y;

        // Determine Mode
        let mode = 'replace';
        if (e.shiftKey) mode = 'add';
        else if (e.ctrlKey || e.altKey) mode = 'subtract';

        if (Math.abs(dx) < 5 && Math.abs(dy) < 5) {
            // Click (Single Pick)
            console.log("[Editor] Click detected. Calling pick()...");
            this.pick(e.clientX, e.clientY, mode);
        } else {
            // Drag (Box Select)
            console.log("[Editor] Drag detected. Calling boxSelect()...");
            this.boxSelect(this.startPos.x, this.startPos.y, e.clientX, e.clientY, mode);
        }
    }

    pick(x, y, mode) {
        // Get Ray from camera
        const ray = this.getRay(x, y);
        const rayOrigin = ray.origin;
        const rayDir = ray.direction;

        // DEBUG LOGGING
        let debugLog = "--- PICK DEBUG ---\n";
        debugLog += `Mouse: (${x}, ${y})\n`;
        debugLog += `Ray Origin: ${rayOrigin.x.toFixed(3)}, ${rayOrigin.y.toFixed(3)}, ${rayOrigin.z.toFixed(3)}\n`;
        debugLog += `Ray Dir:    ${rayDir.x.toFixed(3)}, ${rayDir.y.toFixed(3)}, ${rayDir.z.toFixed(3)}\n`;
        debugLog += `Camera Pos: ${this.camera.position.x.toFixed(3)}, ${this.camera.position.y.toFixed(3)}, ${this.camera.position.z.toFixed(3)}\n`;
        debugLog += `Camera Zoom: ${this.camera.zoom}\n`;

        // Log matrices (first few elements)
        const mw = this.camera.matrixWorld.elements;
        debugLog += `MatWorld: [${mw[0].toFixed(2)}, ${mw[1].toFixed(2)}, ..., ${mw[12].toFixed(2)}, ${mw[13].toFixed(2)}, ${mw[14].toFixed(2)}]\n`;

        let minT = Infinity;
        let pickedId = -1;

        // Loop through all atoms
        for (let i = 0; i < this.system.nAtoms; i++) {
            const ax = this.system.pos[i * 3];
            const ay = this.system.pos[i * 3 + 1];
            const az = this.system.pos[i * 3 + 2];
            const center = new THREE.Vector3(ax, ay, az);

            const radius = 0.5;

            const t = this.raySphere(rayOrigin, rayDir, center, radius);

            debugLog += `Atom ${i}: Pos(${ax.toFixed(3)}, ${ay.toFixed(3)}, ${az.toFixed(3)}) t=${t.toFixed(3)}\n`;

            if (t < minT && t > 0) {
                minT = t;
                pickedId = i;
            }
        }
        debugLog += `Picked ID: ${pickedId}, minT: ${minT}\n`;
        debugLog += "------------------";

        console.log(debugLog);
        window.lastPickDebug = debugLog;
        const debugEl = document.getElementById('debug-output');
        if (debugEl) debugEl.innerText = debugLog;

        if (pickedId !== -1) {
            this.system.select(pickedId, mode);
            this.molRenderer.update();
            if (this.onSelectionChange) this.onSelectionChange();
            window.logger.info(`Selected Atom ${pickedId} (t=${minT.toFixed(2)})`);
            return;
        }

        // If clicked nothing
        if (mode === 'replace') {
            this.system.clearSelection();
            this.molRenderer.update();
            if (this.onSelectionChange) this.onSelectionChange();
        }
    }

    getRay(x, y) {
        // Normalize mouse coordinates
        const ndcX = (x / window.innerWidth) * 2 - 1;
        const ndcY = -(y / window.innerHeight) * 2 + 1;

        const rayOrigin = new THREE.Vector3();
        const rayDir = new THREE.Vector3();

        if (this.camera.isPerspectiveCamera) {
            rayOrigin.setFromMatrixPosition(this.camera.matrixWorld);
            rayDir.set(ndcX, ndcY, 0.5).unproject(this.camera).sub(rayOrigin).normalize();
        } else if (this.camera.isOrthographicCamera) {
            rayOrigin.set(ndcX, ndcY, -1).unproject(this.camera);
            rayDir.set(0, 0, -1).transformDirection(this.camera.matrixWorld);
        }

        return { origin: rayOrigin, direction: rayDir };
    }

    // Custom Ray-Sphere Intersection (Geometric)
    // Returns distance t, or Infinity if no hit
    raySphere(rayOrigin, rayDir, center, radius) {
        const oc = new THREE.Vector3().subVectors(center, rayOrigin);
        const tca = oc.dot(rayDir);
        if (tca < 0) return Infinity; // Behind ray

        const d2 = oc.dot(oc) - tca * tca;
        const r2 = radius * radius;

        if (d2 > r2) return Infinity; // Miss

        const thc = Math.sqrt(r2 - d2);
        const t0 = tca - thc;
        // const t1 = tca + thc; // Far hit

        return t0;
    }

    boxSelect(x1, y1, x2, y2, mode) {
        const minX = Math.min(x1, x2);
        const maxX = Math.max(x1, x2);
        const minY = Math.min(y1, y2);
        const maxY = Math.max(y1, y2);

        const width = window.innerWidth;
        const height = window.innerHeight;

        if (mode === 'replace') this.system.clearSelection();

        let count = 0;
        const vec = new THREE.Vector3();

        for (let i = 0; i < this.system.nAtoms; i++) {
            vec.set(
                this.system.pos[i * 3],
                this.system.pos[i * 3 + 1],
                this.system.pos[i * 3 + 2]
            );

            // Project to screen
            vec.project(this.camera);

            // Map to pixels
            const sx = (vec.x * 0.5 + 0.5) * width;
            const sy = (-(vec.y * 0.5) + 0.5) * height;

            // Check bounds
            if (sx >= minX && sx <= maxX && sy >= minY && sy <= maxY) {
                // Check Z to ensure it's in front of camera? (optional, usually project handles it but check z < 1)
                if (vec.z < 1.0) {
                    if (mode === 'subtract') this.system.selection.delete(i);
                    else this.system.selection.add(i);
                    count++;
                }
            }
        }

        this.system.isDirty = true;
        this.molRenderer.update();
        if (this.onSelectionChange) this.onSelectionChange();
        window.logger.info(`Box Selected ${count} atoms.`);
    }
}
