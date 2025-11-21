class ShortcutManager {
    constructor(editor) {
        this.editor = editor;
        this.shortcuts = {};
        this.init();
    }

    init() {
        window.addEventListener('keydown', this.onKeyDown.bind(this));

        // Define default shortcuts
        this.register('g', () => this.editor.toggleGizmo(), 'Toggle Gizmo');
        this.register('t', () => this.editor.setGizmoMode('translate'), 'Translate Mode');
        this.register('r', () => this.editor.setGizmoMode('rotate'), 'Rotate Mode');
        this.register('s', () => this.editor.setGizmoMode('scale'), 'Scale Mode');
        this.register('Escape', () => this.editor.clearSelection(), 'Clear Selection');
        // this.register('Delete', () => this.editor.deleteSelection(), 'Delete Selection');
        // this.register('Backspace', () => this.editor.deleteSelection(), 'Delete Selection');
    }

    register(key, action, description) {
        this.shortcuts[key.toLowerCase()] = { action, description };
    }

    onKeyDown(e) {
        // Ignore if typing in an input field
        if (e.target.tagName === 'INPUT' || e.target.tagName === 'TEXTAREA') return;

        const key = e.key.toLowerCase();
        if (this.shortcuts[key]) {
            this.shortcuts[key].action();
            // window.logger.info(`Shortcut: ${this.shortcuts[key].description}`);
        }
    }
}
