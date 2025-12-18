export class ShortcutManager {
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
        this.register('Delete', () => this.editor.deleteSelection(), 'Delete Selection');
        this.register('Backspace', () => this.editor.deleteSelection(), 'Delete Selection');
        this.register('b', () => this.editor.recalculateBonds(), 'Recalculate Bonds');
        this.register('a', () => this.editor.addAtom(), 'Add Atom');
    }

    register(key, action, description) {
        this.shortcuts[key.toLowerCase()] = { action, description };
    }

    onKeyDown(e) {
        // Ignore only when user is actively typing/editing text.
        // Using document.activeElement is more robust than event target,
        // because listeners are on window.
        const ae = document.activeElement;
        const tag = ae ? ae.tagName : '';
        const isEditable = (ae && (ae.isContentEditable || tag === 'INPUT' || tag === 'TEXTAREA' || tag === 'SELECT'));
        if (isEditable) return;

        const key = e.key.toLowerCase();

        window.logger.debug(`[ShortcutManager] Key pressed: '${e.key}' (mapped to '${key}')`);

        if (this.shortcuts[key]) {
            window.logger.info(`[ShortcutManager] Executing action: ${this.shortcuts[key].description}`);
            // Prevent browser default actions (esp. Backspace navigation)
            if (key === 'delete' || key === 'backspace') e.preventDefault();
            this.shortcuts[key].action();
        } else {
            window.logger.debug(`[ShortcutManager] No shortcut found for '${key}'`);
        }
    }
}
