"use strict";



export const GUIutils = {

    el: function(parent, tag, props = null, style = null) {
        const el = document.createElement(tag);
        if (props) {
            if (props.className !== undefined) el.className = props.className;
            if (props.text !== undefined) el.textContent = props.text;
            if (props.html !== undefined) el.innerHTML = props.html;
            if (props.type !== undefined) el.type = props.type;
            if (props.value !== undefined) el.value = props.value;
            if (props.placeholder !== undefined) el.placeholder = props.placeholder;
            if (props.checked !== undefined) el.checked = props.checked;
            if (props.min !== undefined) el.min = props.min;
            if (props.max !== undefined) el.max = props.max;
            if (props.step !== undefined) el.step = props.step;
            if (props.attrs) for (const k in props.attrs) el.setAttribute(k, props.attrs[k]);
            if (props.onchange) el.onchange = props.onchange;
            if (props.oninput) el.oninput = props.oninput;
            if (props.onclick) el.onclick = props.onclick;
        }
        if (style) for (const k in style) el.style[k] = style[k];
        if (parent) parent.appendChild(el);
        return el;
    },

    div: function(parent, className = null, style = null) { return this.el(parent, 'div', className ? { className } : null, style); },
    row: function(parent, style = null) { return this.el(parent, 'div', { className: 'gui-row' }, style); },
    span: function(parent, text = '', style = null) { return this.el(parent, 'span', { text }, style); },
    btn: function(parent, text, onclick = null, style = null, className = 'gui-btn') { return this.el(parent, 'button', { className, text, onclick }, style); },
    input: function(parent, opts = null, style = null) { return this.el(parent, 'input', opts, style); },
    num: function(parent, value, opts = null, style = null) {
        const o = opts ? { ...opts } : {};
        o.type = 'number';
        if (value !== undefined) o.value = String(value);
        if (!o.className) o.className = 'gui-input';
        return this.el(parent, 'input', o, style);
    },
    textInput: function(parent, value, opts = null, style = null) {
        const o = opts ? { ...opts } : {};
        o.type = 'text';
        if (value !== undefined) o.value = String(value);
        if (!o.className) o.className = 'gui-input';
        return this.el(parent, 'input', o, style);
    },
    range: function(parent, value, min, max, step, oninput = null, style = null) {
        return this.el(parent, 'input', { type: 'range', className: 'gui-input', value: String(value), min: String(min), max: String(max), step: String(step), oninput }, style);
    },
    selectList: function(parent, items, selected = null, onchange = null, style = null, className = 'gui-select') {
        const sel = this.el(parent, 'select', { className, onchange }, style);
        for (const it of items) {
            const opt = document.createElement('option');
            opt.value = it;
            opt.textContent = it;
            if (selected !== null && it === selected) opt.selected = true;
            sel.appendChild(opt);
        }
        return sel;
    },

    setSelectOptions: function(sel, items, opts = null) {
        sel.innerHTML = '';
        const selectedValue = opts && (opts.selectedValue !== undefined) ? opts.selectedValue : undefined;
        const selectFirst = opts && opts.selectFirst;
        let didSelect = false;
        for (let i = 0; i < items.length; i++) {
            const it = items[i];
            const opt = document.createElement('option');
            opt.value = (it && (it.value !== undefined)) ? it.value : String(it);
            opt.textContent = (it && (it.text !== undefined)) ? it.text : String(it);
            const shouldSelect = (it && it.selected) || (selectedValue !== undefined && opt.value == selectedValue);
            if (shouldSelect && !didSelect) { opt.selected = true; didSelect = true; }
            sel.appendChild(opt);
        }
        if (!didSelect && selectFirst && sel.options.length > 0) sel.options[0].selected = true;
        return sel;
    },
    labelCheck: function(parent, text, checked = false, onchange = null, style = null, className = 'gui-checkbox-label') {
        const lbl = this.el(parent, 'label', { className }, style);
        const chk = this.el(lbl, 'input', { type: 'checkbox', checked, onchange });
        lbl.appendChild(document.createTextNode(text));
        return { label: lbl, input: chk };
    },

    /**
     * Create a container div with flex layout
     */
    box: function(parent, vertical = true) {
        const div = document.createElement('div');
        div.style.display = 'flex';
        div.style.flexDirection = vertical ? 'column' : 'row';
        div.style.gap = '4px';
        if (parent) parent.appendChild(div);
        return div;
    },

    /**
     * Create a label
     */
    label: function(parent, text) {
        const el = document.createElement('span');
        el.textContent = text;
        el.className = 'compact-label'; // Use existing class if available, or style inline
        if (parent) parent.appendChild(el);
        return el;
    },

    /**
     * Create a number input (spin box)
     */
    spinBox: function(parent, value, min, max, step, decimals = 2, callback = null) {
        const input = document.createElement('input');
        input.type = 'number';
        input.value = value;
        if (min !== undefined) input.min = min;
        if (max !== undefined) input.max = max;
        if (step !== undefined) input.step = step;
        input.className = 'compact-input';
        
        // Style for consistency
        input.style.width = '60px';

        if (callback) {
            input.addEventListener('input', () => callback(parseFloat(input.value)));
            input.addEventListener('change', () => callback(parseFloat(input.value)));
        }

        if (parent) parent.appendChild(input);
        return input;
    },

    /**
     * Create a textarea with reasonable defaults
     */
    textArea: function(parent, value = '', opts = null) {
        const ta = document.createElement('textarea');
        ta.value = value;
        ta.className = (opts && opts.className) ? opts.className : 'gui-textarea';
        ta.style.display = (opts && (opts.display !== undefined)) ? opts.display : 'none';
        ta.style.height = (opts && opts.height) ? opts.height : '100px';
        ta.style.width = (opts && opts.width) ? opts.width : '100%';
        ta.style.marginTop = (opts && (opts.marginTop !== undefined)) ? opts.marginTop : '5px';
        ta.style.fontSize = (opts && opts.fontSize) ? opts.fontSize : '0.8em';
        ta.style.fontFamily = (opts && opts.fontFamily) ? opts.fontFamily : 'monospace';
        if (opts && opts.placeholder) ta.placeholder = opts.placeholder;
        if (parent) parent.appendChild(ta);
        return ta;
    },

    /**
     * Create a text input
     */
    text: function(parent, value, callback = null) {
        const input = document.createElement('input');
        input.type = 'text';
        input.value = value;
        input.className = 'compact-input';
        input.style.width = '80px';
        
        if (callback) {
            input.addEventListener('change', () => callback(input.value));
        }

        if (parent) parent.appendChild(input);
        return input;
    },

    /**
     * Create a checkbox
     */
    checkBox: function(parent, text, checked, callback = null) {
        const label = document.createElement('label');
        label.style.display = 'flex';
        label.style.alignItems = 'center';
        label.style.gap = '4px';
        
        const input = document.createElement('input');
        input.type = 'checkbox';
        input.checked = checked;
        
        if (callback) {
            input.addEventListener('change', () => callback(input.checked));
        }

        label.appendChild(input);
        label.appendChild(document.createTextNode(text));

        if (parent) parent.appendChild(label);
        return { label, input };
    },

    /**
     * Create a button
     */
    button: function(parent, text, callback = null) {
        const btn = document.createElement('button');
        btn.textContent = text;
        btn.className = 'small-btn'; // Use existing class
        
        if (callback) {
            btn.addEventListener('click', callback);
        }

        if (parent) parent.appendChild(btn);
        return btn;
    },

    /**
     * Create a select dropdown
     */
    select: function(parent, options, selected, callback = null) {
        const sel = document.createElement('select');
        sel.className = 'compact-select';

        for (const key in options) {
            const opt = document.createElement('option');
            opt.value = options[key]; // Value
            opt.textContent = key;    // Display text
            if (options[key] === selected) opt.selected = true;
            sel.appendChild(opt);
        }

        if (callback) {
            sel.addEventListener('change', () => callback(sel.value));
        }

        if (parent) parent.appendChild(sel);
        return sel;
    },

    /**
     * Create a group box (fieldset-like)
     */
    group: function(parent, title) {
        const div = document.createElement('div');
        div.className = 'panel-section'; // Use existing class for consistent styling
        
        if (title) {
            const h3 = document.createElement('h3');
            h3.textContent = title;
            h3.style.margin = '2px 0 5px 0';
            h3.style.fontSize = '1em';
            div.appendChild(h3);
        }
        
        if (parent) parent.appendChild(div);
        return div;
    },

    /**
     * Data-driven GUI creation
     * @param {HTMLElement} parent - Container element
     * @param {Object} param_specs - Dictionary of parameter specifications
     * @param {Function} on_change - Global callback when any parameter changes (optional)
     * @returns {Object} - Object containing references to widgets and a getValues() function
     */
    create_gui: function(parent, param_specs, on_change = null) {
        const widgets = {};
        const groups = {};
        
        // Helper to get or create group
        const getGroup = (groupName) => {
            if (!groupName) return parent;
            if (!groups[groupName]) {
                groups[groupName] = this.group(parent, groupName);
            }
            return groups[groupName];
        };

        for (const name in param_specs) {
            const spec = param_specs[name];
            const container = getGroup(spec.group);
            
            // Row for label + widget
            const row = this.box(container, false);
            row.style.alignItems = 'center';
            row.style.justifyContent = 'space-between';
            row.style.marginBottom = '2px';

            this.label(row, name + ":");

            let widget;
            const handleChange = (val) => {
                if (on_change) on_change(name, val);
            };

            if (spec.widget === 'double' || spec.widget === 'int') {
                const decimals = spec.widget === 'int' ? 0 : (spec.decimals || 2);
                widget = this.spinBox(row, spec.value, spec.range[0], spec.range[1], spec.step, decimals, handleChange);
            } else if (spec.widget === 'bool') {
                const res = this.checkBox(row, '', spec.value, handleChange);
                widget = res.input;
            } else if (spec.widget === 'text') {
                widget = this.text(row, spec.value, handleChange);
            }

            widgets[name] = widget;
        }

        return {
            widgets: widgets,
            getValues: () => {
                const vals = {};
                for (const name in widgets) {
                    const w = widgets[name];
                    if (w.type === 'checkbox') vals[name] = w.checked;
                    else if (w.type === 'number') vals[name] = parseFloat(w.value);
                    else vals[name] = w.value;
                }
                return vals;
            },
            setValues: (vals) => {
                for (const name in vals) {
                    if (widgets[name]) {
                        const w = widgets[name];
                        if (w.type === 'checkbox') w.checked = vals[name];
                        else w.value = vals[name];
                    }
                }
            }
        };
    }
};

// Optional browser global for legacy code / debugging
if (typeof window !== 'undefined') {
    window.GUIutils = GUIutils;
}
