"""
BaseGUI.py - Common GUI utilities and widget templates for PyQt5 applications

PURPOSE:
This module provides a base class (BaseGUI) with reusable helper methods for creating
common PyQt5 widgets. It is designed to reduce boilerplate code and increase consistency
across GUI applications in the FireCore project.

CODING POLICIES:
- Polymorphic design: Use default arguments to consolidate similar functions
  - Example: spinBox(int_mode=True) instead of separate spinBox/spinBoxInt
  - Example: textEdit(plain=True) instead of separate textEdit/plainTextEdit
  - Example: fileDialog(mode='save'|'open'|'directory') instead of separate functions
- Single-line function calls preferred
- Minimal dependencies: only PyQt5, json, re, os
- Fail loudly: no silent error handling
- Explicit inputs/outputs with default named arguments

USAGE FOR CODE REDUCTION AND REUSE:
Instead of writing verbose widget creation code like:
    btn = QtWidgets.QPushButton("Click Me")
    btn.clicked.connect(callback)
    layout.addWidget(btn)

Use the helper:
    self.button("Click Me", callback, layout=layout)

This pattern applies to buttons, checkboxes, comboboxes, spinboxes, text edits,
file dialogs, labels, and group boxes. By inheriting from BaseGUI, GUI classes
can focus on application logic rather than widget boilerplate.

EXAMPLE REFACTORING:
Before:
    spin = QtWidgets.QSpinBox()
    spin.setRange(0, 100)
    spin.setValue(0)
    spin.setEnabled(False)
    spin.valueChanged.connect(callback)
    layout.addWidget(spin)

After:
    spin = self.spinBox(0, vmin=0, vmax=100, enabled=False, callback=callback, layout=layout, int_mode=True)
"""
from PyQt5 import QtWidgets
from PyQt5 import QtCore, QtWidgets
from PyQt5.QtGui import QFont
from PyQt5 import QtCore, QtGui, QtWidgets
from PyQt5.QtWidgets import QHBoxLayout

import json
import re
import os

# --- helper to extract balanced block ---
def extract_json_block(src, start_pos, open_sym, close_sym, bToText=True ):
    """Extract balanced block (e.g., {}, [], ()) from source string starting at position."""
    #print("extract_json_block(): ", start_pos, open_sym, close_sym, bToText)
    open_idx = src.find(open_sym, start_pos)
    if open_idx == -1: return None
    depth = 0
    for i in range(open_idx, len(src)):
        if   src[i] == open_sym:  depth += 1
        elif src[i] == close_sym: depth -= 1
        if depth == 0:
            if bToText:
                b_start, b_end = open_idx, i
                text = src[b_start+1:b_end].strip('\n\r ') if b_start is not None else ''  
                text = "\n".join([ l.strip() for l in text.split('\n') ])
                print("extract_json_block() extracted text:\n", text)
                return text
            else:
                return (open_idx,i)
    print("extract_json_block() failed to extract block")
    return None

def strip_json_comments(json_str):
    """Remove // and /* */ style comments from JSON string."""
    pattern = r'//.*?$|/\*.*?\*/|\'(?:\\.|[^\\\'])*\'|"(?:\\.|[^\\"])*"'
    return re.sub(
        pattern,
        lambda m: m.group(0) if m.group(0).startswith(('"', "'")) else '',
        json_str,
        flags=re.MULTILINE|re.DOTALL
    )

class BaseGUI(QtWidgets.QMainWindow):
    """Base class containing common GUI utilities and widget templates."""
    
    def __init__(self, title="Application GUI"):
        """Initialize BaseGUI with optional title."""
        super().__init__()
        # set smaller default font
        app = QtWidgets.QApplication.instance()
        if app:
            app.setFont(QFont("Sans", 8))
        self.setAttribute(QtCore.Qt.WA_DeleteOnClose)
        self.setWindowTitle(title)
        self.main_widget = QtWidgets.QWidget(self)
        self.setCentralWidget(self.main_widget)
    
    def button(self, text, callback=None, tooltip=None, layout=None):
        """Create QPushButton with optional callback, tooltip, and auto-add to layout."""
        btn = QtWidgets.QPushButton(text)
        if callback is not None: btn.clicked.connect(callback)
        if tooltip  is not None: btn.setToolTip(tooltip)
        if layout   is not None: layout.addWidget(btn)
        return btn

    def checkBox(self, text, checked=False, callback=None, layout=None):
        """Create QCheckBox with optional callback and auto-add to layout."""
        chk = QtWidgets.QCheckBox(text)
        chk.setChecked(checked)
        if callback is not None: chk.stateChanged.connect(callback)
        if layout   is not None: layout.addWidget(chk)
        return chk

    def comboBox(self, items=None, callback=None, layout=None, pass_index=False):
        """Create QComboBox with optional items, callback, and auto-add to layout."""
        cb = QtWidgets.QComboBox()
        if items is not None:    cb.addItems(items)
        if callback is not None:
            if pass_index:
                cb.currentIndexChanged.connect(callback)
            else:
                cb.currentTextChanged.connect(callback)
        if layout   is not None: layout.addWidget(cb)
        return cb

    def spinBox(self, value=0.0, step=0.1, max_width=80, vmin=-1e9, vmax=1e9, decimals=4, enabled=True, callback=None, layout=None, label=None, int_mode=False):
        """Create QSpinBox (int_mode=True) or QDoubleSpinBox (default) with auto-add to layout."""
        spin = QtWidgets.QSpinBox() if int_mode else QtWidgets.QDoubleSpinBox()
        spin.setDecimals(decimals) if not int_mode else None
        spin.setSingleStep(int(step) if int_mode else step)
        spin.setRange(vmin, vmax)
        spin.setValue(value)
        spin.setMaximumWidth(max_width)
        spin.setEnabled(enabled)
        if callback is not None: spin.valueChanged.connect(callback)
        elif hasattr(self, 'on_param_changed'): spin.valueChanged.connect(self.on_param_changed)
        if layout is not None:
            if label is not None: layout.addRow(label, spin)
            else: layout.addWidget(spin)
        return spin

    def textEdit(self, text="", read_only=False, min_size=None, layout=None, wrap=False, plain=False):
        """Create QTextEdit (default) or QPlainTextEdit (plain=True) with auto-add to layout."""
        txt = QtWidgets.QPlainTextEdit() if plain else QtWidgets.QTextEdit()
        if text: txt.setPlainText(text) if plain else txt.setText(text)
        txt.setReadOnly(read_only)
        if min_size is not None: txt.setMinimumSize(*min_size)
        if not wrap and not plain: txt.setWordWrapMode(QtGui.QTextOption.NoWrap)
        if layout is not None: layout.addWidget(txt)
        return txt

    def fileDialog(self, mode="open", title="Select File", filter_str="All Files (*)", start_dir=None):
        """Open file dialog (mode='open'|'save'|'directory') and return selected path or None."""
        if start_dir is None: start_dir = os.path.expanduser("~")
        if mode == "directory":
            dialog = QtWidgets.QFileDialog()
            dialog.setFileMode(QtWidgets.QFileDialog.Directory)
            dialog.setDirectory(start_dir)
            dialog.setWindowTitle(title)
            if dialog.exec_():
                selected = dialog.selectedFiles()[0]
                return selected if os.path.isdir(selected) else None
            return None
        elif mode == "save":
            fname, _ = QtWidgets.QFileDialog.getSaveFileName(self, title, start_dir, filter_str)
            return fname if fname else None
        else:
            fname, _ = QtWidgets.QFileDialog.getOpenFileName(self, title, start_dir, filter_str)
            return fname if fname else None

    def setSpinBox(self, spin, value=None, vmin=None, vmax=None, enabled=None):
        """Set multiple QSpinBox/QDoubleSpinBox properties at once."""
        if value is not None: spin.setValue(value)
        if vmin is not None or vmax is not None: spin.setRange(vmin if vmin is not None else spin.minimum(), vmax if vmax is not None else spin.maximum())
        if enabled is not None: spin.setEnabled(enabled)
        return spin

    def label(self, text, layout=None, word_wrap=False):
        """Create QLabel with optional word wrap and auto-add to layout."""
        lbl = QtWidgets.QLabel(text)
        lbl.setWordWrap(word_wrap)
        if layout is not None: layout.addWidget(lbl)
        return lbl

    def group(self, title, layout=None):
        """Create QGroupBox with optional layout."""
        grp = QtWidgets.QGroupBox(title)
        if layout is not None: grp.setLayout(layout)
        return grp

    def spin_row(self, defaults, step, layout=None, label=None):
        """Create row of multiple spin boxes in horizontal layout with label."""
        container = QtWidgets.QWidget()
        hbox = QHBoxLayout(container)
        hbox.setContentsMargins(0,0,0,0)
        hbox.setSpacing(2)
        spins = [self.spinBox(d, step, 70, layout=hbox) for d in defaults]
        #for spin in spins: hbox.addWidget(spin)
        self.params_layout.addRow(label, container)
        #if layout: layout.addWidget(container)
        return spins

    # def populate_params_from_json(self, params_dict):
    #     """Create spin boxes from *params_dict*"""
    #     # Should be implemented by child class
    #     raise NotImplementedError("populate_params_from_json must be implemented by child class")

    def populate_params_from_dict(self, params_dict):
        """Create spin boxes from params_dict and populate params_layout."""
        print("---------------\npopulate_params_from_dict()")
        while self.params_layout.rowCount() > 0: self.params_layout.removeRow(0)
        self.param_widgets.clear()
        for name, (typ,defaults,step) in params_dict.items():
            print("name: ", name, "typ: ", typ, "defaults: ", defaults, "step: ", step)
            # Ensure defaults is a list-like for consistent handling
            if not isinstance(defaults, (list, tuple)):
                defaults = [defaults]
            if len(defaults) == 1:
                print("single value: ", name, defaults, step)
                self.param_widgets[name] = self.spinBox(defaults[0], step, layout=self.params_layout, label=name)
            else:
                print("multiple values: ", name, defaults, step)
                self.param_widgets[name] = self.spin_row(defaults, step, layout=self.params_layout, label=name)
        #exit()
        self.update_sim_uniforms()
        print("populate_params_from_dict() DONE")