"""Plugin registry + tab host for VispyMolBrowser — ordered registration and east-side QTabWidget."""
from __future__ import annotations

from typing import Iterable, List, Optional, Sequence

from PyQt5 import QtCore, QtWidgets

from pyBall.GUI.mol_browser_plugins.base import MolBrowserContext, MolBrowserPlugin


class MolBrowserPluginRegistry:
    """Ordered plugin list with optional dynamic registration."""

    def __init__(self, plugins: Optional[Iterable[MolBrowserPlugin]] = None):
        self._plugins: List[MolBrowserPlugin] = []
        if plugins:
            for p in plugins:
                self.register(p)

    def register(self, plugin: MolBrowserPlugin):
        if any(p.plugin_id == plugin.plugin_id for p in self._plugins):
            raise ValueError(f"MolBrowserPluginRegistry: duplicate plugin_id={plugin.plugin_id!r}")
        self._plugins.append(plugin)
        self._plugins.sort(key=lambda p: (-int(p.priority), p.plugin_id))

    @property
    def plugins(self):
        return tuple(self._plugins)

    def filter_directory_entries(self, entries: Sequence[str], ctx: MolBrowserContext) -> List[str]:
        out = list(entries)
        for p in self._plugins:
            out = list(p.filter_directory_entries(out, ctx))
        return out

    def notify_directory_changed(self, ctx: MolBrowserContext, host: 'MolBrowserPluginHost'):
        for p in self._plugins:
            p.on_directory_changed(ctx)
        host.refresh(ctx)

    def notify_molecule_selected(self, ctx: MolBrowserContext, host: 'MolBrowserPluginHost'):
        for p in self._plugins:
            p.on_molecule_selected(ctx)
        host.refresh(ctx)


class MolBrowserPluginHost(QtWidgets.QWidget):
    """Right-hand tab strip — one tab per relevant plugin."""

    def __init__(self, registry: MolBrowserPluginRegistry, parent=None):
        super().__init__(parent)
        self._registry = registry
        self._panels = {}
        self._active_ids = set()
        layout = QtWidgets.QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        self._tabs = QtWidgets.QTabWidget()
        self._tabs.setTabPosition(QtWidgets.QTabWidget.East)
        self._placeholder = QtWidgets.QLabel('(no analysis plugins for this directory)')
        self._placeholder.setAlignment(QtCore.Qt.AlignCenter)
        self._placeholder.setWordWrap(True)
        self._stack = QtWidgets.QStackedWidget()
        self._stack.addWidget(self._placeholder)
        self._stack.addWidget(self._tabs)
        layout.addWidget(self._stack)
        self._ctx = None
        self.hide()

    def refresh(self, ctx: MolBrowserContext):
        self._ctx = ctx
        want = [p for p in self._registry.plugins if p.is_relevant(ctx)]
        want_ids = {p.plugin_id for p in want}
        for pid in list(self._active_ids - want_ids):
            p = next(pp for pp in self._registry.plugins if pp.plugin_id == pid)
            p.on_deactivate()
            w = self._panels.pop(pid, None)
            if w is not None:
                idx = self._tabs.indexOf(w)
                if idx >= 0:
                    self._tabs.removeTab(idx)
                w.deleteLater()
        for p in want:
            if p.plugin_id not in self._panels:
                panel = p.create_panel(self._tabs)
                self._panels[p.plugin_id] = panel
                self._tabs.addTab(panel, p.title)
        self._active_ids = want_ids
        if want:
            self._stack.setCurrentWidget(self._tabs)
            self.show()
        else:
            self._stack.setCurrentWidget(self._placeholder)
            self.hide()
