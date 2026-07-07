"""Mol browser plugin protocol — extend VispyMolBrowser via registered panels + lifecycle hooks.

Plugins receive MolBrowserContext (directory, selection, molecule, VispyMolView) and attach QWidget
panels without modifying the core browser. Registry sorts by priority; host shows one tab per
relevant plugin. Register new plugins in mol_browser_plugins/__init__.py::default_plugin_registry.

Guide: doc/Topics/FTIR_Nanocrystals/Python_Vispy_MolBrowser_Plugins.md
"""
from __future__ import annotations

from abc import ABC, abstractmethod
from dataclasses import dataclass
from typing import TYPE_CHECKING, Optional, Sequence

from PyQt5 import QtWidgets

if TYPE_CHECKING:
    from pyBall.GUI.VispyMolBrowser import MoleculeData, VispyMolView


@dataclass
class MolBrowserContext:
    """Snapshot passed to plugins on directory / selection changes."""
    directory: str
    selected_path: Optional[str]
    molecule: Optional['MoleculeData']
    viewer: 'VispyMolView'

    def pipeline_stages(self):
        from pyBall.io.crystal_npz import find_nanocrystal_pipeline_stages, pipeline_dir_for_molecule_path
        if self.selected_path:
            return find_nanocrystal_pipeline_stages(pipeline_dir_for_molecule_path(self.selected_path))
        if self.directory:
            return find_nanocrystal_pipeline_stages(self.directory)
        return {}


class MolBrowserPlugin(ABC):
    """Base class for VispyMolBrowser side-panel extensions."""

    plugin_id: str = 'base'
    title: str = 'Plugin'
    priority: int = 0

    @abstractmethod
    def is_relevant(self, ctx: MolBrowserContext) -> bool:
        """Return True when this plugin should be visible for the current browser state."""

    def create_panel(self, parent: QtWidgets.QWidget) -> QtWidgets.QWidget:
        raise NotImplementedError(f"{type(self).__name__}.create_panel")

    def on_directory_changed(self, ctx: MolBrowserContext):
        pass

    def on_molecule_selected(self, ctx: MolBrowserContext):
        pass

    def on_deactivate(self):
        pass

    def filter_directory_entries(self, entries: Sequence[str], ctx: MolBrowserContext) -> Sequence[str]:
        return entries
