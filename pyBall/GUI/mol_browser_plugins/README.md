# mol_browser_plugins — VispyMolBrowser extension panels

Extensible side-panel plugins for [`VispyMolBrowser.py`](../VispyMolBrowser.py). Plugins receive a shared `MolBrowserContext` (directory, selection, loaded molecule, 3D viewer) and attach analysis UI without modifying the core browser loop.

- **base.py** — `MolBrowserPlugin` ABC and `MolBrowserContext`; lifecycle hooks (`is_relevant`, `on_molecule_selected`, …)
- **registry.py** — `MolBrowserPluginRegistry` (ordered registration, entry filtering) and `MolBrowserPluginHost` (east tab strip)
- **vibration_spectrum.py** — `VibrationSpectrumPlugin`: zoomable FTIR histogram, mode pick, 3D arrows/animation via `04_hessian.npz` + `05_spectrum.npz`
- **__init__.py** — `default_plugin_registry()` wires built-in plugins; add new plugins here

**User guide:** [`doc/Topics/FTIR_Nanocrystals/Python_Vispy_MolBrowser_Plugins.md`](../../../doc/Topics/FTIR_Nanocrystals/Python_Vispy_MolBrowser_Plugins.md)

**Tests:** `tests/tSiNCs/test_mol_browser_plugins.py`
