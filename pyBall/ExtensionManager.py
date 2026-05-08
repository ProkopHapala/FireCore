"""
ExtensionManager.py - Lazy-loading extension system for KekuleExplorerGUI
=========================================================================
Extensions are loaded on first use. Each extension declares:
  - Python module / class to import
  - Required file-system paths (validated before load)
  - Optional Python-package dependencies (validated before load)
  - A build_ui(window) callable that returns a UIComponents bundle

The manager never silently swallows failures: missing deps / bad paths raise
ExtensionLoadError immediately.  An ExtensionProxy stands in for unloaded
extensions and raises a loud error on any attribute access.
"""

import importlib, os, json

# ---------------------------------------------------------------------------
# Exceptions
# ---------------------------------------------------------------------------
class ExtensionLoadError(RuntimeError):
    pass

class ExtensionNotAvailableError(AttributeError):
    pass

# ---------------------------------------------------------------------------
# Proxy: stands in for an extension that failed to load
# ---------------------------------------------------------------------------
class ExtensionProxy:
    """Raises ExtensionNotAvailableError on every attribute access."""
    def __init__(self, name: str, reason: str):
        object.__setattr__(self, '_name',   name)
        object.__setattr__(self, '_reason', reason)
    def __getattr__(self, attr):
        raise ExtensionNotAvailableError(
            f"Extension '{object.__getattribute__(self,'_name')}' not available: "
            f"{object.__getattribute__(self,'_reason')}"
        )
    def __bool__(self): return False   # easy "if ext:" check

# ---------------------------------------------------------------------------
# UIComponents: what every build_ui() must return
# ---------------------------------------------------------------------------
class UIComponents:
    """Bundle returned by each extension's build_ui(window) function.

    Fields (all optional, default empty):
      panel       – QWidget shown in a CollapsibleSection
      edit_modes  – list of (label: str, callback: callable)
      view_modes  – list of (label: str, callback: callable)
    """
    def __init__(self, panel=None, edit_modes=None, view_modes=None):
        self.panel      = panel
        self.edit_modes = edit_modes or []
        self.view_modes = view_modes or []

# ---------------------------------------------------------------------------
# Extension registry (declarative metadata only)
# ---------------------------------------------------------------------------
#   module       – dotted import path
#   class_name   – None means use the module itself as the instance
#   dependencies – list of importable package names
#   req_paths    – config keys whose values must be existing paths
#   build_ui     – name of a function in the module returning UIComponents
EXTENSION_REGISTRY = {
    'firecore': dict(
        module='pyBall.FireCore', class_name=None,
        dependencies=[], req_paths=['fdata_dir'],
        build_ui='build_fireball_ui',
    ),
    'dftb': dict(
        module='pyBall.dftb_utils', class_name=None,
        dependencies=[], req_paths=[],
        build_ui='build_dftb_ui',
    ),
    'afm': dict(
        module='pyBall.OCL.AFM', class_name='AFMulator',
        dependencies=['pyopencl'], req_paths=['cl_src_dir'],
        build_ui='build_ui',
    ),
    'mmff': dict(
        module='pyBall.MMFF', class_name=None,
        dependencies=[], req_paths=[],
        build_ui='build_ui',
    ),
    'grid': dict(
        module='pyBall.FireballOCL.Grid', class_name='GridProjector',
        dependencies=['pyopencl'], req_paths=['fdata_dir'],
        build_ui='build_ui',
    ),
    'psi4': dict(
        module='pyBall.psi4_utils', class_name=None,
        dependencies=['psi4'], req_paths=[],
        build_ui='build_ui',
    ),
    'pyscf': dict(
        module='pyBall.pyscf_utils', class_name=None,
        dependencies=['pyscf'], req_paths=[],
        build_ui='build_ui',
    ),
    'moldyn': dict(
        module='pyBall.OCL.MolecularDynamics', class_name='MolecularDynamics',
        dependencies=['pyopencl'], req_paths=[],
        build_ui='build_ui',
    ),
    'povray': dict(
        module='pyBall.POVray', class_name=None,
        dependencies=[], req_paths=[],
        build_ui='build_ui',
    ),
}

# ---------------------------------------------------------------------------
# Default configuration (per-extension)
# ---------------------------------------------------------------------------
DEFAULT_CONFIG = {
    'firecore': dict(enabled=True,  fdata_dir='/home/prokop/Fireball/Fdata_HCNOS', verbosity=0),
    'dftb':     dict(enabled=True,  executable='dftb+', workdir='./dftb_workdir'),
    'afm':      dict(enabled=False, cl_src_dir='../../cpp/common_resources/cl'),
    'mmff':     dict(enabled=False, lib_path=''),
    'grid':     dict(enabled=False, fdata_dir=''),
    'psi4':     dict(enabled=False),
    'pyscf':    dict(enabled=False),
    'moldyn':   dict(enabled=False),
    'povray':   dict(enabled=False),
}

# ---------------------------------------------------------------------------
# ExtensionLoader: low-level load / cache
# ---------------------------------------------------------------------------
class ExtensionLoader:
    def __init__(self):
        self._cache  = {}    # name -> module/instance  (loaded OK)
        self._errors = {}    # name -> error string      (load failed)

    def _validate(self, name: str, cfg: dict):
        """Check deps and paths; raise ExtensionLoadError on first failure."""
        meta = EXTENSION_REGISTRY[name]
        for pkg in meta['dependencies']:
            try:
                importlib.import_module(pkg)
            except ImportError:
                raise ExtensionLoadError(f"[{name}] missing Python package: {pkg}")
        for key in meta['req_paths']:
            val = cfg.get(key, '')
            if not val or not os.path.exists(val):
                raise ExtensionLoadError(f"[{name}] required path missing: {key}={val!r}")

    def load(self, name: str, cfg: dict):
        """Return loaded module/instance, raising ExtensionLoadError on failure."""
        if name in self._cache:
            return self._cache[name]
        if name not in EXTENSION_REGISTRY:
            raise ExtensionLoadError(f"Unknown extension: {name!r}")
        self._validate(name, cfg)
        meta = EXTENSION_REGISTRY[name]
        try:
            mod = importlib.import_module(meta['module'])
        except ImportError as e:
            self._errors[name] = str(e)
            raise ExtensionLoadError(f"[{name}] import failed: {e}") from e
        if meta['class_name']:
            obj = getattr(mod, meta['class_name'])
        else:
            obj = mod
        self._cache[name] = obj
        return obj

    def get_or_proxy(self, name: str, cfg: dict):
        """Return extension or an ExtensionProxy if unavailable."""
        try:
            return self.load(name, cfg)
        except ExtensionLoadError as e:
            self._errors[name] = str(e)
            return ExtensionProxy(name, str(e))

    def status(self, name: str) -> str:
        if name in self._cache:  return 'loaded'
        if name in self._errors: return f'error: {self._errors[name]}'
        return 'not_loaded'

# ---------------------------------------------------------------------------
# ExtensionManager: high-level API used by the GUI
# ---------------------------------------------------------------------------
class ExtensionManager:
    """Singleton-like manager.  Instantiate once in KekuleExplorerWindow.__init__."""

    _CONFIG_FILE = os.path.join(os.path.dirname(__file__), 'extensions_config.json')

    def __init__(self):
        self._cfg    = {k: dict(v) for k, v in DEFAULT_CONFIG.items()}
        self._loader = ExtensionLoader()
        self._load_config_file()

    # ---- config ----------------------------------------------------------
    def _load_config_file(self):
        if os.path.isfile(self._CONFIG_FILE):
            try:
                with open(self._CONFIG_FILE) as f:
                    data = json.load(f)
                for name, vals in data.items():
                    if name in self._cfg:
                        self._cfg[name].update(vals)
            except Exception as e:
                print(f"ExtensionManager: could not read config: {e}")

    def save_config(self):
        try:
            with open(self._CONFIG_FILE, 'w') as f:
                json.dump(self._cfg, f, indent=2)
        except Exception as e:
            print(f"ExtensionManager: could not save config: {e}")

    def get_config(self, name: str) -> dict:
        return self._cfg.get(name, {})

    def set_config(self, name: str, key: str, value):
        self._cfg.setdefault(name, {})[key] = value

    # ---- enabled list ----------------------------------------------------
    def enabled_extensions(self):
        return [n for n, c in self._cfg.items() if c.get('enabled', False)]

    # ---- access ----------------------------------------------------------
    def get(self, name: str):
        """Return extension instance or proxy. Never raises."""
        return self._loader.get_or_proxy(name, self._cfg.get(name, {}))

    def require(self, name: str):
        """Return extension instance; raises ExtensionLoadError if unavailable."""
        return self._loader.load(name, self._cfg.get(name, {}))

    def is_loaded(self, name: str) -> bool:
        return self._loader.status(name) == 'loaded'

    def status(self, name: str) -> str:
        return self._loader.status(name)

    # ---- UI builder ------------------------------------------------------
    def build_ui(self, name: str, window) -> UIComponents:
        """Import the extension module and call its build_ui(window) -> UIComponents."""
        ext = self.get(name)
        if not ext:
            return UIComponents()    # proxy → no panel, no modes
        meta = EXTENSION_REGISTRY.get(name, {})
        builder_name = meta.get('build_ui', 'build_ui')
        try:
            # Check if builder is a function in this module (for extensions with UI here)
            if builder_name == 'build_fireball_ui':
                return build_fireball_ui(window)
            elif builder_name == 'build_dftb_ui':
                return build_dftb_ui(window)
            # Otherwise, import from the extension module
            mod = importlib.import_module(meta['module'])
            builder = getattr(mod, builder_name, None)
            if builder is None:
                return UIComponents()
            return builder(window)
        except Exception as e:
            print(f"ExtensionManager: build_ui({name}) failed: {e}")
            return UIComponents()


# ---------------------------------------------------------------------------
# Extension UI builders (GUI-related functions for specific extensions)
# ---------------------------------------------------------------------------

def build_fireball_ui(window):
    """Build Fireball panel for KekuleExplorerGUI.
    Returns ExtensionManager.UIComponents.
    """
    from PyQt5 import QtWidgets
    import numpy as np
    from pyBall import elements

    panel = QtWidgets.QWidget()
    layout = QtWidgets.QVBoxLayout(panel)
    layout.setSpacing(3)
    layout.setContentsMargins(2, 2, 2, 2)

    scf_btn = QtWidgets.QPushButton("Compute SCF")
    scf_btn.clicked.connect(window.compute_orbitals)
    layout.addWidget(scf_btn)

    window.orbital_info_label = QtWidgets.QLabel("Orbitals: Not computed")
    window.orbital_info_label.setWordWrap(True)
    layout.addWidget(window.orbital_info_label)

    row1 = QtWidgets.QHBoxLayout()
    row1.addWidget(QtWidgets.QLabel("Z:"))
    window.z_height_spinbox = window.spinBox(2.0, 0.5, vmin=-10.0, vmax=20.0)
    row1.addWidget(window.z_height_spinbox)
    row1.addWidget(QtWidgets.QLabel("Orb:"))
    window.orbital_spinbox = window.spinBox(0, vmin=0, vmax=999, enabled=False, callback=window.update_orbital_energy_label, int_mode=True)
    row1.addWidget(window.orbital_spinbox)
    layout.addLayout(row1)

    row2 = QtWidgets.QHBoxLayout()
    window.plot_orb_btn = QtWidgets.QPushButton("Plot Orb")
    window.plot_orb_btn.setEnabled(False)
    window.plot_orb_btn.clicked.connect(window.plot_orbital_from_spinbox)
    row2.addWidget(window.plot_orb_btn)
    window.plot_density_btn = QtWidgets.QPushButton("Plot Dens")
    window.plot_density_btn.setEnabled(False)
    window.plot_density_btn.clicked.connect(window.plot_density)
    row2.addWidget(window.plot_density_btn)
    window.plot_delta_btn = QtWidgets.QPushButton("Plot Delta")
    window.plot_delta_btn.setEnabled(False)
    window.plot_delta_btn.clicked.connect(window.plot_delta_rho)
    row2.addWidget(window.plot_delta_btn)
    layout.addLayout(row2)

    fdata_btn = QtWidgets.QPushButton("Set Fdata")
    fdata_btn.clicked.connect(window.set_fdata_path)
    layout.addWidget(fdata_btn)

    view_modes = [
        ('Molecular Orbital', lambda: window.set_view_mode('orbital')),
        ('Density',           lambda: window.set_view_mode('density')),
        ('Delta-Rho',         lambda: window.set_view_mode('delta_rho')),
    ]
    return UIComponents(panel=panel, view_modes=view_modes)


def build_dftb_ui(window):
    """Build DFTB+ panel for KekuleExplorerGUI.
    Returns ExtensionManager.UIComponents.
    """
    from PyQt5 import QtWidgets

    panel = QtWidgets.QWidget()
    layout = QtWidgets.QVBoxLayout(panel)
    layout.setSpacing(3)
    layout.setContentsMargins(2, 2, 2, 2)

    relax_btn = QtWidgets.QPushButton("Relax (DFTB+)")
    relax_btn.clicked.connect(window.run_relaxation)
    layout.addWidget(relax_btn)

    window.dftb_status_label = QtWidgets.QLabel("Status: Ready")
    window.dftb_status_label.setWordWrap(True)
    layout.addWidget(window.dftb_status_label)

    return UIComponents(panel=panel)
