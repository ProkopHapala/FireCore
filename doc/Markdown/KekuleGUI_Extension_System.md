# USER

plese figure out how to make some extensiton manager in @KekuleExplorerGUI.py 
currently it use two external program/libraries DFTB and Fireball, since these progam mey be not installed we must load them lazzyliy, set some paths which may differ on different system. We will face similar issues with all such program (we want to add PPAFM imaging, MMFF/UFF oprimized, grid forcefield and others). We should not tread all such extension separately but make some robust extension system which handles it. We should have some config set of option which extenstion are loaded (just boleand but can be eventually read from some config file like json, no need to do it now), and in that config we can have all the paths described. 

Look how we currently handle DFTB and FireCore, and if we can unify thet in single extension system. Later we also want to add @AFM.py @MMFF.py @Grid.py @POVray.py @psi4_utils.py @pyscf_utils.py @MolecularDynamics.py 

do not program anything yet, rather review the code, figure out what could be specific requirements and how to desing fexible general extension system.

This system should handle safe loading unloading of the extension (extension may be not present on the system, paths may be worng), and I/O of data, which is typically atomic ceometry (atom positions), eventually bdning topology, and some other data (points at which we compute density or moleculer orbitlas, values of function at these sample points, 3D grid, etc.)

Review the code, think about it and make detailed notes and a plan

---

# SWE-1.6

I'll review the code to understand the current extension handling and design a flexible general extension system. Let me start by reading the relevant files.

I need to continue reading the KekuleExplorerGUI.py file to see how it currently handles DFTB and Fireball, and also look at more of the extension files to understand their interfaces.

Now let me look at the KekuleBackend to see how DFTB is handled there:

Let me look at more of the KekuleBackend to understand DFTB handling better:

Based on my review of the code, here are my detailed notes and a plan for designing a flexible extension system:

## Current Extension Handling Analysis

### KekuleExplorerGUI.py Current State

**FireCore (Fireball DFT):**
- **Import strategy:** Direct module-level import with try/except (lines 48-52) - NOT lazy
- **Path configuration:** `fdata_path` stored in QSettings (line 67), defaults to hardcoded path
- **Usage:** [compute_orbitals()](cci:1://file:///home/prokop/git/FireCore/pyBall/KekuleExplorerGUI.py:1026:4-1092:62), [plot_orbital_from_spinbox()](cci:1://file:///home/prokop/git/FireCore/pyBall/KekuleExplorerGUI.py:471:4-501:64), [plot_density()](cci:1://file:///home/prokop/git/FireCore/pyBall/KekuleExplorerGUI.py:1149:4-1171:72), [plot_delta_rho()](cci:1://file:///home/prokop/git/FireCore/pyBall/KekuleExplorerGUI.py:1173:4-1195:74), [set_fdata_path()](cci:1://file:///home/prokop/git/FireCore/pyBall/KekuleExplorerGUI.py:1197:4-1205:104)
- **Data I/O:** 
  - Input: Atomic positions from `self.backend.sys.apos`, element types
  - Output: Eigenvalues, wavefunction coefficients, density values on grids
- **Error handling:** Try/except in [compute_orbitals()](cci:1://file:///home/prokop/git/FireCore/pyBall/KekuleExplorerGUI.py:1026:4-1092:62) with [_raise()](cci:1://file:///home/prokop/git/FireCore/pyBall/KekuleExplorerGUI.py:175:4-195:35) helper

**DFTB+:**
- **Import strategy:** Lazy import inside [KekuleBackend.run_relaxation()](cci:1://file:///home/prokop/git/FireCore/pyBall/KekuleBackend.py:909:4-938:29) (line 912: `from . import dftb_utils`)
- **Path configuration:** None visible (likely handled inside dftb_utils)
- **Usage:** [run_relaxation()](cci:1://file:///home/prokop/git/FireCore/pyBall/KekuleBackend.py:909:4-938:29) delegates to backend
- **Data I/O:** Atomic positions, element types, lattice vectors → relaxed positions, forces, energy

### Extension File Analysis

| Extension | Dependencies | Path Config | Data I/O | Key Requirements |
|-----------|-------------|-------------|----------|-----------------|
| **AFM.py** | PyOpenCL, OpenCL kernels | `cl_src_dir` (defaults to `cpp/common_resources/cl`) | XYZ→AtomicSystem, FF images, scan data | GPU device, ElementTypes.dat |
| **MMFF.py** | C++ library (ctypes) | Build path: `cpp/Build/libs/Molecular` or `cpp/Build-opt/libs/Molecular` | Many spline/grid/force operations | Compiled C++ lib, parameter files |
| **Grid.py** | PyOpenCL, FireCore, FdataParser | `fdata_dir` (Fdata basis files) | Sparse density→3D grid projection | Fdata dir, OpenCL device, basis loading |
| **POVray.py** | None (pure Python) | None | Atomic geometry→POV-Ray format | Minimal dependencies |
| **psi4_utils.py** | psi4, resp (optional) | None (lazy-loaded) | Geometry→QM calc→charges/energy | Lazy import, QM package installation |
| **pyscf_utils.py** | pyscf | None | Geometry→QM calc→energy | QM package installation |
| **MolecularDynamics.py** | PyOpenCL, MMFF, OpenCL kernels | `cpp/common_resources/cl/relax.cl` | Multi-system MD, force fields | GPU, C++ MMFF, OpenCL kernels |

## Specific Requirements Identified

### 1. **Lazy Loading Strategy**
- **Problem:** FireCore imports at module level causes crash if not available
- **Solution:** All extensions should be lazy-loaded on first use
- **Pattern:** Import inside methods or use a lazy import registry

### 2. **Path Configuration**
- **Problem:** Paths are hardcoded or scattered across modules
- **Types of paths needed:**
  - Executable/binary paths (DFTB+, external QM packages)
  - Library paths (C++ MMFF lib, FireCore lib)
  - Data directories (Fdata, ElementTypes.dat, kernel .cl files)
  - Working directories (temporary calculation folders)
- **Solution:** Centralized configuration with system-specific overrides

### 3. **Safe Loading/Unloading**
- **Error scenarios to handle:**
  - Extension not installed (ImportError)
  - Wrong path (FileNotFoundError, OSError)
  - Incompatible version (AttributeError, API mismatch)
  - GPU not available (OpenCL errors)
  - Missing data files
- **Current approach:** Try/except with [_raise()](cci:1://file:///home/prokop/git/FireCore/pyBall/KekuleExplorerGUI.py:175:4-195:35) - good but needs standardization

### 4. **Data I/O Standardization**
- **Common data types across extensions:**
  - **Atomic geometry:** positions (N×3), element types (N), charges (N)
  - **Bonding topology:** bonds (M×2), angles, dihedrals
  - **Grid data:** 3D arrays (nx×ny×nz), origin, spacing
  - **Sample points:** arbitrary points (P×3) for function evaluation
  - **Scalar fields:** values at points or on grids
  - **Force/energy data:** forces (N×3), energies (scalar)
- **Need:** Standardized data containers/formats

### 5. **Extension Lifecycle**
- **States:** Not loaded → Loaded (initialized) → Ready → Error
- **Operations:** Load, configure, validate, execute, cleanup, unload

## Proposed Extension System Design

### Architecture Overview

```
ExtensionManager (singleton)
├── ExtensionRegistry (metadata about all extensions)
├── ExtensionConfig (path/settings storage)
├── ExtensionLoader (lazy loading logic)
└── Extension instances (cached loaded extensions)
```

### Core Components

#### 1. **ExtensionRegistry**
- **Purpose:** Declarative metadata about each extension
- **Structure:**
```python
EXTENSION_REGISTRY = {
    'firecore': {
        'module': 'pyBall.FireCore',
        'class': None,  # Module-level functions
        'required_paths': ['fdata_dir'],
        'optional_paths': [],
        'dependencies': [],  # Python packages
        'provides': ['scf', 'orbitals', 'density', 'forces'],
        'data_io': {
            'input': ['apos', 'atypes'],
            'output': ['eigenvalues', 'wfcoef', 'density']
        }
    },
    'dftb': {
        'module': 'pyBall.dftb_utils',
        'class': None,
        'required_paths': ['dftb_executable'],
        'optional_paths': [],
        'dependencies': [],
        'provides': ['relaxation', 'forces', 'energy'],
        'data_io': {
            'input': ['apos', 'enames', 'lvecs'],
            'output': ['apos_relaxed', 'forces', 'energy']
        }
    },
    'afm': {
        'module': 'pyBall.OCL.AFM',
        'class': 'AFMulator',
        'required_paths': ['cl_src_dir'],
        'optional_paths': ['params_path'],
        'dependencies': ['pyopencl'],
        'provides': ['afm_simulation', 'force_field', 'scan'],
        'data_io': {
            'input': ['xyz_path'],
            'output': ['force_field_images', 'scan_data']
        }
    },
    # ... other extensions
}
```

#### 2. **ExtensionConfig**
- **Purpose:** Centralized path and settings storage
- **Storage:** JSON file + QSettings for GUI integration
- **Structure:**
```json
{
    "firecore": {
        "enabled": true,
        "fdata_dir": "/home/prokop/Fireball/Fdata_HCNOS",
        "verbosity": 0
    },
    "dftb": {
        "enabled": true,
        "executable": "/usr/local/bin/dftb+",
        "workdir": "./dftb_workdir"
    },
    "afm": {
        "enabled": false,
        "cl_src_dir": "../../cpp/common_resources/cl",
        "device_preference": "nvidia"
    },
    "mmff": {
        "enabled": true,
        "lib_path": "../../cpp/Build-opt/libs/Molecular"
    },
    "grid": {
        "enabled": false,
        "fdata_dir": "/home/prokop/Fireball/Fdata_HC_minimal"
    },
    "psi4": {
        "enabled": false
    },
    "pyscf": {
        "enabled": false
    }
}
```

#### 3. **ExtensionLoader**
- **Purpose:** Lazy loading with validation
- **Key methods:**
  - `load_extension(name)` - Import module, validate paths, initialize
  - `unload_extension(name)` - Cleanup, clear cache
  - `is_available(name)` - Check if can be loaded
  - `get_extension(name)` - Get cached instance or load

#### 4. **Extension Base Interface**
- **Purpose:** Standardized API for all extensions
- **Methods:**
```python
class ExtensionBase:
    @classmethod
    def check_dependencies(cls) -> bool:
        """Check if required Python packages are installed."""
        pass
    
    @classmethod
    def check_paths(cls, config: dict) -> tuple[bool, str]:
        """Check if required paths exist. Returns (ok, error_msg)."""
        pass
    
    @classmethod
    def initialize(cls, config: dict):
        """Initialize extension with configuration."""
        pass
    
    def get_capabilities(self) -> list[str]:
        """Return list of provided capabilities (e.g., ['relaxation', 'forces'])."""
        pass
    
    def validate_input(self, data: dict) -> tuple[bool, str]:
        """Validate input data format. Returns (ok, error_msg)."""
        pass
```

### Data I/O Standardization

#### Standard Data Containers
```python
@dataclass
class AtomicGeometry:
    apos: np.ndarray  # (N, 3) positions in Angstrom
    enames: list[str]  # element names
    atypes: np.ndarray  # (N,) atomic numbers
    qs: np.ndarray | None = None  # (N,) charges
    bonds: np.ndarray | None = None  # (M, 2) atom indices

@dataclass
class GridSpec:
    origin: np.ndarray  # (3,)
    dA: np.ndarray  # (3,) spacing vectors
    dB: np.ndarray
    dC: np.ndarray
    ngrid: np.ndarray  # (3,) dimensions

@dataclass
class QMResults:
    energy: float
    forces: np.ndarray  # (N, 3)
    eigenvalues: np.ndarray | None = None
    orbitals: np.ndarray | None = None
    density: np.ndarray | None = None
```

### Integration with KekuleExplorerGUI

#### UI Changes
1. **Extension Settings Dialog:** Central place to configure all extensions
2. **Dynamic Section Creation:** Only create UI sections for enabled extensions
3. **Status Indicators:** Show load status (✓ loaded, ✗ error, ? not configured)

#### Refactored Code Structure
```python
class KekuleExplorerWindow(BaseGUI):
    def __init__(self):
        # ... existing init ...
        self.extension_manager = ExtensionManager()
        self.extension_manager.load_config()
        
    def initUI(self):
        # ... existing sections ...
        # Dynamic extension sections
        for ext_name in self.extension_manager.get_enabled_extensions():
            section = self.create_extension_section(ext_name)
            if section:
                side_layout.addWidget(section)
    
    def create_extension_section(self, ext_name):
        if ext_name == 'firecore':
            return self.create_fireball_section()
        elif ext_name == 'dftb':
            return self.create_dftb_section()
        elif ext_name == 'afm':
            return self.create_afm_section()
        # ...
```

### Implementation Plan

#### Phase 1: Core Infrastructure
1. Create `pyBall/ExtensionManager.py` with:
   - ExtensionRegistry (metadata)
   - ExtensionConfig (JSON storage)
   - ExtensionLoader (lazy loading)
   - ExtensionBase (interface)

2. Create data containers in `pyBall/data_structures.py`:
   - AtomicGeometry
   - GridSpec
   - QMResults

#### Phase 2: Migrate Existing Extensions
1. Refactor FireCore usage:
   - Remove module-level import
   - Use ExtensionManager.lazy_load('firecore')
   - Move fdata_path to ExtensionConfig

2. Refactor DFTB usage:
   - Register in ExtensionRegistry
   - Move path config to ExtensionConfig
   - Keep lazy import pattern

#### Phase 3: Add New Extensions
1. Implement AFM extension wrapper
2. Implement MMFF extension wrapper
3. Implement Grid extension wrapper
4. Implement POVray extension wrapper
5. Implement psi4/pyscf extension wrappers
6. Implement MolecularDynamics extension wrapper

#### Phase 4: GUI Integration
1. Create Extension Settings dialog
2. Dynamic section creation
3. Status indicators
4. Error handling improvements

### Key Design Principles

1. **Fail Loudly:** Following user's preference - throw descriptive errors, never silent fallbacks
2. **Lazy Loading:** Import only when needed
3. **Configuration Separation:** Paths in config, not hardcoded
4. **Standardized Data:** Common data structures for I/O
5. **Extensibility:** Easy to add new extensions via registry
6. **Validation:** Check dependencies and paths before use
7. **Cleanup:** Proper unload/cleanup mechanism

### Open Questions

1. **Should extensions be singletons or per-instance?**
   - Likely singleton for most (GPU context, library handles)
   - Per-instance for multi-system scenarios

2. **How to handle GPU device selection?**
   - Global setting or per-extension?
   - Need OpenCL device enumeration

3. **Working directory management?**
   - Each extension gets its own workdir?
   - Central temp directory?

4. **Version compatibility checking?**
   - How to handle API changes in external packages?

This design provides a robust, flexible foundation for managing all external extensions while maintaining the "fail loudly" philosophy and lazy loading requirements.


---

# USER

OK, some standartization if good, but we should not be to gidi, constrained as that would introduce more friction than solve. We should not preliminary limit what data and what function the KekuleGUI can exchange and call formthe extension. We should alow flexibility to define arbitrary function handling the interface, but use them only when the extension is safely loaded, otherwise call some error/exception

Also we can make some GUI widgets or panels, but these panels will be added (rendered, unfolded) only when extension is loaded. Can we fold/unfold panels in pygt? In web it is possible easily. We are figthing with space in the GUI panel, so foldable panels for extenti which can be closed so they do not take space would be very usefull. Some extension e.g. fact classical forcefield should also add new edti-modes like real time relaxed scan when one can drag molecule around while it is relaxing. Thes should be added to "Edit Mode" selection. Similarly some othe selection boxes (view-mode AFM or Molecular Orbital etc.)

plese desing system aloowing this, but not limiting too much what we can actualy do and what data we can comunicated.

now again plan, do not write any code yet

---

# SWE-1.6

Revised design focusing on flexibility and dynamic UI registration:

## Revised Extension System Design

### Core Philosophy

**Flexible, not prescriptive:**
- Extensions define their own data structures and function signatures
- No enforced standard containers (AtomicGeometry, GridSpec, etc.)
- Extension system provides **infrastructure**, not **constraints**
- Safe loading + configuration + dynamic UI registration = the only requirements

### Architecture Overview

```
ExtensionManager (singleton)
├── ExtensionRegistry (declarative metadata)
├── ExtensionConfig (JSON-based settings)
├── ExtensionLoader (lazy loading with validation)
└── Extension instances (cached, with lifecycle)
```

### ExtensionRegistry (Minimal Metadata)

```python
EXTENSION_REGISTRY = {
    'firecore': {
        'module': 'pyBall.FireCore',
        'class': None,  # Module-level API
        'required_paths': ['fdata_dir'],
        'optional_paths': [],
        'dependencies': [],  # Python packages to check
        'ui_builder': 'build_fireball_ui',  # Function to create UI panel
    },
    'dftb': {
        'module': 'pyBall.dftb_utils',
        'class': None,
        'required_paths': ['dftb_executable'],
        'optional_paths': [],
        'dependencies': [],
        'ui_builder': 'build_dftb_ui',
    },
    'afm': {
        'module': 'pyBall.OCL.AFM',
        'class': 'AFMulator',
        'required_paths': ['cl_src_dir'],
        'optional_paths': ['params_path'],
        'dependencies': ['pyopencl'],
        'ui_builder': 'build_afm_ui',
    },
    'mmff': {
        'module': 'pyBall.MMFF',
        'class': None,
        'required_paths': ['lib_path'],
        'optional_paths': ['params_path'],
        'dependencies': [],
        'ui_builder': 'build_mmff_ui',
    },
    'grid': {
        'module': 'pyBall.FireballOCL.Grid',
        'class': 'GridProjector',
        'required_paths': ['fdata_dir'],
        'optional_paths': [],
        'dependencies': ['pyopencl'],
        'ui_builder': 'build_grid_ui',
    },
    'psi4': {
        'module': 'pyBall.psi4_utils',
        'class': None,
        'required_paths': [],
        'optional_paths': [],
        'dependencies': ['psi4', 'resp'],
        'ui_builder': 'build_psi4_ui',
    },
    'pyscf': {
        'module': 'pyBall.pyscf_utils',
        'class': None,
        'required_paths': [],
        'optional_paths': [],
        'dependencies': ['pyscf'],
        'ui_builder': 'build_pyscf_ui',
    },
    'moldyn': {
        'module': 'pyBall.OCL.MolecularDynamics',
        'class': 'MolecularDynamics',
        'required_paths': [],
        'optional_paths': [],
        'dependencies': ['pyopencl'],
        'ui_builder': 'build_moldyn_ui',
    },
    'povray': {
        'module': 'pyBall.POVray',
        'class': None,
        'required_paths': [],
        'optional_paths': [],
        'dependencies': [],
        'ui_builder': 'build_povray_ui',
    },
}
```

### ExtensionConfig (JSON)

```json
{
    "firecore": {
        "enabled": true,
        "fdata_dir": "/home/prokop/Fireball/Fdata_HCNOS",
        "verbosity": 0
    },
    "dftb": {
        "enabled": true,
        "executable": "/usr/local/bin/dftb+",
        "workdir": "./dftb_workdir"
    },
    "afm": {
        "enabled": false,
        "cl_src_dir": "../../cpp/common_resources/cl",
        "device_preference": "nvidia"
    },
    "mmff": {
        "enabled": true,
        "lib_path": "../../cpp/Build-opt/libs/Molecular"
    },
    "grid": {
        "enabled": false,
        "fdata_dir": "/home/prokop/Fireball/Fdata_HC_minimal"
    },
    "psi4": {
        "enabled": false
    },
    "pyscf": {
        "enabled": false
    },
    "moldyn": {
        "enabled": false
    },
    "povray": {
        "enabled": true
    }
}
```

### ExtensionLoader (Core Loading Logic)

```python
class ExtensionLoader:
    def __init__(self):
        self._loaded = {}  # name -> module/class instance
        self._status = {}  # name -> 'not_loaded', 'loaded', 'error'
        self._errors = {}  # name -> error message
    
    def load(self, name: str, config: dict):
        """Load extension with validation."""
        if name in self._loaded:
            return self._loaded[name]
        
        meta = EXTENSION_REGISTRY.get(name)
        if not meta:
            raise ValueError(f"Unknown extension: {name}")
        
        # Check dependencies
        for dep in meta['dependencies']:
            try:
                importlib.import_module(dep)
            except ImportError:
                self._status[name] = 'error'
                self._errors[name] = f"Missing dependency: {dep}"
                raise ExtensionLoadError(f"Missing dependency: {dep}")
        
        # Check required paths
        for path_key in meta['required_paths']:
            path = config.get(path_key)
            if not path or not os.path.exists(path):
                self._status[name] = 'error'
                self._errors[name] = f"Missing required path: {path_key}={path}"
                raise ExtensionLoadError(f"Missing required path: {path_key}={path}")
        
        # Import module
        try:
            module = importlib.import_module(meta['module'])
        except ImportError as e:
            self._status[name] = 'error'
            self._errors[name] = str(e)
            raise ExtensionLoadError(f"Import failed: {e}")
        
        # Get class or module
        if meta['class']:
            instance = getattr(module, meta['class'])()
        else:
            instance = module
        
        self._loaded[name] = instance
        self._status[name] = 'loaded'
        return instance
    
    def get(self, name: str, config: dict):
        """Get extension, load if needed."""
        try:
            return self.load(name, config)
        except ExtensionLoadError as e:
            # Return a proxy that raises errors on any access
            return ExtensionProxy(name, str(e))
    
    def is_available(self, name: str, config: dict) -> bool:
        """Check if extension can be loaded."""
        try:
            self.load(name, config)
            return True
        except ExtensionLoadError:
            return False

class ExtensionProxy:
    """Proxy object that raises errors on any attribute access."""
    def __init__(self, name: str, error_msg: str):
        self._name = name
        self._error_msg = error_msg
    
    def __getattr__(self, attr):
        raise ExtensionNotLoadedError(
            f"Extension '{self._name}' is not available: {self._error_msg}"
        )
```

### Dynamic UI Registration

#### UI Builder Pattern

Each extension provides a `ui_builder` function that returns UI components:

```python
def build_fireball_ui(parent_window, extension_instance, config):
    """Build Fireball UI panel.
    
    Returns dict with:
    - 'panel': QWidget (the main collapsible panel)
    - 'edit_modes': list of (name, callback) tuples
    - 'view_modes': list of (name, callback) tuples
    - 'selection_boxes': list of (name, widget) tuples
    """
    panel = QWidget()
    layout = QVBoxLayout(panel)
    
    # Add Fireball-specific controls
    scf_btn = QPushButton("Compute SCF")
    scf_btn.clicked.connect(lambda: parent_window.compute_orbitals_fireball())
    layout.addWidget(scf_btn)
    
    return {
        'panel': panel,
        'edit_modes': [],  # Fireball doesn't add edit modes
        'view_modes': [
            ('Molecular Orbital', parent_window.set_view_mode_orbital),
            ('Density', parent_window.set_view_mode_density),
        ],
        'selection_boxes': [],
    }

def build_mmff_ui(parent_window, extension_instance, config):
    """Build MMFF UI panel."""
    panel = QWidget()
    layout = QVBoxLayout(panel)
    
    relax_btn = QPushButton("Real-time Relax")
    relax_btn.clicked.connect(lambda: parent_window.toggle_realtime_relax())
    layout.addWidget(relax_btn)
    
    return {
        'panel': panel,
        'edit_modes': [
            ('Relaxed Scan', parent_window.set_edit_mode_relaxed_scan),
        ],
        'view_modes': [],
        'selection_boxes': [],
    }

def build_afm_ui(parent_window, extension_instance, config):
    """Build AFM UI panel."""
    panel = QWidget()
    layout = QVBoxLayout(panel)
    
    scan_btn = QPushButton("AFM Scan")
    scan_btn.clicked.connect(lambda: parent_window.run_afm_scan())
    layout.addWidget(scan_btn)
    
    return {
        'panel': panel,
        'edit_modes': [],
        'view_modes': [
            ('AFM Image', parent_window.set_view_mode_afm),
        ],
        'selection_boxes': [],
    }
```

#### Foldable Panel Implementation in Qt

Qt has `QCollapsibleBox` or we can implement with `QGroupBox` + toggle button:

```python
class CollapsibleBox(QWidget):
    def __init__(self, title="", parent=None):
        super().__init__(parent)
        self.toggle_btn = QToolButton(title)
        self.toggle_btn.setCheckable(True)
        self.toggle_btn.setChecked(True)
        self.toggle_btn.setStyleSheet("QToolButton { border: none; }")
        self.toggle_btn.setToolButtonStyle(Qt.ToolButtonTextBesideIcon)
        self.toggle_btn.setArrowType(Qt.ArrowType.DownArrow)
        self.toggle_btn.clicked.connect(self.on_clicked)
        
        self.content_area = QWidget()
        self.content_area.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)
        
        layout = QVBoxLayout(self)
        layout.setSpacing(0)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.addWidget(self.toggle_btn)
        layout.addWidget(self.content_area)
        
        self.anim = QPropertyAnimation(self.content_area, b"maximumHeight")
        self.anim.setDuration(200)
        self.anim.setStartValue(0)
        self.anim.setEasingCurve(QEasingCurve.InOutQuint)
    
    def on_clicked(self):
        checked = self.toggle_btn.isChecked()
        self.toggle_btn.setArrowType(
            Qt.ArrowType.DownArrow if not checked else Qt.ArrowType.UpArrow
        )
        if checked:
            self.content_area.setMaximumHeight(16777215)  # "infinite"
            self.anim.setDirection(QPropertyAnimation.Direction.Forward)
        else:
            self.anim.setDirection(QPropertyAnimation.Direction.Backward)
        self.anim.start()
    
    def setContent(self, widget):
        layout = QVBoxLayout(self.content_area)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.addWidget(widget)
```

### Integration with KekuleExplorerGUI

#### Main Window Structure

```python
class KekuleExplorerWindow:
    def __init__(self):
        # ... existing init ...
        self.extension_manager = ExtensionManager()
        self.extension_manager.load_config()
        
        # Extension UI containers
        self.extension_panels = {}  # name -> CollapsibleBox
        self.extension_edit_modes = {}  # name -> list of (name, callback)
        self.extension_view_modes = {}  # name -> list of (name, callback)
        self.extension_selection_boxes = {}  # name -> list of (name, widget)
    
    def initUI(self):
        # ... existing base UI ...
        
        # Build extension panels dynamically
        for ext_name in self.extension_manager.get_enabled_extensions():
            self._add_extension_ui(ext_name)
    
    def _add_extension_ui(self, ext_name):
        """Add extension UI components."""
        config = self.extension_manager.get_config(ext_name)
        
        try:
            ext_instance = self.extension_manager.get(ext_name, config)
            meta = EXTENSION_REGISTRY[ext_name]
            
            # Load UI builder
            ui_module = importlib.import_module(meta['module'])
            ui_builder = getattr(ui_module, meta['ui_builder'])
            ui_components = ui_builder(self, ext_instance, config)
            
            # Add collapsible panel
            panel = ui_components['panel']
            collapsible = CollapsibleBox(ext_name.capitalize())
            collapsible.setContent(panel)
            self.side_layout.addWidget(collapsible)
            self.extension_panels[ext_name] = collapsible
            
            # Register edit modes
            for mode_name, callback in ui_components.get('edit_modes', []):
                self.edit_mode_combo.addItem(mode_name)
                self.extension_edit_modes[ext_name] = ui_components['edit_modes']
            
            # Register view modes
            for mode_name, callback in ui_components.get('view_modes', []):
                self.view_mode_combo.addItem(mode_name)
                self.extension_view_modes[ext_name] = ui_components['view_modes']
            
            # Register selection boxes
            for box_name, widget in ui_components.get('selection_boxes', []):
                # Add to appropriate layout
                self.extension_selection_boxes[ext_name] = ui_components['selection_boxes']
            
        except ExtensionLoadError as e:
            # Extension not available - show disabled placeholder
            disabled_panel = QLabel(f"{ext_name.capitalize()}: Not available\n{e}")
            disabled_panel.setStyleSheet("color: gray; font-style: italic;")
            self.side_layout.addWidget(disabled_panel)
    
    def on_edit_mode_changed(self, mode_name):
        """Handle edit mode selection."""
        # Check if mode belongs to an extension
        for ext_name, modes in self.extension_edit_modes.items():
            for name, callback in modes:
                if name == mode_name:
                    try:
                        callback()
                    except ExtensionNotLoadedError as e:
                        self.statusBar().showMessage(str(e))
                    return
        # ... existing edit mode handling ...
    
    def on_view_mode_changed(self, mode_name):
        """Handle view mode selection."""
        for ext_name, modes in self.extension_view_modes.items():
            for name, callback in modes:
                if name == mode_name:
                    try:
                        callback()
                    except ExtensionNotLoadedError as e:
                        self.statusBar().showMessage(str(e))
                    return
        # ... existing view mode handling ...
```

### Extension Usage Pattern

Extensions can call into the GUI backend directly:

```python
# In pyBall/MMFF.py (ui builder)
def build_mmff_ui(parent_window, mmff_lib, config):
    """Build MMFF UI panel."""
    panel = QWidget()
    layout = QVBoxLayout(panel)
    
    relax_btn = QPushButton("Real-time Relax")
    relax_btn.clicked.connect(lambda: parent_window.toggle_realtime_relax())
    layout.addWidget(relax_btn)
    
    return {
        'panel': panel,
        'edit_modes': [
            ('Relaxed Scan', parent_window.set_edit_mode_relaxed_scan),
        ],
        'view_modes': [],
        'selection_boxes': [],
    }

# In KekuleExplorerGUI
def toggle_realtime_relax(self):
    """Toggle real-time relaxation using MMFF."""
    mmff = self.extension_manager.get('mmff', self.extension_manager.get_config('mmff'))
    
    # Directly use MMFF API - no standardization enforced
    atypes = np.array([elements.ELEMENT_DICT[e][0] for e in self.backend.sys.enames])
    apos = np.array(self.backend.sys.apos, dtype=np.float64)
    
    # Call MMFF functions directly with extension-specific data structures
    mmff.makeFFs(atypes, apos, ...)
    forces = mmff.evalForces(apos, ...)
    
    # Update GUI
    self.backend.sys.apos += forces * dt
    self.refresh_view()
```

### Implementation Plan

#### Phase 1: Core Infrastructure
1. Create `pyBall/ExtensionManager.py`:
   - `ExtensionLoader` class with lazy loading
   - `ExtensionProxy` for error handling
   - `ExtensionManager` singleton wrapper
   - `ExtensionConfig` JSON storage

2. Create `EXTENSION_REGISTRY` dict in `pyBall/ExtensionManager.py`

3. Create `CollapsibleBox` widget in `pyBall/gui_widgets.py`

#### Phase 2: Migrate Existing Extensions
1. FireCore:
   - Add to registry
   - Create `build_fireball_ui()` in FireCore or separate wrapper
   - Remove module-level import from KekuleExplorerGUI
   - Use `extension_manager.get('firecore')` in compute_orbitals()

2. DFTB:
   - Add to registry
   - Create `build_dftb_ui()` in dftb_utils or wrapper
   - Keep lazy import in KekuleBackend but use ExtensionManager for config

#### Phase 3: Add New Extensions
1. For each extension (AFM, MMFF, Grid, POVray, psi4, pyscf, MolecularDynamics):
   - Add to registry
   - Create `build_*_ui()` function
   - Implement UI builder returning panel, edit_modes, view_modes, selection_boxes

#### Phase 4: GUI Integration
1. Refactor KekuleExplorerGUI:
   - Add ExtensionManager instance
   - Dynamic panel building in [initUI()](cci:1://file:///home/prokop/git/FireCore/pyBall/KekuleExplorerGUI.py:74:4-173:56)
   - Dynamic edit/view mode registration
   - Handle ExtensionNotLoadedError in callbacks

2. Add Extension Settings dialog:
   - Enable/disable checkboxes
   - Path configuration fields
   - Validation buttons

### Key Design Principles (Revised)

1. **Flexibility over standardization:** Extensions define their own data structures and APIs
2. **Lazy loading with proxy pattern:** Safe access via ExtensionProxy that raises descriptive errors
3. **Dynamic UI registration:** Extensions contribute panels, edit modes, view modes via ui_builder functions
4. **Foldable panels:** CollapsibleBox widgets to save space
5. **Configuration separation:** JSON config for paths and settings
6. **Fail loudly:** ExtensionProxy raises errors on any access to unavailable extensions
7. **Minimal metadata:** Registry only stores what's needed for loading and UI building

### Open Questions

1. **UI builder location:** Should ui_builder functions be in the extension module itself or in a separate wrapper?
   - Option A: In extension module (e.g., `pyBall/MMFF.py:build_mmff_ui`)
   - Option B: In `pyBall/extensions/mmff_wrapper.py`
   - **Recommendation:** Option A for simple extensions, Option B for complex ones

2. **Extension instance lifecycle:** Should instances be cached or recreated on each use?
   - GPU contexts suggest caching
   - Stateful calculations suggest per-instance
   - **Recommendation:** Cache by default, allow explicit reload

3. **Settings dialog granularity:** Should each extension define its own settings UI or use a generic form?
   - **Recommendation:** Generic form based on required_paths/optional_paths from registry

This design provides maximum flexibility for extensions while still giving them the infrastructure they need for safe loading, configuration, and dynamic UI integration.