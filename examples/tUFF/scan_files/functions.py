import glob
import os
import numpy as np

def generate_confs(apos, nconf, seed, scan_atom, scan_dim, z_range, xy_range, scan_pos):
    """Generates configurations for scanning."""
    base = apos.copy()
    rng = np.random.default_rng(seed)
    confs = np.repeat(base[None, :, :], nconf, axis=0)
    # Don't add random noise when scanning
    # confs += 0.1 * rng.standard_normal(confs.shape)

    z_scan = False
    xy_scan = False

    original_pos = base[scan_atom].copy()

    if scan_dim == 'z':
        z_scan = True
        scan_params = [float(x) for x in z_range.split(',')]
        start, end, step = scan_params
        z_values = np.arange(start, end, step)
        nconf_generated = len(z_values)
        
        # Adjust nconf to match the generated number of points
        if nconf != nconf_generated:
            print(f"Warning: nconf ({nconf}) does not match the number of generated points ({nconf_generated}). Adjusting nconf.")
            nconf = nconf_generated
            confs = np.repeat(base[None, :, :], nconf, axis=0)

        scan_pos_vals = [float(x) for x in scan_pos.split(',')]
        
        for i in range(nconf):
            target_pos = np.array([scan_pos_vals[0], scan_pos_vals[1], z_values[i]])
            displacement = target_pos - original_pos
            confs[i, :, :] += displacement

    elif scan_dim == 'xy':
        xy_scan = True
        scan_params = [float(x) for x in xy_range.split(',')]
        x_start, x_end, x_step, y_start, y_end, y_step = scan_params
        
        x_values = np.arange(x_start, x_end, x_step)
        y_values = np.arange(y_start, y_end, y_step)
        
        X, Y = np.meshgrid(x_values, y_values)
        
        X_flat = X.ravel()
        Y_flat = Y.ravel()
        
        nconf_generated = len(X_flat)
        
        # Adjust nconf to match the generated number of points
        if nconf != nconf_generated:
            print(f"Warning: nconf ({nconf}) does not match the number of generated points ({nconf_generated}). Adjusting nconf.")
            nconf = nconf_generated
            confs = np.repeat(base[None, :, :], nconf, axis=0)

        scan_pos_vals = [float(x) for x in scan_pos.split(',')]
        
        for i in range(nconf):
            target_pos = np.array([X_flat[i], Y_flat[i], scan_pos_vals[0]])
            displacement = target_pos - original_pos
            confs[i, :, :] += displacement

    return confs, z_scan, xy_scan, nconf

# ==================
#  Buffer Specs (mirrored from test_UFF_ocl.py)
# ==================
# Unified buffer specification to collect and print CPU buffers.
# Format:
# 'buffer_name': {
#     'stride': canonical columns (GPU layout),
#     'cpu_stride': optional CPU columns if different,
#     'type': 'int' (topology) or 'float' (parameters)
# }
BUF_SPECS = {
    # --- Topology Buffers (Integers) ---
    'bonAtoms':  {'stride': 2, 'type': 'int'},
    'angAtoms':  {'stride': 4, 'cpu_stride': 3, 'type': 'int'},
    'dihAtoms':  {'stride': 4, 'type': 'int'},
    'invAtoms':  {'stride': 4, 'type': 'int'},
    'neighs':    {'stride': 4, 'type': 'int'},
    'neighBs':   {'stride': 4, 'type': 'int'},
    # --- Parameter Buffers (Floats) ---
    'bonParams': {'stride': 2, 'type': 'float'},
    'angParams': {'stride': 5, 'type': 'float'},
    'dihParams': {'stride': 3, 'type': 'float'},
    'invParams': {'stride': 4, 'type': 'float'},
}

TOPOLOGY_SPECS = {k: v for k, v in BUF_SPECS.items() if v['type'] == 'int'}
PARAMS_SPECS   = {k: v for k, v in BUF_SPECS.items() if v['type'] == 'float'}

UFF_COMPONENTS = ('bonds', 'angles', 'dihedrals', 'inversions')

# ==================
#  Helper Functions
# ==================
def cleanup_xyz_files(patterns):
    """Deletes files matching given glob patterns."""
    print("--- Cleaning up old trajectory files ---")
    for pattern in patterns:
        files = glob.glob(pattern)
        if not files:
            print(f"No files found for pattern: {pattern}")
        for f in files:
            try:
                os.remove(f)
                print(f"Removed: {f}")
            except OSError as e:
                print(f"Error removing file {f}: {e}")

def get_cpu_bufs(uff_obj, specs):
    """Collect CPU buffers from the given UFF object and pad to canonical stride when needed."""
    bufs = {}
    for name, spec in specs.items():
        if hasattr(uff_obj, name):
            buf = getattr(uff_obj, name)
            canonical_stride = spec['stride']
            cpu_stride = spec.get('cpu_stride', canonical_stride)
            if cpu_stride != canonical_stride:
                # e.g. angAtoms have 3 columns on CPU, 4 on GPU
                import numpy as _np
                padded = _np.full((buf.shape[0], canonical_stride), -1, dtype=_np.int32)
                padded[:, :cpu_stride] = buf
                bufs[name] = padded
            else:
                bufs[name] = buf
    return bufs

def build_component_flags(base_components=None, *, enable=None, disable=None, preset='all'):
    preset = (preset or 'all').lower()
    preset_components = {
        'all': set(UFF_COMPONENTS),
        'bonded': set(UFF_COMPONENTS),
        'bonded-only': set(UFF_COMPONENTS),
        'none': set(),
        'grid-only': set(),
        'nonbonded-only': set(),
    }

    if base_components is None:
        base_set = set(preset_components.get(preset, UFF_COMPONENTS))
    else:
        base_set = set(base_components)

    if enable:
        base_set.update(enable)
    if disable:
        base_set.difference_update(disable)

    return {name: 1 if name in base_set else 0 for name in UFF_COMPONENTS}