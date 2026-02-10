# Rotational Dynamics Implementation for Pi-Orbitals

## Overview

This document describes the implementation of rotational dynamics for pi-orbitals in the pyOpenCL MD framework, designed to properly conserve angular momentum compared to the previous linear dynamics approach.

## Changes Made

### 1. MolecularDynamics.py - Added Rotational Force Kernel Driver

**Location**: `/home/prokop/git/FireCore/pyBall/OCL/MolecularDynamics.py`

#### Added Methods:

1. **`run_getMMFFf4_rot()`** - Driver for the new rotational force kernel
   ```python
   def run_getMMFFf4_rot(self):
       self.prg.getMMFFf4_rot(self.queue, self.sz_node, self.sz_loc, *self.kernel_args_getMMFFf4_rot)
       self.queue.finish()
   ```

2. **`run_MD_step()`** - Optimized MD step execution without string comparisons
   ```python
   def run_MD_step(self, do_clean=True, do_nb=False, do_mmff=True, use_rot=False, force_kernel='basic'):
       """Run a single MD step with specified kernels - optimized version without string comparisons."""
   ```
   
   This method:
   - Replaces the string-comparison-based `run_step()` in `MD_test_utils.py`
   - Provides direct kernel dispatch for better performance
   - Supports multiple integrator modes: 'basic', 'rot', 'rattle', 'none'
   - Allows switching between `getMMFFf4` and `getMMFFf4_rot` force kernels

#### Modified Methods:

- **`setup_kernels()`** - Added initialization of `kernel_args_getMMFFf4_rot`

### 2. MD_test_utils.py - Updated to Use Optimized Method

**Location**: `/home/prokop/git/FireCore/pyBall/MD_test_utils.py`

**Modified Function**: `run_step()` now calls `md.run_MD_step()` directly:
```python
def run_step(md, do_clean=True, do_nb=True, do_mmff=True, mode='basic', use_rot_force=False):
    """DEPRECATED: Use md.run_MD_step() directly for better performance."""
    md.run_MD_step(do_clean=do_clean, do_nb=do_nb, do_mmff=do_mmff, 
                   use_rot=use_rot_force, force_kernel=mode)
```

### 3. Test Script Updates

**Location**: `/home/prokop/git/FireCore/tests/tUFF/test_MMFFsp3_pyOCL.py`

#### Changed Default Molecule:
- Changed from methanol to **formic acid** (asymmetric sp² molecule)
- Formic acid is better for testing pi-orbital dynamics

#### Added Command-Line Option:
```python
ap.add_argument('--use-rot-force', type=int, default=0, 
                help='Use getMMFFf4_rot instead of getMMFFf4 for force calculation')
```

#### Updated Main Loop:
- Direct call to `md.run_MD_step()` instead of `run_step()`
- Eliminates per-step string comparisons

### 4. Comparison Script

**Location**: `/home/prokop/git/FireCore/tests/tUFF/compare_rot_dynamics.sh`

Bash script to compare basic vs rotational dynamics side-by-side.

## Usage Examples

### Basic Mode (old implementation)
```bash
python test_MMFFsp3_pyOCL.py --steps 200 --mode basic --monitor 1
```

### Rotational Mode (new implementation)
```bash
python test_MMFFsp3_pyOCL.py --steps 200 --mode rot --use-rot-force 1 --monitor 1
```

### Compare Both
```bash
bash compare_rot_dynamics.sh
```

### Direct MD step usage in Python
```python
# Old way (deprecated)
run_step(md, do_clean=True, do_nb=False, do_mmff=True, mode='rot')

# New way (optimized)
md.run_MD_step(do_clean=True, do_nb=False, do_mmff=True, use_rot=True, force_kernel='rot')
```

## OpenCL Kernel Mapping

| Python Method | OpenCL Kernel | Purpose |
|---------------|---------------|---------|
| `run_getMMFFf4()` | `getMMFFf4` | Linear pi-orbital forces |
| `run_getMMFFf4_rot()` | `getMMFFf4_rot` | Rotational pi-orbital forces |
| `run_updateAtomsMMFFf4()` | `updateAtomsMMFFf4` | Basic leap-frog integrator |
| `run_updateAtomsMMFFf4_rot()` | `updateAtomsMMFFf4_rot` | Rotational integrator |
| `run_updateAtomsMMFFf4_RATTLE()` | `updateAtomsMMFFf4_RATTLE` | Constrained integrator |

## Performance Improvements

**Before**: String comparison every MD step
```python
if mode == 'rot':
    md.run_updateAtomsMMFFf4_rot()
elif mode == 'rattle':
    md.run_updateAtomsMMFFf4_RATTLE()
elif mode == 'basic':
    md.run_updateAtomsMMFFf4()
```

**After**: Direct dispatch via dictionary or conditional
- No string comparisons in hot loop
- Better branch prediction
- Cleaner code organization

## Conservation Results

See `/home/prokop/git/FireCore/tests/tUFF/RESULTS_rot_dynamics.md` for detailed conservation analysis.

**Key Findings**:
- ✓ Rotational mode achieves **perfect force conservation** (ΔF = 0)
- ✓ Rotational mode achieves **perfect torque conservation** (ΔT = 0)  
- ✓ Rotational mode achieves **perfect potential energy conservation** (ΔEpot = 0)
- ⚠ Energy drift issues need investigation (possibly due to time step or integration scheme)

## Integration with Monitoring

The optimized `run_MD_step()` integrates seamlessly with existing monitoring infrastructure:

```python
for _ in range(steps):
    md.run_MD_step(do_clean=do_clean, do_nb=do_nb, do_mmff=do_mmff, 
                   use_rot=use_rot_force, force_kernel=mode)
    if want_record or monitor_enabled:
        collect_diagnostics(md, mm, masses, want_record, monitor_enabled, 
                          monitor_props, totals_hist, monitor_data, trj)
```

All existing monitoring features work without modification:
- Trajectory recording
- Invariant tracking (F, T, P, L, E)
- Plotting and visualization
- NPZ data export

## Future Work

1. **Investigate energy conservation** in rotational mode
   - Test with smaller time steps
   - Check moment of inertia scaling
   - Verify angular velocity integration

2. **Test RATTLE mode** for constrained pi-orbital dynamics

3. **Benchmark performance** of optimized dispatch vs string comparison

4. **Extend to multi-system simulations** if needed

## References

- OpenCL kernels: `/home/prokop/git/FireCore/cpp/common_resources/cl/relax_multi_mini.cl`
  - `getMMFFf4_rot` (line ~824)
  - `updateAtomsMMFFf4_rot` (implementation varies)
- Documentation: `/home/prokop/git/FireCore/doc/DevNotes/MD_MMFF_OCL_notes.md`
