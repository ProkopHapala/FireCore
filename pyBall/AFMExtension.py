"""
AFMExtension.py - AFM simulation extension for KekuleExplorerGUI
===============================================================
Uses tested functions from AFM_utils.py (same as test_full_pipeline.py)
"""

import os
import numpy as np
from PyQt5 import QtWidgets, QtCore
from pyBall.GUI.CollapsibleSection import CollapsibleSection
from pyBall.ExtensionManager import UIComponents


def _get_afm_geometry(window):
    """Convert backend geometry to AFM format."""
    ELEM_Z = {'H': 1, 'C': 6, 'N': 7, 'O': 8, 'P': 15, 'S': 16, 'Br': 35, 'I': 53}
    apos = window.backend.sys.apos.astype(np.float64)
    enames = list(window.backend.sys.enames)
    atomTypes = np.array([ELEM_Z.get(e, 6) for e in enames], dtype=np.int32)
    return apos, atomTypes, enames


def _update_afm_status(window, msg):
    """Update AFM status label and status bar."""
    if hasattr(window, 'afm_status_label'):
        window.afm_status_label.setPlainText(f"Status: {msg}")
    window.statusBar().showMessage(f"AFM: {msg}")
    QtWidgets.QApplication.processEvents()


def run_afm_full_pipeline(window):
    """Run complete AFM pipeline using tested run_afm_from_xyz()."""
    try:
        from pyBall.OCL import AFM_utils as afm_utils
        from pyBall import dftb_utils as du
        import tempfile

        atomPos, atomTypes, enames = _get_afm_geometry(window)
        if len(atomPos) == 0:
            raise ValueError("No atoms in molecule")

        basis = window.afm_basis_combo.currentText()
        step = window.afm_step_spin.value()
        margin = window.afm_margin_spin.value()
        scan_range = window.afm_scan_range_spin.value()
        h_min = window.afm_hmin_spin.value()
        h_max = window.afm_hmax_spin.value()
        h_step = window.afm_hstep_spin.value()

        stm_params = None
        stm_enable = False
        if hasattr(window, 'afm_stm_enable'):
            try:
                stm_enable = bool(window.afm_stm_enable.isChecked())
            except Exception:
                stm_enable = False
        # Convenience: if user selected STM Signal for plotting, auto-enable STM compute
        if (not stm_enable) and hasattr(window, 'afm_component_combo'):
            try:
                stm_enable = (str(window.afm_component_combo.currentText()) == 'STM Signal')
            except Exception:
                pass

        print(f"[AFM] STM enabled={stm_enable}")

        if stm_enable:
            mo_indices = None
            if hasattr(window, 'afm_stm_mo_indices'):
                try:
                    mo_indices = [int(s) for s in window.afm_stm_mo_indices.text().replace(',', ' ').split() if s.strip()]
                except Exception:
                    mo_indices = None
                if mo_indices is not None and len(mo_indices) == 0:
                    mo_indices = None

            lumo_offsets = None
            if mo_indices is None:
                if hasattr(window, 'afm_stm_use_homo_range') and window.afm_stm_use_homo_range.isChecked():
                    off0 = int(window.afm_stm_homo_off0.value())
                    nsel = int(window.afm_stm_homo_n.value())
                    lumo_offsets = list(range(off0, off0 + max(1, nsel)))
                else:
                    try:
                        lumo_offsets = [int(s) for s in window.afm_stm_lumo_offsets.text().replace(',', ' ').split() if s.strip()]
                    except Exception:
                        lumo_offsets = None
                    if (lumo_offsets is None) or (len(lumo_offsets) == 0):
                        lumo_offsets = [1, 2, 3]

            stm_field = 'ldos'
            if hasattr(window, 'afm_stm_field_combo'):
                stm_field = str(window.afm_stm_field_combo.currentText())
            if stm_field != 'ldos':
                if mo_indices is not None:
                    if len(mo_indices) != 1:
                        raise ValueError(f"STM field='{stm_field}' requires exactly 1 MO index, got {mo_indices}")
                else:
                    if len(lumo_offsets) != 1:
                        raise ValueError(f"STM field='{stm_field}' requires exactly 1 MO (set nMOs=1), got lumo_offsets={lumo_offsets}")
            stm_params = {
                'compute': True,
                'lumo_offsets': lumo_offsets,
                'field': stm_field,
                'use_exp_basis': True,
                'exp_beta': float(window.afm_stm_exp_beta.value()),
                'exp_r0': float(window.afm_stm_exp_r0.value()),
                'bond_resolved': bool(window.afm_stm_bond_resolved.isChecked()),
            }
            if mo_indices is not None:
                stm_params['mo_indices'] = mo_indices

        if stm_params is not None:
            _update_afm_status(window, f"STM enabled: field={stm_params.get('field','ldos')}  {'mo_indices='+str(stm_params.get('mo_indices')) if 'mo_indices' in stm_params else 'lumo_offsets='+str(stm_params['lumo_offsets'])}")
        else:
            _update_afm_status(window, "STM disabled")

        _update_afm_status(window, "Running full AFM pipeline (DFTB+ SCF + projection + fields + relax)...")

        # Use tested run_afm_from_xyz with proper CO tip handling
        slako_prefix = du.SK_PATHS.get(basis, basis)
        output_dir = tempfile.mkdtemp(prefix='afm_gui_')
        work_dir = os.path.join(output_dir, 'dftb_work')

        # Write XYZ for run_afm_from_xyz
        from pyBall import atomicUtils as au
        xyz_path = os.path.join(output_dir, 'molecule.xyz')
        au.save_xyz(xyz_path, enames, atomPos)

        results = afm_utils.run_afm_from_xyz(
            xyz_file=xyz_path,
            output_dir=output_dir,
            basis=basis,
            slako_prefix=slako_prefix,
            work_dir=work_dir,
            step=step, margin=margin, z_extra=6.0,
            scan_range=scan_range, scan_step=0.1,
            height_range=(h_min, h_max), height_step=h_step,
            pauli_params={'A': window.afm_pauli_a_spin.value(),
                          'beta': window.afm_pauli_beta_spin.value()},
            vdw_params={'C6_CO': window.afm_vdw_c6_spin.value()},
            relax_params={'K_LAT': window.afm_klat_spin.value()},
            plot_steps=False,
            use_dense_projection=(stm_params is not None),
            ppm_mode=True
            ,stm_params=stm_params
        )

        window._afm_results = results
        # Do not overwrite heights returned by pipeline (it defines STM/AFM z-slices)
        
        # Note: run_afm_pipeline doesn't return densities in intermediates (they were inputs)
        # Store what we do get
        window._afm_potentials = {
            'E_pauli_field': results['intermediates']['E_pauli_field'],
            'E_ES_field': results['intermediates']['E_ES_field'],
            'E_vdw': results['intermediates']['E_vdw'],
            'F_total': results['intermediates']['F_total'],  # (Fx,Fy,Fz,E) from GPU
            'V_ES': results['intermediates']['V_ES'],
            'origin': results['grid_spec']['origin'],
            'step': step,
            'grid_spec': results['grid_spec']
        }

        # STM (optional)
        if stm_params is not None:
            inter = results.get('intermediates', {})
            if 'stm_grid' in inter:
                window._afm_results['stm_grid'] = inter['stm_grid']
            if 'stm_meta' in inter:
                window._afm_results['stm_meta'] = inter['stm_meta']
                if hasattr(window, 'afm_stm_info_label'):
                    try:
                        m = inter['stm_meta']
                        txt = f"nMO={m.get('nmo','?')}  nOcc={m.get('nocc','?')}  HOMO={m.get('homo','?')}  LUMO={m.get('lumo','?')}\nE(HOMO)={m.get('E_homo','?')} eV  E(LUMO)={m.get('E_lumo','?')} eV\nMOs={m.get('mo_list','?')}  mode={m.get('mode','?')}  field={m.get('field','?')}"
                        window.afm_stm_info_label.setText(txt)
                    except Exception:
                        pass
        
        # Also compute and store density for visualization (re-run DFTB quickly)
        _update_afm_status(window, "Computing density for visualization...")
        d = afm_utils.get_density_from_dftb_plus(
            atomPos, atomTypes, basis, slako_prefix, work_dir,
            step=step, margin=margin, z_extra=6.0, verbosity=0
        )
        window._afm_density = d

        nz = results['df'].shape[2]
        # Set z-height to middle of actual returned scan heights for viewing
        heights = window._afm_results.get('heights', None)
        if heights is not None and len(heights) > 0:
            window.afm_z_height_spin.setValue(float(heights[len(heights)//2]))
        else:
            mid_z = (window.afm_hmin_spin.value() + window.afm_hmax_spin.value()) / 2
            window.afm_z_height_spin.setValue(mid_z)

        msg = f"AFM complete: {nz} slices, df range [{results['df'].min():.2f}, {results['df'].max():.2f}] Hz"
        _update_afm_status(window, msg)

    except Exception as e:
        _update_afm_status(window, f"FAILED: {e}")
        raise


def run_afm_density(window):
    """Run density projection only using tested get_density_from_dftb_plus()."""
    try:
        from pyBall.OCL import AFM_utils as afm_utils
        from pyBall import dftb_utils as du

        atomPos, atomTypes, enames = _get_afm_geometry(window)
        if len(atomPos) == 0:
            raise ValueError("No atoms in molecule")

        basis = window.afm_basis_combo.currentText()
        step = window.afm_step_spin.value()
        margin = window.afm_margin_spin.value()
        slako_prefix = du.SK_PATHS.get(basis, basis)

        _update_afm_status(window, "Computing DFTB+ density (this may take a while)...")

        import tempfile
        work_dir = tempfile.mkdtemp(prefix='afm_density_')

        # Use the tested function from AFM_utils
        d = afm_utils.get_density_from_dftb_plus(
            atomPos, atomTypes, basis, slako_prefix, work_dir,
            step=step, margin=margin, z_extra=6.0, verbosity=1
        )

        window._afm_density = d
        window._afm_potentials = None
        window._afm_results = None

        msg = f"Density computed: grid {d['ngrid']}, {len(atomPos)} atoms, step {step:.2f} A"
        _update_afm_status(window, msg)

    except Exception as e:
        _update_afm_status(window, f"FAILED: {e}")
        raise


def run_afm_potentials(window):
    """Compute potentials from density using tested pipeline."""
    try:
        if window._afm_density is None:
            raise ValueError("Density not computed. Run 'Project Density' first.")

        from pyBall.OCL import AFM_utils as afm_utils
        from pyBall.OCL import AFM as afm_mod
        import tempfile

        _update_afm_status(window, "Computing Pauli/ES/vdW potentials (need valid CO tip)...")

        # Get density data
        rho_scf = window._afm_density['rho_scf']
        rho_na = window._afm_density['rho_na']
        rho_diff = window._afm_density['rho_diff']
        V_ES = window._afm_density['V_ES']
        origin = window._afm_density['origin']
        ngrid = window._afm_density['ngrid']
        step = float(window._afm_density['grid_spec']['dA'][0])
        margin = window.afm_margin_spin.value()  # Need margin for CO tip cache key

        # Get proper CO tip using AFM_utils functions (same as run_afm_from_xyz)
        _update_afm_status(window, "Getting CO tip (checking cache or computing)...")
        
        # Setup paths like run_afm_from_xyz does
        _ROOT = os.path.normpath(os.path.join(os.path.dirname(__file__), '..'))
        fdata_dir = os.path.join(_ROOT, 'tests', 'pyFireball', 'Fdata')
        fdata_basis = os.path.join(fdata_dir, 'basis')
        
        target_shape = tuple(ngrid)
        
        # Try to get cached CO tip first
        co_tip_cached = afm_utils._get_cached_co_tip(step, margin, fdata_dir, fdata_basis)
        
        if co_tip_cached is not None:
            _update_afm_status(window, "Using cached CO tip...")
            co_rho_total_raw, co_rho_delta_raw = co_tip_cached
        else:
            _update_afm_status(window, "Computing CO tip on-the-fly (first time, may take a while)...")
            # Compute CO tip using same method as run_afm_from_xyz
            import tempfile
            co_tip_work = tempfile.mkdtemp(prefix='co_tip_')
            co_grid_spec, co_ngrid, co_origin = afm_utils._compute_co_tip_grid(step=step, margin=margin)
            afm_utils._call_compute_co_tip_script(co_tip_work, co_grid_spec, step, 100, fdata_dir, fdata_basis)
            co_rho_total_raw = np.load(os.path.join(co_tip_work, 'co_rho_total.npy'))
            co_rho_delta_raw = np.load(os.path.join(co_tip_work, 'co_rho_delta.npy'))
            # Cache for future
            afm_utils._save_cached_co_tip(co_rho_total_raw, co_rho_delta_raw, step, margin, fdata_dir, fdata_basis)
            _update_afm_status(window, "CO tip computed and cached.")
        
        # Pad and roll to match sample grid
        co_rho_total = afm_utils._pad_and_roll_co_tip(co_rho_total_raw, target_shape)
        co_rho_delta = afm_utils._pad_and_roll_co_tip(co_rho_delta_raw, target_shape)
        
        pauli_A = window.afm_pauli_a_spin.value()
        pauli_beta = window.afm_pauli_beta_spin.value()

        # Compute Pauli (using tested function from AFM.py)
        E_pauli, grads_pauli = afm_mod.compute_pauli_field(
            rho_scf, co_rho_total, step,
            A_pauli=pauli_A, beta_pauli=pauli_beta, tip_rolled=False
        )

        # Compute ES (using tested function)
        E_ES, grads_ES = afm_mod.compute_es_conv_field(
            V_ES, co_rho_delta, step, tip_rolled=False, return_grads=True
        )

        # Compute vdW (using tested function)
        atomPos, atomTypes, _ = _get_afm_geometry(window)
        E_vdw, grads_vdw = afm_mod.compute_dispersion_grid(
            atomPos, atomTypes, origin, step, ngrid,
            C6_CO=window.afm_vdw_c6_spin.value(), use_opencl=True
        )

        window._afm_potentials = {
            'E_pauli_field': E_pauli,
            'grads_pauli': grads_pauli,
            'E_ES_field': E_ES,
            'grads_ES': grads_ES,
            'E_vdw': E_vdw,
            'grads_vdw': grads_vdw,
            'origin': origin,
            'step': step,
            'grid_spec': window._afm_density['grid_spec']
        }

        msg = f"Potentials: Pauli[{E_pauli.min():.2f},{E_pauli.max():.2f}] ES[{E_ES.min():.2f},{E_ES.max():.2f}] vdW[{E_vdw.min():.2f},{E_vdw.max():.2f}]"
        _update_afm_status(window, msg)

    except Exception as e:
        _update_afm_status(window, f"FAILED: {e}")
        raise


def run_afm_relaxation(window):
    """Run probe-particle relaxation using tested compose_and_relax_total()."""
    try:
        if window._afm_potentials is None:
            raise ValueError("Potentials not computed. Run 'Potentials' first.")

        from pyBall.OCL import AFM_utils as afm_utils

        atomPos, _, _ = _get_afm_geometry(window)

        scan_range = window.afm_scan_range_spin.value()
        h_min = window.afm_hmin_spin.value()
        h_max = window.afm_hmax_spin.value()
        h_step = window.afm_hstep_spin.value()

        # Get total gradient - either pre-combined (full pipeline) or compose from parts (staged)
        if 'grads_total' in window._afm_potentials:
            grads_total = window._afm_potentials['grads_total']
        else:
            # Staged mode: compose from individual grads
            grads_total = (window._afm_potentials['grads_pauli'] +
                          window._afm_potentials['grads_ES'] +
                          window._afm_potentials['grads_vdw'])

        origin = window._afm_potentials['origin']
        step = window._afm_potentials['step']

        # Setup scan grid
        x_min, x_max = atomPos[:, 0].min(), atomPos[:, 0].max()
        y_min, y_max = atomPos[:, 1].min(), atomPos[:, 1].max()

        scan_xs = np.linspace(x_min - scan_range, x_max + scan_range,
                              int((x_max - x_min + 2 * scan_range) / 0.1))
        scan_ys = np.linspace(y_min - scan_range, y_max + scan_range,
                              int((y_max - y_min + 2 * scan_range) / 0.1))
        heights = np.arange(h_min, h_max, h_step)

        _update_afm_status(window, f"Running PP relaxation on {len(scan_xs)}x{len(scan_ys)} grid...")

        # Use tested compose_and_relax_total
        df, tip_disp = afm_utils.compose_and_relax_total(
            grads_total, scan_xs, scan_ys, heights,
            origin, step, atomPos,
            K_LAT=window.afm_klat_spin.value(),
            use_gpu_relax=True, ppm_mode=True
        )

        if window._afm_results is None:
            window._afm_results = {}
        window._afm_results['df'] = df
        window._afm_results['tip_disp'] = tip_disp
        window._afm_results['scan_xs'] = scan_xs
        window._afm_results['scan_ys'] = scan_ys
        window._afm_results['heights'] = heights

        nz = df.shape[2]
        # Set z-height to middle of scan range for viewing
        mid_z = (window.afm_hmin_spin.value() + window.afm_hmax_spin.value()) / 2
        window.afm_z_height_spin.setValue(mid_z)

        msg = f"Relaxation complete: df range [{df.min():.2f}, {df.max():.2f}] Hz"
        _update_afm_status(window, msg)

    except Exception as e:
        _update_afm_status(window, f"FAILED: {e}")
        raise


def _get_z_slice(grid_spec, step, z_height):
    """Convert physical z-height to grid index."""
    origin_z = grid_spec['origin'][2]
    iz = int(np.clip(np.round((z_height - origin_z) / step), 0, grid_spec['ngrid'][2] - 1))
    actual_z = origin_z + iz * step
    return iz, actual_z


def plot_afm_slice(window):
    """Plot single z-slice in a GUI window (not to disk)."""
    try:
        import matplotlib
        matplotlib.use('Qt5Agg')
        import matplotlib.pyplot as plt
        from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
        from matplotlib.figure import Figure

        component = window.afm_component_combo.currentText()
        z_height = window.afm_z_height_spin.value()
        auto_limits = window.afm_auto_limits.isChecked()

        # Determine data source and get grid info
        if component == "AFM Image (df)":
            if window._afm_results is None or 'df' not in window._afm_results:
                raise ValueError("No AFM results. Run full pipeline or relaxation first.")
            data_3d = window._afm_results['df']
            # AFM df uses scan heights directly, not grid origin
            heights = window._afm_results.get('heights', [])
            if len(heights) == 0:
                raise ValueError("No heights in AFM results")
            # Find closest height index
            h_idx = np.argmin(np.abs(heights - z_height))
            actual_z = heights[h_idx]
            iz = h_idx
            step = heights[1] - heights[0] if len(heights) > 1 else 0.1
            cmap = 'afmhot'
            symmetric = False
            data_label = "Frequency Shift (Hz)"
            # Extract slice data
            data = data_3d[:, :, iz]

        elif component == "STM Signal":
            if window._afm_results is None or 'stm_grid' not in window._afm_results:
                raise ValueError("No STM results. Enable 'STM' and run full AFM pipeline.")
            data_3d = window._afm_results['stm_grid']
            heights = window._afm_results.get('heights', [])
            if len(heights) == 0:
                raise ValueError("No heights in AFM results")
            h_idx = np.argmin(np.abs(heights - z_height))
            actual_z = heights[h_idx]
            iz = h_idx
            step = heights[1] - heights[0] if len(heights) > 1 else 0.1
            cmap = 'viridis'
            symmetric = False
            data_label = "STM Signal (arb.)"
            data = data_3d[:, :, iz]
            stm_meta = window._afm_results.get('stm_meta', None)
            if stm_meta is not None:
                print(f"[STM Plot] Z={z_height:.2f}A -> iz={iz}, actual_z={actual_z:.2f}A  field={stm_meta.get('field',None)}  MOs={stm_meta.get('mo_list',None)}  HOMO={stm_meta.get('homo',None)} LUMO={stm_meta.get('lumo',None)}")
            
        elif component in ["SCF Density", "Neutral Density", "Delta Density"]:
            if window._afm_density is None:
                raise ValueError("Density not computed. Run 'Project Density' first.")
            density_map = {
                "SCF Density": ("rho_scf", "viridis", False, "SCF Density"),
                "Neutral Density": ("rho_na", "viridis", False, "Neutral Density"),
                "Delta Density": ("rho_diff", "seismic", True, "Delta Density")
            }
            key, cmap, symmetric, data_label = density_map[component]
            data_3d = window._afm_density[key]
            grid_spec = window._afm_density['grid_spec']
            step = float(grid_spec['dA'][0])
            # Get slice at requested z-height
            iz, actual_z = _get_z_slice(grid_spec, step, z_height)
            # Extract slice data
            data = data_3d[:, :, iz]
            
        elif component in ["Pauli Energy", "Electrostatic Energy", "vdW Energy"]:
            if window._afm_potentials is None:
                raise ValueError("Potentials not computed. Run 'Potentials' first.")
            field_map = {
                "Pauli Energy": ("E_pauli_field", "seismic", True, "Pauli Energy (eV)"),
                "Electrostatic Energy": ("E_ES_field", "seismic", True, "ES Energy (eV)"),
                "vdW Energy": ("E_vdw", "seismic", True, "vdW Energy (eV)")
            }
            key, cmap, symmetric, data_label = field_map[component]
            data_3d = window._afm_potentials[key]
            grid_spec = window._afm_potentials['grid_spec']
            step = window._afm_potentials['step']
            # Get slice at requested z-height
            iz, actual_z = _get_z_slice(grid_spec, step, z_height)
            # Extract slice data
            data = data_3d[:, :, iz]
            
        else:  # Total Potential or Total Z-Force
            if window._afm_potentials is None or window._afm_potentials.get('F_total') is None:
                raise ValueError("Force field data not available. Run full pipeline first.")
            F_total = window._afm_potentials['F_total']  # (Fx,Fy,Fz,E)
            print(f"[AFM Plot] F_total shape: {F_total.shape}")
            grid_spec = window._afm_potentials['grid_spec']
            step = window._afm_potentials['step']
            iz, actual_z = _get_z_slice(grid_spec, step, z_height)
            
            if component == "Total Potential":
                # F_total[..., 3] is the energy E
                data_3d = F_total[..., 3]
                cmap = "seismic"
                symmetric = True
                data_label = "Total Potential (eV)"
            else:  # Total Z-Force
                # F_total[..., 2] is Fz, negate so repulsive = positive = red
                data_3d = -F_total[..., 2]
                cmap = "seismic"
                symmetric = True
                data_label = "Total Z-Force (eV/Ang)"
            
            # Extract slice data
            data = data_3d[:, :, iz]
        
        # Data range info
        data_min, data_max = data.min(), data.max()
        data_mean = data.mean()
        
        # Debug output for slice
        print(f"[AFM Plot] Z={z_height:.2f}A -> iz={iz}, actual_z={actual_z:.2f}A, range=[{data_min:.3f},{data_max:.3f}], mean={data_mean:.3f}")

        # Create figure
        fig = Figure(figsize=(7, 6), dpi=100)
        ax = fig.add_subplot(111)

        # Determine vmin/vmax
        if auto_limits:
            if symmetric:
                vmax = np.max(np.abs(data))
                vmin = -vmax
            else:
                vmin, vmax = data_min, data_max
        else:
            vmin = window.afm_vmin_spin.value()
            vmax = window.afm_vmax_spin.value()
        
        print(f"[AFM Plot] colormap: vmin={vmin:.3f}, vmax={vmax:.3f}, auto={auto_limits}, sym={symmetric}")

        im = ax.imshow(data.T, origin='lower', cmap=cmap, vmin=vmin, vmax=vmax)

        # Title with data range
        title = f"{component}\nZ={actual_z:.2f}A (iz={iz}) | Range: [{data_min:.3f}, {data_max:.3f}] Mean: {data_mean:.3f}"
        if component == 'STM Signal':
            stm_meta = window._afm_results.get('stm_meta', None) if (window._afm_results is not None) else None
            if stm_meta is not None:
                title = f"{component}  field={stm_meta.get('field',None)}  MOs={stm_meta.get('mo_list',None)}\nZ={actual_z:.2f}A (iz={iz}) | Range: [{data_min:.3f}, {data_max:.3f}] Mean: {data_mean:.3f}"
        ax.set_title(title, fontsize=10)
        
        fig.colorbar(im, ax=ax, label=data_label)

        # Show/update in window - reuse existing window if available
        if not hasattr(window, '_afm_plot_window') or window._afm_plot_window is None:
            # Create new window
            window._afm_plot_window = QtWidgets.QDialog(window)
            window._afm_plot_window.setWindowTitle("AFM Slice Viewer")
            window._afm_plot_layout = QtWidgets.QVBoxLayout(window._afm_plot_window)
            window._afm_plot_window.resize(700, 600)
            # Handle window close to reset reference
            def on_plot_window_closed():
                window._afm_plot_window = None
            window._afm_plot_window.finished.connect(on_plot_window_closed)
            window._afm_plot_window.show()
        else:
            # Clear existing layout
            while window._afm_plot_layout.count():
                item = window._afm_plot_layout.takeAt(0)
                if item.widget():
                    item.widget().deleteLater()
        
        # Add new canvas
        canvas = FigureCanvas(fig)
        if hasattr(window, 'install_mpl_canvas_screenshot_menu'):
            try:
                window.install_mpl_canvas_screenshot_menu(canvas, fig, default_name=f"{component.replace(' ','_')}.png")
            except Exception:
                pass
        window._afm_plot_layout.addWidget(canvas)
        window._afm_plot_window.setWindowTitle(f"AFM Slice - {component} Z={z_height:.2f}A")
        window._afm_plot_window.show()
        window._afm_plot_window.raise_()
        window._afm_plot_window.activateWindow()

        window.statusBar().showMessage(f"Showing {component} at Z={z_height:.2f}A (range: [{data_min:.3f}, {data_max:.3f}])")

    except Exception as e:
        raise RuntimeError(f"Plot FAILED: {e}")


def plot_afm_diagnostic_panel(window):
    """Plot diagnostic panel with all field components in GUI window."""
    try:
        if window._afm_potentials is None:
            raise ValueError("Potentials not computed. Run 'Potentials' first.")

        import matplotlib
        matplotlib.use('Qt5Agg')
        import matplotlib.pyplot as plt
        from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
        from matplotlib.figure import Figure

        E_pauli = window._afm_potentials['E_pauli_field']
        E_ES = window._afm_potentials['E_ES_field']
        E_vdw = window._afm_potentials['E_vdw']
        E_total = E_pauli + E_ES + E_vdw

        # Create 4-panel figure
        fig = Figure(figsize=(14, 10), dpi=100)

        fields = [
            (E_total, 'Total', 'afmhot', False, "eV"),
            (E_pauli, 'Pauli', 'seismic', True, "eV"),
            (E_ES, 'Electrostatics', 'seismic', True, "eV"),
            (E_vdw, 'vdW', 'seismic', True, "eV")
        ]

        # Get z-height from UI and convert to index
        z_height = window.afm_z_height_spin.value()
        grid_spec = window._afm_potentials['grid_spec']
        step = window._afm_potentials['step']
        iz, actual_z = _get_z_slice(grid_spec, step, z_height)

        for i, (field, name, cmap, sym, unit) in enumerate(fields):
            ax = fig.add_subplot(2, 2, i + 1)
            data = field[:, :, iz]
            data_min, data_max = data.min(), data.max()
            data_mean = data.mean()
            
            if sym:
                vmax = np.max(np.abs(data))
                im = ax.imshow(data.T, origin='lower', cmap=cmap, vmin=-vmax, vmax=vmax)
            else:
                im = ax.imshow(data.T, origin='lower', cmap=cmap)
            
            title = f"{name}\nZ={actual_z:.2f}A | Range: [{data_min:.3f}, {data_max:.3f}] {unit} | Mean: {data_mean:.3f}"
            ax.set_title(title, fontsize=9)
            fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)

        fig.suptitle(f"AFM Energy Components at Z={actual_z:.2f}A", fontsize=12, y=1.02)
        fig.tight_layout()

        # Show in separate window
        canvas = FigureCanvas(fig)
        plot_window = QtWidgets.QDialog(window)
        plot_window.setWindowTitle("AFM Diagnostic Panel")
        layout = QtWidgets.QVBoxLayout(plot_window)
        layout.addWidget(canvas)
        plot_window.resize(1000, 900)
        plot_window.show()

        if not hasattr(window, '_afm_plot_windows'):
            window._afm_plot_windows = []
        window._afm_plot_windows.append(plot_window)

        window.statusBar().showMessage(f"Diagnostic panel shown (iz={iz})")

    except Exception as e:
        raise RuntimeError(f"Diagnostic plot FAILED: {e}")


def build_ui(window):
    """Build AFM panel for KekuleExplorerGUI."""
    panel = QtWidgets.QWidget()
    layout = QtWidgets.QVBoxLayout(panel)
    layout.setSpacing(3)
    layout.setContentsMargins(2, 2, 2, 2)

    # Main pipeline buttons
    full_btn = QtWidgets.QPushButton("Run Full AFM Pipeline")
    full_btn.clicked.connect(lambda: run_afm_full_pipeline(window))
    layout.addWidget(full_btn)

    # Staged computation buttons
    stage_layout = QtWidgets.QHBoxLayout()
    density_btn = QtWidgets.QPushButton("1. Project Density")
    density_btn.clicked.connect(lambda: run_afm_density(window))
    stage_layout.addWidget(density_btn)

    potential_btn = QtWidgets.QPushButton("2. Potentials")
    potential_btn.clicked.connect(lambda: run_afm_potentials(window))
    stage_layout.addWidget(potential_btn)

    relax_btn = QtWidgets.QPushButton("3. PP Relax")
    relax_btn.clicked.connect(lambda: run_afm_relaxation(window))
    stage_layout.addWidget(relax_btn)
    layout.addLayout(stage_layout)

    # Error/Status display (copyable)
    window.afm_status_label = QtWidgets.QPlainTextEdit()
    window.afm_status_label.setPlaceholderText("Status messages will appear here...")
    window.afm_status_label.setMaximumHeight(80)
    window.afm_status_label.setReadOnly(True)
    layout.addWidget(window.afm_status_label)

    # Parameters section
    param_sec = CollapsibleSection("Parameters", collapsed=True, parent=panel)
    param_widget = QtWidgets.QWidget()
    param_layout = QtWidgets.QVBoxLayout(param_widget)
    param_layout.setSpacing(2)
    param_layout.setContentsMargins(0, 0, 0, 0)

    # Density parameters
    density_group = QtWidgets.QGroupBox("Density")
    density_grid = QtWidgets.QGridLayout(density_group)
    density_grid.addWidget(QtWidgets.QLabel("Basis:"), 0, 0)
    window.afm_basis_combo = QtWidgets.QComboBox()
    window.afm_basis_combo.addItems(["mio-1-1", "3ob-3-1"])
    density_grid.addWidget(window.afm_basis_combo, 0, 1)

    density_grid.addWidget(QtWidgets.QLabel("Step:"), 1, 0)
    window.afm_step_spin = QtWidgets.QDoubleSpinBox()
    window.afm_step_spin.setRange(0.05, 0.5)
    window.afm_step_spin.setValue(0.1)
    window.afm_step_spin.setSingleStep(0.05)
    density_grid.addWidget(window.afm_step_spin, 1, 1)

    density_grid.addWidget(QtWidgets.QLabel("Margin:"), 2, 0)
    window.afm_margin_spin = QtWidgets.QDoubleSpinBox()
    window.afm_margin_spin.setRange(2.0, 10.0)
    window.afm_margin_spin.setValue(4.0)
    density_grid.addWidget(window.afm_margin_spin, 2, 1)
    param_layout.addWidget(density_group)

    # Scan parameters
    scan_group = QtWidgets.QGroupBox("Scan")
    scan_grid = QtWidgets.QGridLayout(scan_group)
    scan_grid.addWidget(QtWidgets.QLabel("Range:"), 0, 0)
    window.afm_scan_range_spin = QtWidgets.QDoubleSpinBox()
    window.afm_scan_range_spin.setRange(1.0, 10.0)
    window.afm_scan_range_spin.setValue(3.0)
    scan_grid.addWidget(window.afm_scan_range_spin, 0, 1)

    scan_grid.addWidget(QtWidgets.QLabel("H min:"), 1, 0)
    window.afm_hmin_spin = QtWidgets.QDoubleSpinBox()
    window.afm_hmin_spin.setRange(1.5, 5.0)
    window.afm_hmin_spin.setValue(2.8)
    scan_grid.addWidget(window.afm_hmin_spin, 1, 1)

    scan_grid.addWidget(QtWidgets.QLabel("H max:"), 2, 0)
    window.afm_hmax_spin = QtWidgets.QDoubleSpinBox()
    window.afm_hmax_spin.setRange(3.0, 8.0)
    window.afm_hmax_spin.setValue(3.6)
    scan_grid.addWidget(window.afm_hmax_spin, 2, 1)

    scan_grid.addWidget(QtWidgets.QLabel("H step:"), 3, 0)
    window.afm_hstep_spin = QtWidgets.QDoubleSpinBox()
    window.afm_hstep_spin.setRange(0.05, 0.3)
    window.afm_hstep_spin.setValue(0.1)
    scan_grid.addWidget(window.afm_hstep_spin, 3, 1)
    param_layout.addWidget(scan_group)

    # Physics parameters (using fitted defaults from test_full_pipeline.py)
    physics_group = QtWidgets.QGroupBox("Physics")
    physics_grid = QtWidgets.QGridLayout(physics_group)
    
    # Pauli fitted defaults: mio-1-1: A=787.22, beta=1.2371; 3ob-3-1: A=509.28, beta=1.0586
    physics_grid.addWidget(QtWidgets.QLabel("Pauli A:"), 0, 0)
    window.afm_pauli_a_spin = QtWidgets.QDoubleSpinBox()
    window.afm_pauli_a_spin.setRange(0.1, 2000.0)
    window.afm_pauli_a_spin.setValue(787.22)  # Default for mio-1-1
    window.afm_pauli_a_spin.setDecimals(2)
    physics_grid.addWidget(window.afm_pauli_a_spin, 0, 1)

    physics_grid.addWidget(QtWidgets.QLabel("Beta:"), 1, 0)
    window.afm_pauli_beta_spin = QtWidgets.QDoubleSpinBox()
    window.afm_pauli_beta_spin.setRange(0.5, 3.0)
    window.afm_pauli_beta_spin.setValue(1.2371)  # Default for mio-1-1
    window.afm_pauli_beta_spin.setDecimals(4)
    physics_grid.addWidget(window.afm_pauli_beta_spin, 1, 1)
    
    # Auto-update Pauli defaults when basis changes
    def on_basis_changed(idx):
        basis = window.afm_basis_combo.currentText()
        if basis == 'mio-1-1':
            window.afm_pauli_a_spin.setValue(787.22)
            window.afm_pauli_beta_spin.setValue(1.2371)
        elif basis == '3ob-3-1':
            window.afm_pauli_a_spin.setValue(509.28)
            window.afm_pauli_beta_spin.setValue(1.0586)
    window.afm_basis_combo.currentIndexChanged.connect(on_basis_changed)

    physics_grid.addWidget(QtWidgets.QLabel("C6:"), 2, 0)
    window.afm_vdw_c6_spin = QtWidgets.QDoubleSpinBox()
    window.afm_vdw_c6_spin.setRange(10.0, 100.0)
    window.afm_vdw_c6_spin.setValue(30.0)
    physics_grid.addWidget(window.afm_vdw_c6_spin, 2, 1)

    physics_grid.addWidget(QtWidgets.QLabel("K_LAT:"), 3, 0)
    window.afm_klat_spin = QtWidgets.QDoubleSpinBox()
    window.afm_klat_spin.setRange(0.1, 2.0)
    window.afm_klat_spin.setValue(0.5)
    physics_grid.addWidget(window.afm_klat_spin, 3, 1)
    param_layout.addWidget(physics_group)

    param_sec.setContent(param_widget)
    layout.addWidget(param_sec)

    # Visualization section
    viz_sec = CollapsibleSection("Visualization", collapsed=True, parent=panel)
    viz_widget = QtWidgets.QWidget()
    viz_layout = QtWidgets.QVBoxLayout(viz_widget)
    viz_layout.setSpacing(2)

    window.afm_component_combo = QtWidgets.QComboBox()
    window.afm_component_combo.addItems([
        "AFM Image (df)",           # Default - the actual AFM result
        "STM Signal",               # Optional (requires compute enabled in full pipeline)
        "SCF Density",              # Debug: rho_scf
        "Neutral Density",          # Debug: rho_na
        "Delta Density",            # Debug: rho_diff
        "Pauli Energy",             # Debug: E_pauli
        "Electrostatic Energy",     # Debug: E_es
        "vdW Energy",               # Debug: E_vdw
        "Total Potential",          # Debug: E (from grads_total[...,3])
        "Total Z-Force"             # Debug: Fz (from grads_total[...,2])
    ])
    viz_layout.addWidget(QtWidgets.QLabel("Component:"))
    viz_layout.addWidget(window.afm_component_combo)
    
    # Auto-set z-height when component changes
    def on_component_changed(idx):
        component = window.afm_component_combo.currentText()
        if "Density" in component:
            window.afm_z_height_spin.setValue(1.5)  # 1.5A for density (near molecule)
        else:
            window.afm_z_height_spin.setValue(3.0)  # 3A for AFM/fields (above molecule)
    window.afm_component_combo.currentIndexChanged.connect(on_component_changed)

    # Z-height above molecule (Angstroms)
    z_layout = QtWidgets.QHBoxLayout()
    z_layout.addWidget(QtWidgets.QLabel("Z-height (A):"))
    window.afm_z_height_spin = QtWidgets.QDoubleSpinBox()
    window.afm_z_height_spin.setRange(-2.0, 10.0)
    window.afm_z_height_spin.setValue(3.0)  # Default 3A for AFM
    window.afm_z_height_spin.setSingleStep(0.1)
    window.afm_z_height_spin.setDecimals(2)
    z_layout.addWidget(window.afm_z_height_spin)
    
    # Live update checkbox (enabled by default)
    window.afm_live_update = QtWidgets.QCheckBox("Live")
    window.afm_live_update.setChecked(True)
    window.afm_live_update.setToolTip("Auto-update plot when z-height changes")
    z_layout.addWidget(window.afm_live_update)
    viz_layout.addLayout(z_layout)
    
    # Connect live update
    def on_z_height_changed():
        if window.afm_live_update.isChecked():
            # Check if we have data to plot
            has_data = (window._afm_results is not None and (('df' in window._afm_results) or ('stm_grid' in window._afm_results))) or \
                      window._afm_potentials is not None or window._afm_density is not None
            if has_data:
                try:
                    plot_afm_slice(window)
                except Exception as e:
                    pass  # Silently ignore if can't plot yet
    window.afm_z_height_spin.valueChanged.connect(on_z_height_changed)
    
    # Colormap limits (auto by default, but can override)
    lim_layout = QtWidgets.QHBoxLayout()
    window.afm_auto_limits = QtWidgets.QCheckBox("Auto limits")
    window.afm_auto_limits.setChecked(True)
    lim_layout.addWidget(window.afm_auto_limits)
    
    window.afm_vmin_spin = QtWidgets.QDoubleSpinBox()
    window.afm_vmin_spin.setRange(-1000, 1000)
    window.afm_vmin_spin.setValue(-1.0)
    window.afm_vmin_spin.setEnabled(False)
    window.afm_vmin_spin.setDecimals(3)
    lim_layout.addWidget(QtWidgets.QLabel("vmin:"))
    lim_layout.addWidget(window.afm_vmin_spin)
    
    window.afm_vmax_spin = QtWidgets.QDoubleSpinBox()
    window.afm_vmax_spin.setRange(-1000, 1000)
    window.afm_vmax_spin.setValue(1.0)
    window.afm_vmax_spin.setEnabled(False)
    window.afm_vmax_spin.setDecimals(3)
    lim_layout.addWidget(QtWidgets.QLabel("vmax:"))
    lim_layout.addWidget(window.afm_vmax_spin)
    viz_layout.addLayout(lim_layout)
    
    # Enable/disable manual limits
    def on_auto_changed(state):
        window.afm_vmin_spin.setEnabled(not state)
        window.afm_vmax_spin.setEnabled(not state)
    window.afm_auto_limits.stateChanged.connect(on_auto_changed)

    plot_btn = QtWidgets.QPushButton("Plot Slice (in window)")
    plot_btn.clicked.connect(lambda: plot_afm_slice(window))
    viz_layout.addWidget(plot_btn)

    diag_btn = QtWidgets.QPushButton("Diagnostic Panel (in window)")
    diag_btn.clicked.connect(lambda: plot_afm_diagnostic_panel(window))
    viz_layout.addWidget(diag_btn)

    viz_sec.setContent(viz_widget)
    layout.addWidget(viz_sec)

    # STM section (top-level, separate from Visualization)
    stm_sec = CollapsibleSection("STM", collapsed=True, parent=panel)
    stm_widget = QtWidgets.QWidget()
    stm_grid = QtWidgets.QGridLayout(stm_widget)
    window.afm_stm_enable = QtWidgets.QCheckBox("Compute STM")
    window.afm_stm_enable.setChecked(False)
    stm_grid.addWidget(window.afm_stm_enable, 0, 0, 1, 2)

    window.afm_stm_use_homo_range = QtWidgets.QCheckBox("Use HOMO range")
    window.afm_stm_use_homo_range.setChecked(True)
    stm_grid.addWidget(window.afm_stm_use_homo_range, 1, 0, 1, 2)

    stm_grid.addWidget(QtWidgets.QLabel("HOMO off0:"), 2, 0)
    window.afm_stm_homo_off0 = QtWidgets.QSpinBox()
    window.afm_stm_homo_off0.setRange(-50, 50)
    window.afm_stm_homo_off0.setValue(1)
    stm_grid.addWidget(window.afm_stm_homo_off0, 2, 1)
    stm_grid.addWidget(QtWidgets.QLabel("nMOs:"), 3, 0)
    window.afm_stm_homo_n = QtWidgets.QSpinBox()
    window.afm_stm_homo_n.setRange(1, 50)
    window.afm_stm_homo_n.setValue(3)
    stm_grid.addWidget(window.afm_stm_homo_n, 3, 1)

    stm_grid.addWidget(QtWidgets.QLabel("LUMO offsets:"), 4, 0)
    window.afm_stm_lumo_offsets = QtWidgets.QLineEdit("1 2 3")
    stm_grid.addWidget(window.afm_stm_lumo_offsets, 4, 1)
    stm_grid.addWidget(QtWidgets.QLabel("MO indices:"), 5, 0)
    window.afm_stm_mo_indices = QtWidgets.QLineEdit("")
    stm_grid.addWidget(window.afm_stm_mo_indices, 5, 1)

    stm_grid.addWidget(QtWidgets.QLabel("field:"), 6, 0)
    window.afm_stm_field_combo = QtWidgets.QComboBox()
    window.afm_stm_field_combo.addItems(['ldos', 'psi2', 'psi'])
    window.afm_stm_field_combo.setCurrentText('ldos')
    stm_grid.addWidget(window.afm_stm_field_combo, 6, 1)

    stm_grid.addWidget(QtWidgets.QLabel("exp_beta:"), 7, 0)
    window.afm_stm_exp_beta = QtWidgets.QDoubleSpinBox()
    window.afm_stm_exp_beta.setRange(0.1, 10.0)
    window.afm_stm_exp_beta.setValue(1.0)
    window.afm_stm_exp_beta.setDecimals(3)
    stm_grid.addWidget(window.afm_stm_exp_beta, 7, 1)
    stm_grid.addWidget(QtWidgets.QLabel("exp_r0:"), 8, 0)
    window.afm_stm_exp_r0 = QtWidgets.QDoubleSpinBox()
    window.afm_stm_exp_r0.setRange(0.0, 10.0)
    window.afm_stm_exp_r0.setValue(3.0)
    window.afm_stm_exp_r0.setDecimals(3)
    stm_grid.addWidget(window.afm_stm_exp_r0, 8, 1)
    window.afm_stm_bond_resolved = QtWidgets.QCheckBox("Bond-resolved")
    window.afm_stm_bond_resolved.setChecked(False)
    stm_grid.addWidget(window.afm_stm_bond_resolved, 9, 0, 1, 2)

    window.afm_stm_info_label = QtWidgets.QLabel("")
    window.afm_stm_info_label.setWordWrap(True)
    stm_grid.addWidget(window.afm_stm_info_label, 10, 0, 1, 2)

    def _stm_ui_debug_print(*args):
        try:
            field = str(window.afm_stm_field_combo.currentText()) if hasattr(window, 'afm_stm_field_combo') else 'ldos'
            if window.afm_stm_mo_indices.text().strip() != '':
                print(f"[STM UI] field={field}  mode=mo_indices  mo_indices='{window.afm_stm_mo_indices.text()}'")
            elif window.afm_stm_use_homo_range.isChecked():
                print(f"[STM UI] field={field}  mode=homo_range  off0={int(window.afm_stm_homo_off0.value())}  n={int(window.afm_stm_homo_n.value())}")
            else:
                print(f"[STM UI] field={field}  mode=lumo_offsets  lumo_offsets='{window.afm_stm_lumo_offsets.text()}'")
        except Exception:
            pass

    window.afm_stm_enable.stateChanged.connect(_stm_ui_debug_print)
    window.afm_stm_use_homo_range.stateChanged.connect(_stm_ui_debug_print)
    window.afm_stm_homo_off0.valueChanged.connect(_stm_ui_debug_print)
    window.afm_stm_homo_n.valueChanged.connect(_stm_ui_debug_print)
    window.afm_stm_lumo_offsets.textChanged.connect(_stm_ui_debug_print)
    window.afm_stm_mo_indices.textChanged.connect(_stm_ui_debug_print)
    window.afm_stm_field_combo.currentIndexChanged.connect(_stm_ui_debug_print)
    window.afm_stm_bond_resolved.stateChanged.connect(_stm_ui_debug_print)

    stm_sec.setContent(stm_widget)
    layout.addWidget(stm_sec)

    # AFM state variables
    window._afm_density = None
    window._afm_potentials = None
    window._afm_results = {}

    view_modes = [
        ('AFM Total', lambda: window.set_view_mode('afm_total')),
        ('AFM Pauli', lambda: window.set_view_mode('afm_pauli')),
        ('AFM ES', lambda: window.set_view_mode('afm_es')),
        ('AFM vdW', lambda: window.set_view_mode('afm_vdw')),
    ]
    return UIComponents(panel=panel, view_modes=view_modes)
