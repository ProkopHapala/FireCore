"""
AFMExtension.py - AFM simulation extension for KekuleExplorerGUI
===============================================================
Uses ModularAFMPipeline for staged, cached computation.
Dirty flag system ensures only changed stages are recomputed.
"""

import os
import numpy as np
from PyQt5 import QtWidgets, QtCore
from pyBall.GUI.CollapsibleSection import CollapsibleSection
from pyBall.ExtensionManager import UIComponents


# ============================================================
# Dirty flag management
# ============================================================

class AFMDirtyFlags:
    """
    Tracks which pipeline stages are stale.
    Setting a stage dirty automatically marks all downstream stages dirty too.
    Dependency chain: geometry/basis/step -> S1 -> S2 -> S3 -> S4 -> S5/S6
    """
    _STAGES = ['s1', 's2', 's3', 's4', 's5', 's6']
    _DOWNSTREAM = {'s1': ['s2'], 's2': ['s3'], 's3': ['s4'], 's4': ['s5', 's6'], 's5': [], 's6': []}

    def __init__(self):
        self._flags = {s: True for s in self._STAGES}  # Initially all dirty

    def mark(self, stage):
        """Mark a stage and all its downstream stages dirty."""
        if stage not in self._flags:
            raise KeyError(f"Unknown stage '{stage}'. Valid: {self._STAGES}")
        self._flags[stage] = True
        for ds in self._DOWNSTREAM.get(stage, []):
            self.mark(ds)

    def clean(self, stage):
        self._flags[stage] = False

    def is_dirty(self, stage):
        return self._flags.get(stage, True)

    def mark_geometry_changed(self):
        """Geometry change invalidates entire pipeline."""
        self.mark('s1')

    def mark_density_params_changed(self):
        """Step/margin change invalidates grid projection onwards."""
        self.mark('s2')

    def mark_physics_params_changed(self):
        """Pauli/vdW params change invalidates potentials onwards."""
        self.mark('s3')

    def mark_scan_params_changed(self):
        """Scan range/heights change invalidates relaxation onwards."""
        self.mark('s4')

    def mark_stm_params_changed(self):
        """MO selection change invalidates STM/BR-STM only."""
        self.mark('s5')
        self.mark('s6')

    def status_str(self):
        return "  ".join(f"S{i+1}:{'D' if self._flags[s] else 'C'}" for i, s in enumerate(self._STAGES))


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


def _get_pipeline_params(window):
    """Snapshot current UI parameter values - used for dirty detection."""
    return {
        'basis':      window.afm_basis_combo.currentText(),
        'step':       window.afm_step_spin.value(),
        'margin':     window.afm_margin_spin.value(),
        'z_extra':    6.0,
        'scan_range': window.afm_scan_range_spin.value(),
        'hmin':       window.afm_hmin_spin.value(),
        'hmax':       window.afm_hmax_spin.value(),
        'hstep':      window.afm_hstep_spin.value(),
        'pauli_a':    window.afm_pauli_a_spin.value(),
        'pauli_beta': window.afm_pauli_beta_spin.value(),
        'c6':         window.afm_vdw_c6_spin.value(),
        'klat':       window.afm_klat_spin.value(),
    }


def _ensure_pipeline(window):
    """
    Create or return existing ModularAFMPipeline. 
    Re-creates only if geometry or basis/step parameters changed (S1 dirty).
    Returns the pipeline instance.
    """
    from pyBall.OCL.ModularPipeline import ModularAFMPipeline
    from pyBall import dftb_utils as du
    import tempfile

    atomPos, _, enames = _get_afm_geometry(window)
    params = _get_pipeline_params(window)

    # Check if we need to create/recreate the pipeline
    needs_reinit = (window._afm_pipeline is None)
    if not needs_reinit:
        # Reinit if geometry identity changed (atom count, centroid) or key params
        prev = window._afm_pipeline_params
        reinit_keys = {'basis', 'step', 'margin', 'z_extra', 'scan_range', 'hmin', 'hmax', 'hstep'}
        needs_reinit = any(params[k] != prev.get(k) for k in reinit_keys)
        if not needs_reinit:
            # Check geometry (by atom count + centroid hash)
            prev_geom = window._afm_pipeline_geom_hash
            cur_hash = (len(atomPos), round(float(atomPos[:,0].mean()), 4), round(float(atomPos[:,1].mean()), 4))
            needs_reinit = (prev_geom != cur_hash)

    if needs_reinit:
        _update_afm_status(window, "Initializing modular pipeline...")
        basis = params['basis']
        slako_prefix = du.SK_PATHS.get(basis, basis)
        output_dir = getattr(window, '_afm_output_dir', None)
        if output_dir is None:
            output_dir = tempfile.mkdtemp(prefix='afm_gui_')
            window._afm_output_dir = output_dir

        window._afm_pipeline = ModularAFMPipeline(
            xyz_file=None,  # unused when atomPos/enames injected
            output_dir=output_dir,
            basis=basis, slako_prefix=slako_prefix,
            step=params['step'], margin=params['margin'], z_extra=params['z_extra'],
            scan_range=params['scan_range'], scan_step=0.1,
            height_range=(params['hmin'], params['hmax']), height_step=params['hstep'],
            atomPos=atomPos, enames=enames,
        )
        window._afm_pipeline_params = params.copy()
        window._afm_pipeline_geom_hash = (len(atomPos), round(float(atomPos[:,0].mean()), 4), round(float(atomPos[:,1].mean()), 4))
        window._afm_dirty.mark_geometry_changed()  # Full cascade

    return window._afm_pipeline


def _get_homo_index(window):
    """Return the HOMO index (0-based) from cached eigvals, or None."""
    if window._afm_eigvals is not None:
        homo = int(np.sum(window._afm_eigvals < 0.0)) - 1
        return max(0, homo)
    return None


def _update_homo_label(window):
    """Update the HOMO info label in the STM section after SCF completes."""
    if not hasattr(window, 'afm_homo_label'):
        return
    homo = _get_homo_index(window)
    if homo is None:
        window.afm_homo_label.setText("(run SCF first)")
        return
    nmo = len(window._afm_eigvals)
    lumo = homo + 1
    e_homo = window._afm_eigvals[homo]
    e_lumo = window._afm_eigvals[lumo] if lumo < nmo else float('nan')
    window.afm_homo_label.setText(f"#{homo}  (E={e_homo:.3f} eV,  LUMO #{lumo} E={e_lumo:.3f} eV)")
    # Also bump orbital spin max to nmo-1
    if hasattr(window, 'afm_orbital_spin'):
        window.afm_orbital_spin.setRange(0, nmo - 1)
        window.afm_orbital_spin.setValue(homo)


def _get_stm_params_from_ui(window):
    """Read current STM parameter values from UI widgets.
    MO list is a space/comma-separated list of integers.
    If relative_mo checkbox is checked they are relative to HOMO (0=HOMO, +1=LUMO, -1=HOMO-1).
    Returns absolute mo_indices list (always) and lumo_offsets=None (deprecated).
    """
    raw = window.afm_stm_mo_list.text().replace(',', ' ').split() if hasattr(window, 'afm_stm_mo_list') else []
    try:
        offsets = [int(s) for s in raw if s.strip()]
    except Exception:
        offsets = [0]
    if not offsets:
        offsets = [0]

    relative = hasattr(window, 'afm_stm_relative_mo') and window.afm_stm_relative_mo.isChecked()
    if relative:
        homo = _get_homo_index(window)
        if homo is None:
            raise ValueError("HOMO not determined yet — run Stage 1 (SCF) first.")
        mo_indices = [homo + d for d in offsets]
    else:
        mo_indices = offsets

    return {
        'lumo_offsets': None,
        'mo_indices':   mo_indices,
        'field':        str(window.afm_stm_field_combo.currentText()) if hasattr(window, 'afm_stm_field_combo') else 'ldos',
        'exp_beta':     float(window.afm_stm_exp_beta.value()),
        'exp_r0':       float(window.afm_stm_exp_r0.value()),
        'bond_resolved': bool(window.afm_stm_bond_resolved.isChecked()),
    }


def _ensure_stages_for_component(window, component):
    """Auto-run whichever pipeline stages are needed for `component`, if dirty.
    Returns True if all needed data is ready, raises on failure.
    """
    dirty = window._afm_dirty
    params = _get_pipeline_params(window)

    def _need_s1_to_s4():
        """Run stages 1-4 if any of them are dirty or results missing."""
        pipe = _ensure_pipeline(window)
        if dirty.is_dirty('s1') or window._afm_eigvecs is None:
            _update_afm_status(window, "Auto: Stage 1 SCF...")
            dm_dense, eigvecs, eigvals = pipe.stage1_scf(force_recompute=dirty.is_dirty('s1'))
            dirty.clean('s1')
            window._afm_eigvecs = eigvecs; window._afm_eigvals = eigvals
            _update_homo_label(window)
        else:
            dm_dense = None  # Will be loaded from cache if s2 needs it

        if dirty.is_dirty('s2') or window._afm_density is None:
            if dm_dense is None:
                data = np.load(pipe.cache_stage1, allow_pickle=True)
                dm_dense = data['dm_dense']
            _update_afm_status(window, "Auto: Stage 2 grid...")
            rho_scf, rho_na, rho_diff = pipe.stage2_project(dm_dense, force_recompute=dirty.is_dirty('s2'))
            dirty.clean('s2')
            window._afm_density = {'rho_scf': rho_scf, 'rho_na': rho_na, 'rho_diff': rho_diff,
                                    'origin': pipe.origin, 'ngrid': pipe.ngrid, 'grid_spec': pipe.grid_spec}
        if dirty.is_dirty('s3') or window._afm_potentials is None:
            d = window._afm_density
            _update_afm_status(window, "Auto: Stage 3 potentials...")
            V_ES, E_pauli, E_ES, E_vdw, F_total = pipe.stage3_potentials(
                d['rho_scf'], d['rho_na'], d['rho_diff'], force_recompute=dirty.is_dirty('s3'),
                pauli_params={'A': params['pauli_a'], 'beta': params['pauli_beta']},
                vdw_params={'C6_CO': params['c6']})
            dirty.clean('s3')
            window._afm_potentials = {'V_ES': V_ES, 'E_pauli_field': E_pauli, 'E_ES_field': E_ES,
                                       'E_vdw': E_vdw, 'F_total': F_total,
                                       'origin': pipe.origin, 'step': pipe.step, 'grid_spec': pipe.grid_spec}
        if dirty.is_dirty('s4') or window._afm_results is None or 'df' not in (window._afm_results or {}):
            _update_afm_status(window, "Auto: Stage 4 relax...")
            df, tip_disp, FEs_relax = pipe.stage4_relax(
                window._afm_potentials['F_total'], force_recompute=dirty.is_dirty('s4'),
                relax_params={'K_LAT': params['klat']}, ppm_mode=True)
            dirty.clean('s4')
            window._afm_results = {'df': df, 'tip_disp': tip_disp, 'FEs_relax': FEs_relax,
                                    'heights': pipe.heights, 'scan_xs': pipe.scan_xs, 'scan_ys': pipe.scan_ys}
            mid_z = float(pipe.heights[len(pipe.heights)//2])
            window.afm_z_height_spin.setValue(mid_z)

    if component == "AFM Image (df)":
        _need_s1_to_s4()

    elif component in ("STM Signal", "BR-STM Signal"):
        pipe = _ensure_pipeline(window)
        # S1-S4 must be done first for S5/S6
        _need_s1_to_s4()
        sp = _get_stm_params_from_ui(window)
        need_br = (component == "BR-STM Signal") or sp['bond_resolved']
        if dirty.is_dirty('s5') or 'stm_grid' not in (window._afm_results or {}):
            _update_afm_status(window, f"Auto: Stage 5 STM (MOs={sp['mo_indices']})...")
            stm_grid = pipe.stage5_stm(window._afm_eigvecs, window._afm_eigvals,
                lumo_offsets=sp['lumo_offsets'], mo_indices=sp['mo_indices'],
                field=sp['field'], exp_beta=sp['exp_beta'], exp_r0=sp['exp_r0'])
            window._afm_results['stm_grid'] = stm_grid
            dirty.clean('s5')
        if need_br and (dirty.is_dirty('s6') or 'br_stm_grid' not in (window._afm_results or {})):
            _update_afm_status(window, "Auto: Stage 6 BR-STM...")
            br_stm_grid = pipe.stage6_br_stm(window._afm_eigvecs, window._afm_eigvals, window._afm_results['tip_disp'],
                lumo_offsets=sp['lumo_offsets'], mo_indices=sp['mo_indices'],
                field=sp['field'], exp_beta=sp['exp_beta'], exp_r0=sp['exp_r0'])
            window._afm_results['br_stm_grid'] = br_stm_grid
            dirty.clean('s6')
        if hasattr(window, '_afm_refresh_dirty_label'):
            window._afm_refresh_dirty_label()

    elif component in ("SCF Density", "Neutral Density", "Delta Density"):
        pipe = _ensure_pipeline(window)
        if dirty.is_dirty('s1') or window._afm_eigvecs is None:
            _update_afm_status(window, "Auto: Stage 1 SCF...")
            dm_dense, ev, ev2 = pipe.stage1_scf(force_recompute=dirty.is_dirty('s1'))
            dirty.clean('s1'); window._afm_eigvecs = ev; window._afm_eigvals = ev2
        else:
            dm_dense = np.load(pipe.cache_stage1, allow_pickle=True)['dm_dense']
        if dirty.is_dirty('s2') or window._afm_density is None:
            _update_afm_status(window, "Auto: Stage 2 grid...")
            rho_scf, rho_na, rho_diff = pipe.stage2_project(dm_dense, force_recompute=dirty.is_dirty('s2'))
            dirty.clean('s2')
            window._afm_density = {'rho_scf': rho_scf, 'rho_na': rho_na, 'rho_diff': rho_diff,
                                    'origin': pipe.origin, 'ngrid': pipe.ngrid, 'grid_spec': pipe.grid_spec}

    elif component in ("Pauli Energy", "Electrostatic Energy", "vdW Energy", "Total Potential", "Total Z-Force"):
        _need_s1_to_s4()

    if hasattr(window, '_afm_refresh_dirty_label'):
        window._afm_refresh_dirty_label()


def run_afm_full_pipeline(window):
    """Run complete AFM pipeline (stages 1-4), then STM/BR-STM if enabled.
    Only recomputes stages that are marked dirty."""
    try:
        atomPos, _, enames = _get_afm_geometry(window)
        if len(atomPos) == 0:
            raise ValueError("No atoms in molecule")

        pipe = _ensure_pipeline(window)
        dirty = window._afm_dirty
        params = _get_pipeline_params(window)

        _update_afm_status(window, f"Pipeline [{dirty.status_str()}] - running dirty stages...")

        # --- Stage 1: SCF ---
        _update_afm_status(window, "Stage 1: DFTB+ SCF...")
        dm_dense, eigvecs, eigvals = pipe.stage1_scf(force_recompute=dirty.is_dirty('s1'))
        dirty.clean('s1')
        window._afm_eigvecs = eigvecs
        window._afm_eigvals = eigvals
        _update_homo_label(window)

        # --- Stage 2: Grid projection ---
        _update_afm_status(window, "Stage 2: Grid projection...")
        rho_scf, rho_na, rho_diff = pipe.stage2_project(dm_dense, force_recompute=dirty.is_dirty('s2'))
        dirty.clean('s2')
        window._afm_density = {
            'rho_scf': rho_scf, 'rho_na': rho_na, 'rho_diff': rho_diff,
            'origin': pipe.origin, 'ngrid': pipe.ngrid, 'grid_spec': pipe.grid_spec,
        }

        # --- Stage 3: Potentials ---
        _update_afm_status(window, "Stage 3: FDBM potentials...")
        V_ES, E_pauli_field, E_ES_field, E_vdw, F_total = pipe.stage3_potentials(
            rho_scf, rho_na, rho_diff, force_recompute=dirty.is_dirty('s3'),
            pauli_params={'A': params['pauli_a'], 'beta': params['pauli_beta']},
            vdw_params={'C6_CO': params['c6']},
        )
        dirty.clean('s3')
        window._afm_potentials = {
            'V_ES': V_ES, 'E_pauli_field': E_pauli_field, 'E_ES_field': E_ES_field,
            'E_vdw': E_vdw, 'F_total': F_total,
            'origin': pipe.origin, 'step': pipe.step, 'grid_spec': pipe.grid_spec,
        }

        # --- Stage 4: Relaxation ---
        _update_afm_status(window, "Stage 4: PP relaxation...")
        df, tip_disp, FEs_relax = pipe.stage4_relax(
            F_total, force_recompute=dirty.is_dirty('s4'),
            relax_params={'K_LAT': params['klat']}, ppm_mode=True,
        )
        dirty.clean('s4')
        window._afm_results = {
            'df': df, 'tip_disp': tip_disp, 'FEs_relax': FEs_relax,
            'heights': pipe.heights, 'scan_xs': pipe.scan_xs, 'scan_ys': pipe.scan_ys,
        }

        # --- Stage 5/6: STM / BR-STM (always fast, rerun on request) ---
        stm_enable = hasattr(window, 'afm_stm_enable') and window.afm_stm_enable.isChecked()
        if stm_enable:
            sp = _get_stm_params_from_ui(window)
            _update_afm_status(window, f"Stage 5: STM (field={sp['field']}, lumo_offsets={sp['lumo_offsets']})...")
            stm_grid = pipe.stage5_stm(eigvecs, eigvals,
                lumo_offsets=sp['lumo_offsets'], mo_indices=sp['mo_indices'],
                field=sp['field'], exp_beta=sp['exp_beta'], exp_r0=sp['exp_r0'])
            window._afm_results['stm_grid'] = stm_grid
            dirty.clean('s5')

            if sp['bond_resolved']:
                _update_afm_status(window, "Stage 6: BR-STM...")
                br_stm_grid = pipe.stage6_br_stm(eigvecs, eigvals, tip_disp,
                    lumo_offsets=sp['lumo_offsets'], mo_indices=sp['mo_indices'],
                    field=sp['field'], exp_beta=sp['exp_beta'], exp_r0=sp['exp_r0'])
                window._afm_results['br_stm_grid'] = br_stm_grid
                dirty.clean('s6')

        mid_z = float(pipe.heights[len(pipe.heights)//2])
        window.afm_z_height_spin.setValue(mid_z)

        nz = df.shape[2]
        msg = f"Done [{dirty.status_str()}]  df=[{df.min():.2f},{df.max():.2f}]Hz  {nz} heights"
        _update_afm_status(window, msg)

    except Exception as e:
        _update_afm_status(window, f"FAILED: {e}")
        raise


def run_afm_stage1(window):
    """Run Stage 1 (DFTB+ SCF) only."""
    try:
        atomPos, _, _ = _get_afm_geometry(window)
        if len(atomPos) == 0:
            raise ValueError("No atoms in molecule")
        pipe = _ensure_pipeline(window)
        _update_afm_status(window, "Stage 1: DFTB+ SCF (forced)...")
        dm_dense, eigvecs, eigvals = pipe.stage1_scf(force_recompute=True)
        window._afm_dirty.mark('s1')   # Downstream still dirty
        window._afm_dirty.clean('s1')  # S1 itself is now clean
        window._afm_eigvecs = eigvecs
        window._afm_eigvals = eigvals
        _update_homo_label(window)
        _update_afm_status(window, f"Stage 1 done. [{window._afm_dirty.status_str()}]")
    except Exception as e:
        _update_afm_status(window, f"FAILED: {e}")
        raise


def run_afm_stage2(window):
    """Run Stage 2 (grid projection) only."""
    try:
        pipe = _ensure_pipeline(window)
        if not os.path.exists(pipe.cache_stage1):
            raise ValueError("Stage 1 not computed. Run Stage 1 first.")
        _update_afm_status(window, "Stage 2: Grid projection (forced)...")
        data = np.load(pipe.cache_stage1, allow_pickle=True)
        dm_dense = data['dm_dense']
        rho_scf, rho_na, rho_diff = pipe.stage2_project(dm_dense, force_recompute=True)
        window._afm_dirty.mark('s2')
        window._afm_dirty.clean('s2')
        window._afm_density = {
            'rho_scf': rho_scf, 'rho_na': rho_na, 'rho_diff': rho_diff,
            'origin': pipe.origin, 'ngrid': pipe.ngrid, 'grid_spec': pipe.grid_spec,
        }
        _update_afm_status(window, f"Stage 2 done. [{window._afm_dirty.status_str()}]")
    except Exception as e:
        _update_afm_status(window, f"FAILED: {e}")
        raise


def run_afm_stage3(window):
    """Run Stage 3 (potentials) only."""
    try:
        pipe = _ensure_pipeline(window)
        if not os.path.exists(pipe.cache_stage2):
            raise ValueError("Stage 2 not computed. Run Stage 2 first.")
        _update_afm_status(window, "Stage 3: Potentials (forced)...")
        data = np.load(pipe.cache_stage2)
        params = _get_pipeline_params(window)
        V_ES, E_pauli_field, E_ES_field, E_vdw, F_total = pipe.stage3_potentials(
            data['rho_scf'], data['rho_na'], data['rho_diff'], force_recompute=True,
            pauli_params={'A': params['pauli_a'], 'beta': params['pauli_beta']},
            vdw_params={'C6_CO': params['c6']},
        )
        window._afm_dirty.mark('s3')
        window._afm_dirty.clean('s3')
        window._afm_potentials = {
            'V_ES': V_ES, 'E_pauli_field': E_pauli_field, 'E_ES_field': E_ES_field,
            'E_vdw': E_vdw, 'F_total': F_total,
            'origin': pipe.origin, 'step': pipe.step, 'grid_spec': pipe.grid_spec,
        }
        _update_afm_status(window, f"Stage 3 done. [{window._afm_dirty.status_str()}]")
    except Exception as e:
        _update_afm_status(window, f"FAILED: {e}")
        raise


def run_afm_stage4(window):
    """Run Stage 4 (relaxation) only."""
    try:
        pipe = _ensure_pipeline(window)
        if window._afm_potentials is None or 'F_total' not in window._afm_potentials:
            raise ValueError("Stage 3 not computed. Run Stage 3 first.")
        params = _get_pipeline_params(window)
        _update_afm_status(window, "Stage 4: PP relaxation (forced)...")
        df, tip_disp, FEs_relax = pipe.stage4_relax(
            window._afm_potentials['F_total'], force_recompute=True,
            relax_params={'K_LAT': params['klat']}, ppm_mode=True,
        )
        window._afm_dirty.mark('s4')
        window._afm_dirty.clean('s4')
        window._afm_results = {
            'df': df, 'tip_disp': tip_disp, 'FEs_relax': FEs_relax,
            'heights': pipe.heights, 'scan_xs': pipe.scan_xs, 'scan_ys': pipe.scan_ys,
        }
        mid_z = float(pipe.heights[len(pipe.heights)//2])
        window.afm_z_height_spin.setValue(mid_z)
        _update_afm_status(window, f"Stage 4 done. [{window._afm_dirty.status_str()}]")
    except Exception as e:
        _update_afm_status(window, f"FAILED: {e}")
        raise


def run_stm(window):
    """Run STM/BR-STM independently (stages 5/6 only - very fast from cache)."""
    try:
        pipe = _ensure_pipeline(window)
        if window._afm_eigvecs is None:
            raise ValueError("Eigenvectors not available. Run at least Stage 1 first.")
        if window._afm_results is None or 'tip_disp' not in window._afm_results:
            raise ValueError("Tip displacements not available. Run Stage 4 first.")

        sp = _get_stm_params_from_ui(window)
        eigvecs = window._afm_eigvecs
        eigvals = window._afm_eigvals
        tip_disp = window._afm_results['tip_disp']

        _update_afm_status(window, f"STM (field={sp['field']}, lumo_offsets={sp['lumo_offsets']})...")
        stm_grid = pipe.stage5_stm(eigvecs, eigvals,
            lumo_offsets=sp['lumo_offsets'], mo_indices=sp['mo_indices'],
            field=sp['field'], exp_beta=sp['exp_beta'], exp_r0=sp['exp_r0'])
        window._afm_results['stm_grid'] = stm_grid
        window._afm_dirty.clean('s5')

        if sp['bond_resolved']:
            _update_afm_status(window, "BR-STM...")
            br_stm_grid = pipe.stage6_br_stm(eigvecs, eigvals, tip_disp,
                lumo_offsets=sp['lumo_offsets'], mo_indices=sp['mo_indices'],
                field=sp['field'], exp_beta=sp['exp_beta'], exp_r0=sp['exp_r0'])
            window._afm_results['br_stm_grid'] = br_stm_grid
            window._afm_dirty.clean('s6')

        _update_afm_status(window, f"STM done. [{window._afm_dirty.status_str()}]")
    except Exception as e:
        _update_afm_status(window, f"FAILED: {e}")
        raise


def plot_orbital_map(window):
    """Plot a single molecular orbital with phase (signed psi, not psi^2)."""
    try:
        import matplotlib
        matplotlib.use('Qt5Agg')
        from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
        from matplotlib.figure import Figure

        pipe = _ensure_pipeline(window)
        if window._afm_eigvecs is None:
            raise ValueError("Eigenvectors not available. Run Stage 1 first.")

        mo_idx  = int(window.afm_orbital_spin.value())
        z_height = window.afm_z_height_spin.value()
        eigvecs = window._afm_eigvecs
        eigvals = window._afm_eigvals

        nmo = eigvecs.shape[0]
        if not (0 <= mo_idx < nmo):
            raise ValueError(f"MO index {mo_idx} out of range [0, {nmo-1}]")

        # Find HOMO for labelling
        nocc = nmo // 2  # approximation; real HOMO from eigvals
        homo = int(np.sum(eigvals < 0.0)) - 1 if np.any(eigvals < 0.0) else nocc - 1

        # Sample on a 2D grid at z_height
        origin   = pipe.origin
        ngrid    = pipe.ngrid
        step     = pipe.step
        xs = np.linspace(pipe.scan_xs[0], pipe.scan_xs[-1], len(pipe.scan_xs))
        ys = np.linspace(pipe.scan_ys[0], pipe.scan_ys[-1], len(pipe.scan_ys))
        XX, YY = np.meshgrid(xs, ys, indexing='ij')
        points = np.stack([XX.ravel(), YY.ravel(), np.full(XX.size, z_height)], axis=1).astype(np.float32)

        coeffs = eigvecs[mo_idx].astype(np.float32)
        exp_beta = float(window.afm_stm_exp_beta.value()) if hasattr(window, 'afm_stm_exp_beta') else 1.0
        exp_r0   = float(window.afm_stm_exp_r0.value())  if hasattr(window, 'afm_stm_exp_r0')   else 3.0

        psi = pipe.projector.project_orbital_dense_points_exp(
            points, coeffs, pipe.norb_per_atom, pipe.orb_offsets, pipe.atoms_dict,
            beta=exp_beta, r0=exp_r0
        )
        psi_2d = psi.reshape(len(xs), len(ys))

        # Plot signed wavefunction with seismic colormap (blue=neg, red=pos)
        fig = Figure(figsize=(7, 5.5), dpi=100)
        ax  = fig.add_subplot(111)
        vmax = np.abs(psi_2d).max() or 1.0
        im = ax.imshow(psi_2d.T, origin='lower', cmap='seismic', vmin=-vmax, vmax=vmax,
                       extent=[xs[0], xs[-1], ys[0], ys[-1]], aspect='equal')
        fig.colorbar(im, ax=ax, label='psi (a.u.)')
        rel = mo_idx - homo
        label = f"HOMO{rel:+d}" if rel != 0 else "HOMO"
        ax.set_title(f"MO #{mo_idx} ({label})  E={eigvals[mo_idx]:.3f} eV\nz={z_height:.2f} A")
        ax.set_xlabel('x (A)'); ax.set_ylabel('y (A)')
        _overlay_atoms(ax, window, xs, ys)

        _show_in_plot_window(window, fig, f"Orbital #{mo_idx} ({label})")
        window.statusBar().showMessage(f"Orbital #{mo_idx} ({label}) at z={z_height:.2f} A")
    except Exception as e:
        raise RuntimeError(f"Orbital plot FAILED: {e}")


def _get_z_slice(grid_spec, step, z_height):
    """Convert physical z-height to grid index."""
    origin_z = grid_spec['origin'][2]
    iz = int(np.clip(np.round((z_height - origin_z) / step), 0, grid_spec['ngrid'][2] - 1))
    actual_z = origin_z + iz * step
    return iz, actual_z


def _show_in_plot_window(window, fig, title="AFM Plot"):
    """Show a matplotlib Figure in a reusable Qt dialog window."""
    from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
    if not hasattr(window, '_afm_plot_window') or window._afm_plot_window is None:
        window._afm_plot_window = QtWidgets.QDialog(window)
        window._afm_plot_window.setWindowTitle(title)
        window._afm_plot_layout = QtWidgets.QVBoxLayout(window._afm_plot_window)
        window._afm_plot_window.resize(700, 600)
        def on_closed():
            window._afm_plot_window = None
        window._afm_plot_window.finished.connect(on_closed)
        window._afm_plot_window.show()
    else:
        while window._afm_plot_layout.count():
            item = window._afm_plot_layout.takeAt(0)
            if item.widget():
                item.widget().deleteLater()
    canvas = FigureCanvas(fig)
    if hasattr(window, 'install_mpl_canvas_screenshot_menu'):
        try:
            window.install_mpl_canvas_screenshot_menu(canvas, fig, default_name=f"{title.replace(' ','_')}.png")
        except Exception:
            pass
    window._afm_plot_layout.addWidget(canvas)
    window._afm_plot_window.setWindowTitle(title)
    window._afm_plot_window.show()
    window._afm_plot_window.raise_()
    window._afm_plot_window.activateWindow()


def _overlay_atoms(ax, window, xs, ys):
    """Overlay atom positions as small dots if checkbox is enabled."""
    if not (hasattr(window, 'afm_show_atoms') and window.afm_show_atoms.isChecked()):
        return
    if not hasattr(window, 'backend') or window.backend.sys is None:
        return
    apos = window.backend.sys.apos
    enames = window.backend.sys.enames
    ELEM_COLOR = {'H': 'white', 'C': 'gray', 'N': 'blue', 'O': 'red', 'S': 'yellow'}
    for i, (pos, e) in enumerate(zip(apos, enames)):
        if xs[0] <= pos[0] <= xs[-1] and ys[0] <= pos[1] <= ys[-1]:
            c = ELEM_COLOR.get(e, 'magenta')
            ax.plot(pos[0], pos[1], '.', color=c, markersize=4, markeredgecolor='k', markeredgewidth=0.3)


def plot_afm_slice(window):
    """Plot single z-slice in a GUI window (not to disk). Auto-runs needed stages."""
    try:
        import matplotlib
        matplotlib.use('Qt5Agg')
        from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
        from matplotlib.figure import Figure

        component = window.afm_component_combo.currentText()
        z_height = window.afm_z_height_spin.value()
        auto_limits = window.afm_auto_limits.isChecked()

        # Auto-run any needed pipeline stages
        _ensure_stages_for_component(window, component)

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

        elif component in ("STM Signal", "BR-STM Signal"):
            grid_key = 'stm_grid' if component == "STM Signal" else 'br_stm_grid'
            if window._afm_results is None or grid_key not in window._afm_results:
                raise ValueError(f"No {component} results. Enable STM{'+ Bond-resolved' if component=='BR-STM Signal' else ''} and run pipeline.")
            data_3d = window._afm_results[grid_key]
            heights = window._afm_results.get('heights', [])
            if len(heights) == 0:
                raise ValueError("No heights in AFM results")
            h_idx = np.argmin(np.abs(heights - z_height))
            actual_z = heights[h_idx]
            iz = h_idx
            step = heights[1] - heights[0] if len(heights) > 1 else 0.1
            cmap = 'viridis'
            symmetric = False
            data_label = f"{component} (arb.)"
            data = data_3d[:, :, iz]
            
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
            # Get slice at requested z-height (may be zero if in vacuum)
            iz, actual_z = _get_z_slice(grid_spec, step, z_height)
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

        # Build physical extent for atom overlay alignment
        nx, ny = data.shape
        if window._afm_results is not None and 'scan_xs' in window._afm_results:
            xs = window._afm_results['scan_xs']
            ys = window._afm_results['scan_ys']
            extent = [xs[0], xs[-1], ys[0], ys[-1]]
        elif window._afm_potentials is not None:
            gs = window._afm_potentials['grid_spec']
            xs = np.linspace(gs['origin'][0], gs['origin'][0] + nx * window._afm_potentials['step'], nx)
            ys = np.linspace(gs['origin'][1], gs['origin'][1] + ny * window._afm_potentials['step'], ny)
            extent = [xs[0], xs[-1], ys[0], ys[-1]]
        elif window._afm_density is not None:
            gs = window._afm_density['grid_spec']
            st = float(gs['dA'][0])
            xs = np.linspace(gs['origin'][0], gs['origin'][0] + nx * st, nx)
            ys = np.linspace(gs['origin'][1], gs['origin'][1] + ny * st, ny)
            extent = [xs[0], xs[-1], ys[0], ys[-1]]
        else:
            xs = np.arange(nx); ys = np.arange(ny)
            extent = None

        im = ax.imshow(data.T, origin='lower', cmap=cmap, vmin=vmin, vmax=vmax,
                       extent=extent, aspect='equal')
        ax.set_title(f"{component}\nZ={actual_z:.2f}A (iz={iz}) | [{data_min:.3f}, {data_max:.3f}]", fontsize=10)
        ax.set_xlabel('x (A)'); ax.set_ylabel('y (A)')
        fig.colorbar(im, ax=ax, label=data_label)
        _overlay_atoms(ax, window, xs, ys)

        _show_in_plot_window(window, fig, f"AFM Slice - {component} Z={z_height:.2f}A")
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

    # --- State variables ---
    window._afm_density   = None
    window._afm_potentials = None
    window._afm_results   = None
    window._afm_eigvecs   = None
    window._afm_eigvals   = None
    window._afm_pipeline  = None
    window._afm_pipeline_params   = {}
    window._afm_pipeline_geom_hash = None
    window._afm_output_dir = None
    window._afm_dirty = AFMDirtyFlags()

    # --- Dirty-flag status label (always visible) ---
    window.afm_dirty_label = QtWidgets.QLabel("Cache: [all dirty]")
    window.afm_dirty_label.setStyleSheet("color: #cc6600; font-size: 10px;")
    layout.addWidget(window.afm_dirty_label)

    def _refresh_dirty_label():
        window.afm_dirty_label.setText(f"Cache: [{window._afm_dirty.status_str()}]")
    window._afm_refresh_dirty_label = _refresh_dirty_label

    # --- Main pipeline button ---
    full_btn = QtWidgets.QPushButton("Run Full AFM Pipeline (smart)")
    full_btn.setToolTip("Runs only dirty stages. Change geometry or params to force recompute.")
    full_btn.clicked.connect(lambda: (run_afm_full_pipeline(window), _refresh_dirty_label()))
    layout.addWidget(full_btn)

    # --- Individual stage buttons ---
    stage_layout = QtWidgets.QHBoxLayout()
    s1_btn = QtWidgets.QPushButton("S1: SCF")
    s1_btn.setToolTip("DFTB+ SCF - density matrix and eigenvectors")
    s1_btn.clicked.connect(lambda: (run_afm_stage1(window), _refresh_dirty_label()))
    stage_layout.addWidget(s1_btn)

    s2_btn = QtWidgets.QPushButton("S2: Grid")
    s2_btn.setToolTip("Project density onto real-space grid")
    s2_btn.clicked.connect(lambda: (run_afm_stage2(window), _refresh_dirty_label()))
    stage_layout.addWidget(s2_btn)

    s3_btn = QtWidgets.QPushButton("S3: Pots")
    s3_btn.setToolTip("Compute Pauli/ES/vdW FDBM potentials")
    s3_btn.clicked.connect(lambda: (run_afm_stage3(window), _refresh_dirty_label()))
    stage_layout.addWidget(s3_btn)

    s4_btn = QtWidgets.QPushButton("S4: Relax")
    s4_btn.setToolTip("Probe-particle relaxation -> AFM df + tip_disp")
    s4_btn.clicked.connect(lambda: (run_afm_stage4(window), _refresh_dirty_label()))
    stage_layout.addWidget(s4_btn)
    layout.addLayout(stage_layout)

    # STM/orbital row
    stm_row = QtWidgets.QHBoxLayout()
    stm_run_btn = QtWidgets.QPushButton("Run STM/BR-STM")
    stm_run_btn.setToolTip("Run stages 5/6 only (fast from cached S1+S4 data)")
    stm_run_btn.clicked.connect(lambda: (run_stm(window), _refresh_dirty_label()))
    stm_row.addWidget(stm_run_btn)

    orb_btn = QtWidgets.QPushButton("Plot Orbital")
    orb_btn.setToolTip("Plot selected MO with phase (needs S1)")
    orb_btn.clicked.connect(lambda: plot_orbital_map(window))
    stm_row.addWidget(orb_btn)
    layout.addLayout(stm_row)

    # --- Status display ---
    window.afm_status_label = QtWidgets.QPlainTextEdit()
    window.afm_status_label.setPlaceholderText("Status messages will appear here...")
    window.afm_status_label.setMaximumHeight(70)
    window.afm_status_label.setReadOnly(True)
    layout.addWidget(window.afm_status_label)

    # --- Parameters section ---
    param_sec = CollapsibleSection("Parameters", collapsed=True, parent=panel)
    param_widget = QtWidgets.QWidget()
    param_layout = QtWidgets.QVBoxLayout(param_widget)
    param_layout.setSpacing(2)
    param_layout.setContentsMargins(0, 0, 0, 0)

    density_group = QtWidgets.QGroupBox("Density / Grid")
    density_grid = QtWidgets.QGridLayout(density_group)
    density_grid.addWidget(QtWidgets.QLabel("Basis:"), 0, 0)
    window.afm_basis_combo = QtWidgets.QComboBox()
    window.afm_basis_combo.addItems(["mio-1-1", "3ob-3-1"])
    density_grid.addWidget(window.afm_basis_combo, 0, 1)
    density_grid.addWidget(QtWidgets.QLabel("Step:"), 1, 0)
    window.afm_step_spin = QtWidgets.QDoubleSpinBox()
    window.afm_step_spin.setRange(0.05, 0.5); window.afm_step_spin.setValue(0.1); window.afm_step_spin.setSingleStep(0.05)
    density_grid.addWidget(window.afm_step_spin, 1, 1)
    density_grid.addWidget(QtWidgets.QLabel("Margin:"), 2, 0)
    window.afm_margin_spin = QtWidgets.QDoubleSpinBox()
    window.afm_margin_spin.setRange(2.0, 10.0); window.afm_margin_spin.setValue(4.0)
    density_grid.addWidget(window.afm_margin_spin, 2, 1)
    param_layout.addWidget(density_group)

    scan_group = QtWidgets.QGroupBox("Scan")
    scan_grid = QtWidgets.QGridLayout(scan_group)
    scan_grid.addWidget(QtWidgets.QLabel("Range:"), 0, 0)
    window.afm_scan_range_spin = QtWidgets.QDoubleSpinBox()
    window.afm_scan_range_spin.setRange(1.0, 10.0); window.afm_scan_range_spin.setValue(3.0)
    scan_grid.addWidget(window.afm_scan_range_spin, 0, 1)
    scan_grid.addWidget(QtWidgets.QLabel("H min:"), 1, 0)
    window.afm_hmin_spin = QtWidgets.QDoubleSpinBox()
    window.afm_hmin_spin.setRange(1.5, 5.0); window.afm_hmin_spin.setValue(2.8)
    scan_grid.addWidget(window.afm_hmin_spin, 1, 1)
    scan_grid.addWidget(QtWidgets.QLabel("H max:"), 2, 0)
    window.afm_hmax_spin = QtWidgets.QDoubleSpinBox()
    window.afm_hmax_spin.setRange(3.0, 8.0); window.afm_hmax_spin.setValue(3.6)
    scan_grid.addWidget(window.afm_hmax_spin, 2, 1)
    scan_grid.addWidget(QtWidgets.QLabel("H step:"), 3, 0)
    window.afm_hstep_spin = QtWidgets.QDoubleSpinBox()
    window.afm_hstep_spin.setRange(0.05, 0.3); window.afm_hstep_spin.setValue(0.1)
    scan_grid.addWidget(window.afm_hstep_spin, 3, 1)
    param_layout.addWidget(scan_group)

    physics_group = QtWidgets.QGroupBox("Physics")
    physics_grid = QtWidgets.QGridLayout(physics_group)
    physics_grid.addWidget(QtWidgets.QLabel("Pauli A:"), 0, 0)
    window.afm_pauli_a_spin = QtWidgets.QDoubleSpinBox()
    window.afm_pauli_a_spin.setRange(0.1, 2000.0); window.afm_pauli_a_spin.setValue(787.22); window.afm_pauli_a_spin.setDecimals(2)
    physics_grid.addWidget(window.afm_pauli_a_spin, 0, 1)
    physics_grid.addWidget(QtWidgets.QLabel("Beta:"), 1, 0)
    window.afm_pauli_beta_spin = QtWidgets.QDoubleSpinBox()
    window.afm_pauli_beta_spin.setRange(0.5, 3.0); window.afm_pauli_beta_spin.setValue(1.2371); window.afm_pauli_beta_spin.setDecimals(4)
    physics_grid.addWidget(window.afm_pauli_beta_spin, 1, 1)
    physics_grid.addWidget(QtWidgets.QLabel("C6:"), 2, 0)
    window.afm_vdw_c6_spin = QtWidgets.QDoubleSpinBox()
    window.afm_vdw_c6_spin.setRange(10.0, 100.0); window.afm_vdw_c6_spin.setValue(30.0)
    physics_grid.addWidget(window.afm_vdw_c6_spin, 2, 1)
    physics_grid.addWidget(QtWidgets.QLabel("K_LAT:"), 3, 0)
    window.afm_klat_spin = QtWidgets.QDoubleSpinBox()
    window.afm_klat_spin.setRange(0.1, 2.0); window.afm_klat_spin.setValue(0.5)
    physics_grid.addWidget(window.afm_klat_spin, 3, 1)
    param_layout.addWidget(physics_group)

    def on_basis_changed(idx):
        basis = window.afm_basis_combo.currentText()
        window._afm_dirty.mark_geometry_changed()
        _refresh_dirty_label()
        if basis == 'mio-1-1':
            window.afm_pauli_a_spin.setValue(787.22); window.afm_pauli_beta_spin.setValue(1.2371)
        elif basis == '3ob-3-1':
            window.afm_pauli_a_spin.setValue(509.28); window.afm_pauli_beta_spin.setValue(1.0586)
    window.afm_basis_combo.currentIndexChanged.connect(on_basis_changed)

    def _mark_s2(): window._afm_dirty.mark_density_params_changed(); _refresh_dirty_label()
    def _mark_s3(): window._afm_dirty.mark_physics_params_changed(); _refresh_dirty_label()
    def _mark_s4(): window._afm_dirty.mark_scan_params_changed(); _refresh_dirty_label()
    def _mark_s56(): window._afm_dirty.mark_stm_params_changed(); _refresh_dirty_label()

    window.afm_step_spin.valueChanged.connect(_mark_s2)
    window.afm_margin_spin.valueChanged.connect(_mark_s2)
    window.afm_pauli_a_spin.valueChanged.connect(_mark_s3)
    window.afm_pauli_beta_spin.valueChanged.connect(_mark_s3)
    window.afm_vdw_c6_spin.valueChanged.connect(_mark_s3)
    window.afm_klat_spin.valueChanged.connect(_mark_s4)
    window.afm_scan_range_spin.valueChanged.connect(_mark_s4)
    window.afm_hmin_spin.valueChanged.connect(_mark_s4)
    window.afm_hmax_spin.valueChanged.connect(_mark_s4)
    window.afm_hstep_spin.valueChanged.connect(_mark_s4)

    param_sec.setContent(param_widget)
    layout.addWidget(param_sec)

    # --- Visualization section ---
    viz_sec = CollapsibleSection("Visualization", collapsed=True, parent=panel)
    viz_widget = QtWidgets.QWidget()
    viz_layout = QtWidgets.QVBoxLayout(viz_widget)
    viz_layout.setSpacing(2)

    window.afm_component_combo = QtWidgets.QComboBox()
    window.afm_component_combo.addItems([
        "AFM Image (df)", "STM Signal", "BR-STM Signal",
        "SCF Density", "Neutral Density", "Delta Density",
        "Pauli Energy", "Electrostatic Energy", "vdW Energy",
        "Total Potential", "Total Z-Force",
    ])
    viz_layout.addWidget(QtWidgets.QLabel("Component:"))
    viz_layout.addWidget(window.afm_component_combo)

    z_layout = QtWidgets.QHBoxLayout()
    z_layout.addWidget(QtWidgets.QLabel("Z-height (A):"))
    window.afm_z_height_spin = QtWidgets.QDoubleSpinBox()
    window.afm_z_height_spin.setRange(-20.0, 20.0); window.afm_z_height_spin.setValue(3.0)
    window.afm_z_height_spin.setSingleStep(0.1); window.afm_z_height_spin.setDecimals(2)
    z_layout.addWidget(window.afm_z_height_spin)
    window.afm_live_update = QtWidgets.QCheckBox("Live")
    window.afm_live_update.setChecked(True)
    z_layout.addWidget(window.afm_live_update)
    viz_layout.addLayout(z_layout)

    def on_z_height_changed():
        if window.afm_live_update.isChecked():
            has_data = (window._afm_results is not None) or (window._afm_potentials is not None) or (window._afm_density is not None)
            if has_data:
                try: plot_afm_slice(window)
                except Exception: pass
    window.afm_z_height_spin.valueChanged.connect(on_z_height_changed)

    lim_layout = QtWidgets.QHBoxLayout()
    window.afm_auto_limits = QtWidgets.QCheckBox("Auto limits")
    window.afm_auto_limits.setChecked(True)
    lim_layout.addWidget(window.afm_auto_limits)
    window.afm_vmin_spin = QtWidgets.QDoubleSpinBox()
    window.afm_vmin_spin.setRange(-1000, 1000); window.afm_vmin_spin.setValue(-1.0)
    window.afm_vmin_spin.setEnabled(False); window.afm_vmin_spin.setDecimals(3)
    lim_layout.addWidget(QtWidgets.QLabel("vmin:")); lim_layout.addWidget(window.afm_vmin_spin)
    window.afm_vmax_spin = QtWidgets.QDoubleSpinBox()
    window.afm_vmax_spin.setRange(-1000, 1000); window.afm_vmax_spin.setValue(1.0)
    window.afm_vmax_spin.setEnabled(False); window.afm_vmax_spin.setDecimals(3)
    lim_layout.addWidget(QtWidgets.QLabel("vmax:")); lim_layout.addWidget(window.afm_vmax_spin)
    viz_layout.addLayout(lim_layout)
    window.afm_auto_limits.stateChanged.connect(lambda s: (window.afm_vmin_spin.setEnabled(not s), window.afm_vmax_spin.setEnabled(not s)))

    window.afm_show_atoms = QtWidgets.QCheckBox("Overlay atom positions")
    window.afm_show_atoms.setChecked(True)
    window.afm_show_atoms.setToolTip("Show atom positions as colored dots on AFM/STM/orbital plots")
    viz_layout.addWidget(window.afm_show_atoms)

    plot_btn = QtWidgets.QPushButton("Plot Slice")
    plot_btn.clicked.connect(lambda: plot_afm_slice(window))
    viz_layout.addWidget(plot_btn)
    diag_btn = QtWidgets.QPushButton("Diagnostic Panel")
    diag_btn.clicked.connect(lambda: plot_afm_diagnostic_panel(window))
    viz_layout.addWidget(diag_btn)
    viz_sec.setContent(viz_widget)
    layout.addWidget(viz_sec)

    # --- STM / Orbital section ---
    stm_sec = CollapsibleSection("STM / Orbitals", collapsed=True, parent=panel)
    stm_widget = QtWidgets.QWidget()
    stm_grid = QtWidgets.QGridLayout(stm_widget)
    row = 0

    window.afm_stm_enable = QtWidgets.QCheckBox("Compute STM in full pipeline")
    stm_grid.addWidget(window.afm_stm_enable, row, 0, 1, 2); row += 1

    window.afm_stm_bond_resolved = QtWidgets.QCheckBox("Bond-resolved (BR-STM)")
    stm_grid.addWidget(window.afm_stm_bond_resolved, row, 0, 1, 2); row += 1

    # HOMO reference (read-only info label, updated after Stage 1)
    stm_grid.addWidget(QtWidgets.QLabel("HOMO iMO:"), row, 0)
    window.afm_homo_label = QtWidgets.QLabel("(run SCF first)")
    window.afm_homo_label.setStyleSheet("font-weight: bold; color: #006600;")
    stm_grid.addWidget(window.afm_homo_label, row, 1); row += 1

    # MO list: space/comma-separated integers
    stm_grid.addWidget(QtWidgets.QLabel("MO list:"), row, 0)
    window.afm_stm_mo_list = QtWidgets.QLineEdit("-1 0 1")
    window.afm_stm_mo_list.setToolTip("Space/comma separated. Relative to HOMO if checkbox below is ticked (0=HOMO, +1=LUMO, -1=HOMO-1). Absolute otherwise.")
    stm_grid.addWidget(window.afm_stm_mo_list, row, 1); row += 1

    window.afm_stm_relative_mo = QtWidgets.QCheckBox("Relative to HOMO")
    window.afm_stm_relative_mo.setChecked(True)
    window.afm_stm_relative_mo.setToolTip("If checked, MO list is relative to HOMO index. 0=HOMO, +1=LUMO, -1=HOMO-1 etc.")
    stm_grid.addWidget(window.afm_stm_relative_mo, row, 0, 1, 2); row += 1

    stm_grid.addWidget(QtWidgets.QLabel("field:"), row, 0)
    window.afm_stm_field_combo = QtWidgets.QComboBox()
    window.afm_stm_field_combo.addItems(['ldos', 'psi2', 'psi'])
    stm_grid.addWidget(window.afm_stm_field_combo, row, 1); row += 1

    stm_grid.addWidget(QtWidgets.QLabel("exp_beta:"), row, 0)
    window.afm_stm_exp_beta = QtWidgets.QDoubleSpinBox()
    window.afm_stm_exp_beta.setRange(0.1, 10.0); window.afm_stm_exp_beta.setValue(1.0); window.afm_stm_exp_beta.setDecimals(3)
    stm_grid.addWidget(window.afm_stm_exp_beta, row, 1); row += 1

    stm_grid.addWidget(QtWidgets.QLabel("exp_r0:"), row, 0)
    window.afm_stm_exp_r0 = QtWidgets.QDoubleSpinBox()
    window.afm_stm_exp_r0.setRange(0.0, 10.0); window.afm_stm_exp_r0.setValue(3.0); window.afm_stm_exp_r0.setDecimals(3)
    stm_grid.addWidget(window.afm_stm_exp_r0, row, 1); row += 1

    # Orbital map (single MO with phase)
    stm_grid.addWidget(QtWidgets.QLabel("── Orbital Map ──"), row, 0, 1, 2); row += 1
    stm_grid.addWidget(QtWidgets.QLabel("iMO (abs):"), row, 0)
    window.afm_orbital_spin = QtWidgets.QSpinBox()
    window.afm_orbital_spin.setRange(0, 999); window.afm_orbital_spin.setValue(0)
    window.afm_orbital_spin.setToolTip("Absolute MO index (0=lowest). HOMO shown above.")
    stm_grid.addWidget(window.afm_orbital_spin, row, 1); row += 1

    # Mark STM dirty when relevant params change
    for w in [window.afm_stm_exp_beta, window.afm_stm_exp_r0]:
        w.valueChanged.connect(_mark_s56)
    window.afm_stm_mo_list.textChanged.connect(_mark_s56)
    window.afm_stm_relative_mo.stateChanged.connect(_mark_s56)
    window.afm_stm_field_combo.currentIndexChanged.connect(_mark_s56)
    window.afm_stm_bond_resolved.stateChanged.connect(_mark_s56)

    stm_sec.setContent(stm_widget)
    layout.addWidget(stm_sec)

    # --- Geometry change hook ---
    # Connect to KekuleExplorerGUI geometry-change signal if available
    if hasattr(window, 'sig_geometry_changed'):
        def _on_geom_changed():
            window._afm_dirty.mark_geometry_changed()
            _refresh_dirty_label()
        window.sig_geometry_changed.connect(_on_geom_changed)

    view_modes = []
    return UIComponents(panel=panel, view_modes=view_modes)
