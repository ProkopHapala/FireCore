import numpy as np
import pyopencl as cl
import pyopencl.array as cl_array
import time
import os

from pyBall.AtomicSystem import AtomicSystem
from pyBall.OCL.MMFF import MMFF as MMFF_pyocl
from pyBall.OCL.MolecularDynamics import MolecularDynamics
import matplotlib.pyplot as plt

class ManipulationPathOpt:
    def __init__(self, mol, tip_handle_idx, n_pop, n_steps, enable_nonbond=False, tip_k=1.0):
        self.mol = mol
        self.tip_handle_idx = tip_handle_idx
        self.n_pop = n_pop
        self.n_steps = n_steps
        self.n_rep = n_pop * n_steps
        self.tip_k = tip_k
        
        # Build MMFF
        self.mm = MMFF_pyocl(bTorsion=False, verbosity=0)
        self.mm.capping_atoms = set() # treat all atoms as nodes
        self.mm.reorder_nodes_first = False
        self.mm.toMMFFsp3_loc(mol=mol, atom_types=self.mm.atom_types, bRealloc=True, bEPairs=False, bUFF=False)

        qs = getattr(mol, 'qs', None)
        if qs is not None:
            qs = np.asarray(qs, dtype=np.float32)
            if qs.shape[0] == self.mm.natoms:
                self.mm.REQs[:, 2] = qs
        self.mm.make_back_neighs(b_cap_neighs=False)
        if enable_nonbond and (self.mm.excl is None):
            self.mm.excl = self.mm._make_excl_1_2_3(self.mm.neighs, neighCell=self.mm.neighCell, npbc=self.mm.npbc, EXCL_MAX=16)

        # Setup MD
        self.md = MolecularDynamics(nloc=32, debug_build_options='-DDBG_UFF=0', enable_nonbond=enable_nonbond)
        # IMPORTANT: relax_multi.cl updateAtomsMMFFf4 reads sysneighs/sysbonds as int[nSystems*nMaxSysNeighs]
        # but MolecularDynamics allocates these buffers only for nMaxSysNeighs==1.
        # Therefore we must keep nDOFs.w (nMaxSysNeighs) == 1, otherwise kernel reads out-of-bounds and crashes the queue.
        self.md.perBatch = 1
        self.md.realloc(self.mm, nSystems=self.n_rep)

        self.md.setup_kernels()
        # Pack identical mmff into all systems
        for iSys in range(self.n_rep):
            self.md.pack_system(iSys, self.mm)

        # Use kernels compiled in MolecularDynamics (relax_multi.cl)
        self.prog = self.md.prg

        # Buffers for PathDiffusion
        self.cl_tip_positions = None
        self.cl_apos_initial = None
        self.cl_max_tip_force = cl_array.zeros(self.md.queue, self.n_rep, dtype=np.float32)
        self.cl_gap_sizes = cl_array.zeros(self.md.queue, self.n_rep, dtype=np.float32)
        self.morse_params = np.array([1.0, 2.0, 1.5, 0.0], dtype=np.float32) # De, r0, alpha, padding
        self.morse_params_host = self.morse_params.copy()

        self.use_tip_constraint = False
        self.tip_constraint_K = 0.0

        self.gridff_enabled = False

    def set_substrate_gridff(self, bspline_PLQ, grid_p0=(-2.0, -2.0, 0.0), grid_step=(0.1, 0.1, 0.1), use_texture=False, r_damp=0.0, alpha_morse=1.5):
        bspline_PLQ = np.asarray(bspline_PLQ)
        if bspline_PLQ.ndim != 4:
            raise ValueError(f"GridFF bspline_PLQ must have shape (nx,ny,nz,nch); got {bspline_PLQ.shape}")
        if bspline_PLQ.shape[3] != 4:
            raise ValueError(f"GridFF bspline_PLQ last dim must be 4 (float4 per grid point); got {bspline_PLQ.shape}")

        grid_shape = tuple(int(x) for x in bspline_PLQ.shape[:3])
        self.md.initGridFF(
            grid_shape=grid_shape,
            bspline_data=np.ascontiguousarray(bspline_PLQ.astype(np.float32)),
            grid_p0=tuple(float(x) for x in grid_p0),
            grid_step=tuple(float(x) for x in grid_step),
            use_texture=bool(use_texture),
            r_damp=float(r_damp),
            alpha_morse=float(alpha_morse),
            bKernels=True,
        )
        self.gridff_enabled = True

    def set_tip_trajectory(self, tip_positions):
        """
        tip_positions: array of shape (n_pop, n_steps, 3)
        """
        assert tip_positions.shape == (self.n_pop, self.n_steps, 3)
        tip_flat = tip_positions.reshape(self.n_rep, 3)
        
        # We will use getTipMorse instead of constr, but we still need the positions on GPU
        tip_flat_4 = np.zeros((self.n_rep, 4), dtype=np.float32)
        tip_flat_4[:, :3] = tip_flat
        self.cl_tip_positions = cl_array.to_device(self.md.queue, tip_flat_4)

    def set_tip_pos_constraints(self, K=1000.0):
        """Hard positional constraint of `tip_handle_idx` to the per-replica tip positions.

        Uses buffers `constr` and `constrK` consumed by `updateAtomsMMFFf4`.

        NOTE: If enabled, you usually want to disable `_run_tip_morse()` because otherwise
        you're applying two different tip models at once.
        """
        if self.cl_tip_positions is None:
            raise RuntimeError("set_tip_pos_constraints(): call set_tip_trajectory() first")
        K = float(K)
        self.use_tip_constraint = True
        self.tip_constraint_K = K

        # updateAtomsMMFFf4 indexes constraints by iaa=iG+iS*natoms (only atoms have constraints)
        natoms = self.mm.natoms
        constr = np.zeros((self.n_rep, natoms, 4), dtype=np.float32)
        constrK = np.zeros((self.n_rep, natoms, 4), dtype=np.float32)

        tip = self.cl_tip_positions.get()  # (n_rep,4)
        ia = int(self.tip_handle_idx)
        constr[:, ia, :3] = tip[:, :3]
        constr[:, ia, 3] = 1.0
        constrK[:, ia, 0] = K
        constrK[:, ia, 1] = K
        constrK[:, ia, 2] = K

        self.md.toGPU('constr', constr)
        self.md.toGPU('constrK', constrK)

    def init_replica_states(self, initial_pos, tip_positions, bShift=True):
        # Upload tip trajectory first (used for shifts)
        self.set_tip_trajectory(tip_positions)

        # IMPORTANT: bonded forces depend on node positions (nnode) in the MMFF layout.
        # The OpenCL helper kernel init_trajectory_shift initializes ONLY natoms and leaves nodes undefined,
        # which corrupts bonded forces and can destroy H2O geometry. Therefore we build full nvecs init on host.

        # Start from the already-packed reference system0 geometry (includes atoms + nodes)
        apos0_full = np.empty((self.mm.nvecs, 4), dtype=np.float32)
        cl.enqueue_copy(self.md.queue, apos0_full, self.md.buffer_dict['apos'], device_offset=0)
        self.md.queue.finish()

        # Overwrite atom positions from provided initial_pos (keep node layout as in MMFF packing)
        apos0_full[:self.mm.natoms, :3] = np.asarray(initial_pos, dtype=np.float32)

        # Replicate to all systems, optionally apply rigid shift by tip displacement
        apos_all = np.zeros((self.n_rep, self.mm.nvecs, 4), dtype=np.float32)
        apos_all[:, :, :] = apos0_full[None, :, :]

        if bShift:
            tip = self.cl_tip_positions.get().astype(np.float32)  # (n_rep,4)
            tip0 = tip[(np.arange(self.n_rep, dtype=np.int32) // np.int32(self.n_steps)) * np.int32(self.n_steps), :]
            shifts = (tip[:, :3] - tip0[:, :3]).astype(np.float32)  # (n_rep,3)
            # IMPORTANT: only atoms are translated. The "node" vectors (pi-orbitals) stored in apos[natoms:natoms+nnode]
            # represent directions and must NOT be translated, otherwise bonded forces become garbage.
            apos_all[:, :self.mm.natoms, :3] += shifts[:, None, :]

        self.md.toGPU('apos', apos_all)

        # If hard tip constraints are enabled, refresh them now (cl_tip_positions was just re-uploaded).
        if self.use_tip_constraint and (self.tip_constraint_K > 0.0):
            self.set_tip_pos_constraints(K=self.tip_constraint_K)

    def set_causal_tethers(self, K_band, L_allowed):
        sysneighs = np.zeros(self.n_rep, dtype=np.int32)
        sysbonds = np.zeros((self.n_rep, 4), dtype=np.float32)
        
        for ipop in range(self.n_pop):
            for istep in range(self.n_steps):
                irep = ipop * self.n_steps + istep
                if istep == 0:
                    sysneighs[irep] = irep # no pull for the first step
                    sysbonds[irep] = [0.0, 1e5, 0.0, 0.0]
                else:
                    sysneighs[irep] = irep - 1
                    # Ktens must be negative to pull towards neighbor (pj), because kernel uses (Lmax - l) * Ktens
                    sysbonds[irep] = [0.0, L_allowed, 0.0, -K_band] # Lmin, Lmax, Kpress, Ktens
                    
        self.md.toGPU('sysneighs', sysneighs)
        self.md.toGPU('sysbonds', sysbonds)

    def set_morse_params(self, De, r0, alpha):
        self.morse_params_host[0] = De
        self.morse_params_host[1] = r0
        self.morse_params_host[2] = alpha

    def _run_tip_morse(self):
        self.prog.getTipMorse(
            self.md.queue, (self.n_rep,), None,
            np.int32(self.mm.natoms), np.int32(self.tip_handle_idx), np.int32(self.mm.nvecs),
            self.md.buffer_dict['apos'],
            self.md.buffer_dict['aforce'],
            self.cl_tip_positions.data,
            np.float32(self.morse_params_host[0]),
            np.float32(self.morse_params_host[1]),
            np.float32(self.morse_params_host[2]),
            self.cl_max_tip_force.data
        )

    def _reduce_trajectory_gaps_gpu(self):
        self.prog.reduce_trajectory_gaps(
            self.md.queue, (self.n_rep,), None,
            np.int32(self.mm.natoms), np.int32(self.n_steps), np.int32(self.mm.nvecs),
            self.md.buffer_dict['apos'],
            self.cl_gap_sizes.data
        )
        self.md.queue.finish()
        gaps_flat = self.cl_gap_sizes.get()
        return gaps_flat.reshape((self.n_pop, self.n_steps))

    def run_relax(self, n_steps_md, K_band_sched, L_allowed, dt=0.01, damp=0.95, Flimit=10.0, out_dir="out", bTrace=False, trace_interval=10):
        # Populate MDparams for ALL replicas
        # NOTE: updateAtomsMMFFf4 uses MDpars.x (dt) and MDpars.z (friction, velocity damping: ve.xyz *= MDpars.z)
        md_params_arr = np.zeros((self.n_rep, 4), dtype=np.float32)
        # Keep semantics consistent with tests/tMMFF/test_ditetraceno_surface.py and test_relaxed_scan_tip.py
        # where MDparams = [dt, 1e6, vfac, 0] and vfac = 1-damp (damp is a small number like 0.01)
        vfac = 1.0 - float(damp)
        md_params_arr[:, :] = [dt, 1e6, vfac, 0.0]
        self.md.toGPU('MDparams', md_params_arr)
        
        epochs = len(K_band_sched)
        steps_per_epoch = n_steps_md // epochs
        
        os.makedirs(out_dir, exist_ok=True)
        
        history_apos = []
        history_gaps = []
        trace_apos = []

        def _geom_check(apos_rep, istep=None, irep_hint=None, hard_fail=True):
            # apos_rep: (n_rep,natoms,3)
            if apos_rep.shape[1] < 3:
                return
            rO  = apos_rep[:, 0, :]
            rH1 = apos_rep[:, 1, :]
            rH2 = apos_rep[:, 2, :]
            OH1 = np.linalg.norm(rH1 - rO, axis=1)
            OH2 = np.linalg.norm(rH2 - rO, axis=1)
            HH  = np.linalg.norm(rH2 - rH1, axis=1)
            # hard_fail=False: warn only (used mid-annealing when band forces can transiently stretch bonds)
            # hard_fail=True:  raise if epoch-end geometry is unphysical
            bad = (HH < 0.7) | (OH1 < 0.7) | (OH2 < 0.7) | (OH1 > 1.6) | (OH2 > 1.6)
            if np.any(bad):
                bad_reps = np.where(bad)[0]
                i = int(bad_reps[0])  # first bad replica (not argmin HH)
                if irep_hint is not None:
                    i = int(irep_hint)
                print(f"{'ERROR' if hard_fail else 'WARN'} _geom_check bad reps: {bad_reps.tolist()} step={istep}")
                print(f"  rep={i}: HH={HH[i]:.4f} OH1={OH1[i]:.4f} OH2={OH2[i]:.4f}")
                if not hard_fail:
                    return
                # Dump offending replica for post-mortem
                try:
                    f_rep = self.download_forces()[i, :, :].copy()
                except Exception:
                    f_rep = None
                np.save(os.path.join(out_dir, f"BAD_step_{istep}_rep_{i}_apos.npy"), apos_rep[i])
                if f_rep is not None:
                    np.save(os.path.join(out_dir, f"BAD_step_{istep}_rep_{i}_force.npy"), f_rep)
                xyzp = os.path.join(out_dir, f"BAD_step_{istep}_rep_{i}.xyz")
                with open(xyzp, 'w') as fout:
                    fout.write(f"{self.mm.natoms}\n")
                    fout.write(f"BAD step={istep} rep={i} HH={HH[i]:.6f} OH1={OH1[i]:.6f} OH2={OH2[i]:.6f}\n")
                    for ia in range(self.mm.natoms):
                        name = self.mol.enames[ia] if hasattr(self.mol, 'enames') else 'X'
                        r = apos_rep[i, ia, :]
                        fout.write(f"{name} {r[0]:.6f} {r[1]:.6f} {r[2]:.6f}\n")
                raise RuntimeError(f"Unphysical H2O geometry detected at step={istep} rep={i}: HH={HH[i]:.4f} OH1={OH1[i]:.4f} OH2={OH2[i]:.4f}")

        def _forces_fmax():
            f = self.download_forces()[:, :, :3]  # (n_rep,natoms,3)
            fm = np.linalg.norm(f, axis=2)        # (n_rep,natoms)
            return fm.max(axis=1)                 # (n_rep,)
        
        # Convergence parameters (match sequential scan style)
        Fconv = 1e-4
        check_interval = 1

        for ep in range(epochs):
            K_band = K_band_sched[ep]
            self.set_causal_tethers(K_band, L_allowed)
            self.cl_max_tip_force.fill(0.0)

            # Integrate until converged or hit step budget for this epoch
            for step in range(steps_per_epoch):
                if self.md.kernel_args_cleanForceMMFFf4 is not None:
                    self.md.run_cleanForceMMFFf4()
                else:
                    self.md.toGPU('aforce', np.zeros((self.n_rep, self.mm.nvecs, 4), dtype=np.float32))

                # Non-bonded / GridFF first (per user request)
                if getattr(self.md, 'has_gridff', False):
                    if getattr(self.md, 'use_texture', False) and (getattr(self.md, 'kernel_args_getNonBond_GridFF_Bspline_tex', None) is not None):
                        self.md.run_getNonBond_GridFF_Bspline_tex()
                    elif getattr(self.md, 'kernel_args_getNonBond_GridFF_Bspline_ex2', None) is not None:
                        self.md.run_getNonBond_GridFF_Bspline_ex2()
                    elif getattr(self.md, 'kernel_args_getNonBond_GridFF_Bspline', None) is not None:
                        self.md.run_getNonBond_GridFF_Bspline()
                    else:
                        raise RuntimeError("GridFF enabled but getNonBond_GridFF_Bspline kernel args are not initialized")
                else:
                    if self.md.enable_nonbond:
                        if self.md.kernel_args_getNonBond_ex2 is not None:
                            self.md.run_getNonBond_ex2()
                        elif self.md.kernel_args_getNonBond is not None:
                            self.md.run_getNonBond()

                # Bonded after non-bonded
                if self.md.kernel_args_getMMFFf4 is not None:
                    self.md.run_getMMFFf4()
                    self.md.queue.finish()
                    
                # Tip constraint
                if not self.use_tip_constraint:
                    self._run_tip_morse()
                    self.md.queue.finish()
                    
                # Update (integrator & tethers)
                if self.md.kernel_args_updateAtomsMMFFf4 is not None:
                    self.md.run_updateAtomsMMFFf4()
                    self.md.queue.finish()

                if bTrace and (step % trace_interval == 0):
                    trace_apos.append(self.download_states())

                if (step % check_interval) == 0:
                    apos_now = self.download_states()
                    _geom_check(apos_now, istep=step, hard_fail=False)
                    fmax = _forces_fmax()
                    if np.all(fmax < Fconv):
                        break

            # End of epoch: download states, save xyz, record history
            apos_ep = self.download_states()
            _geom_check(apos_ep, istep='epoch_end')
            gaps_ep = self._reduce_trajectory_gaps_gpu()
            
            history_apos.append(apos_ep.copy())
            history_gaps.append(gaps_ep.copy())
            
            self.save_xyz_movie(apos_ep, os.path.join(out_dir, f"traj_K_{K_band:.2f}.xyz"))

        return history_apos, history_gaps, trace_apos

    def save_xyz_movie(self, apos, filename):
        """
        apos: (n_rep, natoms, 3)
        """
        with open(filename, 'w') as f:
            for iSys in range(self.n_rep):
                f.write(f"{self.mm.natoms}\n")
                f.write(f"Replica {iSys} pop {iSys // self.n_steps} step {iSys % self.n_steps}\n")
                for ia in range(self.mm.natoms):
                    name = self.mol.enames[ia] if hasattr(self.mol, 'enames') else "C"
                    r = apos[iSys, ia, :]
                    f.write(f"{name} {r[0]:.5f} {r[1]:.5f} {r[2]:.5f}\n")

    def plot_trajectories(self, history_apos, history_gaps, K_band_sched, plot_atom_idx=0, out_dir="out"):
        epochs = len(K_band_sched)
        
        # Plot 1: Trajectory of selected atom
        plt.figure(figsize=(8, 6))
        for ep in range(epochs):
            K_band = K_band_sched[ep]
            # Get pop=0, steps
            apos_pop = history_apos[ep].reshape((self.n_pop, self.n_steps, self.mm.natoms, 3))
            traj_x = apos_pop[0, :, plot_atom_idx, 0]
            traj_y = apos_pop[0, :, plot_atom_idx, 1]
            plt.plot(traj_x, traj_y, marker='o', label=f'K_band = {K_band:.2f}')
            
        plt.title(f"Trajectory of Atom {plot_atom_idx}")
        plt.xlabel("X coordinate (A)")
        plt.ylabel("Y coordinate (A)")
        plt.legend()
        plt.grid(True)
        plt.savefig(os.path.join(out_dir, "atom_trajectory.png"))
        plt.close()
        
        # Plot 2: Gap sizes (penalty) along the path
        plt.figure(figsize=(8, 6))
        for ep in range(epochs):
            K_band = K_band_sched[ep]
            gaps_pop0 = history_gaps[ep][0, :] # shape: (n_steps)
            plt.plot(range(self.n_steps), gaps_pop0, marker='x', label=f'K_band = {K_band:.2f}')
            
        plt.title("Inter-replica Gaps (Penalty) vs Step")
        plt.xlabel("Step index")
        plt.ylabel("Max Atom Gap (A)")
        plt.legend()
        plt.grid(True)
        plt.savefig(os.path.join(out_dir, "gap_penalties.png"))
        plt.close()

    def download_states(self):
        apos_buf = np.empty((self.n_rep, self.mm.nvecs, 4), dtype=np.float32)
        cl.enqueue_copy(self.md.queue, apos_buf, self.md.buffer_dict['apos'])
        self.md.queue.finish()
        # Return only atom positions (natoms,3)
        return apos_buf[:, :self.mm.natoms, :3]

    def download_forces(self):
        f_buf = np.empty((self.n_rep, self.mm.nvecs, 4), dtype=np.float32)
        cl.enqueue_copy(self.md.queue, f_buf, self.md.buffer_dict['aforce'])
        self.md.queue.finish()
        return f_buf[:, :self.mm.natoms, :]

    def compute_gaps(self, apos):
        """
        Fall-back CPU gap compute if needed.
        """
        apos_pop = apos.reshape(self.n_pop, self.n_steps, self.mm.natoms, 3)
        diff = apos_pop[:, 1:, :, :] - apos_pop[:, :-1, :, :]
        dist = np.linalg.norm(diff, axis=-1) 
        gaps = np.max(dist, axis=-1) 
        return gaps

    def get_fine_points(self, gaps, tip_positions, gap_threshold=0.5, budget_fine=100):
        """
        Calculates dynamic n-section refinement for each population string based on gaps.
        Returns a list of refined tip positions (if gaps > threshold) for each trajectory.
        """
        refined_tips_per_pop = []
        for p in range(self.n_pop):
            torn = np.where(gaps[p] > gap_threshold)[0]
            if len(torn) == 0:
                refined_tips_per_pop.append(None)
                continue
            
            pts_per_interval = budget_fine // len(torn)
            fine_points = []
            for t_idx in torn:
                p0 = tip_positions[p, t_idx]
                p1 = tip_positions[p, t_idx + 1]
                # interpolate strictly inside the interval
                interp = np.linspace(p0, p1, pts_per_interval + 2)[1:-1]
                fine_points.append(interp)
                
            refined_tips_per_pop.append(np.vstack(fine_points))
            
        return refined_tips_per_pop

