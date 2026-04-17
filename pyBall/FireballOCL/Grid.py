import numpy as np
import pyopencl as cl
import pyopencl.cltypes
import os
import pyopencl.array as cl_array
import time
from ..OCL.OpenCLBase import OpenCLBase
from .FdataParser import FdataParser

class GridProjector(OpenCLBase):
    """
    Host class for projecting sparse density matrices to a real-space grid using OpenCL.
    """
    def __init__(self, fdata_dir, ctx=None, queue=None, nloc=32, debug_early_exit=False, debug_clear_only=False, debug_return0=False, debug_read_task=False, debug_read_grid=False, verbosity=0):
        super().__init__(nloc=nloc)
        self.fdata_dir = fdata_dir
        self.parser = FdataParser(fdata_dir)
        self.debug_early_exit = bool(debug_early_exit)
        self.debug_clear_only = bool(debug_clear_only)
        self.debug_return0 = bool(debug_return0)
        self.debug_read_task = bool(debug_read_task)
        self.debug_read_grid = bool(debug_read_grid)
        self.verbosity = int(verbosity)
        if not hasattr(self.parser, "species_info"):
            try:
                self.parser.parse_info()
            except Exception as e:
                # Keep going; parse_info will be invoked lazily later if needed
                pass
        if ctx:
            self.ctx = ctx
            self.queue = queue if queue else cl.CommandQueue(self.ctx)
        self.task_dtype = [
            ('x', 'i4'), ('y', 'i4'), ('z', 'i4'), ('w', 'i4'),
            ('na', 'i4'), ('nj', 'i4'), ('pad1', 'i4'), ('pad2', 'i4')
        ]
        self.task_dtype_np = np.dtype(self.task_dtype)
        self._load_kernels()
        self.basis_data = {}

    def load_basis(self, species_nz):
        """Loads radial basis functions for given species."""
        missing = []
        for nz in species_nz:
            if nz in self.basis_data: continue
            wfs = self.parser.find_wf(nz)
            if len(wfs)==0:
                missing.append(nz); continue
            self.basis_data[nz] = [self.parser.read_wf(f) for f in wfs]
        if missing:
            raise RuntimeError(f"No .wf files found for species {missing} under {self.fdata_dir}; ensure Fdata dir has *.ZZ.wf")
        
        # Prepare for GPU: pack into a single buffer
        # Resample all wavefunctions to a common uniform grid (finest dr, largest rcutoff)
        all_nz = sorted(self.basis_data.keys())
        if len(all_nz)==0:
            raise RuntimeError("load_basis called with empty species list (species_nz).")
        max_shells = max(len(v) for v in self.basis_data.values())
        if max_shells==0:
            raise RuntimeError(f"No wavefunctions loaded for species {all_nz}")
        # Find finest dr and largest rcutoff across all shells/species
        # IMPORTANT: .wf files store rcutoff in Bohr. Convert to Angstrom via abohr.
        # See Fortran read_wf.f90 line 208: drr_wf = abohr * rcutoffwf / (mesh-1)
        ABOHR = 0.529177       # Bohr -> Angstrom
        all_dr = []
        all_rc_ang = []
        for nz in all_nz:
            for wf in self.basis_data[nz]:
                wf_dr_ang = ABOHR * wf['rcutoff'] / (wf['mesh'] - 1)
                all_dr.append(wf_dr_ang)
                all_rc_ang.append(ABOHR * wf['rcutoff'])
        dr = min(all_dr)                   # finest spacing in Angstrom
        rc_max_ang = max(all_rc_ang)       # largest cutoff in Angstrom
        n_nodes = int(np.ceil(rc_max_ang / dr)) + 1
        if self.verbosity > 0: print(f"[DEBUG] load_basis: common grid dr={dr:.6f} Å  rc_max={rc_max_ang:.3f} Å  n_nodes={n_nodes}")
        
        packed_basis = np.zeros((len(all_nz), max_shells, n_nodes), dtype=np.float32)
        for i, nz in enumerate(all_nz):
            for ish, wf in enumerate(self.basis_data[nz]):
                wf_mesh = wf['mesh']
                wf_rc_bohr = wf['rcutoff']
                wf_dr_ang  = ABOHR * wf_rc_bohr / (wf_mesh - 1)
                wf_data    = wf['data']
                # Resample from wf's own grid (Angstrom) to common grid via linear interpolation
                r_common = np.arange(n_nodes) * dr
                r_orig   = np.arange(wf_mesh) * wf_dr_ang
                resampled = np.interp(r_common, r_orig, wf_data, right=0.0)
                # Verify normalization: ∫ R² r² dr should be ~1.0 with correct Angstrom grid
                S_rad = np.trapz(resampled**2 * r_common**2, r_common)
                packed_basis[i, ish, :] = resampled.astype(np.float32)
                if self.verbosity > 0: print(f"[DEBUG]   species {nz} shell {ish} (l={wf.get('l','?')}): mesh={wf_mesh} rc={wf_rc_bohr:.3f} Bohr = {ABOHR*wf_rc_bohr:.3f} Å  S_rad={S_rad:.6f} (should be ~1.0)")
        
        self.d_basis = cl.Buffer(self.ctx, cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR, hostbuf=packed_basis)
        self.basis_meta = {'n_species': len(all_nz), 'max_shells': max_shells, 'n_nodes': n_nodes, 'dr': dr, 'nz_map': {nz: i for i, nz in enumerate(all_nz)}}
        return packed_basis

    def _load_kernels(self):
        cl_path = os.path.join(os.path.dirname(__file__), "cl/Grid.cl")
        # Ensure the directory and file exist
        os.makedirs(os.path.dirname(cl_path), exist_ok=True)
        if not os.path.exists(cl_path):
            with open(cl_path, "w") as f:
                f.write("// Grid projection kernels\n")
        
        # We might want to pass some constants to the kernel during build
        build_opts = []
        if self.debug_early_exit:
            build_opts.append("-DDEBUG_EARLY_EXIT=1")
        if self.debug_clear_only:
            build_opts.append("-DDEBUG_CLEAR_ONLY=1")
        if self.debug_return0:
            build_opts.append("-DDEBUG_RETURN0=1")
        if self.debug_read_task:
            build_opts.append("-DDEBUG_READ_TASK=1")
        if self.debug_read_grid:
            build_opts.append("-DDEBUG_READ_GRID=1")
        if hasattr(self, 'basis_meta') and ('max_shells' in self.basis_meta):
            build_opts.append(f"-DMAX_ORBS={self.basis_meta['max_shells'] * 9}") # conservative upper bound

        self.load_program(kernel_path=cl_path, build_options=build_opts if len(build_opts)>0 else None)

    def check_overlap_sphere_aabb(self, center, radius, box_min, box_max):
        """ Fast AABB-Sphere collision: Find closest point in box to sphere center """
        closest_p = np.clip(center, box_min, box_max)
        distance_sq = np.sum((center - closest_p)**2)
        return distance_sq < (radius**2)

    def build_tasks_gpu(self, atoms, grid_spec, block_res=8, nMaxAtom=64):
        """
        GPU-based task building using OpenCL kernels.
        Pseudocode:
        1) count_atoms_per_block: for each atom, find overlapping blocks (via floor-index range + sphere/AABB), atomic_inc block_counts[b].
        2) fill_task_atoms: for each atom, again walk overlapping blocks, atomic_inc block_offsets[b], write atom id into task_atoms_raw[b][slot] if slot < nMaxAtom.
        3) On host: read block_counts, derive mask, check max_count<=nMaxAtom, compute task_offsets = prefix over (mask).
        4) compact_tasks: for each block with count>0, write TaskData(x,y,z,na,nj=-1) at task_offsets[b], copy task_atoms_raw[b] into compacted task_atoms_out.
        5) Host copies tasks_np/task_atoms_np back; optional host sort by na desc.
        Note: compaction is only at block level (drop empty blocks); task_atoms remains padded to nMaxAtom per task (holes stay).
        """
        nx, ny, nz = grid_spec['ngrid'][:3]
        n_blocks_xyz = np.array([nx // block_res, ny // block_res, nz // block_res], dtype=np.int32)
        n_blocks_total = int(np.prod(n_blocks_xyz))
        natoms = len(atoms['pos'])

        # 1. Prepare AtomData buffer
        atom_data = np.zeros(natoms, dtype=[
            ('pos_rcut', 'f4', 4),
            ('type', 'i4'),
            ('i0orb', 'i4'),
            ('norb', 'i4'),
            ('pad', 'i4')
        ])
        for i in range(natoms):
            atom_data[i]['pos_rcut'][:3] = atoms['pos'][i]
            atom_data[i]['pos_rcut'][3]  = atoms['Rcut'][i]
            atom_data[i]['type'] = atoms['type'][i]
            atom_data[i]['norb'] = 4
            atom_data[i]['i0orb'] = 0
            
        # DEBUG: print first atom
        if natoms > 0 and self.verbosity > 0:
            print(f"[DEBUG] atom_data[0]: pos_rcut={atom_data[0]['pos_rcut']} type={atom_data[0]['type']}")

        mf = cl.mem_flags
        d_grid  = cl.Buffer(self.ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=self.grid_to_np(grid_spec))
        d_atoms = cl.Buffer(self.ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=atom_data)
        
        T0 = time.perf_counter_ns()
        # 2. Kernel 1: Count atoms per block
        d_block_counts = cl.Buffer(self.ctx, mf.READ_WRITE, n_blocks_total * 4)
        cl.enqueue_fill_buffer(self.queue, d_block_counts, np.int32(0), 0, n_blocks_total * 4)
        self.prg.count_atoms_per_block(
            self.queue, (natoms,), None,
            d_grid, np.int32(natoms), d_atoms, np.int32(block_res),
            np.int32(n_blocks_xyz[0]), np.int32(n_blocks_xyz[1]), np.int32(n_blocks_xyz[2]),
            d_block_counts
        )
        self.queue.finish()
        T1 = time.perf_counter_ns()
        print(f"[TIME] count_atoms_per_block {(T1-T0)*1e-6:.3f} [ms]")

        T0 = time.perf_counter_ns()
        # 3. Kernel 2: Fill task_atoms
        d_task_atoms_raw = cl.Buffer(self.ctx, mf.READ_WRITE, n_blocks_total * nMaxAtom * 4)
        cl.enqueue_fill_buffer(self.queue, d_task_atoms_raw, np.int32(-1), 0, n_blocks_total * nMaxAtom * 4)
        # We need a secondary counter for atomic increments during filling
        d_block_fill_counts = cl.Buffer(self.ctx, mf.READ_WRITE, n_blocks_total * 4)
        cl.enqueue_fill_buffer(self.queue, d_block_fill_counts, np.int32(0), 0, n_blocks_total * 4)
        self.prg.fill_task_atoms(
            self.queue, (natoms,), None,
            d_grid, np.int32(natoms), d_atoms, np.int32(block_res),
            np.int32(n_blocks_xyz[0]), np.int32(n_blocks_xyz[1]), np.int32(n_blocks_xyz[2]),
            d_block_fill_counts, d_task_atoms_raw, np.int32(nMaxAtom)
        )
        # 4. Compact tasks
        # Read back counts to host to identify non-empty blocks and compute stats
        h_block_counts = np.empty(n_blocks_total, dtype=np.int32)
        cl.enqueue_copy(self.queue, h_block_counts, d_block_counts)
        self.queue.finish()
        T1 = time.perf_counter_ns()
        print(f"[TIME] count_atoms_per_block.compact_tasks {(T1-T0)*1e-6:.3f} [ms]")

        mask = h_block_counts > 0
        n_tasks = np.sum(mask)
        
        # Stats
        max_count    = h_block_counts.max() if n_blocks_total > 0 else 0
        empty_blocks = np.sum(h_block_counts == 0)
        one_blocks   = np.sum(h_block_counts == 1)
        multi_blocks = n_blocks_total - empty_blocks - one_blocks
        print(f"[DEBUG GPU] block atom stats: na_max={max_count}, nbloks: empty={empty_blocks}, one={one_blocks}, multi={multi_blocks}")
        self.last_block_atom_counts = h_block_counts

        if max_count > nMaxAtom:
             raise RuntimeError(f"GPU build_tasks: block has {max_count} atoms > nMaxAtom={nMaxAtom}")

        # tasks_np must have the correct structured dtype even when empty
        self.task_dtype_np = np.dtype(self.task_dtype)
        if n_tasks == 0:
            return np.zeros(0, dtype=self.task_dtype_np), np.zeros((0, nMaxAtom), dtype=np.int32)



        # Compute task offsets for compaction
        h_task_offsets = np.zeros(n_blocks_total, dtype=np.int32)
        h_task_offsets[mask] = np.arange(n_tasks, dtype=np.int32)
        d_task_offsets   = cl.Buffer(self.ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=h_task_offsets)
        d_tasks_out      = cl.Buffer(self.ctx, mf.READ_WRITE, n_tasks * 32) # TaskData size is 32 bytes
        d_task_atoms_out = cl.Buffer(self.ctx, mf.READ_WRITE, n_tasks * nMaxAtom * 4)

        T0 = time.perf_counter_ns()
        self.prg.compact_tasks(
            self.queue, (int(n_blocks_xyz[0]), int(n_blocks_xyz[1]), int(n_blocks_xyz[2])), None,
            np.int32(n_blocks_xyz[0]), np.int32(n_blocks_xyz[1]), np.int32(n_blocks_xyz[2]),
            d_block_counts, d_task_offsets, d_task_atoms_raw,
            d_tasks_out, d_task_atoms_out, np.int32(nMaxAtom)
        )
        # 5. Read back results
        tasks_np      = np.empty(n_tasks, dtype=self.task_dtype_np)
        task_atoms_np = np.empty((n_tasks, nMaxAtom), dtype=np.int32)
        cl.enqueue_copy(self.queue, tasks_np,      d_tasks_out     )
        cl.enqueue_copy(self.queue, task_atoms_np, d_task_atoms_out)
        self.queue.finish()
        T1 = time.perf_counter_ns()
        print(f"[TIME] compact_tasks + readback {(T1-T0)*1e-6:.3f} [ms]")

        # Optional: sorting by na (descending) on host
        idx = np.argsort(tasks_np['na'])[::-1]
        tasks_np = tasks_np[idx]
        task_atoms_np = task_atoms_np[idx]
        
        return tasks_np, task_atoms_np

    def build_tasks(self, atoms, grid_spec, block_res=8, nMaxAtom=64):
        """
        Partition the grid into tasks (active blocks).
        """
        nx, ny, nz = grid_spec['ngrid'][:3]
        n_blocks = (nx // block_res, ny // block_res, nz // block_res)
        
        tasks = []
        atom_pos = atoms['pos']
        atom_Rcut = atoms['Rcut']
        natoms = len(atom_pos)
        
        origin = np.array(grid_spec['origin'][:3])
        dA = np.array(grid_spec['dA'][:3])
        dB = np.array(grid_spec['dB'][:3])
        dC = np.array(grid_spec['dC'][:3])

        block_counts = []
        max_count = 0
        empty_blocks = 0
        one_blocks = 0

        for fix in range(n_blocks[0]):
            for fiy in range(n_blocks[1]):
                for fiz in range(n_blocks[2]):
                    block_min = origin    + np.array([fix*block_res*dA[0], fiy*block_res*dB[1], fiz*block_res*dC[2]])
                    block_max = block_min + np.array([block_res*dA[0], block_res*dB[1], block_res*dC[2]])

                    atoms_in_block = []
                    for ia in range(natoms):
                        if self.check_overlap_sphere_aabb(atom_pos[ia], atom_Rcut[ia], block_min, block_max):
                            atoms_in_block.append(ia)

                    block_counts.append(len(atoms_in_block))
                    if len(atoms_in_block) == 0:
                        empty_blocks += 1
                        continue
                    if len(atoms_in_block) == 1:
                        one_blocks += 1
                    if len(atoms_in_block) > max_count:
                        max_count = len(atoms_in_block)
                    if len(atoms_in_block) > nMaxAtom:
                        raise RuntimeError(f"Block ({fix},{fiy},{fiz}) has {len(atoms_in_block)} atoms > nMaxAtom={nMaxAtom}")

                    # We want ONE task per voxel block to avoid atomic adds.
                    # We assume up to nMaxAtom (64) fits.
                    tasks.append({
                        'block_idx': (fix, fiy, fiz),
                        'na': min(len(atoms_in_block), nMaxAtom),
                        'nj': -1,
                        'atoms': atoms_in_block[:nMaxAtom]
                    })

        # Sort tasks by workload (na)
        tasks.sort(key=lambda x: x['na'], reverse=True)

        multi_blocks = len(block_counts) - empty_blocks - one_blocks
        if self.verbosity > 0: print(f"[DEBUG] block atom stats: na_max={max_count}, nbloks: empty={empty_blocks}, one={one_blocks}, multi={multi_blocks}")
        self.last_block_atom_counts = np.array(block_counts, dtype=np.int32)

        tasks_np = np.zeros(len(tasks), dtype=self.task_dtype_np)
        
        task_atoms_np = np.zeros((len(tasks), nMaxAtom), dtype=np.int32)
        
        for i, t in enumerate(tasks):
            tasks_np[i]['x'], tasks_np[i]['y'], tasks_np[i]['z'] = t['block_idx']
            tasks_np[i]['na'] = t['na']
            tasks_np[i]['nj'] = t['nj']
            task_atoms_np[i, :t['na']] = t['atoms']

        if self.verbosity > 0: print(f"[DEBUG] build_tasks finished: n_tasks={len(tasks)}")
        return tasks_np, task_atoms_np

    def grid_to_np(self, grid_spec):
        """Convert grid spec dictionary to numpy struct for GPU."""
        grid_spec_np = np.zeros(1, dtype=[
            ('origin', 'f4', 4),
            ('dA', 'f4', 4),
            ('dB', 'f4', 4),
            ('dC', 'f4', 4),
            ('ngrid', 'i4', 4)
        ])
        grid_spec_np[0]['origin'][:3] = grid_spec['origin']
        grid_spec_np[0]['dA'][:3] = grid_spec['dA']
        grid_spec_np[0]['dB'][:3] = grid_spec['dB']
        grid_spec_np[0]['dC'][:3] = grid_spec['dC']
        grid_spec_np[0]['ngrid'][:3] = grid_spec['ngrid']
        
        # DEBUG: print grid_spec_np values
        if self.verbosity > 0: print(f"[DEBUG] grid_spec_np: origin={grid_spec_np[0]['origin']} dA={grid_spec_np[0]['dA']} dB={grid_spec_np[0]['dB']} dC={grid_spec_np[0]['dC']} ngrid={grid_spec_np[0]['ngrid']}")
        
        return grid_spec_np

    def project(self, rho, neighs, atoms, grid_spec, tasks=None, nMaxAtom=64, use_gpu_tasks=False, use_tiled=True):
        """
        Main entry point for density projection using the tiled kernel.
        """
        if tasks is None:
            T0 = time.perf_counter_ns()
            if use_gpu_tasks:
                tasks_np, task_atoms_np = self.build_tasks_gpu(atoms, grid_spec, nMaxAtom=nMaxAtom)
            else:
                tasks_np, task_atoms_np = self.build_tasks(atoms, grid_spec, nMaxAtom=nMaxAtom)
            T1 = time.perf_counter_ns()
            if self.verbosity > 0: print(f"[TIME] build_tasks finished in {(T1-T0)*1e-6:.3f} [ms]")
        else:
            tasks_np, task_atoms_np = tasks

        n_tasks = len(tasks_np)
        ngrid_in = grid_spec['ngrid']
        if self.verbosity > 0: print(f"[DEBUG] grid_spec['ngrid']={ngrid_in} type={type(ngrid_in)}")
        nx, ny, nz = [int(x) for x in ngrid_in[:3]]
        if self.verbosity > 0: print(f"[DEBUG] derived grid dims nx,ny,nz=({nx},{ny},{nz})")

        # Prepare other buffers
        natoms = len(atoms['pos'])

        # DEBUG/ASSERT: validate task_atoms indices for active entries
        if n_tasks > 0:
            na_arr = tasks_np['na'].astype(np.int32)
            bad = []
            for it in range(n_tasks):
                na = int(na_arr[it])
                if na <= 0: continue
                idxs = task_atoms_np[it, :na]
                if (idxs < 0).any() or (idxs >= natoms).any():
                    bad.append((it, na, int(idxs.min()), int(idxs.max())))
                    if len(bad) >= 5:
                        break
            if bad:
                raise RuntimeError(f"GridProjector.project(): invalid atom index in task_atoms for tasks={bad} natoms={natoms}")

        atom_data = np.zeros(natoms, dtype=[
            ('pos_rcut', 'f4', 4),
            ('type', 'i4'),
            ('i0orb', 'i4'),
            ('norb', 'i4'),
            ('pad', 'i4')
        ])
        if (not hasattr(self, 'basis_meta')) or ('nz_map' not in self.basis_meta):
            raise RuntimeError('GridProjector.project(): basis_meta.nz_map missing; call load_basis(species_nz) before project().')
        if ('n_species' not in self.basis_meta) or ('max_shells' not in self.basis_meta) or ('n_nodes' not in self.basis_meta):
            raise RuntimeError(f"GridProjector.project(): basis_meta incomplete keys={list(self.basis_meta.keys())}")
        for i in range(natoms):
            atom_data[i]['pos_rcut'][:3] = atoms['pos'][i]
            atom_data[i]['pos_rcut'][3]  = atoms['Rcut'][i]
            # IMPORTANT: kernel expects a compact species index into packed basis_data, not atomic Z
            Z = int(atoms['type'][i])
            try:
                atom_data[i]['type'] = int(self.basis_meta['nz_map'][Z])
            except Exception as e:
                raise RuntimeError(f"GridProjector.project(): species nz={Z} not in loaded basis nz_map keys={list(self.basis_meta['nz_map'].keys())}") from e
            atom_data[i]['norb'] = 4 # Default for C, H with s,p
            atom_data[i]['i0orb'] = 0

        # DEBUG/ASSERT: mapped species indices must be in-range for packed basis buffer
        it_min = int(atom_data['type'].min()) if natoms > 0 else -1
        it_max = int(atom_data['type'].max()) if natoms > 0 else -1
        if self.verbosity > 0: print(f"[DEBUG] basis_meta: n_species={self.basis_meta['n_species']} max_shells={self.basis_meta['max_shells']} n_nodes={self.basis_meta['n_nodes']} dr={self.basis_meta['dr']:.6f}")
        if self.verbosity > 0: print(f"[DEBUG] atom_data.type range=[{it_min},{it_max}] unique={sorted(set(atom_data['type'].tolist()))}")
        if it_min < 0 or it_max >= int(self.basis_meta['n_species']):
            raise RuntimeError(f"GridProjector.project(): atom_data.type out of range [0,{self.basis_meta['n_species']-1}] got range=[{it_min},{it_max}]")

        # 2. Buffers
        mf = cl.mem_flags
        
        # DEBUG: check tasks_np size and dtype
        if self.verbosity > 0: print(f"[DEBUG] tasks_np: len={len(tasks_np)} itemsize={tasks_np.dtype.itemsize} nbytes={tasks_np.nbytes}")
        
        d_grid = cl.Buffer(self.ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=self.grid_to_np(grid_spec))
        
        if len(tasks_np) > 0:
            d_tasks = cl.Buffer(self.ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=tasks_np)
        else:
            # Fallback for empty buffer to avoid INVALID_BUFFER_SIZE
            d_tasks = cl.Buffer(self.ctx, mf.READ_ONLY, size=32) 
            
        d_atoms = cl.Buffer(self.ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=atom_data)
        if len(task_atoms_np) > 0:
            d_task_atoms = cl.Buffer(self.ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=task_atoms_np)
        else:
            d_task_atoms = cl.Buffer(self.ctx, mf.READ_ONLY, size=nMaxAtom * 4)
        
        rho32 = rho.astype(np.float32)
        d_rho = cl.Buffer(self.ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=rho32)
        d_neigh_j = cl.Buffer(self.ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=neighs.neigh_j.astype(np.int32))
        
        # species_info placeholder
        species_info = np.zeros((10, 4), dtype=np.int32)
        d_species_info = cl.Buffer(self.ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=species_info)

        out_nbytes = int(nx) * int(ny) * int(nz) * 4
        if self.verbosity > 0: print(f"[DEBUG] allocating d_out: nx,ny,nz=({nx},{ny},{nz}) out_nbytes={out_nbytes}")
        d_out = cl.Buffer(self.ctx, mf.WRITE_ONLY, out_nbytes)
        cl.enqueue_fill_buffer(self.queue, d_out, np.float32(0), 0, out_nbytes)

        # 3. Kernel launch
        ls = (32,)  # local size
        gs = (n_tasks * ls[0],)
        
        # d_basis placeholder
        if not hasattr(self, 'd_basis'):
             self.d_basis = cl.Buffer(self.ctx, mf.READ_ONLY, size=4)
             self.basis_meta = {'n_nodes': 0, 'dr': 0.0, 'max_shells': 0}

        if use_tiled:
            if self.verbosity > 0: print(f"[DEBUG] project_tiled: gs={gs}, ls={ls}, n_tasks={n_tasks}")
        else:
            if self.verbosity > 0: print(f"[DEBUG] project (non-tiled): gs={gs}, ls={ls}, n_tasks={n_tasks}")

        T0_ns = time.perf_counter_ns()
        if use_tiled:
            self.prg.project_density_sparse_tiled(
                self.queue, gs, ls,
                d_grid,
                np.int32(n_tasks),
                d_tasks, d_atoms, d_task_atoms,
                d_rho, 
                d_neigh_j,
                self.d_basis,
                d_species_info,
                np.int32(self.basis_meta['n_nodes']),
                np.float32(self.basis_meta['dr']),
                np.int32(self.basis_meta['max_shells']),
                np.int32(rho.shape[1]), # neigh_max
                np.int32(rho.shape[2]), # numorb_max
                np.int32(nMaxAtom),
                d_out
            )
        else:
            self.prg.project_density_sparse(
                self.queue, gs, ls,
                d_grid,
                np.int32(n_tasks),
                d_tasks, d_atoms, d_task_atoms,
                d_rho, 
                d_neigh_j,
                self.d_basis,
                d_species_info,
                np.int32(self.basis_meta['n_nodes']),
                np.float32(self.basis_meta['dr']),
                np.int32(self.basis_meta['max_shells']),
                np.int32(rho.shape[1]), # neigh_max
                np.int32(rho.shape[2]), # numorb_max
                np.int32(nMaxAtom),
                d_out
            )
        self.queue.finish()
        dt_ns = time.perf_counter_ns() - T0_ns
        if self.verbosity > 0: print(f"[TIME] project_tiled finished in {dt_ns*1e-6:.9f} [ms]")

        if self.verbosity > 0: print(f"[DEBUG] allocating host res: shape=({nx},{ny},{nz}) nbytes={int(nx)*int(ny)*int(nz)*4}")
        res = np.empty((int(nx), int(ny), int(nz)), dtype=np.float32)
        cl.enqueue_copy(self.queue, res, d_out)
        self.queue.finish()

        return res

    def project_orbital(self, coeffs, norb_per, atoms_dict, grid_spec, nMaxAtom=64, _debug_Fortran_order=False):
        """
        Project a single molecular orbital onto a 3D grid using the orbital projection kernel.

        Computes ψ(r) = Σ_i C_i φ_i(r) (signed wavefunction, not density)

        Args:
            coeffs: (natoms, numorb_max) MO coefficients. If _debug_Fortran_order=True, expects
                    Fortran order [s, py, pz, px]. Otherwise expects OpenCL order [s, px, py, pz].
            norb_per: (natoms,) number of orbitals per atom
            atoms_dict: dict with 'pos', 'Rcut', 'type'
            grid_spec: dict with 'origin', 'dA', 'dB', 'dC', 'ngrid'
            nMaxAtom: max atoms per task
            _debug_Fortran_order: If True, coeffs are in Fortran order and will be remapped

        Returns:
            psi: (nx, ny, nz) signed wavefunction
        """
        import numpy as np
        import time

        # Build tasks
        tasks_np, task_atoms_np = self.build_tasks(atoms_dict, grid_spec, nMaxAtom=64, block_res=8)
        if self.verbosity > 0: print(f"[DEBUG] project_orbital: n_tasks={len(tasks_np)}")

        # Prepare coefficient buffer with remapping from Fortran to OpenCL order
        natoms = len(atoms_dict['pos'])
        numorb_max = max(norb_per)
        coeffs_flat = np.zeros(natoms * numorb_max, dtype=np.float32)
        
        # Remapping: Fortran [s, py, pz, px] -> OpenCL [s, px, py, pz]
        _ORT_SPP_TO_STD = np.array([0, 3, 1, 2], dtype=np.int32)
        
        for ia in range(natoms):
            i0 = ia * numorb_max
            no = norb_per[ia]
            if _debug_Fortran_order and no == 4:
                # Remap p-orbitals from Fortran to OpenCL order
                coeffs_flat[i0:i0+no] = coeffs[ia, :no][_ORT_SPP_TO_STD]
            else:
                coeffs_flat[i0:i0+no] = coeffs[ia, :no]

        # Prepare atom data
        atom_data = np.zeros(natoms, dtype=[
            ('pos_rcut', 'f4', 4),
            ('type', 'i4'),
            ('i0orb', 'i4'),
            ('norb', 'i4'),
            ('pad', 'i4')
        ])
        for ia in range(natoms):
            atom_data[ia]['pos_rcut'][:3] = atoms_dict['pos'][ia]
            atom_data[ia]['pos_rcut'][3] = atoms_dict['Rcut'][ia]
            Z = int(atoms_dict['type'][ia])
            try:
                atom_data[ia]['type'] = int(self.basis_meta['nz_map'][Z])
            except Exception as e:
                raise RuntimeError(f"GridProjector.project_orbital(): species nz={Z} not in loaded basis nz_map keys={list(self.basis_meta['nz_map'].keys())}") from e
            atom_data[ia]['norb'] = norb_per[ia]
            atom_data[ia]['i0orb'] = ia * numorb_max

        # Grid spec
        d_grid = cl.Buffer(self.ctx, cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR,
                         hostbuf=self.grid_to_np(grid_spec))

        # Tasks
        if len(tasks_np) > 0:
            d_tasks = cl.Buffer(self.ctx, cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR,
                             hostbuf=tasks_np)
        else:
            d_tasks = cl.Buffer(self.ctx, cl.mem_flags.READ_ONLY, size=32)

        d_atoms = cl.Buffer(self.ctx, cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR,
                          hostbuf=atom_data)

        if len(task_atoms_np) > 0:
            d_task_atoms = cl.Buffer(self.ctx, cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR,
                                    hostbuf=task_atoms_np)
        else:
            d_task_atoms = cl.Buffer(self.ctx, cl.mem_flags.READ_ONLY, size=nMaxAtom * 4)

        d_coeffs = cl.Buffer(self.ctx, cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR,
                           hostbuf=coeffs_flat.astype(np.float32))

        # Output grid
        nx, ny, nz = grid_spec['ngrid'][:3]
        out_nbytes = int(nx) * int(ny) * int(nz) * 4
        d_out = cl.Buffer(self.ctx, cl.mem_flags.WRITE_ONLY, out_nbytes)
        cl.enqueue_fill_buffer(self.queue, d_out, np.float32(0), 0, out_nbytes)

        # Load kernel
        self._load_kernels()
        n_tasks = len(tasks_np)

        # Launch kernel
        ls = (32,)
        gs = (n_tasks * ls[0],)

        T0 = time.perf_counter_ns()
        self.prg.project_orbital(
            self.queue, gs, ls,
            d_grid, np.int32(n_tasks),
            d_tasks, d_atoms, d_task_atoms,
            d_coeffs,
            self.d_basis,
            np.int32(self.basis_meta['n_nodes']),
            np.float32(self.basis_meta['dr']),
            np.int32(self.basis_meta['max_shells']),
            np.int32(numorb_max),
            np.int32(nMaxAtom),
            d_out
        )
        self.queue.finish()
        T1 = time.perf_counter_ns()
        if self.verbosity > 0: print(f"[TIME] project_orbital {(T1-T0)*1e-6:.3f} [ms]")

        # Read back
        res = np.empty((int(nx), int(ny), int(nz)), dtype=np.float32)
        cl.enqueue_copy(self.queue, res, d_out)
        self.queue.finish()

        return res
