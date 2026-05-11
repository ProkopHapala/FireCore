import numpy as np
import pyopencl as cl
import pyopencl.cltypes
import os
import pyopencl.array as cl_array
import time
from ..OCL.OpenCLBase import OpenCLBase

class GridProjector(OpenCLBase):
    """
    Host class for projecting sparse density matrices to a real-space grid using OpenCL.
    """
    def __init__(self, fdata_dir, ctx=None, queue=None, nloc=32, debug_early_exit=False, debug_clear_only=False, debug_return0=False, debug_read_task=False, debug_read_grid=False, verbosity=0):
        super().__init__(nloc=nloc)
        self.fdata_dir = fdata_dir
        self.debug_early_exit = bool(debug_early_exit)
        self.debug_clear_only = bool(debug_clear_only)
        self.debug_return0 = bool(debug_return0)
        self.debug_read_task = bool(debug_read_task)
        self.debug_read_grid = bool(debug_read_grid)
        self.verbosity = int(verbosity)
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

    def load_basis_sto(self, species_list, dr=None, rc_max=None):
        """
        Load STO (Slater-type orbital) basis for DFTB+ into the same GPU buffer as load_basis().
        
        Evaluates STO analytically on a uniform grid: R_l(r) = sum_i c_i * r^(l+pow-1) * exp(-alpha_i*r)
        Packs identically to load_basis() as float2(value, d2) per node.
        Always uses max_shells=2 (s and p) so the kernel is uniform (H p-shell = zeros).
        
        Args:
            species_list: list of dicts with keys:
                'atomic_number': int
                'orbitals': list of {'l', 'exponents', 'coefficients', 'cutoff'}
                'resolution': float (Å, grid spacing hint)
            dr: override common grid spacing (Å); if None, uses min resolution/2
            rc_max: override max cutoff (Å); if None, uses max across all shells
        """
        from .DFTBplusParser import _spline_d2_uniform, compute_sto_radial

        all_nz = sorted(set(sp['atomic_number'] for sp in species_list))
        max_shells = 2  # Fixed: always s (shell 0) and p (shell 1), pad missing with zeros

        # Determine common grid
        all_rc = []
        all_res = []
        for sp in species_list:
            for orb in sp['orbitals']:
                all_rc.append(orb['cutoff'])
            all_res.append(sp.get('resolution', 0.15))
        if rc_max is None:
            rc_max = max(all_rc)
        if dr is None:
            dr = min(all_res) / 2.0

        n_nodes = int(np.ceil(rc_max / dr)) + 2

        if self.verbosity > 0:
            print(f"[GridProjector STO] Common grid: dr={dr:.6f} Å, rc_max={rc_max:.3f} Å, n_nodes={n_nodes}")

        packed_basis = np.zeros((len(all_nz), max_shells, n_nodes, 2), dtype=np.float32)

        # Map atomic_number -> index in all_nz
        nz_map = {nz: i for i, nz in enumerate(all_nz)}

        # Build species dict by atomic_number for easy lookup
        sp_by_nz = {sp['atomic_number']: sp for sp in species_list}

        r = np.arange(n_nodes) * dr

        for nz in all_nz:
            i_spec = nz_map[nz]
            sp = sp_by_nz[nz]
            # Sort orbitals by l so shell 0 = s, shell 1 = p
            orbs_by_l = {}
            for orb in sp['orbitals']:
                l = orb['l']
                if l not in orbs_by_l:
                    orbs_by_l[l] = orb

            for ish in range(max_shells):
                l = ish  # shell index = angular momentum (0=s, 1=p)
                if l not in orbs_by_l:
                    # No orbital for this l (e.g. H has no p) — leave as zeros
                    continue
                orb = orbs_by_l[l]
                aa = np.asarray(orb['coefficients'], dtype=np.float64)
                alpha = np.asarray(orb['exponents'], dtype=np.float64)

                # Evaluate STO on common grid
                vals = compute_sto_radial(r, aa, alpha, l).astype(np.float32)
                d2 = _spline_d2_uniform(vals.astype(np.float64), dr).astype(np.float32)

                packed_basis[i_spec, ish, :, 0] = vals
                packed_basis[i_spec, ish, :, 1] = d2

                if self.verbosity > 0:
                    nAlpha = len(alpha)
                    nPow = aa.shape[0] if aa.ndim > 1 else 1
                    print(f"[GridProjector STO]   Z={nz} shell {ish} (l={l}): nAlpha={nAlpha}, nPow={nPow}, cutoff={orb['cutoff']:.2f}")

        self.d_basis = cl.Buffer(self.ctx, cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR, hostbuf=packed_basis)
        self.basis_meta = {
            'n_species': len(all_nz),
            'max_shells': max_shells,
            'n_nodes': n_nodes,
            'dr': dr,
            'nz_map': nz_map,
        }
        self.species_nz = all_nz
        return packed_basis

    def _load_kernels(self):
        cl_path = os.path.join(os.path.dirname(__file__), "cl/Grid_dftb.cl")
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

        print(f"[GridProjector] Loading program from {cl_path}...")
        self.load_program(kernel_path=cl_path, build_options=build_opts if len(build_opts)>0 else None)
        print(f"[GridProjector] Program load complete")

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
        n_blocks = (
            (int(nx) + block_res - 1) // block_res,
            (int(ny) + block_res - 1) // block_res,
            (int(nz) + block_res - 1) // block_res,
        )
        
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

    def project_density(self, rho, neighs, atoms, grid_spec, tasks=None, nMaxAtom=64, use_gpu_tasks=False, use_tiled=True):
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


    def project_orbital_points(self, points, coeffs, norb_per, atoms_dict, _debug_Fortran_order=False):
        """Evaluate a single orbital at arbitrary points (debugging parity with Fortran orb2points).

        This avoids any grid sampling / slicing ambiguity.

        Args:
            points: (n_points,3) float32/float64 positions in Angstrom
            coeffs: (natoms,4) coefficients. Default expects [px,py,pz,s].
                    If _debug_Fortran_order=True and atom has 4 orbitals, expects [s,py,pz,px].
            norb_per: (natoms,) number of orbitals per atom (1 or 4 for H/O in H2O)
            atoms_dict: dict with 'pos','Rcut','type'
        Returns:
            psi: (n_points,) float32
        """
        import numpy as np
        points = np.asarray(points, dtype=np.float32)
        if points.ndim != 2 or points.shape[1] != 3:
            raise ValueError(f"project_orbital_points: points must be (n,3), got {points.shape}")
        natoms = len(atoms_dict['pos'])

        # Pack coefficients exactly like project_orbital()
        numorb_max = 4
        coeffs_flat = np.zeros(natoms * numorb_max, dtype=np.float32)
        _ORT_SPP_TO_OCL = np.array([3, 1, 2, 0], dtype=np.int32)  # [s,py,pz,px] -> [px,py,pz,s]
        for ia in range(natoms):
            i0 = ia * numorb_max
            no = int(norb_per[ia])
            if _debug_Fortran_order and (no == 4):
                coeffs_flat[i0:i0+4] = coeffs[ia, :4][_ORT_SPP_TO_OCL]
            else:
                coeffs_flat[i0:i0+4] = coeffs[ia, :4]

        # AtomData
        atom_data = np.zeros(natoms, dtype=[
            ('pos_rcut', 'f4', 4),
            ('type', 'i4'),
            ('i0orb', 'i4'),
            ('norb', 'i4'),
            ('pad', 'i4')
        ])
        for ia in range(natoms):
            atom_data[ia]['pos_rcut'][:3] = atoms_dict['pos'][ia]
            atom_data[ia]['pos_rcut'][3]  = atoms_dict['Rcut'][ia]
            Z = int(atoms_dict['type'][ia])
            if Z not in self.basis_meta['nz_map']:
                raise RuntimeError(f"project_orbital_points: species nz={Z} not loaded; loaded={list(self.basis_meta['nz_map'].keys())}")
            atom_data[ia]['type'] = int(self.basis_meta['nz_map'][Z])
            atom_data[ia]['norb'] = int(norb_per[ia])
            atom_data[ia]['i0orb'] = ia * numorb_max

        # Buffers
        mf = cl.mem_flags
        d_points = cl.Buffer(self.ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=np.c_[points, np.zeros((len(points),1),np.float32)].astype(np.float32))
        d_atoms  = cl.Buffer(self.ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=atom_data)
        d_coeffs = cl.Buffer(self.ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=coeffs_flat)
        d_out    = cl.Buffer(self.ctx, mf.WRITE_ONLY, size=len(points)*4)

        # Build program (ensures new kernel is available)
        self._load_kernels()

        # Launch
        gs = (int(len(points)),)
        ls = None
        self.prg.project_orbital_points(
            self.queue, gs, ls,
            np.int32(len(points)),
            d_points,
            d_atoms,
            np.int32(natoms),
            d_coeffs,
            self.d_basis,
            np.int32(self.basis_meta['n_nodes']),
            np.float32(self.basis_meta['dr']),
            np.int32(self.basis_meta['max_shells']),
            d_out
        )
        self.queue.finish()

        out = np.empty(len(points), dtype=np.float32)
        cl.enqueue_copy(self.queue, out, d_out)
        self.queue.finish()
        return out

    def prepare_orbital_points_projection(self, points, atoms_dict, nMaxAtom=64):
        """
        One-time setup for batched orbital projection at arbitrary points.
        Pre-uploads all static GPU buffers (points, atoms).
        """
        mf = cl.mem_flags
        t0 = time.perf_counter_ns()

        points = np.asarray(points, dtype=np.float32)
        n_points = len(points)
        natoms   = len(atoms_dict['pos'])
        numorb_max = 4

        atom_data = np.zeros(natoms, dtype=[
            ('pos_rcut', 'f4', 4), ('type', 'i4'), ('i0orb', 'i4'), ('norb', 'i4'), ('pad', 'i4')
        ])
        for ia in range(natoms):
            atom_data[ia]['pos_rcut'][:3] = atoms_dict['pos'][ia]
            atom_data[ia]['pos_rcut'][3]  = atoms_dict['Rcut'][ia]
            Z = int(atoms_dict['type'][ia])
            atom_data[ia]['type']  = int(self.basis_meta['nz_map'][Z])
            atom_data[ia]['norb']  = int(4) # numorb_max
            atom_data[ia]['i0orb'] = ia * numorb_max

        # Upload static buffers
        # points are expanded to float4 for alignment
        points_f4 = np.zeros((n_points, 4), dtype=np.float32)
        points_f4[:, :3] = points
        d_points = cl.Buffer(self.ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=points_f4)
        d_atoms  = cl.Buffer(self.ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=atom_data)

        # Mutable/Reused buffers
        d_coeffs = cl.Buffer(self.ctx, mf.READ_ONLY,  size=natoms * numorb_max * 4)
        d_out    = cl.Buffer(self.ctx, mf.WRITE_ONLY, size=n_points * 4)

        self._load_kernels()
        t1 = time.perf_counter_ns()
        print(f"[prepare_orbital_points_projection] Setup: {(t1-t0)*1e-6:.1f} ms, n_points={n_points}")
        return dict(
            n_points=n_points, natoms=natoms, numorb_max=numorb_max,
            d_points=d_points, d_atoms=d_atoms, d_coeffs=d_coeffs, d_out=d_out
        )

    def project_orbital_points_prepped(self, coeffs_flat, ctx):
        """
        Project orbital at pre-loaded points.
        """
        cl.enqueue_copy(self.queue, ctx['d_coeffs'], coeffs_flat.astype(np.float32))
        
        gs = (ctx['n_points'],)
        ls = None
        self.prg.project_orbital_points(
            self.queue, gs, ls,
            np.int32(ctx['n_points']),
            ctx['d_points'],
            ctx['d_atoms'],
            np.int32(ctx['natoms']),
            ctx['d_coeffs'],
            self.d_basis,
            np.int32(self.basis_meta['n_nodes']),
            np.float32(self.basis_meta['dr']),
            np.int32(self.basis_meta['max_shells']),
            ctx['d_out']
        )
        self.queue.finish()
        out = np.empty(ctx['n_points'], dtype=np.float32)
        cl.enqueue_copy(self.queue, out, ctx['d_out'])
        self.queue.finish()
        return out




    def mo_overlap_points_exp_sk(
            self,
            tip_centers,
            tip_pos_rel,
            smp_pos,
            coeffs_tip,
            coeffs_smp,
            tip_quat=None,
            beta=1.0,
            r0=3.0,
            rcut=8.0,
        ):
        """GPU orbital-overlap scan map for molecular tip vs sample using exp+SK.

        Each work-item computes one scan pixel corresponding to one tip-center position.

        Args:
            tip_centers: (npts,3) float32/64 tip center positions
            tip_pos_rel: (ntip_atoms,3) float32/64 tip atom positions relative to tip center
            smp_pos:     (nsmp_atoms,3) float32/64 sample atom positions
            coeffs_tip:  (ntip_atoms,4) float32 coeffs [px,py,pz,s]
            coeffs_smp:  (nsmp_atoms,4) float32 coeffs [px,py,pz,s]
        Returns:
            t: (npts,) float32 signed amplitude
            I: (npts,) float32 intensity t^2
        """
        import numpy as np
        import pyopencl as cl

        tip_centers = np.asarray(tip_centers, dtype=np.float32)
        if tip_quat is None:
            tip_quat = np.zeros((len(tip_centers), 4), dtype=np.float32)
            tip_quat[:, 3] = 1.0
        tip_quat = np.asarray(tip_quat, dtype=np.float32)
        assert tip_quat.shape[0] == tip_centers.shape[0]
        assert tip_quat.shape[1] == 4
        tip_pos_rel = np.asarray(tip_pos_rel, dtype=np.float32)
        smp_pos = np.asarray(smp_pos, dtype=np.float32)
        coeffs_tip = np.asarray(coeffs_tip, dtype=np.float32)
        coeffs_smp = np.asarray(coeffs_smp, dtype=np.float32)

        if tip_centers.ndim != 2 or tip_centers.shape[1] != 3:
            raise ValueError(f"mo_overlap_points_exp_sk: tip_centers must be (n,3), got {tip_centers.shape}")
        if tip_pos_rel.ndim != 2 or tip_pos_rel.shape[1] != 3:
            raise ValueError(f"mo_overlap_points_exp_sk: tip_pos_rel must be (n,3), got {tip_pos_rel.shape}")
        if smp_pos.ndim != 2 or smp_pos.shape[1] != 3:
            raise ValueError(f"mo_overlap_points_exp_sk: smp_pos must be (n,3), got {smp_pos.shape}")
        if coeffs_tip.shape != (len(tip_pos_rel), 4):
            raise ValueError(f"mo_overlap_points_exp_sk: coeffs_tip must be (ntip,4), got {coeffs_tip.shape}")
        if coeffs_smp.shape != (len(smp_pos), 4):
            raise ValueError(f"mo_overlap_points_exp_sk: coeffs_smp must be (nsmp,4), got {coeffs_smp.shape}")

        # Pack to float4 buffers
        tip_centers4 = np.c_[tip_centers, np.zeros((len(tip_centers), 1), dtype=np.float32)].astype(np.float32)
        tip_pos_rel4 = np.c_[tip_pos_rel, np.zeros((len(tip_pos_rel), 1), dtype=np.float32)].astype(np.float32)
        smp_pos4 = np.c_[smp_pos, np.zeros((len(smp_pos), 1), dtype=np.float32)].astype(np.float32)

        mf = cl.mem_flags
        d_tip_centers = cl.Buffer(self.ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=tip_centers4)
        d_tip_quat    = cl.Buffer(self.ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=tip_quat)
        d_tip_pos_rel = cl.Buffer(self.ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=tip_pos_rel4)
        d_smp_pos     = cl.Buffer(self.ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=smp_pos4)
        d_ct = cl.Buffer(self.ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=coeffs_tip.astype(np.float32))
        d_cs = cl.Buffer(self.ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=coeffs_smp.astype(np.float32))
        d_out_t = cl.Buffer(self.ctx, mf.WRITE_ONLY, size=len(tip_centers) * 4)
        d_out_I = cl.Buffer(self.ctx, mf.WRITE_ONLY, size=len(tip_centers) * 4)

        self._load_kernels()

        gs = (int(len(tip_centers4)),)
        ls = None
        self.prg.mo_overlap_points_exp_sk(
            self.queue, gs, ls,
            np.int32(len(tip_centers)),
            d_tip_centers,
            d_tip_quat,
            d_tip_pos_rel,
            d_smp_pos,
            np.int32(len(tip_pos_rel)),
            np.int32(len(smp_pos)),
            d_ct,
            d_cs,
            np.float32(beta),
            np.float32(r0),
            np.float32(rcut),
            d_out_t,
            d_out_I
        )
        self.queue.finish()

        out_t = np.empty(len(tip_centers), dtype=np.float32)
        out_I = np.empty(len(tip_centers), dtype=np.float32)
        cl.enqueue_copy(self.queue, out_t, d_out_t)
        cl.enqueue_copy(self.queue, out_I, d_out_I)
        self.queue.finish()
        return out_t, out_I


    def stm_dyson_wg_scan(
            self,
            tip_centers,
            tip_pos_rel,
            smp_pos,
            GT_global=None,
            GS_global=None,
            uT_source=None,
            beta=1.0,
            r0=3.0,
            rcut=8.0,
            local_size=32,
        ):
        import numpy as np
        import pyopencl as cl

        tip_centers = np.asarray(tip_centers, dtype=np.float32)
        tip_pos_rel = np.asarray(tip_pos_rel, dtype=np.float32)
        smp_pos     = np.asarray(smp_pos, dtype=np.float32)

        if tip_centers.ndim != 2 or tip_centers.shape[1] != 3:
            raise ValueError(f"stm_dyson_wg_scan: tip_centers must be (n,3), got {tip_centers.shape}")
        if tip_pos_rel.ndim != 2 or tip_pos_rel.shape[1] != 3:
            raise ValueError(f"stm_dyson_wg_scan: tip_pos_rel must be (n,3), got {tip_pos_rel.shape}")
        if smp_pos.ndim != 2 or smp_pos.shape[1] != 3:
            raise ValueError(f"stm_dyson_wg_scan: smp_pos must be (n,3), got {smp_pos.shape}")

        n_pixels   = int(tip_centers.shape[0])
        ntip_atoms = int(tip_pos_rel.shape[0])
        nsmp_atoms = int(smp_pos.shape[0])
        nt = 4 * ntip_atoms
        ns = 4 * nsmp_atoms

        if GT_global is None:
            GT_global = np.eye(nt, dtype=np.complex64)
        else:
            GT_global = np.asarray(GT_global)
        if GS_global is None:
            GS_global = np.eye(ns, dtype=np.complex64)
        else:
            GS_global = np.asarray(GS_global)

        if GT_global.shape != (nt, nt):
            raise ValueError(f"stm_dyson_wg_scan: GT_global must be ({nt},{nt}), got {GT_global.shape}")
        if GS_global.shape != (ns, ns):
            raise ValueError(f"stm_dyson_wg_scan: GS_global must be ({ns},{ns}), got {GS_global.shape}")

        if uT_source is None:
            uT_source = np.zeros(nt, dtype=np.complex64)
            uT_source[3] = 1.0 + 0.0j
        else:
            uT_source = np.asarray(uT_source)
        if uT_source.shape != (nt,):
            raise ValueError(f"stm_dyson_wg_scan: uT_source must be ({nt},), got {uT_source.shape}")

        tip_centers4 = np.c_[tip_centers, np.zeros((n_pixels, 1), dtype=np.float32)].astype(np.float32)
        tip_pos_rel4 = np.c_[tip_pos_rel, np.zeros((ntip_atoms, 1), dtype=np.float32)].astype(np.float32)
        smp_pos4     = np.c_[smp_pos,     np.zeros((nsmp_atoms, 1), dtype=np.float32)].astype(np.float32)

        # Pack complex64 -> float2
        GT_f2 = np.asarray(GT_global, dtype=np.complex64).view(np.float32).reshape(nt*nt, 2)
        GS_f2 = np.asarray(GS_global, dtype=np.complex64).view(np.float32).reshape(ns*ns, 2)
        uT_f2 = np.asarray(uT_source, dtype=np.complex64).view(np.float32).reshape(nt, 2)

        mf = cl.mem_flags
        d_tip_centers = cl.Buffer(self.ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=tip_centers4)
        d_tip_pos_rel = cl.Buffer(self.ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=tip_pos_rel4)
        d_smp_pos     = cl.Buffer(self.ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=smp_pos4)
        d_GT          = cl.Buffer(self.ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=GT_f2.astype(np.float32))
        d_GS          = cl.Buffer(self.ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=GS_f2.astype(np.float32))
        d_uT          = cl.Buffer(self.ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=uT_f2.astype(np.float32))
        d_out         = cl.Buffer(self.ctx, mf.WRITE_ONLY, size=n_pixels * 4)

        self._load_kernels()

        ls = (int(local_size),)
        gs = (int(n_pixels) * int(local_size),)
        self.prg.solve_stm_dyson_wg(
            self.queue, gs, ls,
            np.int32(n_pixels),
            d_tip_centers,
            d_tip_pos_rel,
            d_smp_pos,
            np.int32(ntip_atoms),
            np.int32(nsmp_atoms),
            d_GT,
            d_GS,
            d_uT,
            np.float32(beta),
            np.float32(r0),
            np.float32(rcut),
            d_out,
        )
        self.queue.finish()

        out = np.empty(n_pixels, dtype=np.float32)
        cl.enqueue_copy(self.queue, out, d_out)
        self.queue.finish()
        return out

    def stm_gf_dyson_2mol_mo_scan(
            self,
            tip_centers,
            tip_pos_rel,
            smp_pos,
            GT_global,
            GS_global,
            c_tip,
            c_smp,
            tip_norb_per,
            smp_norb_per,
            beta=1.0,
            r0=3.0,
            rcut=8.0,
        ):
        """GPU GF-Dyson MO scan for 2-molecule STM (work-item per pixel).

        Math: amp(p) = c_tip^H · GT · M_ts(p) · GS · c_smp

        Precomputes on CPU:
            v_S = GS @ c_smp   (remapped Fortran→OCL [px,py,pz,s] order)
            u_T = c_tip^H @ GT (remapped Fortran→OCL [px,py,pz,s] order)

        GPU kernel computes M_ts via simplified exponential SK hopping and
        accumulates amp = Σ_{it,is} u_T[it] * V_{it,is} * v_S[is].

        Args:
            tip_centers:  (n_pixels, 3) float32 tip center positions in Å
            tip_pos_rel:  (ntip_atoms, 3) float32 tip atom positions relative to center
            smp_pos:      (nsmp_atoms, 3) float32 sample atom positions in Å
            GT_global:    (tip_norb_fort, tip_norb_fort) complex Green's fcn (Fortran conv)
            GS_global:    (smp_norb_fort, smp_norb_fort) complex Green's fcn (Fortran conv)
            c_tip:        (tip_norb_fort,) complex MO coefficients (Fortran conv)
            c_smp:        (smp_norb_fort,) complex MO coefficients (Fortran conv)
            tip_norb_per: (ntip_atoms,) int32 orbitals per tip atom
            smp_norb_per: (nsmp_atoms,) int32 orbitals per sample atom
            beta, r0, rcut: exponential SK parameters

        Returns:
            out: (n_pixels,) float32 current intensity |amp|^2
        """
        import numpy as np
        import pyopencl as cl

        tip_centers = np.asarray(tip_centers, dtype=np.float32)
        tip_pos_rel = np.asarray(tip_pos_rel, dtype=np.float32)
        smp_pos     = np.asarray(smp_pos, dtype=np.float32)
        GT_global   = np.asarray(GT_global, dtype=np.complex128)
        GS_global   = np.asarray(GS_global, dtype=np.complex128)
        c_tip       = np.asarray(c_tip, dtype=np.complex128)
        c_smp       = np.asarray(c_smp, dtype=np.complex128)
        tip_norb_per = np.asarray(tip_norb_per, dtype=np.int32)
        smp_norb_per = np.asarray(smp_norb_per, dtype=np.int32)

        n_pixels    = int(tip_centers.shape[0])
        ntip_atoms  = int(tip_pos_rel.shape[0])
        nsmp_atoms  = int(smp_pos.shape[0])

        # --- Remap Fortran→OCL utility ---
        # Fortran per-atom: [s, py, pz, px] (Ortega)
        # OpenCL Grid:      [px, py, pz, s]  (Cartesian, padded to 4)
        _PERM_F2O = np.array([3, 1, 2, 0], dtype=np.int32)  # Fort[r] → OCL[perm[r]]

        def _remap_vec_fortran_to_ocl(v_fort, norb_per):
            natoms = len(norb_per)
            v_ocl = np.zeros(natoms * 4, dtype=np.complex128)
            starts = np.zeros(natoms + 1, dtype=np.int32)
            starts[1:] = np.cumsum(norb_per)
            for ia in range(natoms):
                no = int(norb_per[ia])
                i0f = int(starts[ia])
                i0o = ia * 4
                if no == 1:
                    v_ocl[i0o + 3] = v_fort[i0f]  # s → OCL slot 3
                elif no == 4:
                    for k in range(4):
                        v_ocl[i0o + k] = v_fort[i0f + int(_PERM_F2O[k])]
                else:
                    v_ocl[i0o:i0o + no] = v_fort[i0f:i0f + no]
            return v_ocl

        # Precompute vectors in Fortran convention
        v_S_fort = GS_global @ c_smp       # (smp_norb_fort,)
        u_T_fort = np.conj(c_tip) @ GT_global  # (tip_norb_fort,)  row vector result

        # Remap to OCL [px,py,pz,s] order
        v_S_ocl = _remap_vec_fortran_to_ocl(v_S_fort, smp_norb_per)
        u_T_ocl = _remap_vec_fortran_to_ocl(u_T_fort, tip_norb_per)

        # Build orb2atom in OpenCL convention (padded to 4 per atom)
        tip_norb_ocl = ntip_atoms * 4
        smp_norb_ocl = nsmp_atoms * 4
        tip_orb2atom_ocl = np.repeat(np.arange(ntip_atoms, dtype=np.int32), 4)
        smp_orb2atom_ocl = np.repeat(np.arange(nsmp_atoms, dtype=np.int32), 4)

        # Pack float4 buffers
        tip_centers4 = np.c_[tip_centers, np.zeros((n_pixels, 1), dtype=np.float32)].astype(np.float32)
        tip_pos_rel4 = np.c_[tip_pos_rel, np.zeros((ntip_atoms, 1), dtype=np.float32)].astype(np.float32)
        smp_pos4     = np.c_[smp_pos,     np.zeros((nsmp_atoms, 1), dtype=np.float32)].astype(np.float32)

        # Pack complex → float2
        uT_f2 = np.asarray(u_T_ocl, dtype=np.complex64).view(np.float32).reshape(tip_norb_ocl, 2)
        vS_f2 = np.asarray(v_S_ocl, dtype=np.complex64).view(np.float32).reshape(smp_norb_ocl, 2)

        mf = cl.mem_flags
        d_tip_centers = cl.Buffer(self.ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=tip_centers4)
        d_tip_pos_rel = cl.Buffer(self.ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=tip_pos_rel4)
        d_smp_pos     = cl.Buffer(self.ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=smp_pos4)
        d_tip_o2a     = cl.Buffer(self.ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=tip_orb2atom_ocl.astype(np.int32))
        d_smp_o2a     = cl.Buffer(self.ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=smp_orb2atom_ocl.astype(np.int32))
        d_uT = cl.Buffer(self.ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=uT_f2.astype(np.float32))
        d_vS = cl.Buffer(self.ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=vS_f2.astype(np.float32))
        d_out = cl.Buffer(self.ctx, mf.WRITE_ONLY, size=n_pixels * 4)

        self._load_kernels()

        gs = (int(n_pixels),)
        ls = None
        self.prg.stm_gf_dyson_2mol_mo_scan(
            self.queue, gs, ls,
            np.int32(n_pixels),
            d_tip_centers,
            d_tip_pos_rel,
            d_smp_pos,
            d_tip_o2a,
            d_smp_o2a,
            d_uT,
            d_vS,
            np.int32(ntip_atoms),
            np.int32(nsmp_atoms),
            np.int32(tip_norb_ocl),
            np.int32(smp_norb_ocl),
            np.float32(float(beta)),
            np.float32(float(r0)),
            np.float32(float(rcut)),
            d_out,
        )
        self.queue.finish()

        out = np.empty(n_pixels, dtype=np.float32)
        cl.enqueue_copy(self.queue, out, d_out)
        self.queue.finish()
        return out

    def mo_overlap_points_exp_sk_2mol(
            self,
            tip_centers,
            tip_pos_rel,
            smp_pos,
            coeffs_tip,
            coeffs_smp,
            tip_quat=None,
            beta=1.0,
            r0=3.0,
            rcut=8.0,
        ):
        """Same as mo_overlap_points_exp_sk, but calls an explicit two-molecule kernel entrypoint.

        This is intended for workflows where the tip and sample are different molecules.
        The math is identical; we keep a separate kernel name to avoid breaking existing
        call sites and make scripts self-documenting.
        """
        import numpy as np
        import pyopencl as cl

        tip_centers = np.asarray(tip_centers, dtype=np.float32)
        if tip_quat is None:
            tip_quat = np.zeros((len(tip_centers), 4), dtype=np.float32)
            tip_quat[:, 3] = 1.0
        tip_quat = np.asarray(tip_quat, dtype=np.float32)
        assert tip_quat.shape[0] == tip_centers.shape[0]
        assert tip_quat.shape[1] == 4
        tip_pos_rel = np.asarray(tip_pos_rel, dtype=np.float32)
        smp_pos = np.asarray(smp_pos, dtype=np.float32)
        coeffs_tip = np.asarray(coeffs_tip, dtype=np.float32)
        coeffs_smp = np.asarray(coeffs_smp, dtype=np.float32)

        if tip_centers.ndim != 2 or tip_centers.shape[1] != 3:
            raise ValueError(f"mo_overlap_points_exp_sk_2mol: tip_centers must be (n,3), got {tip_centers.shape}")
        if tip_pos_rel.ndim != 2 or tip_pos_rel.shape[1] != 3:
            raise ValueError(f"mo_overlap_points_exp_sk_2mol: tip_pos_rel must be (n,3), got {tip_pos_rel.shape}")
        if smp_pos.ndim != 2 or smp_pos.shape[1] != 3:
            raise ValueError(f"mo_overlap_points_exp_sk_2mol: smp_pos must be (n,3), got {smp_pos.shape}")
        if coeffs_tip.shape != (len(tip_pos_rel), 4):
            raise ValueError(f"mo_overlap_points_exp_sk_2mol: coeffs_tip must be (ntip,4), got {coeffs_tip.shape}")
        if coeffs_smp.shape != (len(smp_pos), 4):
            raise ValueError(f"mo_overlap_points_exp_sk_2mol: coeffs_smp must be (nsmp,4), got {coeffs_smp.shape}")

        tip_centers4 = np.c_[tip_centers, np.zeros((len(tip_centers), 1), dtype=np.float32)].astype(np.float32)
        tip_pos_rel4 = np.c_[tip_pos_rel, np.zeros((len(tip_pos_rel), 1), dtype=np.float32)].astype(np.float32)
        smp_pos4 = np.c_[smp_pos, np.zeros((len(smp_pos), 1), dtype=np.float32)].astype(np.float32)

        mf = cl.mem_flags
        d_tip_centers = cl.Buffer(self.ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=tip_centers4)
        d_tip_quat    = cl.Buffer(self.ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=tip_quat)
        d_tip_pos_rel = cl.Buffer(self.ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=tip_pos_rel4)
        d_smp_pos     = cl.Buffer(self.ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=smp_pos4)
        d_ct = cl.Buffer(self.ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=coeffs_tip.astype(np.float32))
        d_cs = cl.Buffer(self.ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=coeffs_smp.astype(np.float32))
        d_out_t = cl.Buffer(self.ctx, mf.WRITE_ONLY, size=len(tip_centers) * 4)
        d_out_I = cl.Buffer(self.ctx, mf.WRITE_ONLY, size=len(tip_centers) * 4)

        self._load_kernels()

        gs = (int(len(tip_centers4)),)
        ls = None
        self.prg.mo_overlap_points_exp_sk_2mol(
            self.queue, gs, ls,
            np.int32(len(tip_centers)),
            d_tip_centers,
            d_tip_quat,
            d_tip_pos_rel,
            d_smp_pos,
            np.int32(len(tip_pos_rel)),
            np.int32(len(smp_pos)),
            d_ct,
            d_cs,
            np.float32(beta),
            np.float32(r0),
            np.float32(rcut),
            d_out_t,
            d_out_I
        )
        self.queue.finish()

        out_t = np.empty(len(tip_centers), dtype=np.float32)
        out_I = np.empty(len(tip_centers), dtype=np.float32)
        cl.enqueue_copy(self.queue, out_t, d_out_t)
        cl.enqueue_copy(self.queue, out_I, d_out_I)
        self.queue.finish()
        return out_t, out_I

    def project_orbital_points_exp(self, points, coeffs, norb_per, atoms_dict, beta=1.0, r0=3.0, _debug_Fortran_order=False):
        """Evaluate a single orbital at arbitrary points using exponential radial decay.

        Uses OpenCL kernel `project_orbital_points_exp` from `cl/Grid.cl`.

        Args:
            points: (n_points,3) float32/float64 positions in Angstrom
            coeffs: (natoms,4) coefficients in [px,py,pz,s] order (or Fortran order if _debug_Fortran_order=True)
            norb_per: (natoms,) number of orbitals per atom
            atoms_dict: dict with 'pos','Rcut','type'
            beta, r0: exp(-beta*(r-r0)) parameters
        Returns:
            psi: (n_points,) float32
        """
        import numpy as np
        points = np.asarray(points, dtype=np.float32)
        if points.ndim != 2 or points.shape[1] != 3:
            raise ValueError(f"project_orbital_points_exp: points must be (n,3), got {points.shape}")
        natoms = len(atoms_dict['pos'])

        # Pack coefficients exactly like project_orbital_points
        numorb_max = 4
        coeffs_flat = np.zeros(natoms * numorb_max, dtype=np.float32)
        _ORT_SPP_TO_OCL = np.array([3, 1, 2, 0], dtype=np.int32)  # [s,py,pz,px] -> [px,py,pz,s]
        for ia in range(natoms):
            i0 = ia * numorb_max
            no = int(norb_per[ia])
            if _debug_Fortran_order and (no == 4):
                coeffs_flat[i0:i0+4] = coeffs[ia, :4][_ORT_SPP_TO_OCL]
            else:
                coeffs_flat[i0:i0+4] = coeffs[ia, :4]

        atom_data = np.zeros(natoms, dtype=[
            ('pos_rcut', 'f4', 4),
            ('type', 'i4'),
            ('i0orb', 'i4'),
            ('norb', 'i4'),
            ('pad', 'i4')
        ])
        for ia in range(natoms):
            atom_data[ia]['pos_rcut'][:3] = atoms_dict['pos'][ia]
            atom_data[ia]['pos_rcut'][3]  = atoms_dict['Rcut'][ia]
            Z = int(atoms_dict['type'][ia])
            if Z not in self.basis_meta['nz_map']:
                raise RuntimeError(f"project_orbital_points_exp: species nz={Z} not loaded; loaded={list(self.basis_meta['nz_map'].keys())}")
            atom_data[ia]['type'] = int(self.basis_meta['nz_map'][Z])
            atom_data[ia]['norb'] = int(norb_per[ia])
            atom_data[ia]['i0orb'] = ia * numorb_max

        mf = cl.mem_flags
        d_points = cl.Buffer(self.ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=np.c_[points, np.zeros((len(points),1),np.float32)].astype(np.float32))
        d_atoms  = cl.Buffer(self.ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=atom_data)
        d_coeffs = cl.Buffer(self.ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=coeffs_flat)
        d_out    = cl.Buffer(self.ctx, mf.WRITE_ONLY, size=len(points)*4)

        self._load_kernels()

        gs = (int(len(points)),)
        ls = None
        self.prg.project_orbital_points_exp(
            self.queue, gs, ls,
            np.int32(len(points)),
            d_points,
            d_atoms,
            np.int32(natoms),
            d_coeffs,
            np.float32(beta),
            np.float32(r0),
            d_out
        )
        self.queue.finish()

        out = np.empty(len(points), dtype=np.float32)
        cl.enqueue_copy(self.queue, out, d_out)
        self.queue.finish()
        return out

    def response_amplitude_exp(
        self, points, atoms_dict_s, norb_per_s, starts_s,
        v, G0, E, eta, E_tip=0.0,
        beta=1.0, r0=3.0, A_ss=-1.0, A_sp=-1.0, rcut=20.0
    ):
        """GPU-accelerated response amplitude map via OpenCL kernel `response_amplitude_exp`.

        Precompute on CPU:
            A_ss = (E + i*eta) * S_s - H_s
            G0 = inv(A_ss)            # complex (ns, ns)
            v  = C_MO^T @ G0          # complex (ns,)

        GPU kernel builds coupling a_st = (E+iη)S_ts - H_ts per grid point
        and computes resp = |v·a_st^H|^2 / |(E+iη-E_tip) - a_st·G0·a_st^H|^2.

        Args:
            points:      (npts, 3) float32 tip positions
            atoms_dict_s: dict with 'pos', 'Rcut', 'type' for sample atoms
            norb_per_s:  (natoms_s,) orbital counts
            starts_s:    (natoms_s+1,) orbital offsets (cumsum)
            v:           (ns,) complex64/128 precomputed v = C^T G0
            G0:          (ns, ns) complex64/128 precomputed Green's function
            E, eta:      float energy and broadening
            E_tip:       float tip onsite energy
            beta, r0, A_ss, A_sp, rcut: exponential SK parameters

        Returns:
            resp: (npts,) float32 response amplitudes
        """
        import numpy as np
        import pyopencl as cl

        points = np.asarray(points, dtype=np.float32)
        if points.ndim != 2 or points.shape[1] != 3:
            raise ValueError(f"points must be (n,3), got {points.shape}")
        npts = len(points)
        natoms_s = len(atoms_dict_s['pos'])
        ns = int(starts_s[-1])
        if ns > 256:
            raise ValueError(f"ns={ns} > 256 kernel private array limit")

        v = np.asarray(v, dtype=np.complex64)
        G0 = np.asarray(G0, dtype=np.complex64)
        starts_s = np.asarray(starts_s, dtype=np.int32)

        atom_data = np.zeros(natoms_s, dtype=[
            ('pos_rcut', 'f4', 4),
            ('type', 'i4'),
            ('i0orb', 'i4'),
            ('norb', 'i4'),
            ('pad', 'i4')
        ])
        for ia in range(natoms_s):
            atom_data[ia]['pos_rcut'][:3] = atoms_dict_s['pos'][ia]
            atom_data[ia]['pos_rcut'][3] = atoms_dict_s['Rcut'][ia]
            Z = int(atoms_dict_s['type'][ia])
            if Z not in self.basis_meta['nz_map']:
                raise RuntimeError(f"response_amplitude_exp: species nz={Z} not loaded")
            atom_data[ia]['type'] = int(self.basis_meta['nz_map'][Z])
            atom_data[ia]['norb'] = int(norb_per_s[ia])
            atom_data[ia]['i0orb'] = int(starts_s[ia])

        mf = cl.mem_flags
        d_points = cl.Buffer(self.ctx, mf.READ_ONLY | mf.COPY_HOST_PTR,
                                hostbuf=np.c_[points, np.zeros((npts, 1), np.float32)].astype(np.float32))
        d_atoms = cl.Buffer(self.ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=atom_data)
        d_starts = cl.Buffer(self.ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=starts_s)
        d_vre = cl.Buffer(self.ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=v.real.astype(np.float32))
        d_vim = cl.Buffer(self.ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=v.imag.astype(np.float32))
        d_G0re = cl.Buffer(self.ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=G0.real.astype(np.float32).ravel())
        d_G0im = cl.Buffer(self.ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=G0.imag.astype(np.float32).ravel())
        d_out = cl.Buffer(self.ctx, mf.WRITE_ONLY, size=npts * 4)

        self._load_kernels()

        gs = (int(npts),)
        ls = None
        self.prg.response_amplitude_exp(
            self.queue, gs, ls,
            np.int32(npts),
            d_points,
            np.int32(natoms_s),
            d_atoms,
            d_starts,
            np.int32(ns),
            d_vre, d_vim,
            d_G0re, d_G0im,
            np.float32(float(E)),
            np.float32(float(eta)),
            np.float32(float(E_tip)),
            np.float32(float(beta)),
            np.float32(float(r0)),
            np.float32(float(A_ss)),
            np.float32(float(A_sp)),
            np.float32(float(rcut)),
            d_out
        )
        self.queue.finish()

        out = np.empty(npts, dtype=np.float32)
        cl.enqueue_copy(self.queue, out, d_out)
        self.queue.finish()
        return out

    def project_orbital(self, coeffs, norb_per, atoms_dict, grid_spec, nMaxAtom=64, _debug_Fortran_order=False):
        """
        Project a single molecular orbital onto a 3D grid using the orbital projection kernel.

        Computes ψ(r) = Σ_i C_i φ_i(r) (signed wavefunction, not density)

        Args:
            coeffs: (natoms, 4) MO coefficients.
                    By default expects FireballOCL convention [px, py, pz, s].
                    If _debug_Fortran_order=True, expects Fortran order [s, py, pz, px] for sp3 atoms.
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
        numorb_max = 4
        coeffs_flat = np.zeros(natoms * numorb_max, dtype=np.float32)

        # Remapping: Fortran [s, py, pz, px] -> FireballOCL [px, py, pz, s]
        _ORT_SPP_TO_OCL = np.array([3, 1, 2, 0], dtype=np.int32)

        for ia in range(natoms):
            i0 = ia * numorb_max
            no = int(norb_per[ia])
            if _debug_Fortran_order and (no == 4):
                coeffs_flat[i0:i0+4] = coeffs[ia, :4][_ORT_SPP_TO_OCL]
            else:
                # Expect coeffs already in [px,py,pz,s]
                # IMPORTANT: even for H (no==1) we still need the s coefficient in slot 3.
                coeffs_flat[i0:i0+4] = coeffs[ia, :4]

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
        d_grid = cl.Buffer(self.ctx, cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR, hostbuf=self.grid_to_np(grid_spec))

        # Tasks
        if len(tasks_np) > 0:
            d_tasks = cl.Buffer(self.ctx, cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR,  hostbuf=tasks_np)
        else:
            d_tasks = cl.Buffer(self.ctx, cl.mem_flags.READ_ONLY, size=32)

        d_atoms = cl.Buffer(self.ctx, cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR, hostbuf=atom_data)

        if len(task_atoms_np) > 0:
            d_task_atoms = cl.Buffer(self.ctx, cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR,  hostbuf=task_atoms_np)
        else:
            d_task_atoms = cl.Buffer(self.ctx, cl.mem_flags.READ_ONLY, size=nMaxAtom * 4)

        d_coeffs = cl.Buffer(self.ctx, cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR, hostbuf=coeffs_flat.astype(np.float32))

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

    def prepare_orbital_projection(self, atoms_dict, grid_spec, nMaxAtom=64):
        """
        One-time setup for batched orbital projection.
        Pre-uploads all static GPU buffers (tasks, atom_data, grid_spec).
        Returns a context dict to pass to project_orbital_prepped() in a loop.
        Eliminates ~6 buffer allocations + build_tasks per orbital call.
        """
        mf = cl.mem_flags
        t0 = time.perf_counter_ns()

        tasks_np, task_atoms_np = self.build_tasks(atoms_dict, grid_spec, nMaxAtom=nMaxAtom)
        n_tasks = len(tasks_np)
        natoms   = len(atoms_dict['pos'])
        numorb_max = 4
        nx, ny, nz = [int(x) for x in grid_spec['ngrid'][:3]]
        out_nbytes    = nx * ny * nz * 4
        coeffs_nbytes = natoms * numorb_max * 4

        atom_data = np.zeros(natoms, dtype=[
            ('pos_rcut', 'f4', 4), ('type', 'i4'), ('i0orb', 'i4'), ('norb', 'i4'), ('pad', 'i4')
        ])
        for ia in range(natoms):
            atom_data[ia]['pos_rcut'][:3] = atoms_dict['pos'][ia]
            atom_data[ia]['pos_rcut'][3]  = atoms_dict['Rcut'][ia]
            Z = int(atoms_dict['type'][ia])
            atom_data[ia]['type']  = int(self.basis_meta['nz_map'][Z])
            atom_data[ia]['norb']  = numorb_max
            atom_data[ia]['i0orb'] = ia * numorb_max

        # Upload static buffers once
        d_grid       = cl.Buffer(self.ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=self.grid_to_np(grid_spec))
        d_tasks      = cl.Buffer(self.ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=tasks_np)      if n_tasks > 0            else cl.Buffer(self.ctx, mf.READ_ONLY, size=32)
        d_atoms      = cl.Buffer(self.ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=atom_data)
        d_task_atoms = cl.Buffer(self.ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=task_atoms_np) if len(task_atoms_np) > 0 else cl.Buffer(self.ctx, mf.READ_ONLY, size=nMaxAtom * 4)

        # Mutable buffers: reused every orbital (only coeffs changes)
        d_coeffs = cl.Buffer(self.ctx, mf.READ_ONLY,  size=coeffs_nbytes)
        d_out    = cl.Buffer(self.ctx, mf.WRITE_ONLY, size=out_nbytes)

        self._load_kernels()
        t1 = time.perf_counter_ns()
        print(f"[prepare_orbital_projection] Setup: {(t1-t0)*1e-6:.1f} ms, n_tasks={n_tasks}")
        return dict(
            n_tasks=n_tasks, natoms=natoms, numorb_max=numorb_max,
            nx=nx, ny=ny, nz=nz, out_nbytes=out_nbytes, nMaxAtom=nMaxAtom,
            d_grid=d_grid, d_tasks=d_tasks, d_atoms=d_atoms,
            d_task_atoms=d_task_atoms, d_coeffs=d_coeffs, d_out=d_out,
        )

    def project_orbital_prepped(self, coeffs_flat, ctx):
        """
        Project a single orbital using pre-built GPU buffers from prepare_orbital_projection().
        Only uploads coeffs_flat per call. Returns (nx,ny,nz) psi grid (float32).
        """
        cl.enqueue_copy(self.queue, ctx['d_coeffs'], coeffs_flat.astype(np.float32))
        cl.enqueue_fill_buffer(self.queue, ctx['d_out'], np.float32(0), 0, ctx['out_nbytes'])

        ls = (32,)
        gs = (ctx['n_tasks'] * ls[0],)
        self.prg.project_orbital(
            self.queue, gs, ls,
            ctx['d_grid'], np.int32(ctx['n_tasks']),
            ctx['d_tasks'], ctx['d_atoms'], ctx['d_task_atoms'],
            ctx['d_coeffs'],
            self.d_basis,
            np.int32(self.basis_meta['n_nodes']),
            np.float32(self.basis_meta['dr']),
            np.int32(self.basis_meta['max_shells']),
            np.int32(ctx['numorb_max']),
            np.int32(ctx['nMaxAtom']),
            ctx['d_out']
        )
        self.queue.finish()
        res = np.empty((ctx['nx'], ctx['ny'], ctx['nz']), dtype=np.float32)
        cl.enqueue_copy(self.queue, res, ctx['d_out'])
        self.queue.finish()
        return res

# ================================================================
# High-level evaluation functions
# ================================================================

def evaluate_mos_on_points(projector, mo_indices, points, evecs, natoms, species_per_atom,
                          species_names, species_list_ang, norb_per_atom, atoms_dict):
    """
    Evaluate multiple MOs at arbitrary points using OpenCL.
    
    Args:
        projector: GridProjector instance (already configured with load_basis_sto)
        mo_indices: list of 0-indexed MO indices
        points: array [npoints, 3] in Angstrom
        evecs: (nstates, norb) eigenvector array
        natoms: number of atoms
        species_per_atom: (natoms,) 0-based species indices
        species_names: list of species names
        species_list_ang: list from parse_basis_hsd_ang
        norb_per_atom: (natoms,) number of orbitals per atom
        atoms_dict: dict with 'pos', 'Rcut', 'type' for atoms
    
    Returns:
        point_vals: list of arrays, each [npoints]
    """
    from .DFTBplusParser import evec_to_kernel_coeffs
    
    point_vals = []
    for imo in mo_indices:
        coeffs = evec_to_kernel_coeffs(evecs[imo], natoms, species_per_atom,
                                      species_names, species_list_ang)
        psi = projector.project_orbital_points(points, coeffs, norb_per_atom, atoms_dict)
        point_vals.append(psi)
    
    return point_vals


def setup_gridprojector_from_dftb(dftb_data, species_list_ang, ctx=None, queue=None, verbosity=0):
    """
    Configure GridProjector instance from parsed DFTB+ data.
    
    Args:
        dftb_data: dict with keys:
            - coords_bohr: (natoms, 3) array (Bohr)
            - species_per_atom: (natoms,) 0-based species indices
            - species_names: list of species names
        species_list_ang: list from parse_basis_hsd_ang (Å units)
        ctx: OpenCL context (optional)
        queue: OpenCL command queue (optional)
        verbosity: verbosity level
    
    Returns:
        projector: configured GridProjector instance
        atoms_dict: dict with 'pos', 'Rcut', 'type' for projection
    """
    projector = GridProjector(fdata_dir=None, ctx=ctx, queue=queue, verbosity=verbosity)
    projector.load_basis_sto(species_list_ang)
    
    # Build atoms dict
    coords_ang = dftb_data['coords_bohr'] * 0.5291772109  # Bohr -> Angstrom
    sp_by_name = {sp['name']: sp for sp in species_list_ang}
    
    natoms = len(coords_ang)
    atomic_numbers = np.array([sp_by_name[dftb_data['species_names'][si]]['atomic_number'] for si in dftb_data['species_per_atom']])
    cutoffs = np.array([sp_by_name[dftb_data['species_names'][si]]['orbitals'][0]['cutoff'] for si in dftb_data['species_per_atom']])
    
    atoms_dict = {
        'pos': coords_ang,
        'Rcut': cutoffs,
        'type': atomic_numbers
    }
    
    return projector, atoms_dict


BOHR2ANG = 0.5291772109
B3_FACTOR = 1.0 / (BOHR2ANG**3)   # orbital-normalized-in-Bohr → density in Ang^-3


def project_dftb_density(geo, evecs, projector, atoms_dict, grid_spec, basis, B3_FACTOR=B3_FACTOR):
    """
    Project SCF electron density onto a real-space grid.
    Occupied MOs are accumulated as sum_i occ_i * |psi_i|^2.

    Returns rho_grid (nx,ny,nz) float32 in e/Å³.
    """
    from .DFTBplusParser import precompute_coeff_gather
    import time
    natoms = geo['natoms']
    occs   = geo['occupations'][:, 0, 0]
    occ_idx = np.where(occs > 1e-6)[0]
    print(f"[project_dftb_density] {len(occ_idx)} occupied states")

    src_idx, dst_idx = precompute_coeff_gather(natoms, geo['species_per_atom'], geo['species_names'], basis)
    proj_ctx = projector.prepare_orbital_projection(atoms_dict, grid_spec)

    nx, ny, nz = proj_ctx['nx'], proj_ctx['ny'], proj_ctx['nz']
    rho_grid    = np.zeros((nx, ny, nz), dtype=np.float32)
    coeffs_flat = np.zeros(natoms * 4, dtype=np.float32)

    t0 = time.time()
    for i in occ_idx:
        coeffs_flat[:] = 0.0
        coeffs_flat[dst_idx] = evecs[i][src_idx]
        psi = projector.project_orbital_prepped(coeffs_flat, proj_ctx)
        rho_grid += occs[i] * (psi ** 2)
    n = len(occ_idx)
    print(f"[project_dftb_density] {time.time()-t0:.2f}s  ({(time.time()-t0)/n*1e3:.1f} ms/orbital)")
    return (rho_grid * B3_FACTOR).astype(np.float32)


def project_neutral_density(geo, projector, atoms_dict, grid_spec, basis, B3_FACTOR=B3_FACTOR):
    """
    Project superposition of neutral atom densities onto a real-space grid.
    Uses reference valence occupations per (l, Z).

    Returns rho_na_grid (nx,ny,nz) float32 in e/Å³.
    """
    from .DFTBplusParser import precompute_coeff_gather
    import time
    OCC_NA = {1: {0: 1.0}, 6: {0: 2.0, 1: 2/3}, 7: {0: 2.0, 1: 1.0}, 8: {0: 2.0, 1: 4/3}}

    natoms          = geo['natoms']
    species_per_atom = geo['species_per_atom']
    species_names   = geo['species_names']
    sp_by_name      = {sp['name']: sp for sp in basis}

    src_idx, dst_idx = precompute_coeff_gather(natoms, species_per_atom, species_names, basis)
    proj_ctx = projector.prepare_orbital_projection(atoms_dict, grid_spec)

    nx, ny, nz  = proj_ctx['nx'], proj_ctx['ny'], proj_ctx['nz']
    rho_na      = np.zeros((nx, ny, nz), dtype=np.float32)
    coeffs_flat = np.zeros(natoms * 4, dtype=np.float32)

    t0 = time.time()
    orb_offset = 0
    for ia in range(natoms):
        Z  = int(atoms_dict['type'][ia])
        sp = sp_by_name[species_names[species_per_atom[ia]]]
        occ_by_l = OCC_NA.get(Z, {0: 2.0, 1: 2/3})
        for orb in sp['orbitals']:
            l  = orb['l']
            nm = 2 * l + 1
            f_na = occ_by_l.get(l, 0.0)
            if f_na > 0:
                for m in range(nm):
                    coeffs_flat[:] = 0.0
                    coeffs_flat[dst_idx[orb_offset + m]] = 1.0
                    psi = projector.project_orbital_prepped(coeffs_flat, proj_ctx)
                    rho_na += f_na * (psi ** 2)
            orb_offset += nm
    print(f"[project_neutral_density] {time.time()-t0:.2f}s")
    return (rho_na * B3_FACTOR).astype(np.float32)