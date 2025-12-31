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
    def __init__(self, fdata_dir, ctx=None, queue=None, nloc=32):
        super().__init__(nloc=nloc)
        self.fdata_dir = fdata_dir
        self.parser = FdataParser(fdata_dir)
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
        # For simplicity, assume same mesh and dr for all wfs initially
        all_nz = sorted(self.basis_data.keys())
        if len(all_nz)==0:
            raise RuntimeError("load_basis called with empty species list (species_nz).")
        max_shells = max(len(v) for v in self.basis_data.values())
        if max_shells==0:
            raise RuntimeError(f"No wavefunctions loaded for species {all_nz}")
        first_wf = self.basis_data[all_nz[0]][0]
        n_nodes = first_wf['mesh']
        dr = first_wf['rcutoff'] / (n_nodes - 1)
        
        packed_basis = np.zeros((len(all_nz), max_shells, n_nodes), dtype=np.float32)
        for i, nz in enumerate(all_nz):
            for ish, wf in enumerate(self.basis_data[nz]):
                packed_basis[i, ish, :len(wf['data'])] = wf['data']
        
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
        if hasattr(self, 'basis_meta'):
            build_opts.append(f"-DMAX_ORBS={self.basis_meta['max_shells'] * 9}") # conservative upper bound
        
        self.load_program(kernel_path=cl_path)

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
        if natoms > 0:
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
        print(f"[DEBUG] block atom stats: na_max={max_count}, nbloks: empty={empty_blocks}, one={one_blocks}, multi={multi_blocks}")
        self.last_block_atom_counts = np.array(block_counts, dtype=np.int32)

        tasks_np = np.zeros(len(tasks), dtype=self.task_dtype_np)
        
        task_atoms_np = np.zeros((len(tasks), nMaxAtom), dtype=np.int32)
        
        for i, t in enumerate(tasks):
            tasks_np[i]['x'], tasks_np[i]['y'], tasks_np[i]['z'] = t['block_idx']
            tasks_np[i]['na'] = t['na']
            tasks_np[i]['nj'] = t['nj']
            task_atoms_np[i, :t['na']] = t['atoms']

        print(f"[DEBUG] build_tasks finished: n_tasks={len(tasks)}")
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
        print(f"[DEBUG] grid_spec_np: origin={grid_spec_np[0]['origin']} dA={grid_spec_np[0]['dA']} dB={grid_spec_np[0]['dB']} dC={grid_spec_np[0]['dC']} ngrid={grid_spec_np[0]['ngrid']}")
        
        return grid_spec_np

    def project(self, rho, neighs, atoms, grid_spec, tasks=None, nMaxAtom=64, use_gpu_tasks=False):
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
            print(f"[TIME] build_tasks finished in {(T1-T0)*1e-6:.3f} [ms]")
        else:
            tasks_np, task_atoms_np = tasks

        n_tasks = len(tasks_np)
        nx, ny, nz = grid_spec['ngrid'][:3]
        
        # Prepare other buffers
        natoms = len(atoms['pos'])
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
            atom_data[i]['norb'] = 4 # Default for C, H with s,p
            atom_data[i]['i0orb'] = 0

        # 2. Buffers
        mf = cl.mem_flags
        
        # DEBUG: check tasks_np size and dtype
        print(f"[DEBUG] tasks_np: len={len(tasks_np)} itemsize={tasks_np.dtype.itemsize} nbytes={tasks_np.nbytes}")
        
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
        
        d_out = cl.Buffer(self.ctx, mf.WRITE_ONLY, nx * ny * nz * 4)
        cl.enqueue_fill_buffer(self.queue, d_out, np.float32(0), 0, nx * ny * nz * 4)

        # 3. Kernel launch
        ls = (32,)  # local size
        gs = (n_tasks * ls[0],)
        
        # d_basis placeholder
        if not hasattr(self, 'd_basis'):
             self.d_basis = cl.Buffer(self.ctx, mf.READ_ONLY, size=4)
             self.basis_meta = {'n_nodes': 0, 'dr': 0.0, 'max_shells': 0}

        print(f"[DEBUG] project_tiled: gs={gs}, ls={ls}, n_tasks={n_tasks}")

        T0_ns = time.perf_counter_ns()
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
        self.queue.finish()
        dt_ns = time.perf_counter_ns() - T0_ns
        print(f"[TIME] project_tiled finished in {dt_ns*1e-6:.9f} [ms]")

        res = np.empty((nx, ny, nz), dtype=np.float32)
        cl.enqueue_copy(self.queue, res, d_out)
        self.queue.finish()

        return res
