import numpy as np
import pyopencl as cl
import pyopencl.cltypes
import os
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

    def build_tasks(self, atoms, fb_dims, neighs, grid_spec, block_res=8, nMaxAtom=16):
        """
        Partition the grid into tasks (active blocks).
        grid_spec: {'origin': [x,y,z], 'dA': [x,y,z], 'dB': [x,y,z], 'dC': [x,y,z], 'ngrid': [nx,ny,nz]}
        """
        # neigh_max = neighs.neigh_j.shape[1]
        n_atoms_input = len(atoms['pos'])
        print(f"[DEBUG] build_tasks: n_atoms={n_atoms_input}")
        
        origin = np.array(grid_spec['origin'])
        ngrid = np.array(grid_spec['ngrid'])
        dA = np.array(grid_spec['dA'])
        dB = np.array(grid_spec['dB'])
        dC = np.array(grid_spec['dC'])
        
        # Grid steps
        step = np.array([np.linalg.norm(dA), np.linalg.norm(dB), np.linalg.norm(dC)])
        fine_side = step * block_res
        
        apos = np.array(atoms['pos'])
        ar = np.array(atoms['Rcut'])

        # 1. Map atoms to fine blocks once
        active_fine_blocks = {} # (fix, fiy, fiz) -> list of atoms
        for i in range(n_atoms_input):
            p, r = apos[i], ar[i]
            f_idx_min = np.floor((p - r - origin) / fine_side).astype(int)
            f_idx_max = np.floor((p + r - origin) / fine_side).astype(int)
            
            for fix in range(max(0, f_idx_min[0]), min(int(np.ceil(ngrid[0]/block_res)), f_idx_max[0] + 1)):
                for fiy in range(max(0, f_idx_min[1]), min(int(np.ceil(ngrid[1]/block_res)), f_idx_max[1] + 1)):
                    for fiz in range(max(0, f_idx_min[2]), min(int(np.ceil(ngrid[2]/block_res)), f_idx_max[2] + 1)):
                        # AABB check
                        b_min = origin + np.array([fix, fiy, fiz]) * fine_side
                        b_max = b_min + fine_side
                        if self.check_overlap_sphere_aabb(p, r, b_min, b_max):
                            active_fine_blocks.setdefault((fix, fiy, fiz), []).append(i)

        tasks = []
        task_atoms_list = []
        
        S = nMaxAtom // 2
        for (fix, fiy, fiz), overlapping_atoms in active_fine_blocks.items():
            if len(overlapping_atoms) == 0:
                continue
            
            # Split atoms in block into chunks of size S
            chunks = [overlapping_atoms[i:i + S] for i in range(0, len(overlapping_atoms), S)]
            n_chunks = len(chunks)
            
            # Tile interactions
            covered_diags = set()
            for i in range(n_chunks):
                for j in range(i, n_chunks):
                    if i == j:
                        if i in covered_diags: continue
                        # Diagonal task: try to combine C_i and C_{i+1}
                        atoms_to_add = list(chunks[i])
                        if i + 1 < n_chunks:
                            atoms_to_add.extend(chunks[i+1])
                            covered_diags.add(i+1)
                        
                        tasks.append({
                            'block_idx': (fix, fiy, fiz),
                            'na': len(atoms_to_add),
                            'nj': -1,
                            'atoms': atoms_to_add
                        })
                    else:
                        if j == i + 1 and (i+1) in covered_diags:
                            continue
                        # Off-diagonal task for C_i and C_j
                        atoms_to_add = list(chunks[i]) + list(chunks[j])
                        tasks.append({
                            'block_idx': (fix, fiy, fiz),
                            'na': len(atoms_to_add),
                            'nj': len(chunks[i]),
                            'atoms': atoms_to_add
                        })

        # Sort tasks by workload (na)
        tasks.sort(key=lambda x: x['na'], reverse=True)

        tasks_np = np.zeros(len(tasks), dtype=[
            ('block_idx', 'i4', 4),
            ('na', 'i4'),
            ('nj', 'i4'),
            ('pad1', 'i4'),
            ('pad2', 'i4')
        ])
        
        task_atoms_np = np.zeros((len(tasks), nMaxAtom), dtype=np.int32)
        
        for i, t in enumerate(tasks):
            tasks_np[i]['block_idx'][:3] = t['block_idx']
            tasks_np[i]['na'] = t['na']
            tasks_np[i]['nj'] = t['nj']
            task_atoms_np[i, :t['na']] = t['atoms']

        for i in range(len(tasks_np)):
            print(f"Task {i}: block_idx={tasks_np[i]['block_idx']} na: {tasks_np[i]['na']:<5} nj: {tasks_np[i]['nj']:<5}")    

        print(f"[DEBUG] build_tasks finished: n_tasks={len(tasks)}")
        return tasks_np, task_atoms_np

    def project(self, rho, neighs, atoms, grid_spec, tasks=None, nMaxAtom=16):
        """
        Main entry point for density projection.
        rho: np.ndarray shape (natoms, neigh_max, numorb_max, numorb_max)
        neighs: FireballData object containing neighn, neigh_j, etc.
        atoms: dict with 'pos', 'Rcut', 'type' (Z)
        grid_spec: grid parameters
        """
        if tasks is None:
            tasks_np, task_atoms = self.build_tasks(atoms, None, neighs, grid_spec, nMaxAtom=nMaxAtom)
        else:
            tasks_np, task_atoms = tasks

        n_tasks = len(tasks_np)
        if n_tasks == 0:
            return np.zeros(grid_spec['ngrid'], dtype=np.float32)

        # DEBUG: print buffer sizes to track resource usage
        nx, ny, nz = grid_spec['ngrid']
        ngrid_total = int(nx) * int(ny) * int(nz)
        print(f"[DEBUG] project: n_tasks={n_tasks}")
        print(f"[DEBUG] grid_spec: ngrid={grid_spec['ngrid']}, origin={grid_spec['origin']}")
        print(f"[DEBUG] rho shape={rho.shape}, dtype={rho.dtype}, bytes={rho.nbytes}")
        print(f"[DEBUG] out buffer bytes={ngrid_total*4}")

        # 1. Prepare atom and species data for GPU
        natoms_sys = len(atoms['pos'])
        atom_data = np.zeros(natoms_sys, dtype=[
            ('pos_rcut', 'f4', 4),
            ('type', 'i4'),
            ('i0orb', 'i4'),
            ('norb', 'i4'),
            ('pad', 'i4')
        ])
        
        nz_map = self.basis_meta['nz_map']
        all_nz_basis = sorted(nz_map.keys())
        
        # Species info for orbital loops (packed as int4 for OpenCL)
        species_info_flat = np.zeros(len(all_nz_basis) * 4, dtype=np.int32)
        for i, nz in enumerate(all_nz_basis):
            info = self.parser.species_info[nz]
            species_info_flat[i*4] = info['nssh']
            ls_vals = info['lssh']
            for j in range(min(3, len(ls_vals))):
                species_info_flat[i*4 + 1 + j] = ls_vals[j]

        # Precompute authoritative norb per species from parser (sum of (2l+1) per shell)
        species_norb = {}
        for nz in all_nz_basis:
            info = self.parser.species_info[nz]
            species_norb[nz] = sum(2 * l + 1 for l in info['lssh'])

        for i in range(natoms_sys):
            atom_data[i]['pos_rcut'][:3] = atoms['pos'][i]
            atom_data[i]['pos_rcut'][3] = atoms['Rcut'][i]
            nz = atoms['type'][i]
            atom_data[i]['type'] = nz_map[nz]
            iatyp_z = int(neighs.iatyp[i])
            norb_val = None
            # Prefer num_orb if it is sized to cover Z
            if (iatyp_z - 1) >= 0 and (iatyp_z - 1) < len(neighs.num_orb) and neighs.num_orb[iatyp_z - 1] > 0:
                norb_val = int(neighs.num_orb[iatyp_z - 1])
            elif iatyp_z in species_norb:
                norb_val = species_norb[iatyp_z]
            else:
                raise RuntimeError(f"Missing orbital count for atom {i} Z={iatyp_z}: num_orb len={len(neighs.num_orb)}, species_norb keys={list(species_norb.keys())}")
            atom_data[i]['norb'] = norb_val
            atom_data[i]['i0orb'] = neighs.degelec[i]

        # 2. Buffers
        d_tasks = cl.Buffer(self.ctx, cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR, hostbuf=tasks_np)
        d_atoms = cl.Buffer(self.ctx, cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR, hostbuf=atom_data)
        d_task_atoms = cl.Buffer(self.ctx, cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR, hostbuf=task_atoms)
        d_species_info = cl.Buffer(self.ctx, cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR, hostbuf=species_info_flat)
        
        rho32 = rho.astype(np.float32)
        d_rho = cl.Buffer(self.ctx, cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR, hostbuf=rho32)
        
        d_out = cl.Buffer(self.ctx, cl.mem_flags.WRITE_ONLY, nx * ny * nz * 4)
        cl.enqueue_fill_buffer(self.queue, d_out, np.float32(0), 0, nx * ny * nz * 4)

        # 3. Kernel launch (1D). Each task covers 8*8*8 voxels (512). Flatten tasks*voxels.
        vox_per_task = 512
        total_work = n_tasks * vox_per_task
        ls = (16,)  # 1D local size
        gs = (n_tasks * ls[0],)
        
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
        d_grid = cl.Buffer(self.ctx, cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR, hostbuf=grid_spec_np)

        d_neigh_j = cl.Buffer(self.ctx, cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR, hostbuf=neighs.neigh_j.astype(np.int32))

        # Always create a fresh queue for this context (avoid stale/invalid queue issues)
        self.queue = cl.CommandQueue(self.ctx)

        # Debug device/work sizes to diagnose resource issues
        dev = self.ctx.devices[0]
        print(f"[DEBUG] device={dev.name}, max_work_group={dev.max_work_group_size}, local_mem={dev.local_mem_size}")
        print(f"[DEBUG] gs={gs}, ls={ls}, n_tasks={n_tasks}")

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

        res = np.empty((nx, ny, nz), dtype=np.float32)
        cl.enqueue_copy(self.queue, res, d_out)
        self.queue.finish()
        return res
