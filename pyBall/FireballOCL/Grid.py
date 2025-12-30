import numpy as np
import pyopencl as cl
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

    def build_tasks(self, atoms, fb_dims, neighs, grid_spec, block_res=8):
        """
        Partition the grid into tasks (active blocks).
        grid_spec: {'origin': [x,y,z], 'dA': [x,y,z], 'dB': [x,y,z], 'dC': [x,y,z], 'ngrid': [nx,ny,nz]}
        """
        neigh_max = getattr(neighs, "neigh_max", None)
        if neigh_max is None:
            # try to infer from arrays
            if hasattr(neighs, "neigh_j"):
                neigh_max = neighs.neigh_j.shape[1] if neighs.neigh_j.ndim==2 else neighs.neigh_j.shape[-1]
            else:
                raise AttributeError("neighs.neigh_max missing and cannot infer from neigh_j")
        n_atoms_input = len(atoms['pos'])
        print(f"[DEBUG] build_tasks: n_atoms_input={n_atoms_input}, neigh_max={neigh_max}")
        print(f"[DEBUG] atoms pos min={np.min(atoms['pos'],axis=0)}, max={np.max(atoms['pos'],axis=0)}")
        if hasattr(neighs, "neighn"):
            print(f"[DEBUG] neighn stats: min={np.min(neighs.neighn)}, max={np.max(neighs.neighn)}, mean={np.mean(neighs.neighn)}")
        if hasattr(neighs, "neigh_j"):
            nj = neighs.neigh_j
            print(f"[DEBUG] neigh_j shape={nj.shape}, sample row0={nj[0,:min(8,nj.shape[1])]}")
            # count valid neighbor ids ( >0 )
            valid = nj>0
            print(f"[DEBUG] neigh_j valid count={valid.sum()}, max_val={nj.max()}")
        # 1. Setup Macro-grid
        max_rcut = np.max(atoms['Rcut'])
        macro_res = 2.0 * max_rcut + 0.1
        
        origin = np.array(grid_spec['origin'])
        ngrid = np.array(grid_spec['ngrid'])
        dA = np.array(grid_spec['dA'])
        dB = np.array(grid_spec['dB'])
        dC = np.array(grid_spec['dC'])
        
        # Bounding box of atoms
        apos = np.array(atoms['pos'])
        ar = np.array(atoms['Rcut'])
        
        # Macro cells mapping
        macro_grid = {} # (ix, iy, iz) -> [atom_indices]
        for i, (p, r) in enumerate(zip(apos, ar)):
            idx_min = np.floor((p - r - origin) / macro_res).astype(int)
            idx_max = np.floor((p + r - origin) / macro_res).astype(int)
            
            for ix in range(max(0, idx_min[0]), min(int(np.ceil(ngrid[0]*dA[0]/macro_res)), idx_max[0] + 1)):
                for iy in range(max(0, idx_min[1]), min(int(np.ceil(ngrid[1]*dB[1]/macro_res)), idx_max[1] + 1)):
                    for iz in range(max(0, idx_min[2]), min(int(np.ceil(ngrid[2]*dC[2]/macro_res)), idx_max[2] + 1)):
                        macro_grid.setdefault((ix, iy, iz), []).append(i)
        
        tasks = []
        active_atoms_list = []
        active_pairs_list = []
        
        # Precompute atom-atom adjacency based on cutoffs
        n_atoms = len(atoms['pos'])
        adj = np.zeros((n_atoms, n_atoms), dtype=bool)
        for i in range(n_atoms):
            for j in range(i, n_atoms):
                dist = np.linalg.norm(apos[i] - apos[j])
                if dist < (ar[i] + ar[j]):
                    adj[i, j] = adj[j, i] = True

        # Fine blocks side lengths in Angstroem
        fine_side = np.array([dA[0], dB[1], dC[2]]) * block_res

        for (mix, miy, miz), m_atoms in macro_grid.items():
            # Range of fine blocks in this macro cell
            f_min = np.floor(np.array([mix, miy, miz]) * macro_res / fine_side).astype(int)
            f_max = np.floor((np.array([mix, miy, miz]) + 1) * macro_res / fine_side).astype(int)
            
            for fix in range(f_min[0], f_max[0] + 1):
                for fiy in range(f_min[1], f_max[1] + 1):
                    for fiz in range(f_min[2], f_max[2] + 1):
                        if fix * block_res >= ngrid[0] or fiy * block_res >= ngrid[1] or fiz * block_res >= ngrid[2]: continue
                        
                        b_min = origin + np.array([fix, fiy, fiz]) * fine_side
                        b_max = b_min + fine_side
                        
                        overlapping_atoms = []
                        for i in m_atoms:
                            if self.check_overlap_sphere_aabb(apos[i], ar[i], b_min, b_max):
                                overlapping_atoms.append(i)
                        
                        if not overlapping_atoms: continue
                        
                        # Count active pairs
                        overlapping_pairs = [] # list of (i, j, ineigh_ij)
                        na = len(overlapping_atoms)
                        for idx_a in range(na):
                            iatom = overlapping_atoms[idx_a]
                            for idx_b in range(na):
                                jatom = overlapping_atoms[idx_b]
                                if not adj[iatom, jatom]: continue
                                
                                # Find ineigh such that neigh_j[iatom, ineigh] == jatom + 1
                                # We can precompute this mapping for efficiency
                                # For now, simple search
                                ineigh_ij = -1
                                for k in range(neigh_max):
                                    if neighs.neigh_j[iatom, k] == jatom + 1:
                                        ineigh_ij = k
                                        break
                                
                                if ineigh_ij >= 0:
                                    overlapping_pairs.append((iatom, jatom, ineigh_ij))
                        
                        if not overlapping_pairs: continue
                        
                        task = {
                            'block_idx': (fix, fiy, fiz),
                            'atom_start': len(active_atoms_list),
                            'n_atoms': len(overlapping_atoms),
                            'pair_start': len(active_pairs_list),
                            'n_pairs': len(overlapping_pairs)
                        }
                        tasks.append(task)
                        active_atoms_list.extend(overlapping_atoms)
                        # Store as (iatom, jatom, ineigh_ij, pad)
                        for p in overlapping_pairs:
                            active_pairs_list.append((p[0], p[1], p[2], 0))

        # Convert to numpy arrays for GPU
        tasks_np = np.zeros(len(tasks), dtype=[
            ('block_idx', 'i4', 4),
            ('atom_start', 'i4'),
            ('n_atoms', 'i4'),
            ('pair_start', 'i4'),
            ('n_pairs', 'i4')
        ])
        for i, t in enumerate(tasks):
            tasks_np[i]['block_idx'][:3] = t['block_idx']
            tasks_np[i]['atom_start'] = t['atom_start']
            tasks_np[i]['n_atoms'] = t['n_atoms']
            tasks_np[i]['pair_start'] = t['pair_start']
            tasks_np[i]['n_pairs'] = t['n_pairs']
        
        # DEBUG summary
        uniq_atoms = len(set(active_atoms_list))
        total_atoms_refs = len(active_atoms_list)
        total_pairs = len(active_pairs_list)
        max_pairs = max((t['n_pairs'] for t in tasks), default=0)
        max_atoms = max((t['n_atoms'] for t in tasks), default=0)
        print(f"[DEBUG] tasks: n_tasks={len(tasks)}, uniq_atoms={uniq_atoms}, atom_refs={total_atoms_refs}, pairs={total_pairs}, max_atoms/task={max_atoms}, max_pairs/task={max_pairs}")
        return tasks_np, np.array(active_atoms_list, dtype=np.int32), np.array(active_pairs_list, dtype=np.int32)

    def project(self, rho, neighs, atoms, grid_spec, tasks=None):
        """
        Main entry point for density projection.
        rho: np.ndarray shape (natoms, neigh_max, numorb_max, numorb_max)
        neighs: FireballData object containing neighn, neigh_j, etc.
        atoms: dict with 'pos', 'Rcut', 'type' (Z)
        grid_spec: grid parameters
        """
        if tasks is None:
            tasks_np, active_atoms, active_pairs = self.build_tasks(atoms, None, neighs, grid_spec)
        else:
            tasks_np, active_atoms, active_pairs = tasks

        n_tasks = len(tasks_np)
        if n_tasks == 0:
            return np.zeros(grid_spec['ngrid'], dtype=np.float32)

        # DEBUG: print buffer sizes to track resource usage
        nx, ny, nz = grid_spec['ngrid']
        ngrid_total = int(nx) * int(ny) * int(nz)
        print(f"[DEBUG] project: n_tasks={n_tasks}, active_atoms={len(active_atoms)}, active_pairs={len(active_pairs)}")
        print(f"[DEBUG] grid_spec: ngrid={grid_spec['ngrid']}, origin={grid_spec['origin']}, dA={grid_spec['dA']}, dB={grid_spec['dB']}, dC={grid_spec['dC']}")
        print(f"[DEBUG] rho shape={rho.shape}, dtype={rho.dtype}, bytes={rho.nbytes}")
        print(f"[DEBUG] out buffer bytes={ngrid_total*4}")

        # 1. Prepare atom and species data for GPU
        natoms = len(atoms['pos'])
        atom_data = np.zeros(natoms, dtype=[
            ('pos_rcut', 'f4', 4),
            ('type', 'i4'),
            ('i0orb', 'i4'),
            ('norb', 'i4')
        ])
        
        nz_map = self.basis_meta['nz_map']
        
        # Species info for orbital loops
        # nzx in neighs maps fb_idx -> Z
        # We need mapping fb_idx -> basis_idx
        fb_idx_to_basis_idx = {}
        for fb_idx, z in enumerate(neighs.nzx):
            if z in nz_map:
                fb_idx_to_basis_idx[fb_idx] = nz_map[z]

        all_nz_basis = sorted(nz_map.keys())
        species_info_np = np.zeros(len(all_nz_basis), dtype=[
            ('nssh', 'i4'),
            ('lssh', 'i4', 4)
        ])
        for i, nz in enumerate(all_nz_basis):
            info = self.parser.species_info[nz]
            species_info_np[i]['nssh'] = info['nssh']
            species_info_np[i]['lssh'][:len(info['lssh'])] = info['lssh']

        Z_to_fb_idx = {z: i for i, z in enumerate(neighs.nzx)}

        for i in range(natoms):
            atom_data[i]['pos_rcut'][:3] = atoms['pos'][i]
            atom_data[i]['pos_rcut'][3] = atoms['Rcut'][i]
            nz = atoms['type'][i]
            atom_data[i]['type'] = nz_map[nz]
            
            fb_idx = Z_to_fb_idx[nz]
            atom_data[i]['norb'] = 4 if nz == 6 else 1 # TODO: use info from FireballData
            atom_data[i]['i0orb'] = neighs.degelec[i]

        # 2. Buffers
        d_tasks = cl.Buffer(self.ctx, cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR, hostbuf=tasks_np)
        d_atoms = cl.Buffer(self.ctx, cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR, hostbuf=atom_data)
        d_active_atoms = cl.Buffer(self.ctx, cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR, hostbuf=active_atoms)
        d_active_pairs = cl.Buffer(self.ctx, cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR, hostbuf=active_pairs.astype(np.int32))
        d_species_info = cl.Buffer(self.ctx, cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR, hostbuf=species_info_np)
        
        rho32 = rho.astype(np.float32)
        d_rho = cl.Buffer(self.ctx, cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR, hostbuf=rho32)
        
        # d_neigh_j not strictly needed if we use precomputed active_pairs with ineigh
        
        d_out = cl.Buffer(self.ctx, cl.mem_flags.WRITE_ONLY, nx * ny * nz * 4)
        cl.enqueue_fill_buffer(self.queue, d_out, np.float32(0), 0, nx * ny * nz * 4)

        # 3. Kernel launch
        gs = (n_tasks * 8, 8, 8)
        ls = (8, 8, 8)
        
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

        # Dummy neigh_j buffer if not used by refined kernel
        d_dummy = cl.Buffer(self.ctx, cl.mem_flags.READ_ONLY, 4)

        self.prg.project_density_sparse(
            self.queue, gs, ls,
            grid_spec_np[0],
            np.int32(n_tasks),
            d_tasks, d_atoms, d_active_atoms, d_active_pairs,
            d_rho, 
            d_dummy, # neigh_j placeholder
            self.d_basis,
            d_species_info,
            np.int32(self.basis_meta['n_nodes']),
            np.float32(self.basis_meta['dr']),
            np.int32(self.basis_meta['max_shells']),
            np.int32(rho.shape[1]), # neigh_max
            np.int32(rho.shape[2]), # numorb_max
            d_out
        )

        res = np.empty((nx, ny, nz), dtype=np.float32)
        cl.enqueue_copy(self.queue, res, d_out)
        return res
