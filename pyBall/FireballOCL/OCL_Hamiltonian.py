import numpy as np
import pyopencl as cl
import os
from .FdataParser import FdataParser

class OCL_Hamiltonian:
    def __init__(self, fdata_dir, ctx=None, queue=None):
        self.fdata_dir = fdata_dir
        self.parser = FdataParser(fdata_dir)
        if not ctx:
            os.environ['PYOPENCL_CTX'] = '0' # Try to force platform 0
        self.ctx = ctx if ctx else cl.create_some_context(interactive=False)
        self.queue = queue if queue else cl.CommandQueue(self.ctx)
        self.kernels = {}
        self._load_kernels()

    def _load_kernels(self):
        kernel_path = os.path.join(os.path.dirname(__file__), "cl/hamiltonian.cl")
        with open(kernel_path, 'r') as f:
            self.kernel_src = f.read()
        self.prg = cl.Program(self.ctx, self.kernel_src).build()

    def prepare_splines(self, species_nz):
        """Prepares spline data for GPU."""
        d2c, _ = self.parser.load_species_data(species_nz)
        
        numz_max = 0
        n_nz_max = 0
        for k, v in d2c.items():
            numz_max = max(numz_max, v['numz'])
            n_nz_max = max(n_nz_max, v['num_nonzero'])
            
        self.numz_max = numz_max
        self.n_nz_max = n_nz_max
        self.species_pair_map = {}
        
        n_pairs = len(d2c)
        # Pack into [n_pairs, numz_max, n_nz_max, 4]
        spline_data = np.zeros((n_pairs, numz_max, n_nz_max, 4), dtype=np.float32)
        
        sorted_keys = sorted(d2c.keys())
        h_grids = np.zeros(n_pairs, dtype=np.float32)
        for i, k in enumerate(sorted_keys):
            v = d2c[k]
            self.species_pair_map[k] = i
            h_grids[i] = v['zmax'] / (v['numz'] - 1)
            for j in range(v['num_nonzero']):
                spline = self.parser.build_spline_1d(v['data'][:, j], v['zmax'])
                spline_data[i, :v['numz'], j, 0] = spline[0, :]
                spline_data[i, :v['numz'], j, 1] = spline[1, :]
                spline_data[i, :v['numz'], j, 2] = spline[2, :]
                spline_data[i, :v['numz'], j, 3] = spline[3, :]
                
        self.d_splines = cl.Buffer(self.ctx, cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR, hostbuf=spline_data)
        self.d_h_grids = cl.Buffer(self.ctx, cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR, hostbuf=h_grids)

    def prepare_data_3c(self, species_nz):
        """Prepares 3-center data for GPU."""
        _, d3c = self.parser.load_species_data(species_nz)
        
        if not d3c: return
        
        # d3c: {(root, it, isorp, nz1, nz2, nz3): data}
        # data: [numy, numx, n_nz]
        
        numx_max = 0
        numy_max = 0
        n_nz_max = 0
        for k, v in d3c.items():
            numy_max = max(numy_max, v['numy'])
            numx_max = max(numx_max, v['numx'])
            n_nz_max = max(n_nz_max, v['num_nonzero'])
            
        self.numy_3c = numy_max
        self.numx_3c = numx_max
        self.n_nz_3c_max = n_nz_max
        
        # ntheta is always 5 for den3/bcna
        self.ntheta_3c = 5
        
        # We need to map (root, nz1, nz2, nz3) to an index
        self.species_triplet_map = {}
        triplets = sorted(list(set((k[0], k[3], k[4], k[5]) for k in d3c.keys())))
        for i, t in enumerate(triplets):
            self.species_triplet_map[t] = i
            
        n_triplets = len(triplets)
        # Pack into [n_triplets, ntheta, numy, numx, n_nz]
        # For simplicity, we use float (no splines for 3-center in Fireball, just grid)
        data_3c = np.zeros((n_triplets, self.ntheta_3c, numy_max, numx_max, n_nz_max), dtype=np.float32)
        h_grids_3c = np.zeros((n_triplets, 2), dtype=np.float32) # [hx, hy]
        
        for k, v in d3c.items():
            t = (k[0], k[3], k[4], k[5])
            i_triplet = self.species_triplet_map[t]
            i_theta = k[1] - 1 # itheta=1..5
            
            data_3c[i_triplet, i_theta, :v['numy'], :v['numx'], :v['num_nonzero']] = v['data']
            h_grids_3c[i_triplet, 0] = v['xmax'] / (v['numx'] - 1)
            h_grids_3c[i_triplet, 1] = v['ymax'] / (v['numy'] - 1)
            
        self.d_data_3c = cl.Buffer(self.ctx, cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR, hostbuf=data_3c)
        self.d_h_grids_3c = cl.Buffer(self.ctx, cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR, hostbuf=h_grids_3c)
        
    def assemble_2c(self, ratoms, neighbors, pair_types):
        """
        Runs the assemble_2c kernel.
        ratoms: [natoms, 3]
        neighbors: [n_pairs, 2]
        pair_types: [n_pairs]
        """
        n_pairs = len(neighbors)
        # Pad ratoms to float4
        ratoms4 = np.zeros((ratoms.shape[0], 4), dtype=np.float32)
        ratoms4[:, :3] = ratoms
        
        d_ratoms = cl.Buffer(self.ctx, cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR, hostbuf=ratoms4)
        d_neighs = cl.Buffer(self.ctx, cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR, hostbuf=neighbors.astype(np.int32))
        d_types = cl.Buffer(self.ctx, cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR, hostbuf=pair_types.astype(np.int32))
        
        # Output: [n_pairs, 4, 4] float blocks
        d_blocks = cl.Buffer(self.ctx, cl.mem_flags.WRITE_ONLY, n_pairs * 16 * 4) 
        
        self.prg.assemble_2c(self.queue, (n_pairs,), None, 
                           np.int32(n_pairs), np.int32(self.n_nz_max), np.int32(self.numz_max),
                           d_ratoms, d_neighs, self.d_splines, d_types, self.d_h_grids, d_blocks)
        
        blocks = np.zeros((n_pairs, 4, 4), dtype=np.float32)
        cl.enqueue_copy(self.queue, blocks, d_blocks)
        return blocks

    def assemble_full(self, ratoms, species, neighbors):
        """
        Assembles full H and S matrices.
        species: list of nuclear charges for each atom [natoms]
        neighbors: list of (i, j) for 2-center
        """
        n_atoms = len(species)
        n_pairs = len(neighbors)
        
        # 1. 2-center components
        pairs_S = []
        pairs_T = []
        pairs_Vna = []
        for idx, (i, j) in enumerate(neighbors):
            nz1, nz2 = species[i], species[j]
            tS = self.species_pair_map.get(('overlap', nz1, nz2))
            tT = self.species_pair_map.get(('kinetic', nz1, nz2))
            tV = self.species_pair_map.get(('vna', nz1, nz2))
            if tS is not None: pairs_S.append((i, j, tS, idx))
            if tT is not None: pairs_T.append((i, j, tT, idx))
            if tV is not None: pairs_Vna.append((i, j, tV, idx))
            
        # Helper to run assembly and map back to neighbor-indexed blocks
        def run_2c(pairs):
            if not pairs: return np.zeros((n_pairs, 4, 4), dtype=np.float32)
            res = self.assemble_2c(ratoms, np.array([p[:2] for p in pairs]), np.array([p[2] for p in pairs]))
            full_res = np.zeros((n_pairs, 4, 4), dtype=np.float32)
            for idx, p in enumerate(pairs):
                full_res[p[3]] = res[idx]
            return full_res

        S_blocks = run_2c(pairs_S)
        T_blocks = run_2c(pairs_T)
        Vna_blocks = run_2c(pairs_Vna)
        
        H_blocks = T_blocks + Vna_blocks
        
        return H_blocks, S_blocks

if __name__ == "__main__":
    import os
    fdata_dir = "/home/prokophapala/git/FireCore/tests/pyFireball/Fdata"
    ocl = OCL_Hamiltonian(fdata_dir)
    
    # H2 molecule test
    species_nz = [1, 1]
    ocl.prepare_splines(species_nz)
    
    ratoms = np.array([[0, 0, 0], [0, 0, 0.74]], dtype=np.float32)
    neighbors = [(0, 1), (1, 0)]
    
    H, S = ocl.assemble_full(ratoms, species_nz, neighbors)
    print("H2 Overlap S[0,1]:\n", S[0,0,0])
    print("H2 Kinetic T[0,1]:\n", H[0,0,0])
