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
        build_opts = []
        if os.environ.get("DEBUG_INTERP", "0") == "1":
            build_opts.append("-DDEBUG_INTERP=1")
        opts = " ".join(build_opts) if build_opts else None
        if opts:
            self.prg = cl.Program(self.ctx, self.kernel_src).build(options=[opts])
        else:
            self.prg = cl.Program(self.ctx, self.kernel_src).build()

    def _resolve_pair_type(self, root, nz1, nz2):
        """
        Resolve pair type with vna fallback. DEBUG prints to track missing data.
        """
        key = (root, nz1, nz2)
        pt = self.species_pair_map.get(key)
        # DEBUG: if vna missing, try atom/ontop variants
        if pt is None and root == 'vna':
            for alt in ('vna_atom_00', 'vna_ontopl_00', 'vna_ontopr_00'):
                pt = self.species_pair_map.get((alt, nz1, nz2))
                if pt is not None:
                    print(f"[DEBUG] vna missing for ({nz1},{nz2}), using {alt}")
                    return pt
        return pt

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

        # Alias bare vna to neutral-atom tables if only those exist (Fortran interaction 4 uses vna_atom)
        for (root, nz1, nz2), idx in list(self.species_pair_map.items()):
            if root.startswith('vna_atom_'):
                bare = ('vna', nz1, nz2)
                if bare not in self.species_pair_map:
                    self.species_pair_map[bare] = idx
                
        self.d_splines = cl.Buffer(self.ctx, cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR, hostbuf=spline_data)
        self.d_h_grids = cl.Buffer(self.ctx, cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR, hostbuf=h_grids)

    def _is_s_only(self, nz):
        """Returns True if species nz has only s-shell (l=0) as per info.dat."""
        if not hasattr(self.parser, 'species_info'):
            self.parser.parse_info()
        info = self.parser.species_info.get(nz, {})
        lssh = info.get('lssh', [])
        return len(lssh) == 1 and lssh[0] == 0

    def prepare_data_3c(self, species_nz):
        """Prepares 3-center data for GPU."""
        _, d3c = self.parser.load_species_data(species_nz)
        
        if not d3c: return
        
        # d3c: {(root, it, isorp, nz1, nz2, nz3): data}
        # FdataParser.read_3c returns data in file/read order:
        #   data: [numy, numx, n_nz] (y,x) which matches OpenCL buffer packing.
        
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

    def assemble_3c(self, ratoms, triplets):
        """
        Runs the assemble_3c kernel.
        ratoms: [natoms, 3]
        triplets: [n_triplets, 4] (i, j, k, type_idx)
        """
        n_triplets = len(triplets)
        ratoms4 = np.zeros((ratoms.shape[0], 4), dtype=np.float32)
        ratoms4[:, :3] = ratoms
        
        d_ratoms = cl.Buffer(self.ctx, cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR, hostbuf=ratoms4)
        d_triplets = cl.Buffer(self.ctx, cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR, hostbuf=triplets.astype(np.int32))
        
        # Output: [n_triplets, n_nz_3c_max]
        d_results = cl.Buffer(self.ctx, cl.mem_flags.WRITE_ONLY, n_triplets * self.n_nz_3c_max * 4)
        
        self.prg.assemble_3c(self.queue, (n_triplets,), None,
                           np.int32(n_triplets), np.int32(self.n_nz_3c_max),
                           np.int32(self.numx_3c), np.int32(self.numy_3c),
                           d_ratoms, d_triplets, self.d_data_3c, self.d_h_grids_3c, d_results)
        
        results = np.zeros((n_triplets, self.n_nz_3c_max), dtype=np.float32)
        cl.enqueue_copy(self.queue, results, d_results)
        return results

    def scanHamPiece2c(self, root, nz1, nz2, dR, applyRotation=True):
        """Probes a 2-center piece for verification."""
        ratoms = np.array([[0,0,0], dR], dtype=np.float32)
        neighbors = np.array([[0, 1]], dtype=np.int32)
        pair_type = self._resolve_pair_type(root, nz1, nz2)
        if pair_type is None:
            print(f"[WARN] scanHamPiece2c missing pair_type for {root} ({nz1},{nz2}); returning None")
            return None
        
        if applyRotation:
            blocks = self.assemble_2c(ratoms, neighbors, np.array([pair_type]))
        else:
            # Move along Z to avoid rotation
            r = np.linalg.norm(dR)
            ratoms_z = np.array([[0,0,0], [0,0,r]], dtype=np.float32)
            blocks = self.assemble_2c(ratoms_z, neighbors, np.array([pair_type]))
        return blocks[0]

    def scanHamPiece2c_batch(self, root, nz1, nz2, dRs, applyRotation=True):
        """Batch probe of 2-center pieces; one work-item per point."""
        pair_type = self._resolve_pair_type(root, nz1, nz2)
        if pair_type is None:
            print(f"[WARN] scanHamPiece2c_batch missing pair_type for {root} ({nz1},{nz2}); returning None")
            return None
        dRs_np = np.ascontiguousarray(dRs, dtype=np.float32)
        npoints = dRs_np.shape[0]
        # outputs [npoints,4,4]
        d_in = cl.Buffer(self.ctx, cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR, hostbuf=dRs_np)
        d_out = cl.Buffer(self.ctx, cl.mem_flags.WRITE_ONLY, npoints * 16 * 4)
        self.prg.scan_2c_points(
            self.queue, (npoints,), None,
            np.int32(npoints), np.int32(self.n_nz_max), np.int32(self.numz_max),
            np.int32(pair_type), np.int32(int(applyRotation)),
            d_in, self.d_splines, self.d_h_grids, d_out
        )
        blocks = np.zeros((npoints, 4, 4), dtype=np.float32)
        cl.enqueue_copy(self.queue, blocks, d_out)
        return blocks

    def scanHamPiece3c(self, root, nz1, nz2, nz3, dRj, dRk, applyRotation=True):
        """Probes a 3-center piece for verification."""
        # Basis 1 at (0,0,0), Basis 2 at dRj, Potential at dRk
        ratoms = np.array([[0,0,0], dRj, dRk], dtype=np.float32)
        triplets = np.array([[0, 1, 2, 0]], dtype=np.int32) # triplet index 0
        
        type_idx = self.species_triplet_map.get((root, nz1, nz2, nz3))
        if type_idx is None: return None
        triplets[0, 3] = type_idx
        
        # Currently assemble_3c only returns hlist (Legendre coeffs sum)
        # and doesn't rotate. Verification should match this or we add rotation.
        # FireCore's scanHamPiece3c returns bcnax (rotated).
        # To match, we need orbital mapping and rotation in OpenCL.
        # For now, let's just return the hlist result.
        results = self.assemble_3c(ratoms, triplets)

        # Special-case recovery for s-only species (e.g., H-H s-s vna): map hlist[0] to (0,0)
        if self._is_s_only(nz1) and self._is_s_only(nz2):
            mat = np.zeros((1,1), dtype=np.float32)
            mat[0,0] = results[0,0]
            return mat

        return results[0]

    def scanHamPiece3c_batch(self, root, nz1, nz2, nz3, dRjs, dRks, applyRotation=True):
        """Batch probe of 3-center pieces; one work-item per point."""
        type_idx = self.species_triplet_map.get((root, nz1, nz2, nz3))
        if type_idx is None: return None
        dRjs_np = np.ascontiguousarray(dRjs, dtype=np.float32)
        dRks_np = np.ascontiguousarray(dRks, dtype=np.float32)
        npoints = dRjs_np.shape[0]
        d_j = cl.Buffer(self.ctx, cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR, hostbuf=dRjs_np)
        d_k = cl.Buffer(self.ctx, cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR, hostbuf=dRks_np)
        d_out = cl.Buffer(self.ctx, cl.mem_flags.WRITE_ONLY, npoints * self.n_nz_3c_max * 4)
        self.prg.scan_3c_points(
            self.queue, (npoints,), None,
            np.int32(npoints), np.int32(self.n_nz_3c_max),
            np.int32(self.numx_3c), np.int32(self.numy_3c),
            np.int32(type_idx),
            d_j, d_k,
            self.d_data_3c, self.d_h_grids_3c, d_out
        )
        results = np.zeros((npoints, self.n_nz_3c_max), dtype=np.float32)
        cl.enqueue_copy(self.queue, results, d_out)
        return results

    def scanHamPiece3c_raw_batch(self, root, nz1, nz2, nz3, dRjs, dRks):
        """OpenCL analog of firecore_scanHamPiece3c_raw_batch (no rotation/Legendre)."""
        type_idx = self.species_triplet_map.get((root, nz1, nz2, nz3))
        if type_idx is None:
            print(f"[WARN] scanHamPiece3c_raw_batch missing triplet {root} ({nz1},{nz2},{nz3})")
            return None
        dRjs_np = np.ascontiguousarray(dRjs, dtype=np.float32)
        dRks_np = np.ascontiguousarray(dRks, dtype=np.float32)
        npoints = dRjs_np.shape[0]
        d_j = cl.Buffer(self.ctx, cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR, hostbuf=dRjs_np)
        d_k = cl.Buffer(self.ctx, cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR, hostbuf=dRks_np)
        out_size = npoints * 5 * self.n_nz_3c_max * 4
        d_out = cl.Buffer(self.ctx, cl.mem_flags.WRITE_ONLY, out_size)
        self.prg.scan_3c_raw_points(
            self.queue, (npoints,), None,
            np.int32(npoints), np.int32(self.n_nz_3c_max),
            np.int32(self.numx_3c), np.int32(self.numy_3c),
            np.int32(type_idx),
            d_j, d_k,
            self.d_data_3c, self.d_h_grids_3c, d_out
        )
        results = np.zeros((npoints, 5, self.n_nz_3c_max), dtype=np.float32)
        cl.enqueue_copy(self.queue, results, d_out)
        return results

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

        EQ2 = 14.39975
        S_blocks = run_2c(pairs_S)
        T_blocks = run_2c(pairs_T)
        
        # Aggregate all potential Vna 2-center components
        Vna_total = np.zeros((n_pairs, 4, 4), dtype=np.float32)
        vna_roots = ['vna', 'vna_atom_00', 'vna_ontopl_00', 'vna_ontopr_00']
        for root in vna_roots:
            pairs_v = []
            for idx, (i, j) in enumerate(neighbors):
                t = self.species_pair_map.get((root, species[i], species[j]))
                if t is not None:
                    pairs_v.append((i, j, t, idx))
            if pairs_v:
                Vna_total += run_2c(pairs_v)
        
        Vna_blocks = Vna_total * EQ2
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
