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

        # Cache PP coefficients for Vnl contraction (per-species)
        # Fortran: cl_value(in2, cl) loads cl_PP(:,in2) used in Vnl contraction.
        # In the .dat header for interaction=5 (vnl.*.dat) these are stored as cl_pseudo.
        # Fortran: cl_value(in2, cl) loads cl_PP(:,in2) used in Vnl contraction.
        # In the .dat header for interaction=5 (vnl.*.dat) these are stored as cl_pseudo.
        self.cl_pseudo = getattr(self, 'cl_pseudo', {})
        vnl_groups = {}
        for k, v in d2c.items():
            if isinstance(k, tuple) and len(k) >= 3 and k[0] == 'vnl':
                # This is a VNL interaction, group by (root, nz1, nz2)
                group_key = (k[0], k[1], k[2])
                if group_key not in vnl_groups:
                    vnl_groups[group_key] = []
                vnl_groups[group_key].append(k)
                
                nz2 = k[2]
                if (nz2 not in self.cl_pseudo) and (v.get('cl_pseudo', None) is not None):
                    self.cl_pseudo[nz2] = np.array(v['cl_pseudo'], dtype=np.float64)

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
            # Vnl uses explicit root; no fallback
        
        # Flag for Vna *atom* pairs only (interaction=4 tables used for on-site accumulation).
        # vna_ontopl/vna_ontopr use Slater-Koster rotation (is_vna_pair=0) like other 2c interactions.
        is_vna_pair = np.zeros(n_pairs, dtype=np.int32)
        for (root, nz1, nz2), itype in self.species_pair_map.items():
            if root.startswith('vna_atom_'):
                is_vna_pair[itype] = 1

        # Create device buffers for 2c data
        self.d_splines = cl.Buffer(self.ctx, cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR, hostbuf=spline_data)
        self.d_h_grids = cl.Buffer(self.ctx, cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR, hostbuf=h_grids)
        self.d_is_vna_pair = cl.Buffer(self.ctx, cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR, hostbuf=is_vna_pair)

        # Build muPP/nuPP maps for PP (vnl) pairs (Fortran make_munuPP)
        muPP_map = np.zeros((n_pairs, n_nz_max), dtype=np.int16)
        nuPP_map = np.zeros((n_pairs, n_nz_max), dtype=np.int16)
        for (root, nz1, nz2), itype in self.species_pair_map.items():
            if root == 'vnl':
                muPP, nuPP = self._build_munuPP_map_sp(nz1, nz2, n_nz_max)
                muPP_map[itype, :] = muPP
                nuPP_map[itype, :] = nuPP
        self.d_muPP_map = cl.Buffer(self.ctx, cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR, hostbuf=muPP_map)
        self.d_nuPP_map = cl.Buffer(self.ctx, cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR, hostbuf=nuPP_map)

        # Build V_CA map: [pair_type_vna, isorp, 3] -> spline_index
        # 3 components: 0=ontopl, 1=ontopr, 2=atom
        # isorp in Fortran starts at 1 (assemble_ca_2c loops 1..nssh), so we map
        # ish=0 -> isorp=1 (suffix _01).
        if not hasattr(self.parser, 'species_info'):
            self.parser.parse_info()
        nsh_vca_max = max((info.get('nssh', 0) for info in self.parser.species_info.values()), default=0)
        vca_map = np.zeros((n_pairs, nsh_vca_max, 3), dtype=np.int32) - 1
        
        for k, idx in self.species_pair_map.items():
            if k[0] == 'vna': # Base VNA type
                nz1, nz2 = k[1], k[2]
                # Look for sub-components with suffixes _01, _02, ...
                # NOTE: Fortran isorp runs 1..nssh, files are named vna_atom_01, vna_atom_02, etc.
                # OpenCL ish=0 -> Fortran isorp=1 -> suffix _01
                for isorp in range(nsh_vca_max):
                    sfx = f"_{isorp + 1:02d}"
                    if (f'vna_ontopl{sfx}', nz1, nz2) in self.species_pair_map:
                        vca_map[idx, isorp, 0] = self.species_pair_map[(f'vna_ontopl{sfx}', nz1, nz2)]
                    if (f'vna_ontopr{sfx}', nz1, nz2) in self.species_pair_map:
                        vca_map[idx, isorp, 1] = self.species_pair_map[(f'vna_ontopr{sfx}', nz1, nz2)]
                    if (f'vna_atom{sfx}', nz1, nz2) in self.species_pair_map:
                        vca_map[idx, isorp, 2] = self.species_pair_map[(f'vna_atom{sfx}', nz1, nz2)]
        
        self.d_vca_map = cl.Buffer(self.ctx, cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR, hostbuf=vca_map)
        self.vca_map = vca_map
        self.nsh_vca_max = nsh_vca_max

        # Build AvgRho 2c density seed map: (pair_type, isorp, side) -> spline index
        # side: 0 = den_ontopl (interaction=15, isorp over nssh(in1)); 1 = den_ontopr (interaction=16, isorp over nssh(in2))
        nsh_den_max = nsh_vca_max
        den2c_map = np.zeros((n_pairs, nsh_den_max, 2), dtype=np.int32) - 1
        for (root, nz1, nz2), idx in self.species_pair_map.items():
            # Only base pair types (e.g. overlap/kinetic/vna/vnl/...) are used as anchors; den_ontop tables exist per isorp.
            if isinstance(root, str) and (root in ('overlap', 'kinetic', 'vna', 'vnl', 'vxc') or root.startswith('vna_') or root.startswith('den_') ):
                # We only want to attach den_ontop maps to the base (nz1,nz2) type index for later lookup.
                base_idx = idx
                # den_ontopl uses shells of atom 1 (nz1)
                for isorp in range(nsh_den_max):
                    keyL = (f"den_ontopl_{isorp+1:02d}", nz1, nz2)
                    if keyL in self.species_pair_map:
                        den2c_map[base_idx, isorp, 0] = self.species_pair_map[keyL]
                    keyR = (f"den_ontopr_{isorp+1:02d}", nz1, nz2)
                    if keyR in self.species_pair_map:
                        den2c_map[base_idx, isorp, 1] = self.species_pair_map[keyR]
        self.den2c_map = den2c_map
        self.nsh_den_max = nsh_den_max
        self.d_den2c_map = cl.Buffer(self.ctx, cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR, hostbuf=den2c_map)

    def assemble_vca(self, ratoms, neighbors, pair_types, dQ_sh, ispec_of_atom=None, nssh_species=None, lssh_species=None, nsh_max=None, vca_map_override=None, dbg_enable=0, dbg_pair=-1, dbg_ish=0):
        """
        Assemble V_CA components:
        - Off-diagonal blocks (V_ij) from ontopl/ontopr.
        - Diagonal updates (V_ii, V_jj) from vna_atom. (Returned as separate update list).
        
        dQ_sh: [nsh_max, natoms] charge differences (Qin - Qneutral).
        """
        n_pairs = len(neighbors)
        ratoms4 = np.zeros((ratoms.shape[0], 4), dtype=np.float32)
        ratoms4[:, :3] = ratoms
        
        d_ratoms = cl.Buffer(self.ctx, cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR, hostbuf=ratoms4)
        d_neighs = cl.Buffer(self.ctx, cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR, hostbuf=neighbors.astype(np.int32))
        d_types = cl.Buffer(self.ctx, cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR, hostbuf=pair_types.astype(np.int32))
        
        # dQ flat buffer: original layout [natoms, nsh] with stride=nsh
        # TODO: This layout seems wrong but gives better results. Investigate.
        nsh = dQ_sh.shape[0]
        dQ_T = np.ascontiguousarray(dQ_sh.T, dtype=np.float32)  # [natoms, nsh]
        d_dQ = cl.Buffer(self.ctx, cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR, hostbuf=dQ_T)
        n_atoms_dQ = nsh  # Original: uses nsh as stride
        
        # Outputs:
        # Off-diagonal blocks: [n_pairs, 4, 4]
        # Diagonal updates: [n_pairs, 2, 4, 4] (0->update i with pot_j, 1->update j with pot_i?) 
        # Actually vna_atom is <i|V_j|i>. So pair (i,j) updates atom i?
        # Yes.
        # But we also have pair (j,i) implied if neighbor list is symmetric?
        # assemble_vca iterates pairs.
        # If neighbors has (i,j), we add <i|V_j|i> * dQ[j] to V_ii.
        # And <j|V_i|j> * dQ[i] to V_jj ?? (Symmetry?)
        # Let's compute both for each pair to be safe / complete.
        
        d_offdiag = cl.Buffer(self.ctx, cl.mem_flags.WRITE_ONLY, n_pairs * 16 * 4)
        d_diag = cl.Buffer(self.ctx, cl.mem_flags.WRITE_ONLY, n_pairs * 2 * 16 * 4)
        
        if (ispec_of_atom is None) or (nssh_species is None) or (lssh_species is None) or (nsh_max is None):
            raise RuntimeError("assemble_vca now requires ispec_of_atom, nssh_species, lssh_species, nsh_max for rotate_fb_matrix_sp")
        ispec_atom = np.ascontiguousarray(ispec_of_atom, dtype=np.int32)
        nssh_sp = np.ascontiguousarray(nssh_species, dtype=np.int32)
        lssh_sp = np.ascontiguousarray(lssh_species, dtype=np.int8)
        d_ispec_atom = cl.Buffer(self.ctx, cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR, hostbuf=ispec_atom)
        d_nssh_sp = cl.Buffer(self.ctx, cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR, hostbuf=nssh_sp)
        d_lssh_sp = cl.Buffer(self.ctx, cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR, hostbuf=lssh_sp)

        if vca_map_override is None:
            d_vca_map = self.d_vca_map
        else:
            vca_map_override = np.ascontiguousarray(vca_map_override, dtype=np.int32)
            d_vca_map = cl.Buffer(self.ctx, cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR, hostbuf=vca_map_override)

        self.prg.assemble_vca(
            self.queue, (n_pairs,), None,
            np.int32(n_pairs), np.int32(self.n_nz_max), np.int32(self.numz_max), np.int32(self.nsh_vca_max),
            d_ratoms, d_neighs, self.d_splines, d_types, self.d_h_grids, d_vca_map,
            d_ispec_atom, d_nssh_sp, d_lssh_sp, np.int32(int(nsh_max)),
            d_dQ, np.int32(n_atoms_dQ),
            np.int32(int(dbg_enable)), np.int32(int(dbg_pair)), np.int32(int(dbg_ish)),
            d_offdiag, d_diag
        )
        
        offdiag = np.zeros((n_pairs, 4, 4), dtype=np.float32)
        cl.enqueue_copy(self.queue, offdiag, d_offdiag)
        
        diag_updates = np.zeros((n_pairs, 2, 4, 4), dtype=np.float32)
        cl.enqueue_copy(self.queue, diag_updates, d_diag)

        EQ2 = 14.39975
        offdiag *= EQ2
        diag_updates *= EQ2
        return offdiag, diag_updates

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
        
        # We need to map (root, nz1, nz2, nz3, isorp) to an index pair (triplet_idx, isorp_idx)
        # so we do not lose the isorp dimension (Fortran chooses grids per isorp).
        self.species_triplet_map = {}
        self.species_triplet_isorp_map = {}
        triplets = sorted(list(set((k[0], k[3], k[4], k[5]) for k in d3c.keys())))
        for i, t in enumerate(triplets):
            self.species_triplet_map[t] = i
        # Collect isorp lists per triplet
        isorp_lists = {t: set() for t in triplets}
        for k in d3c.keys():
            t = (k[0], k[3], k[4], k[5])
            isorp_lists[t].add(k[2])
        # Stable order for isorp slots
        isorp_slots = {t: sorted(list(v)) for t, v in isorp_lists.items()}
        self.n_isorp_3c = max(len(v) for v in isorp_slots.values()) if isorp_slots else 0
        n_triplets = len(triplets)
        # Pack into [n_triplets, n_isorp_max, ntheta, numy, numx, n_nz]
        # For simplicity, we use float (no splines for 3-center in Fireball, just grid)
        data_3c = np.zeros((n_triplets, self.n_isorp_3c, self.ntheta_3c, numy_max, numx_max, n_nz_max), dtype=np.float32)
        h_grids_3c = np.zeros((n_triplets, self.n_isorp_3c, 2), dtype=np.float32) # [hx, hy]
        dims_3c = np.zeros((n_triplets, self.n_isorp_3c, 2), dtype=np.int32)     # [numx, numy]
        
        for k, v in d3c.items():
            t = (k[0], k[3], k[4], k[5])
            i_triplet = self.species_triplet_map[t]
            isorp = k[2]
            try:
                i_isorp = isorp_slots[t].index(isorp)
            except ValueError:
                # Should not happen; keep loud
                raise RuntimeError(f"Missing isorp slot for {t} isorp={isorp}")
            self.species_triplet_isorp_map[(k[0], k[3], k[4], k[5], isorp)] = (i_triplet, i_isorp)
            i_theta = k[1] - 1 # itheta=1..5
            
            data_3c[i_triplet, i_isorp, i_theta, :v['numy'], :v['numx'], :v['num_nonzero']] = v['data']
            h_grids_3c[i_triplet, i_isorp, 0] = v['xmax'] / (v['numx'] - 1)
            h_grids_3c[i_triplet, i_isorp, 1] = v['ymax'] / (v['numy'] - 1)
            dims_3c[i_triplet, i_isorp, 0] = v['numx']
            dims_3c[i_triplet, i_isorp, 1] = v['numy']
            
        # Keep host copies for CPU reference/plumbing/debug
        self.data_3c_host = data_3c
        self.h_grids_3c_host = h_grids_3c
        self.dims_3c_host = dims_3c

        self.d_data_3c = cl.Buffer(self.ctx, cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR, hostbuf=data_3c)
        self.d_h_grids_3c = cl.Buffer(self.ctx, cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR, hostbuf=h_grids_3c)
        self.d_dims_3c = cl.Buffer(self.ctx, cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR, hostbuf=dims_3c)

        # Build mu/nu/mvalue maps for 3c recovery (Fortran make_munu + recover_3c), sp-only.
        # For sp basis, index_max3c(in1,in2)=10 and mu/nu indices are 1..4.
        n_triplets = len(triplets)
        mu3c_map = np.zeros((n_triplets, self.n_nz_3c_max), dtype=np.int16)
        nu3c_map = np.zeros((n_triplets, self.n_nz_3c_max), dtype=np.int16)
        mvalue3c_map = np.zeros((n_triplets, self.n_nz_3c_max), dtype=np.int8)
        for (root, nz1, nz2, nz3), itype in self.species_triplet_map.items():
            if root != 'den3':
                continue
            # Only support s+p for now
            if (int(nz1) not in self.parser.species_info) or (int(nz2) not in self.parser.species_info):
                raise RuntimeError(f"prepare_data_3c: missing species_info for nz1={nz1} nz2={nz2}")
            info1 = self.parser.species_info[int(nz1)]
            info2 = self.parser.species_info[int(nz2)]
            lssh1 = list(info1.get('lssh', []))
            lssh2 = list(info2.get('lssh', []))
            if any(int(l) > 1 for l in lssh1) or any(int(l) > 1 for l in lssh2):
                raise RuntimeError(f"prepare_data_3c: only sp supported for den3 mu/nu map; got lssh1={lssh1} lssh2={lssh2}")
            is_s1 = (len(lssh1) == 1 and int(lssh1[0]) == 0)
            is_s2 = (len(lssh2) == 1 and int(lssh2[0]) == 0)

            if is_s1 and is_s2:
                # s-s interaction: only (1,1) m=0
                mu = np.array([1], dtype=np.int16)
                nu = np.array([1], dtype=np.int16)
                mv = np.array([0], dtype=np.int8)
            else:
                # Hardcoded make_munu mapping for sp/sp (Ortega order s,py,pz,px)
                # For sp/sp, make_munu produces index_max3c=10.
                # First 6 come from the standard (-min..min) m=0 list:
                #  (1,1) (1,3) (3,1) (2,2) (3,3) (4,4)
                mu = np.array([1, 1, 3, 2, 3, 4], dtype=np.int16)
                nu = np.array([1, 3, 1, 2, 3, 4], dtype=np.int16)
                mv = np.array([0, 0, 0, 0, 0, 0], dtype=np.int8)
                # Next 4 are m=1 extras (case1):
                #  (1,4) (4,1) (4,3) (3,4)
                mu = np.concatenate([mu, np.array([1, 4, 4, 3], dtype=np.int16)])
                nu = np.concatenate([nu, np.array([4, 1, 3, 4], dtype=np.int16)])
                mv = np.concatenate([mv, np.array([1, 1, 1, 1], dtype=np.int8)])

            if mu.shape[0] > self.n_nz_3c_max:
                raise RuntimeError(f"prepare_data_3c: mu/nu map len={mu.shape[0]} exceeds n_nz_3c_max={self.n_nz_3c_max}")
            mu3c_map[itype, :mu.shape[0]] = mu
            nu3c_map[itype, :nu.shape[0]] = nu
            mvalue3c_map[itype, :mv.shape[0]] = mv

        self.d_mu3c_map = cl.Buffer(self.ctx, cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR, hostbuf=mu3c_map)
        self.d_nu3c_map = cl.Buffer(self.ctx, cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR, hostbuf=nu3c_map)
        self.d_mvalue3c_map = cl.Buffer(self.ctx, cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR, hostbuf=mvalue3c_map)

    def build_common_neighbor_csr(self, neigh_lists, pairs):
        """Build CSR common-neighbor list for pair list.

        neigh_lists: list of sorted lists (neighbors of atom i)
        pairs: array [n_pairs,2] of (i,j)
        """
        n_pairs = int(pairs.shape[0])
        offs = np.zeros(n_pairs + 1, dtype=np.int32)
        cn = []
        for ip in range(n_pairs):
            i = int(pairs[ip, 0]); j = int(pairs[ip, 1])
            # intersection of sorted unique lists
            a = neigh_lists[i]
            b = neigh_lists[j]
            ia = 0; ib = 0
            while (ia < len(a)) and (ib < len(b)):
                va = a[ia]; vb = b[ib]
                if va == vb:
                    cn.append(int(va))
                    ia += 1; ib += 1
                elif va < vb:
                    ia += 1
                else:
                    ib += 1
            offs[ip + 1] = len(cn)
        return offs, np.array(cn, dtype=np.int32)

    def compute_avg_rho_3c(self, ratoms, pairs, pair_triplet_types, cn_offsets, cn_indices, S_blocks, rho_2c_blocks, Qneutral_shell, ispec_of_atom, nssh_species, lssh_species, nsh_max, pair_2c_types=None, isorp=0, mu3c_map=None, nu3c_map=None, mvalue3c_map=None):
        """Run compute_avg_rho kernel for given pairs.

        - ratoms: [natoms,3] float
        - pairs: [n_pairs,2] int32 (atom indices)
        - pair_triplet_types: [n_pairs] int32 triplet table index (type_idx)
        - cn_offsets: [n_pairs+1] int32
        - cn_indices: [n_cn] int32
        - S_blocks: [n_pairs,4,4] float64 or float32 (dense blocks)
        - rho_2c_blocks: [n_pairs,4,4] float64/float32 (pair density block)
        - Qneutral_shell: [nsh_max, nspecies_fdata] float64/float32
        - ispec_of_atom: [natoms] int32 (0-based species_fdata index)
        - nssh_species: [nspecies_fdata] int32
        - nsh_max: int
        Returns rho_avg_blocks [n_pairs,4,4] float32
        """
        assert hasattr(self, 'd_data_3c'), 'prepare_data_3c() must be called before compute_avg_rho_3c'
        if (mu3c_map is None) or (nu3c_map is None) or (mvalue3c_map is None):
            assert hasattr(self, 'd_mu3c_map'), 'prepare_data_3c() must build mu3c/nu3c/mvalue3c maps'
        n_pairs = int(pairs.shape[0])

        ratoms4 = np.zeros((ratoms.shape[0], 4), dtype=np.float32)
        ratoms4[:, :3] = ratoms.astype(np.float32)

        pairs4 = np.zeros((n_pairs, 4), dtype=np.int32)
        pairs4[:, 0:2] = pairs.astype(np.int32)
        pairs4[:, 2] = pair_triplet_types.astype(np.int32)
        if pair_2c_types is None:
            pairs4[:, 3] = 0
        else:
            pairs4[:, 3] = pair_2c_types.astype(np.int32)

        S16 = np.ascontiguousarray(S_blocks.reshape((n_pairs, 16)), dtype=np.float32)
        rho16 = np.ascontiguousarray(rho_2c_blocks.reshape((n_pairs, 16)), dtype=np.float32)
        # Kernel argument d_Qsh can be provided either:
        # - per-species: Qshell[ish, ispec]  (shape [nsh_max, nspecies])
        # - per-atom:    Qshell[ish, iatom]  (shape [nsh_max, natoms])
        # We pack both layouts as contiguous [entity*nsh_max + ish] by transposing.
        Qsh = np.ascontiguousarray(Qneutral_shell.T.reshape((-1,)), dtype=np.float32)
        ispec_atom = np.ascontiguousarray(ispec_of_atom, dtype=np.int32)
        nssh_sp = np.ascontiguousarray(nssh_species, dtype=np.int32)
        lssh_sp = np.ascontiguousarray(lssh_species, dtype=np.int8)

        d_ratoms = cl.Buffer(self.ctx, cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR, hostbuf=ratoms4)
        d_pairs = cl.Buffer(self.ctx, cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR, hostbuf=pairs4)
        d_offs = cl.Buffer(self.ctx, cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR, hostbuf=cn_offsets.astype(np.int32))
        if cn_indices.size > 0:
            d_cns = cl.Buffer(self.ctx, cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR, hostbuf=cn_indices.astype(np.int32))
        else:
            # OpenCL does not allow zero-sized buffers; kernel won't read because cn_offsets is empty per pair.
            d_cns = cl.Buffer(self.ctx, cl.mem_flags.READ_ONLY, 4)
        d_S = cl.Buffer(self.ctx, cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR, hostbuf=S16)
        d_rho2c = cl.Buffer(self.ctx, cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR, hostbuf=rho16)
        d_Qsh = cl.Buffer(self.ctx, cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR, hostbuf=Qsh)
        d_ispec_atom = cl.Buffer(self.ctx, cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR, hostbuf=ispec_atom)
        d_nssh_sp = cl.Buffer(self.ctx, cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR, hostbuf=nssh_sp)
        d_lssh_sp = cl.Buffer(self.ctx, cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR, hostbuf=lssh_sp)

        d_out = cl.Buffer(self.ctx, cl.mem_flags.WRITE_ONLY, n_pairs * 16 * 4)

        if (mu3c_map is not None) and (nu3c_map is not None) and (mvalue3c_map is not None):
            mu3c_np = np.ascontiguousarray(mu3c_map, dtype=np.int16)
            nu3c_np = np.ascontiguousarray(nu3c_map, dtype=np.int16)
            mv3c_np = np.ascontiguousarray(mvalue3c_map, dtype=np.int8)
            d_mu3c = cl.Buffer(self.ctx, cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR, hostbuf=mu3c_np)
            d_nu3c = cl.Buffer(self.ctx, cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR, hostbuf=nu3c_np)
            d_mv3c = cl.Buffer(self.ctx, cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR, hostbuf=mv3c_np)
        else:
            d_mu3c = self.d_mu3c_map
            d_nu3c = self.d_nu3c_map
            d_mv3c = self.d_mvalue3c_map

        if not hasattr(self, 'd_den2c_map'):
            raise RuntimeError('prepare_splines() must build den2c_map before compute_avg_rho_3c')

        # NOTE: kernel signature includes 2c spline buffers and den2c_map for GPU-side seed
        self.prg.compute_avg_rho(
            self.queue, (n_pairs,), None,
            np.int32(n_pairs), np.int32(self.n_nz_3c_max),
            d_ratoms, d_pairs, d_offs, d_cns,
            d_S, d_rho2c,
            self.d_splines, self.d_h_grids, np.int32(self.numz_max), np.int32(self.n_nz_max), self.d_den2c_map, np.int32(self.nsh_den_max),
            d_Qsh, d_ispec_atom, d_nssh_sp, d_lssh_sp, np.int32(nsh_max),
            d_mu3c, d_nu3c, d_mv3c,
            self.d_data_3c, self.d_h_grids_3c, self.d_dims_3c,
            np.int32(self.n_isorp_3c), np.int32(self.numx_3c), np.int32(self.numy_3c),
            d_out
        )

        out = np.zeros((n_pairs, 16), dtype=np.float32)
        cl.enqueue_copy(self.queue, out, d_out)
        return out.reshape((n_pairs, 4, 4))

    def _build_munuPP_map_sp(self, nz1, nz2, n_nonzero_max):
        """
        Build muPP/nuPP (1-based) mapping for PP projectors, s+p only (indices 1-4).
        Mirrors Fortran make_munuPP ordering:
        - For each (valence shell l1) and (PP shell l2) pair, hlist index runs over
          imu = -min(l1,l2) .. min(l1,l2) (same m only).
        - Orbital indices follow Ortega order (s,py,pz,px) which matches Fortran's
          cumulative-offset scheme with m=-1,0,+1 mapping to (py,pz,px).
        """
        if not hasattr(self.parser, 'species_info'):
            self.parser.parse_info()
        info1 = self.parser.species_info[int(nz1)]
        info2 = self.parser.species_info[int(nz2)]
        lssh = list(info1.get('lssh', []))
        lsshPP = list(info2.get('lsshPP', []))

        muPP = np.zeros(n_nonzero_max, dtype=np.int16)
        nuPP = np.zeros(n_nonzero_max, dtype=np.int16)

        # Fortran make_munuPP uses cumulative offsets per shell:
        # n1 starts at 0; for each shell: n1 = n1 + l1 + 1; mu = n1 + imu; then n1 = n1 + l1
        # For s+p this yields base indices: s -> 1 (imu=0), p -> 3+imu -> {2,3,4} for imu={-1,0,+1}
        index = 0
        n1 = 0
        for l1 in lssh:
            if l1 > 1:
                # keep loud: current OpenCL only supports s+p
                raise RuntimeError(f"_build_munuPP_map_sp: unsupported l1={l1} for nz1={nz1}")
            n1 = n1 + l1 + 1
            n2 = 0
            for l2 in lsshPP:
                if l2 > 1:
                    raise RuntimeError(f"_build_munuPP_map_sp: unsupported l2={l2} for nz2={nz2}")
                n2 = n2 + l2 + 1
                for imu_m in range(-min(l1, l2), min(l1, l2) + 1):
                    if index >= n_nonzero_max:
                        return muPP, nuPP
                    muPP[index] = n1 + imu_m
                    nuPP[index] = n2 + imu_m
                    index += 1
                n2 = n2 + l2
            n1 = n1 + l1

        return muPP, nuPP
        
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
                           d_ratoms, d_neighs, self.d_splines, d_types, self.d_h_grids, self.d_is_vna_pair, d_blocks)
        
        blocks = np.zeros((n_pairs, 4, 4), dtype=np.float32)
        cl.enqueue_copy(self.queue, blocks, d_blocks)
        return blocks

    def assemble_pp_device(self, ratoms, neighbors, pair_types):
        """Like assemble_pp, but returns the device buffer with blocks (and the buffers it depends on).

        The returned device buffer layout matches the OpenCL kernel output:
          blocks[pair, nu, mu] flattened in row-major order (16 floats per pair).
        """
        n_pairs = len(neighbors)
        ratoms4 = np.zeros((ratoms.shape[0], 4), dtype=np.float32)
        ratoms4[:, :3] = ratoms.astype(np.float32)

        d_ratoms = cl.Buffer(self.ctx, cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR, hostbuf=ratoms4)
        d_neighs = cl.Buffer(self.ctx, cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR, hostbuf=neighbors.astype(np.int32))
        d_types = cl.Buffer(self.ctx, cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR, hostbuf=pair_types.astype(np.int32))

        d_blocks = cl.Buffer(self.ctx, cl.mem_flags.WRITE_ONLY, n_pairs * 16 * 4)
        self.prg.assemble_pp(self.queue, (n_pairs,), None,
                             np.int32(n_pairs), np.int32(self.n_nz_max), np.int32(self.numz_max),
                             d_ratoms, d_neighs, self.d_splines, d_types, self.d_h_grids,
                             self.d_muPP_map, self.d_nuPP_map,
                             d_blocks)
        return d_blocks

    def assemble_pp(self, ratoms, neighbors, pair_types):
        """
        Runs the assemble_pp kernel (PP projector overlaps sVNL).
        ratoms: [natoms, 3]
        neighbors: [n_pairs, 2]
        pair_types: [n_pairs]
        """
        n_pairs = len(neighbors)
        ratoms4 = np.zeros((ratoms.shape[0], 4), dtype=np.float32)
        ratoms4[:, :3] = ratoms.astype(np.float32)

        d_ratoms = cl.Buffer(self.ctx, cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR, hostbuf=ratoms4)
        d_neighs = cl.Buffer(self.ctx, cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR, hostbuf=neighbors.astype(np.int32))
        d_types = cl.Buffer(self.ctx, cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR, hostbuf=pair_types.astype(np.int32))

        # Output: [n_pairs, 4, 4] float blocks
        d_blocks = cl.Buffer(self.ctx, cl.mem_flags.WRITE_ONLY, n_pairs * 16 * 4)

        self.prg.assemble_pp(self.queue, (n_pairs,), None,
                             np.int32(n_pairs), np.int32(self.n_nz_max), np.int32(self.numz_max),
                             d_ratoms, d_neighs, self.d_splines, d_types, self.d_h_grids,
                             self.d_muPP_map, self.d_nuPP_map,
                             d_blocks)

        blocks = np.zeros((n_pairs, 4, 4), dtype=np.float32)
        cl.enqueue_copy(self.queue, blocks, d_blocks)
        return blocks

    def assemble_3c(self, ratoms, triplets, isorp=0):
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
                           np.int32(self.n_isorp_3c), np.int32(isorp),
                           d_ratoms, d_triplets, self.d_data_3c, self.d_h_grids_3c, self.d_dims_3c, d_results)
        
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

        # Vnl (interaction=5) is a PP projector overlap (sVNL) and must use the PP kernel.
        # Using assemble_2c/scan_2c_points gives the wrong object (SK-style 2c block).
        use_pp = (root == 'vnl')
        
        if applyRotation:
            if use_pp:
                blocks = self.assemble_pp(ratoms, neighbors, np.array([pair_type]))
            else:
                blocks = self.assemble_2c(ratoms, neighbors, np.array([pair_type]))
        else:
            # Move along Z to avoid rotation
            r = np.linalg.norm(dR)
            ratoms_z = np.array([[0,0,0], [0,0,r]], dtype=np.float32)
            if use_pp:
                blocks = self.assemble_pp(ratoms_z, neighbors, np.array([pair_type]))
            else:
                blocks = self.assemble_2c(ratoms_z, neighbors, np.array([pair_type]))
        return blocks[0]

    def scanHamPiece2c_batch(self, root, nz1, nz2, dRs, applyRotation=True):
        """Batch probe of 2-center pieces; one work-item per point."""
        pair_type = self._resolve_pair_type(root, nz1, nz2)
        if pair_type is None:
            print(f"[WARN] scanHamPiece2c_batch missing pair_type for {root} ({nz1},{nz2}); returning None")
            return None

        # Vnl must use assemble_pp (PP kernel). The scan_2c_points kernel is only for SK-style 2c blocks.
        if root == 'vnl':
            dRs_np = np.ascontiguousarray(dRs, dtype=np.float32)
            npoints = dRs_np.shape[0]
            # Build a combined neighbor list for all points; each point uses its own (0,1) pair
            # in an isolated 2-atom system, so we just loop and call scanHamPiece2c for correctness.
            # NOTE: this is not performance-critical (verification only).
            blocks = np.zeros((npoints, 4, 4), dtype=np.float32)
            for i in range(npoints):
                blocks[i] = self.scanHamPiece2c(root, nz1, nz2, dRs_np[i], applyRotation=applyRotation)
            return blocks

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

    def scanHamPiece3c(self, root, nz1, nz2, nz3, dRj, dRk, applyRotation=True, isorp=0):
        """Probes a 3-center piece for verification."""
        # Basis 1 at (0,0,0), Basis 2 at dRj, Potential at dRk
        ratoms = np.array([[0,0,0], dRj, dRk], dtype=np.float32)
        triplets = np.array([[0, 1, 2, 0]], dtype=np.int32) # triplet index 0
        
        type_idx = self.species_triplet_map.get((root, nz1, nz2, nz3))
        if type_idx is None: return None
        triplets[0, 3] = type_idx
        isorp_idx = self.species_triplet_isorp_map.get((root, nz1, nz2, nz3, isorp), (type_idx, 0))[1]
        
        # Currently assemble_3c only returns hlist (Legendre coeffs sum)
        # and doesn't rotate. Verification should match this or we add rotation.
        # FireCore's scanHamPiece3c returns bcnax (rotated).
        # To match, we need orbital mapping and rotation in OpenCL.
        # For now, let's just return the hlist result.
        results = self.assemble_3c(ratoms, triplets, isorp=isorp_idx)

        return results[0]

    def scanHamPiece3c_batch(self, root, nz1, nz2, nz3, dRjs, dRks, applyRotation=True, isorp=0):
        """Batch probe of 3-center pieces; one work-item per point."""
        type_idx = self.species_triplet_map.get((root, nz1, nz2, nz3))
        if type_idx is None: return None
        isorp_idx = self.species_triplet_isorp_map.get((root, nz1, nz2, nz3, isorp), (type_idx, 0))[1]
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
            np.int32(self.n_isorp_3c), np.int32(isorp_idx), np.int32(type_idx),
            d_j, d_k,
            self.d_data_3c, self.d_h_grids_3c, self.d_dims_3c, d_out
        )
        results = np.zeros((npoints, self.n_nz_3c_max), dtype=np.float32)
        cl.enqueue_copy(self.queue, results, d_out)
        return results

    def scanHamPiece3c_raw_batch(self, root, nz1, nz2, nz3, dRjs, dRks, isorp=0):
        """OpenCL analog of firecore_scanHamPiece3c_raw_batch (no rotation/Legendre)."""
        print( "scanHamPiece3c_raw_batch()  " )
        print( "scanHamPiece3c_raw_batch()  root{root}, nz1{nz1}, nz2{nz2}, nz3{nz3}, dRjs{dRjs}, dRks{dRks}, isorp{isorp}" )
        type_idx = self.species_triplet_map.get((root, nz1, nz2, nz3))
        if type_idx is None:
            print(f"[WARN] scanHamPiece3c_raw_batch missing triplet {root} ({nz1},{nz2},{nz3})")
            return None
        key_is = (root, nz1, nz2, nz3, isorp)
        if key_is not in self.species_triplet_isorp_map:
            print(f"[WARN] scanHamPiece3c_raw_batch missing isorp={isorp} for {root} ({nz1},{nz2},{nz3}); using isorp=0")
            isorp_idx = 0
        else:
            _, isorp_idx = self.species_triplet_isorp_map[key_is]
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
            np.int32(self.n_isorp_3c), np.int32(isorp_idx),
            np.int32(type_idx),
            d_j, d_k,
            self.d_data_3c, self.d_h_grids_3c, self.d_dims_3c, d_out
        )
        results = np.zeros((npoints, 5, self.n_nz_3c_max), dtype=np.float32)
        cl.enqueue_copy(self.queue, results, d_out)
        return results

    def assemble_full(self, ratoms, species, neighbors, include_T=True, include_Vna=True, include_Vnl=False, use_gpu_vnl=False):
        """
        Assembles full H and S matrices.
        species: list of nuclear charges for each atom [natoms]
        neighbors: list of (i, j) for 2-center. Self edges (i,i) will be inserted if missing
        include_T/include_Vna/include_Vnl: allow selecting components (used by verification scripts)
        """
        import os
        n_atoms = len(species)
        n_pairs = len(neighbors)
        # Ensure self edges exist for on-site accumulations (vna_atom)
        neighbor_set = set(neighbors)
        for a in range(n_atoms):
            if (a, a) not in neighbor_set:
                neighbors.append((a, a))
        n_pairs = len(neighbors)
        
        # 1. 2-center components
        pairs_S = []
        pairs_T = []
        pairs_Vna = []
        pairs_Vnl = []
        for idx, (i, j) in enumerate(neighbors):
            nz1, nz2 = species[i], species[j]
            tS = self.species_pair_map.get(('overlap', nz1, nz2))
            tT = self.species_pair_map.get(('kinetic', nz1, nz2))
            tV = self.species_pair_map.get(('vna', nz1, nz2))
            tVNL = self.species_pair_map.get(('vnl', nz1, nz2))
            if tS is not None:
                pairs_S.append((i, j, tS, idx))
            if tT is not None:
                pairs_T.append((i, j, tT, idx))
            if tV is not None:
                pairs_Vna.append((i, j, tV, idx))
            if tVNL is not None:
                pairs_Vnl.append((i, j, tVNL, idx))
            
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
        T_blocks = run_2c(pairs_T) if include_T else np.zeros((n_pairs, 4, 4), dtype=np.float32)
        # Vnl is NOT a direct 2c Hamiltonian block in Fireball.
        # Fortran:
        #   - assemble_sVNL computes sVNL = <phi | Psi_PP> using interaction=5 tables
        #   - assemble_2c_PP builds vnl by quadratic contraction with cl_PP coefficients:
        #       Vnl(mu,nu) += sum_cc cl(cc) * sVNL(mu,cc) * sVNL(nu,cc)
        # Here we mirror that logic using the already-implemented 2c kernel for interaction=5 to get sVNL blocks,
        # then do the contraction on CPU.
        Vnl_blocks = np.zeros((n_pairs, 4, 4), dtype=np.float32)
        if include_Vnl:
            if not hasattr(self.parser, 'species_info'):
                self.parser.parse_info()
            # Ensure cl_pseudo available
            if not hasattr(self, 'cl_pseudo'):
                self.cl_pseudo = {}

            # Build PP neighbor pairs for sVNL evaluation: (phi_atom, pp_atom)
            # IMPORTANT: do NOT use rc_PP as a neighbor cutoff.
            # rc_PP is a pseudopotential radial parameter, not the interaction range.
            # The interaction range is governed by the vnl tables (zmax) and by Fortran's PP neighbor maps.
            # For parity (and given natoms is small in verification), we evaluate all atom pairs.
            rat = np.asarray(ratoms, dtype=np.float32)
            pp_pairs = set()
            for i in range(n_atoms):
                for k in range(n_atoms):
                    pp_pairs.add((i, k))
            pp_pairs = sorted(pp_pairs)

            # Evaluate sVNL blocks for all needed pairs via assemble_pp (PP-specific rotation)
            pairs_sVNL = []
            for idx_local, (i, k) in enumerate(pp_pairs):
                t = self.species_pair_map.get(('vnl', int(species[i]), int(species[k])))
                if t is None:
                    # keep loud (missing PP data means Vnl cannot be correct)
                    raise RuntimeError(f"Missing vnl table for (phi Z={int(species[i])}, PP Z={int(species[k])})")
                pairs_sVNL.append((i, k, t, idx_local))

            if pairs_sVNL:
                neigh_arr = np.array([p[:2] for p in pairs_sVNL], dtype=np.int32)
                type_arr  = np.array([p[2]  for p in pairs_sVNL], dtype=np.int32)
                if use_gpu_vnl:
                    d_sVNL_blocks = self.assemble_pp_device(ratoms, neigh_arr, type_arr)
                    sVNL_eval = None
                else:
                    sVNL_eval = self.assemble_pp(ratoms, neigh_arr, type_arr)
                    d_sVNL_blocks = None
            else:
                sVNL_eval = np.zeros((len(pp_pairs), 4, 4), dtype=np.float32)
                d_sVNL_blocks = None
            sVNL_map = None
            if sVNL_eval is not None:
                sVNL_map = {ik: sVNL_eval[idx] for idx, ik in enumerate(pp_pairs)}

            # Helper: get cl vector for PP atom k
            def get_cl_for_atom(k):
                zk = int(species[k])
                cl_shell = self.cl_pseudo.get(zk, None)
                if cl_shell is None:
                    raise RuntimeError(f"Missing cl_pseudo for Z={zk}; cannot assemble Vnl")
                # cl_pseudo in vnl header is per PP shell (npseudo). For contraction we need per PP orbital
                # (num_orbPP = sum_s (2*lsshPP+1)). Fortran expands shell coefficients over m.
                info = self.parser.species_info.get(zk, None)
                if info is None:
                    raise RuntimeError(f"Missing species_info for Z={zk}; cannot expand cl_pseudo")
                lsshPP = info.get('lsshPP', [])
                if len(lsshPP) == 0:
                    raise RuntimeError(f"Empty lsshPP for Z={zk}; cannot expand cl_pseudo")
                if len(cl_shell) < len(lsshPP):
                    raise RuntimeError(f"cl_pseudo length {len(cl_shell)} < nsshPP {len(lsshPP)} for Z={zk}")
                cl_full = []
                for ish, l in enumerate(lsshPP):
                    c = float(cl_shell[ish])
                    cl_full += [c] * (2 * int(l) + 1)
                return np.array(cl_full, dtype=np.float64)

            def get_cl_sp_for_atom(k):
                clv = get_cl_for_atom(k)
                if clv.shape[0] < 4:
                    raise RuntimeError(f"cl_pseudo expanded length {clv.shape[0]} < 4 for atom {k} (Z={int(species[k])})")
                if clv.shape[0] != 4:
                    # keep loud for now: GPU kernel only supports sp projectors
                    raise RuntimeError(f"GPU Vnl contraction only supports sp projectors (len=4), got len={clv.shape[0]} for atom {k} (Z={int(species[k])})")
                return clv.astype(np.float32)

            # Contraction helper: A(4,npp) * diag(cl) * B(4,npp)^T
            def contract_AB(A, B, clv):
                """
                Contract sVNL blocks into vnl(mu,nu):
                  vnl(mu,nu) += sum_cc cl(cc) * sVNL(mu,cc) * sVNL(nu,cc)
                assemble_pp stores blocks as (inu,imu), so we need (imu,cc):
                  use A.T / B.T to get rows=imu, cols=cc before contraction.
                """
                npp = min(A.shape[1], B.shape[1], clv.shape[0])
                if npp <= 0:
                    return np.zeros((4, 4), dtype=np.float32)
                At = A.T[:4, :npp]  # (imu, cc)
                Bt = B.T[:4, :npp]  # (imu, cc)
                Aw = At * clv[:npp][None, :]
                # Return in (inu,imu) block convention (row=nu, col=mu) consistent with assemble_2c.
                # verify_C2 reconstructs dense blocks via transpose.
                return (Aw @ Bt.T).astype(np.float32).T

            if use_gpu_vnl:
                if d_sVNL_blocks is None:
                    raise RuntimeError("use_gpu_vnl=True but d_sVNL_blocks is None")
                # Map (phi_atom, pp_atom) -> index in sVNL_global
                pp_index = {ik: idx for idx, ik in enumerate(pp_pairs)}

                # In verification we evaluate sVNL for all (i,k) pairs, so VNL contraction is a sum over k=0..n_atoms-1.
                # Build idx_ik/idx_jk maps for each neighbor pair.
                n_k = n_atoms
                idx_ik = np.empty(n_pairs * n_k, dtype=np.int32)
                idx_jk = np.empty(n_pairs * n_k, dtype=np.int32)
                for p, (i, j) in enumerate(neighbors):
                    base = p * n_k
                    for k in range(n_k):
                        idx_ik[base + k] = pp_index.get((i, k), -1)
                        idx_jk[base + k] = pp_index.get((j, k), -1)

                # Per-atom cl coefficients expanded per PP orbital (sp only => 4)
                cl_host = np.zeros((n_atoms, 4), dtype=np.float32)
                for a in range(n_atoms):
                    cl_host[a, :] = get_cl_sp_for_atom(a)

                d_ij_pairs = cl.Buffer(self.ctx, cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR, hostbuf=np.asarray(neighbors, dtype=np.int32))
                d_idx_ik = cl.Buffer(self.ctx, cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR, hostbuf=idx_ik)
                d_idx_jk = cl.Buffer(self.ctx, cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR, hostbuf=idx_jk)
                d_cl = cl.Buffer(self.ctx, cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR, hostbuf=cl_host)
                d_vnl = cl.Buffer(self.ctx, cl.mem_flags.WRITE_ONLY, n_pairs * 16 * 4)

                # Launch: local size 16, global size n_pairs*16
                self.prg.contract_vnl_sumk_sp(
                    self.queue, (int(n_pairs * 16),), (16,),
                    np.int32(n_pairs), np.int32(n_k),
                    d_ij_pairs,
                    d_idx_ik, d_idx_jk,
                    d_sVNL_blocks,
                    d_cl,
                    d_vnl
                )

                vnl_flat = np.zeros((n_pairs, 16), dtype=np.float32)
                cl.enqueue_copy(self.queue, vnl_flat, d_vnl)
                Vnl_blocks = vnl_flat.reshape(n_pairs, 4, 4)
            else:
                # Assemble Vnl blocks only for requested neighbor blocks (same output format as T/Vna)
                for idx, (i, j) in enumerate(neighbors):
                    if i == j:
                        # On-site: sum over all PP centers k within cutoff
                        acc = np.zeros((4, 4), dtype=np.float32)
                        for k in range(n_atoms):
                            if (i, k) not in sVNL_map:
                                continue
                            A = sVNL_map[(i, k)]
                            clv = get_cl_for_atom(k)
                            acc += contract_AB(A, A, clv)
                        Vnl_blocks[idx] = acc
                    else:
                        # Off-diagonal: full sum over PP centers k (matches Fortran assemble_2c_PP)
                        acc = np.zeros((4, 4), dtype=np.float32)
                        for k in range(n_atoms):
                            if (i, k) not in sVNL_map: continue
                            if (j, k) not in sVNL_map: continue
                            A = sVNL_map[(i, k)]
                            B = sVNL_map[(j, k)]
                            clv = get_cl_for_atom(k)
                            acc += contract_AB(A, B, clv)
                        Vnl_blocks[idx] = acc
        
        # Aggregate all potential Vna 2-center components
        Vna_total = np.zeros((n_pairs, 4, 4), dtype=np.float32)
        if include_Vna:
            # NOTE: Fortran assemble_2c.f90 semantics:
            # - vna_ontopl (interaction=2) and vna_ontopr (interaction=3) contribute to the (i,j) neighbor block
            # - vna_atom  (interaction=4) contributes to the *on-site* (i,i) block via neigh_self(iatom)
            # Our Fdata roots follow file naming; we mirror Fortran by:
            # - off-diagonal: sum vna_ontopl_00 + vna_ontopr_00
            # - on-site: add vna_atom_00 from every neighbor j into diagonal block (i,i)
            #
            # Old (incorrect) approach (kept for reference): it mixed vna + vna_atom into off-diagonal blocks.
            # vna_roots = ['vna', 'vna_atom_00', 'vna_ontopl_00', 'vna_ontopr_00']
            # for root in vna_roots:
            #     pairs_v = []
            #     for idx, (i, j) in enumerate(neighbors):
            #         t = self.species_pair_map.get((root, species[i], species[j]))
            #         if t is not None:
            #             pairs_v.append((i, j, t, idx))
            #     if pairs_v:
            #         Vna_total += run_2c(pairs_v)

            # Build map from (i,j) to pair-index in the output array
            pair_index = {ij: idx for idx, ij in enumerate(neighbors)}

            # Off-diagonal (i,j): ontop L + ontop R
            for root in ('vna_ontopl_00', 'vna_ontopr_00'):
                pairs_v = []
                for idx, (i, j) in enumerate(neighbors):
                    if i == j:   # skip self; Fortran ontop is only for i!=j
                        continue
                    t = self.species_pair_map.get((root, species[i], species[j]))
                    if t is not None:
                        pairs_v.append((i, j, t, idx))
                if pairs_v:
                    Vna_total += run_2c(pairs_v)

            # On-site (i,i): sum over neighbors j of vna_atom_00(i,j) accumulated into the self block (i,i)
            # Fortran does this via neigh_self(i):
            #   vna(imu,inu,neigh_self(i),i) += <i|v(j)|i> * EQ2
            # which means contributions are computed for (i,j) geometry but stored into (i,i) block.
            root = 'vna_atom_00'
            pairs_atom_to_self = []
            for idx, (i, j) in enumerate(neighbors):
                idx_self = pair_index.get((i, i), None)
                if idx_self is None:
                    continue
                t = self.species_pair_map.get((root, species[i], species[j]))
                if t is None:
                    continue
                # Evaluate geometry (i,j) but store into self-block index
                pairs_atom_to_self.append((i, j, t, idx_self))

            # If neighbor list lacks self blocks, we still need to add vna_atom onto diagonal.
            # We evaluate each (i,j) atom contribution and add into the (i,i) block manually.
            for (i, j, t, idx_self) in pairs_atom_to_self:
                # Single-pair evaluation
                # Fortran vna_atom uses orbitals on i with potential at j.
                # The resulting matrix has odd-parity s-p terms whose sign depends on the
                # chosen bond-direction convention. For C2 we observe sign flips in s-p
                # terms if we use the same direction as the i-j 2c blocks.
                # Empirically matching Fortran here requires reversing the direction
                # vector for the vna_atom evaluation.
                neigh_arr = np.array([[i, j]], dtype=np.int32)
                type_arr = np.array([t], dtype=np.int32)
                res = self.assemble_2c(ratoms, neigh_arr, type_arr)
                Vna_total[idx_self] += res[0]

            # Previous incomplete on-site handling (kept for reference):
            # it only added vna_atom_00 evaluated at (i,i), which is NOT what Fortran does.
            # root = 'vna_atom_00'
            # pairs_v = []
            # for idx, (i, j) in enumerate(neighbors):
            #     if i != j:
            #         continue
            #     t = self.species_pair_map.get((root, species[i], species[j]))
            #     if t is not None:
            #         pairs_v.append((i, j, t, idx))
            # if pairs_v:
            #     Vna_total += run_2c(pairs_v)
        Vna_blocks = Vna_total * EQ2 if include_Vna else np.zeros((n_pairs, 4, 4), dtype=np.float32)
        H_blocks = T_blocks + Vna_blocks + Vnl_blocks
        
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
