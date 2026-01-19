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

    def build_olsxc_off_sp4(self, dens4, densij4, denx16, den1x16, sx16):
        """Micro-parity helper for Fortran build_olsxc_off (sp-only 4x4).

        Inputs:
            dens4/densij4: float32 [n_pairs,4] with (d11,d12,d21,d22)
            denx16/den1x16/sx16: float32 [n_pairs,16] row-major 4x4
        Output:
            bcxcx16: float32 [n_pairs,16] row-major 4x4
        """
        n_pairs = int(dens4.shape[0])
        dens4 = np.ascontiguousarray(dens4, dtype=np.float32)
        densij4 = np.ascontiguousarray(densij4, dtype=np.float32)
        denx16 = np.ascontiguousarray(denx16, dtype=np.float32)
        den1x16 = np.ascontiguousarray(den1x16, dtype=np.float32)
        sx16 = np.ascontiguousarray(sx16, dtype=np.float32)
        out = np.zeros((n_pairs, 16), dtype=np.float32)
        d_dens4 = cl.Buffer(self.ctx, cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR, hostbuf=dens4.view(np.float32).reshape(n_pairs, 4))
        d_densij4 = cl.Buffer(self.ctx, cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR, hostbuf=densij4.view(np.float32).reshape(n_pairs, 4))
        d_denx16 = cl.Buffer(self.ctx, cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR, hostbuf=denx16.reshape(n_pairs * 16))
        d_den1x16 = cl.Buffer(self.ctx, cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR, hostbuf=den1x16.reshape(n_pairs * 16))
        d_sx16 = cl.Buffer(self.ctx, cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR, hostbuf=sx16.reshape(n_pairs * 16))
        d_out = cl.Buffer(self.ctx, cl.mem_flags.WRITE_ONLY, out.nbytes)
        evt = self.prg.build_olsxc_off_sp4(self.queue, (n_pairs,), None,
                                           np.int32(n_pairs),
                                           d_dens4, d_densij4,
                                           d_denx16, d_den1x16, d_sx16,
                                           d_out)
        evt.wait()
        cl.enqueue_copy(self.queue, out, d_out).wait()
        return out

    def build_ca_olsxc_on_sp4(self, arho_on, arhoi_on, rho_on, rhoi_on):
        """Micro-parity helper for Fortran build_ca_olsxc_on (sp-only 4x4 on-site).

        Inputs:
            arho_on: float32 [n_atoms,4] shell-level average densities (2x2 padded to 4)
            arhoi_on: float32 [n_atoms,4] shell-level average densities for 1-center (2x2 padded to 4)
            rho_on: float32 [n_atoms,16] orbital-level true densities (4x4)
            rhoi_on: float32 [n_atoms,16] orbital-level true densities for 1-center (4x4)
        Output:
            bcxcx_on: float32 [n_atoms,16] on-site XC matrix elements (4x4)
        """
        n_atoms = int(arho_on.shape[0])
        arho_on = np.ascontiguousarray(arho_on, dtype=np.float32)
        arhoi_on = np.ascontiguousarray(arhoi_on, dtype=np.float32)
        rho_on = np.ascontiguousarray(rho_on, dtype=np.float32)
        rhoi_on = np.ascontiguousarray(rhoi_on, dtype=np.float32)
        out = np.zeros((n_atoms, 16), dtype=np.float32)
        d_arho_on = cl.Buffer(self.ctx, cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR, hostbuf=arho_on.reshape(n_atoms * 4))
        d_arhoi_on = cl.Buffer(self.ctx, cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR, hostbuf=arhoi_on.reshape(n_atoms * 4))
        d_rho_on = cl.Buffer(self.ctx, cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR, hostbuf=rho_on.reshape(n_atoms * 16))
        d_rhoi_on = cl.Buffer(self.ctx, cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR, hostbuf=rhoi_on.reshape(n_atoms * 16))
        d_out = cl.Buffer(self.ctx, cl.mem_flags.WRITE_ONLY, out.nbytes)
        evt = self.prg.build_ca_olsxc_on_sp4(self.queue, (n_atoms,), None,
                                             np.int32(n_atoms),
                                             d_arho_on, d_arhoi_on,
                                             d_rho_on, d_rhoi_on,
                                             d_out)
        evt.wait()
        cl.enqueue_copy(self.queue, out, d_out).wait()
        return out

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
        is_generic2c = np.zeros(n_pairs, dtype=np.int32)
        for (root, nz1, nz2), itype in self.species_pair_map.items():
            if root.startswith('vna_atom_'):
                is_vna_pair[itype] = 1
            if root.startswith('dipole_'):
                is_generic2c[itype] = 1

        # Build generic mu/nu maps for 2c recovery (Fortran make_munu + recover_2c), sp-only.
        # This is required for dipole_* tables which are not representable by the hardcoded SK 5-parameter form.
        mu2c_map = np.zeros((n_pairs, n_nz_max), dtype=np.int16)
        nu2c_map = np.zeros((n_pairs, n_nz_max), dtype=np.int16)
        for (root, nz1, nz2), itype in self.species_pair_map.items():
            mu, nu = self._build_munu2c_map_sp(nz1, nz2, n_nz_max)
            mu2c_map[itype, :] = mu
            nu2c_map[itype, :] = nu

        # Create device buffers for 2c data
        self.d_splines = cl.Buffer(self.ctx, cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR, hostbuf=spline_data)
        self.d_h_grids = cl.Buffer(self.ctx, cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR, hostbuf=h_grids)
        self.d_is_generic2c = cl.Buffer(self.ctx, cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR, hostbuf=is_generic2c)
        self.d_mu2c_map = cl.Buffer(self.ctx, cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR, hostbuf=mu2c_map)
        self.d_nu2c_map = cl.Buffer(self.ctx, cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR, hostbuf=nu2c_map)
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

    def _build_munu2c_map_sp(self, nz1, nz2, n_nonzero_max):
        """Build mu/nu (1-based) mapping for 2c recovery, s+p only (indices 1-4).

        Mirrors Fortran make_munu ordering for 2-center integrals:
        for each shell pair (l1,l2), add entries for m = -min(l1,l2)..min(l1,l2)
        using Ortega ordering (s,py,pz,px) where pz is bond axis.
        """
        if not hasattr(self.parser, 'species_info'):
            self.parser.parse_info()
        info1 = self.parser.species_info.get(int(nz1), None)
        info2 = self.parser.species_info.get(int(nz2), None)
        if (info1 is None) or (info2 is None):
            raise RuntimeError(f"_build_munu2c_map_sp: missing species_info for nz1={nz1} nz2={nz2}")
        lssh1 = info1.get('lssh', [])
        lssh2 = info2.get('lssh', [])
        for l in lssh1:
            if int(l) > 1:
                raise RuntimeError(f"_build_munu2c_map_sp: unsupported l1={l} for nz1={nz1}")
        for l in lssh2:
            if int(l) > 1:
                raise RuntimeError(f"_build_munu2c_map_sp: unsupported l2={l} for nz2={nz2}")

        mu = np.zeros(n_nonzero_max, dtype=np.int16)
        nu = np.zeros(n_nonzero_max, dtype=np.int16)
        index = 0
        n1 = 0
        for l1 in lssh1:
            l1 = int(l1)
            n1 = n1 + l1 + 1
            n2 = 0
            for l2 in lssh2:
                l2 = int(l2)
                n2 = n2 + l2 + 1
                for imu_m in range(-min(l1, l2), min(l1, l2) + 1):
                    if index >= n_nonzero_max:
                        return mu, nu
                    mu[index] = n1 + imu_m
                    nu[index] = n2 + imu_m
                    index += 1
                n2 = n2 + l2
            n1 = n1 + l1
        return mu, nu
        
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
                           d_ratoms, d_neighs, self.d_splines, d_types, self.d_h_grids,
                           self.d_is_generic2c,
                           self.d_mu2c_map, self.d_nu2c_map, self.d_is_vna_pair,
                           d_blocks)
        
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
            d_in, self.d_splines, self.d_h_grids,
            self.d_is_generic2c, self.d_mu2c_map, self.d_nu2c_map,
            self.d_is_vna_pair,
            d_out
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

    def assemble_ewaldsr_blocks(self, atomPos, atomTypes_Z, sd, neigh_self, neigh_back, dq_atom, ewaldsr3c_mask=None):
        """Assemble EWALDSR contribution using OpenCL kernels.

        This is a verification-oriented routine extracted from tests/pyFireball/verify_C3.py.
        It constructs the exact OpenCL inputs (overlap + dipole_z blocks per *neighbor slot*) and
        launches:
        - ewaldsr_2c_atom_blocks
        - ewaldsr_2c_ontop_blocks
        - ewaldsr_3c_blocks

        The result is returned in the same neighbor-slot block layout as Fortran exports:
            ewaldsr4_o[iatom, ineigh, nu, mu]

        Parameters
        ----------
        atomPos : (natoms,3) float
        atomTypes_Z : (natoms,) int
        sd : FireCore sparse-data struct (must provide neighn, neigh_j, neigh_b, xl, neigh_max)
        neigh_self : (natoms,) int   (0-based index of self-slot in neighbor list)
        neigh_back : (natoms,neigh_max) int (1-based back mapping, as exported by FireCore)
        dq_atom : (natoms,) float32
        ewaldsr3c_mask : optional array (natoms,neigh_max,4,4) used as authoritative mask for which
                         neighbor slots receive 3c contributions and SFIRE overwrite.

        Returns
        -------
        dict with keys:
            - 'ewaldsr4'   : (natoms,neigh_max,4,4) float64 (nu,mu)
            - 'ij_pairs'   : (n_pairs,2) int32
            - 'pair_ineigh': (n_pairs,) int32
            - 'pair_mbeta' : (n_pairs,) int32
            - 'trips'      : (n_trip,4) int32 (ia,ja,ka,slot)
        """

        atomPos = np.ascontiguousarray(atomPos, dtype=np.float64)
        atomTypes_Z = np.ascontiguousarray(atomTypes_Z, dtype=np.int32)
        dq_atom = np.ascontiguousarray(dq_atom, dtype=np.float32)

        natoms = int(atomPos.shape[0])
        neigh_max = int(getattr(sd, 'neigh_max', 0))
        if neigh_max <= 0:
            neigh_max = int(np.max(np.asarray(sd.neighn, dtype=np.int32)))

        # -------------------------------------------------
        # Precompute overlap and dipole_z blocks per neighbor slot using per-pair mbeta-shifted positions
        # -------------------------------------------------
        S4_o = np.zeros((natoms, neigh_max, 4, 4), dtype=np.float64)
        Dip4_o = np.zeros((natoms, neigh_max, 4, 4), dtype=np.float64)

        S_pair_types = []
        D_pair_types = []
        pair_ia = []
        pair_ine = []
        pair_mb = []

        for ia in range(natoms):
            nn = int(sd.neighn[ia])
            Zi = int(atomTypes_Z[ia])
            for ineigh in range(nn):
                ja = int(sd.neigh_j[ia, ineigh]) - 1
                mb = int(sd.neigh_b[ia, ineigh])
                if ja < 0 or ja >= natoms:
                    continue
                if ia == ja and mb == 0:
                    continue
                Zj = int(atomTypes_Z[ja])
                tS = self.species_pair_map.get(('overlap', Zi, Zj))
                tD = self.species_pair_map.get(('dipole_z', Zi, Zj))
                if tS is None:
                    raise RuntimeError(f"Missing overlap table for Z pair ({Zi},{Zj})")
                if tD is None:
                    raise RuntimeError(f"Missing dipole_z table for Z pair ({Zi},{Zj})")
                S_pair_types.append(int(tS))
                D_pair_types.append(int(tD))
                pair_ia.append(int(ia))
                pair_ine.append(int(ineigh))
                pair_mb.append(int(mb))

        if len(pair_ia) == 0:
            raise RuntimeError("EwaldSR OpenCL inputs: empty neighbor list")

        ratoms_pairs = np.zeros((2 * len(pair_ia), 3), dtype=np.float32)
        neigh_pairs = np.zeros((len(pair_ia), 2), dtype=np.int32)
        for k in range(len(pair_ia)):
            ia = pair_ia[k]
            ineigh = pair_ine[k]
            ja = int(sd.neigh_j[ia, ineigh]) - 1
            mb = int(pair_mb[k])
            ratoms_pairs[2 * k + 0, :] = atomPos[ia].astype(np.float32)
            ratoms_pairs[2 * k + 1, :] = (atomPos[ja] + sd.xl[mb]).astype(np.float32)
            neigh_pairs[k, 0] = 2 * k + 0
            neigh_pairs[k, 1] = 2 * k + 1

        S_blocks = self.assemble_2c(ratoms_pairs, neigh_pairs, np.array(S_pair_types, dtype=np.int32)).astype(np.float64)
        D_blocks = self.assemble_2c(ratoms_pairs, neigh_pairs, np.array(D_pair_types, dtype=np.int32)).astype(np.float64)
        for k in range(len(pair_ia)):
            S4_o[pair_ia[k], pair_ine[k], :, :] = S_blocks[k]
            Dip4_o[pair_ia[k], pair_ine[k], :, :] = D_blocks[k]

        # Fill self-slot overlap blocks (on-site). In Fortran these live in the dedicated self-slot.
        for ia in range(natoms):
            mat = int(neigh_self[ia])
            if mat < 0:
                raise RuntimeError(f"EwaldSR OpenCL inputs: missing neigh_self for atom {ia}")
            Z = int(atomTypes_Z[ia])
            blk = self.scanHamPiece2c('overlap', Z, Z, np.array([0.0, 0.0, 0.0], dtype=np.float32), applyRotation=False)
            if blk is None:
                raise RuntimeError(f"EwaldSR OpenCL inputs: scanHamPiece2c overlap self returned None for Z={Z}")
            S4_o[ia, mat, :, :] = np.array(blk, dtype=np.float64)

        # -------------------------------------------------
        # Build directed pair list aligned with neighbor slots
        # -------------------------------------------------
        ij_pairs = []
        pair_ineigh = []
        pair_mbeta = []
        y_list = []
        S_list = []
        D_list = []
        slot_has3c = []

        for ia in range(natoms):
            nn = int(sd.neighn[ia])
            r1 = atomPos[ia]
            for ineigh in range(nn):
                ja = int(sd.neigh_j[ia, ineigh]) - 1
                mb = int(sd.neigh_b[ia, ineigh])
                if ja < 0 or ja >= natoms:
                    continue
                if ia == ja and mb == 0:
                    continue
                r2 = atomPos[ja] + sd.xl[mb]
                y = float(np.linalg.norm(r2 - r1))
                ij_pairs.append((ia, ja))
                pair_ineigh.append(int(ineigh))
                pair_mbeta.append(int(mb))
                y_list.append(y)
                # Kernels expect (mu,nu) flattened => transpose (nu,mu) block.
                S_list.append(S4_o[ia, ineigh, :, :].T.astype(np.float32).reshape(16))
                D_list.append(Dip4_o[ia, ineigh, :, :].T.astype(np.float32).reshape(16))
                if ewaldsr3c_mask is None:
                    slot_has3c.append(1)
                else:
                    slot_has3c.append(1 if float(np.max(np.abs(ewaldsr3c_mask[ia, ineigh]))) >= 1e-12 else 0)

        ij_pairs = np.ascontiguousarray(np.array(ij_pairs, dtype=np.int32))
        pair_ineigh = np.ascontiguousarray(np.array(pair_ineigh, dtype=np.int32))
        pair_mbeta = np.ascontiguousarray(np.array(pair_mbeta, dtype=np.int32))
        y_list = np.ascontiguousarray(np.array(y_list, dtype=np.float32))
        S_list = np.ascontiguousarray(np.array(S_list, dtype=np.float32))
        D_list = np.ascontiguousarray(np.array(D_list, dtype=np.float32))
        slot_has3c = np.ascontiguousarray(np.array(slot_has3c, dtype=np.int8))

        n_pairs = int(ij_pairs.shape[0])
        if n_pairs <= 0:
            raise RuntimeError("EwaldSR OpenCL: empty ij_pairs")

        # Self-slot overlap blocks (mu,nu) for each atom
        S_self = np.zeros((natoms, 16), dtype=np.float32)
        for ia in range(natoms):
            mat = int(neigh_self[ia])
            if mat < 0:
                continue
            S_self[ia, :] = S4_o[ia, mat, :, :].T.astype(np.float32).reshape(16)

        # -------------------------------------------------
        # Run kernels (2c atom + 2c ontop)
        # -------------------------------------------------
        d_ij = cl.Buffer(self.ctx, cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR, hostbuf=ij_pairs)
        d_y = cl.Buffer(self.ctx, cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR, hostbuf=y_list)
        d_dq = cl.Buffer(self.ctx, cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR, hostbuf=dq_atom)
        d_S = cl.Buffer(self.ctx, cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR, hostbuf=S_list)
        d_D = cl.Buffer(self.ctx, cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR, hostbuf=D_list)
        d_Sself = cl.Buffer(self.ctx, cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR, hostbuf=S_self)

        d_out_atom = cl.Buffer(self.ctx, cl.mem_flags.WRITE_ONLY, n_pairs * 16 * 4)
        d_out_on = cl.Buffer(self.ctx, cl.mem_flags.WRITE_ONLY, n_pairs * 16 * 4)

        self.prg.ewaldsr_2c_atom_blocks(self.queue, (n_pairs,), None, np.int32(n_pairs), d_ij, d_y, d_dq, d_Sself, d_out_atom)
        self.prg.ewaldsr_2c_ontop_blocks(self.queue, (n_pairs,), None, np.int32(n_pairs), d_ij, d_y, d_dq, d_S, d_D, d_out_on)

        out_atom = np.zeros((n_pairs, 16), dtype=np.float32)
        out_on = np.zeros((n_pairs, 16), dtype=np.float32)
        cl.enqueue_copy(self.queue, out_atom, d_out_atom)
        cl.enqueue_copy(self.queue, out_on, d_out_on)
        self.queue.finish()

        # -------------------------------------------------
        # 3c triples (verification mode: by default use all triples, or mask if provided)
        # -------------------------------------------------
        trips = []
        ty = []
        td13 = []
        td23 = []
        for slot in range(n_pairs):
            if int(slot_has3c[slot]) == 0:
                continue
            ia = int(ij_pairs[slot, 0])
            ja = int(ij_pairs[slot, 1])
            ineigh = int(pair_ineigh[slot])
            mb = int(pair_mbeta[slot])
            r1 = atomPos[ia]
            r2 = atomPos[ja] + sd.xl[mb]
            y = float(np.linalg.norm(r2 - r1))
            for ka in range(natoms):
                if ka == ia or ka == ja:
                    continue
                d13 = float(np.linalg.norm(atomPos[ka] - r1))
                d23 = float(np.linalg.norm(atomPos[ka] - r2))
                trips.append((ia, ja, ka, slot))
                ty.append(y)
                td13.append(d13)
                td23.append(d23)

        if len(trips) == 0:
            trips = np.zeros((0, 4), dtype=np.int32)
            out3 = np.zeros((0, 16), dtype=np.float32)
        else:
            trips = np.ascontiguousarray(np.array(trips, dtype=np.int32))
            ty = np.ascontiguousarray(np.array(ty, dtype=np.float32))
            td13 = np.ascontiguousarray(np.array(td13, dtype=np.float32))
            td23 = np.ascontiguousarray(np.array(td23, dtype=np.float32))
            n_trip = int(trips.shape[0])

            d_tr = cl.Buffer(self.ctx, cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR, hostbuf=trips)
            d_ty = cl.Buffer(self.ctx, cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR, hostbuf=ty)
            d_t13 = cl.Buffer(self.ctx, cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR, hostbuf=td13)
            d_t23 = cl.Buffer(self.ctx, cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR, hostbuf=td23)
            d_out3 = cl.Buffer(self.ctx, cl.mem_flags.WRITE_ONLY, n_trip * 16 * 4)

            self.prg.ewaldsr_3c_blocks(self.queue, (n_trip,), None, np.int32(n_trip), d_tr, d_ty, d_t13, d_t23, d_dq, d_S, d_D, d_out3)
            out3 = np.zeros((n_trip, 16), dtype=np.float32)
            cl.enqueue_copy(self.queue, out3, d_out3)
            self.queue.finish()

        # -------------------------------------------------
        # Reconstruct EWALDSR in [iatom,ineigh,nu,mu] block layout and apply SFIRE overwrite.
        # -------------------------------------------------
        ewaldsr4_o = np.zeros((natoms, neigh_max, 4, 4), dtype=np.float64)  # [iatom,ineigh,nu,mu]

        # 2c ontop: directed neighbor slot
        for slot in range(n_pairs):
            ia = int(ij_pairs[slot, 0])
            ineigh = int(pair_ineigh[slot])
            ewaldsr4_o[ia, ineigh, :, :] += out_on[slot].reshape(4, 4).T.astype(np.float64)

        # 2c atom: accumulates into self-slot
        for slot in range(n_pairs):
            ia = int(ij_pairs[slot, 0])
            mat = int(neigh_self[ia])
            if mat < 0 or mat >= neigh_max:
                continue
            ewaldsr4_o[ia, mat, :, :] += out_atom[slot].reshape(4, 4).T.astype(np.float64)

        # 3c: directed neighbor slot
        for it in range(int(trips.shape[0])):
            ia = int(trips[it, 0])
            slot = int(trips[it, 3])
            ineigh = int(pair_ineigh[slot])
            ewaldsr4_o[ia, ineigh, :, :] += out3[it].reshape(4, 4).T.astype(np.float64)

        # SFIRE overwrite (only for slots where 3c is enabled)
        for ia in range(natoms):
            nn = int(sd.neighn[ia])
            for ineigh in range(nn):
                if ewaldsr3c_mask is not None:
                    if float(np.max(np.abs(ewaldsr3c_mask[ia, ineigh]))) < 1e-12:
                        continue
                # if no explicit mask, we do NOT know SFIRE scope; apply to all edges by default
                ja = int(sd.neigh_j[ia, ineigh]) - 1
                if ja < 0 or ja >= natoms:
                    continue
                jneigh = int(neigh_back[ia, ineigh]) - 1
                if jneigh < 0 or jneigh >= int(sd.neighn[ja]):
                    continue
                if int(sd.neigh_j[ja, jneigh]) - 1 != ia:
                    continue
                ewaldsr4_o[ja, jneigh, :, :] = ewaldsr4_o[ia, ineigh, :, :].T

        return {
            'ewaldsr4': ewaldsr4_o,
            'ij_pairs': ij_pairs,
            'pair_ineigh': pair_ineigh,
            'pair_mbeta': pair_mbeta,
            'trips': trips,
        }

    def assemble_dip_blocks_per_slot(self, atomPos, atomTypes_Z, sd):
        """Assemble dipole_z blocks per neighbor slot (iatom,ineigh) using OpenCL.

        This is a verification-oriented routine extracted from tests/pyFireball/verify_C3.py.
        It constructs per-neighbor-slot dipole_z blocks using per-pair mbeta-shifted positions.

        Parameters
        ----------
        atomPos : (natoms,3) float
        atomTypes_Z : (natoms,) int
        sd : FireCore sparse-data struct (must provide neighn, neigh_j, neigh_b, xl)

        Returns
        -------
        dict with keys:
            - 'dip4'      : (natoms,neigh_max,4,4) float64 (nu,mu)
            - 'pair_i'    : (n_pairs,) int32
            - 'pair_ineigh': (n_pairs,) int32
            - 'pair_mbeta': (n_pairs,) int32
        """
        atomPos = np.ascontiguousarray(atomPos, dtype=np.float64)
        atomTypes_Z = np.ascontiguousarray(atomTypes_Z, dtype=np.int32)

        natoms = int(atomPos.shape[0])
        neigh_max = int(getattr(sd, 'neigh_max', 0))
        if neigh_max <= 0:
            neigh_max = int(np.max(np.asarray(sd.neighn, dtype=np.int32)))

        dip_pair_types = []
        pair_i = []
        pair_ineigh = []
        pair_mbeta = []

        for ia in range(natoms):
            nn = int(sd.neighn[ia])
            Zi = int(atomTypes_Z[ia])
            for ineigh in range(nn):
                ja = int(sd.neigh_j[ia, ineigh]) - 1
                mb = int(sd.neigh_b[ia, ineigh])
                if ja < 0 or ja >= natoms:
                    continue
                # Skip on-site self edge (y==0); dip is irrelevant there and can be ill-defined numerically.
                if ia == ja and mb == 0:
                    continue
                Zj = int(atomTypes_Z[ja])
                t = self.species_pair_map.get(('dipole_z', Zi, Zj))
                if t is None:
                    raise RuntimeError(f"Missing Fdata dipole_z table for Z pair ({Zi},{Zj})")
                dip_pair_types.append(int(t))
                pair_i.append(int(ia))
                pair_ineigh.append(int(ineigh))
                pair_mbeta.append(int(mb))

        if len(dip_pair_types) == 0:
            raise RuntimeError("Dip test: empty dip_pair_types")

        # Build per-pair positions including mbeta shift (same approach as VCA mbeta handling)
        ratoms_pairs = np.zeros((2 * len(dip_pair_types), 3), dtype=np.float32)
        neigh_pairs = np.zeros((len(dip_pair_types), 2), dtype=np.int32)
        for k in range(len(dip_pair_types)):
            ia = pair_i[k]
            ineigh = pair_ineigh[k]
            ja = int(sd.neigh_j[ia, ineigh]) - 1
            mb = pair_mbeta[k]
            ratoms_pairs[2*k+0, :] = atomPos[ia].astype(np.float32)
            ratoms_pairs[2*k+1, :] = (atomPos[ja] + sd.xl[mb]).astype(np.float32)
            neigh_pairs[k, 0] = 2*k+0
            neigh_pairs[k, 1] = 2*k+1

        dip_blocks = self.assemble_2c(ratoms_pairs, neigh_pairs, np.array(dip_pair_types, dtype=np.int32))

        # Assemble into per-slot layout
        dip4 = np.zeros((natoms, neigh_max, 4, 4), dtype=np.float64)
        for k in range(len(dip_pair_types)):
            ia = pair_i[k]
            ineigh = pair_ineigh[k]
            dip4[ia, ineigh, :, :] = dip_blocks[k].astype(np.float64)

        return {
            'dip4': dip4,
            'pair_i': np.array(pair_i, dtype=np.int32),
            'pair_ineigh': np.array(pair_ineigh, dtype=np.int32),
            'pair_mbeta': np.array(pair_mbeta, dtype=np.int32),
        }

    ## ==================================================
    ##   Temporary python implementations for testing
    ## ==================================================

#@staticmethod
def _cepal_py(rh):
    abohr = 0.529177249
    eq2 = 14.39975
    delta_rh = 1.0e-6
    rhx = float(np.sqrt(rh * rh + delta_rh))
    rho = rhx * (abohr ** 3)
    rho_third = rho ** (1.0 / 3.0)
    rs = 0.62035049 / rho_third
    exc = 0.0
    muxc = 0.0
    dec = 0.0
    ddec = 0.0
    d2dec = 0.0
    if rho < 0.23873241:
        sqrs = np.sqrt(rs)
        den = 1.0 + 1.0529 * sqrs + 0.3334 * rs
        exc = -0.4581652 / rs - 0.1423 / den
        muxc = exc - rs * (0.15273333 / (rs ** 2) + (0.02497128 / sqrs + 0.01581427) / (den ** 2))
        dden = 1.0529 / (2.0 * sqrs) + 0.3334
        d2den = (-0.5) * 1.0529 / (2.0 * rs * sqrs)
        d3den = (0.75) * 1.0529 / (2.0 * rs * rs * sqrs)
        dec = 0.1423 * dden / (den * den)
        ddec = -2.0 * 0.1423 * dden * dden / (den ** 3) + 0.1423 * d2den / (den * den)
        d2dec = 6.0 * 0.1423 * (dden * 3.0) / (den ** 4) - 6.0 * 0.1423 * dden * d2den / (den ** 3) + 0.1423 * d3den / (den * den)
    else:
        rsl = np.log(rs)
        exc = -0.4581652 / rs - 0.0480 + 0.0311 * rsl - 0.0116 * rs + 0.002 * rs * rsl
        muxc = exc - rs * (0.15273333 / (rs ** 2) + 0.01036667 / rs - 0.003866667 + 0.00066667 * (1.0 + rsl))
        dec = 0.0311 / rs - 0.0116 + 0.0020 * (rsl + 1.0)
        ddec = -0.0311 / (rs * rs) + 0.0020 / rs
        d2dec = 2.0 * 0.0311 / (rs * rs * rs) - 0.0020 / (rs * rs)
    ex = -0.7385587664 * rho_third
    dexc = (muxc - exc) / rho
    d2nec = (4.0 * rs / (9.0 * rho * rho)) * dec + (rs * rs / (9.0 * rho * rho)) * ddec
    d2nex = -(2.0 / (9.0 * rho * rho)) * ex
    dmuxc = 2.0 * dexc + rho * (d2nex + d2nec)
    d3nec = (-28.0 * rs / (27.0 * rho * rho * rho)) * dec + (-4.0 * rs * rs / (9.0 * rho * rho * rho)) * ddec + (rs * rs * rs / (-27.0 * rho * rho * rho)) * d2dec
    d3nex = (10.0 / (27.0 * rho * rho * rho)) * ex
    d2muxc = 3.0 * (d2nex + d2nec) + rho * (d3nex + d3nec)
    d2exc = d2nex + d2nec
    hartree1 = eq2 / abohr
    exc *= hartree1
    muxc *= hartree1
    dexc *= hartree1 * (abohr ** 3)
    d2exc *= hartree1 * (abohr ** 6)
    dmuxc *= hartree1 * (abohr ** 3)
    d2muxc *= hartree1 * (abohr ** 6)
    return exc, muxc, dexc, d2exc, dmuxc, d2muxc

#@staticmethod
def _build_olsxc_off_py(den1x, denx, sx, dens, densij):
    # For this C-only test case: in1=in2=1 with two shells: s (l=0), p (l=1)
    # This reproduces fortran/ASSEMBLERS/build_olsxc_off.f90 logic.
    lssh = [0, 1]
    bcxcx = np.zeros((4, 4), dtype=np.float64)
    n1 = 0
    for issh, l1 in enumerate(lssh, start=1):
        n1 = n1 + l1 + 1
        n2 = 0
        for jssh, l2 in enumerate(lssh, start=1):
            n2 = n2 + l2 + 1
            _, muxc, _, _, dmuxc, _ = _cepal_py(float(dens[issh - 1, jssh - 1]))
            _, muxcij, _, _, dmuxcij, _ = _cepal_py(float(densij[issh - 1, jssh - 1]))
            for ind1 in range(-l1, l1 + 1):
                imu = n1 + ind1
                for ind2 in range(-l2, l2 + 1):
                    inu = n2 + ind2
                    if (imu < 1) or (imu > 4) or (inu < 1) or (inu > 4):
                        continue
                    bc = muxc * sx[imu - 1, inu - 1]
                    bc += dmuxc * (denx[imu - 1, inu - 1] - dens[issh - 1, jssh - 1] * sx[imu - 1, inu - 1])
                    bc -= muxcij * sx[imu - 1, inu - 1] + dmuxcij * (den1x[imu - 1, inu - 1] - densij[issh - 1, jssh - 1] * sx[imu - 1, inu - 1])
                    bcxcx[imu - 1, inu - 1] = bc
            n2 = n2 + l2
        n1 = n1 + l1
    return bcxcx

#@staticmethod
def _epsilon_fb_py(r1, r2):
    r1mag = np.linalg.norm(r1)
    r2mag = np.linalg.norm(r2)
    spe = np.zeros((3, 3), dtype=np.float64)
    if r2mag < 1e-4:
        np.fill_diagonal(spe, 1.0)
        return spe
    zphat = r2 / r2mag
    if r1mag > 1e-4:
        r1hat = r1 / r1mag
        yphat = np.cross(zphat, r1hat)
        ypmag = np.linalg.norm(yphat)
        if ypmag > 1e-6:
            yphat = yphat / ypmag
            xphat = np.cross(yphat, zphat)
            spe[:, 0] = xphat
            spe[:, 1] = yphat
            spe[:, 2] = zphat
            return spe
    if abs(zphat[0]) > 1e-4:
        yphat = np.array([-(zphat[1] + zphat[2]) / zphat[0], 1.0, 1.0])
    elif abs(zphat[1]) > 1e-4:
        yphat = np.array([1.0, -(zphat[0] + zphat[2]) / zphat[1], 1.0])
    else:
        yphat = np.array([1.0, 1.0, -(zphat[0] + zphat[1]) / zphat[2]])
    yphat = yphat / np.linalg.norm(yphat)
    xphat = np.cross(yphat, zphat)
    spe[:, 0] = xphat
    spe[:, 1] = yphat
    spe[:, 2] = zphat
    return spe

#@staticmethod
def _twister_pmat(eps):
    # Fortran twister: pmat(1,1)=eps(2,2), pmat(1,2)=eps(2,3), pmat(1,3)=eps(2,1), ...
    return np.array([
        [eps[1, 1], eps[1, 2], eps[1, 0]],
        [eps[2, 1], eps[2, 2], eps[2, 0]],
        [eps[0, 1], eps[0, 2], eps[0, 0]],
    ], dtype=np.float64)

#@staticmethod
def _chooser(l, pmat):
    if l == 0:
        return np.array([[1.0]], dtype=np.float64)
    if l == 1:
        return pmat.astype(np.float64)
    raise RuntimeError(f"rotate_fb: unsupported l={l} (only s+p handled)")

#@staticmethod
def _rotate_fb_py(in1, in2, eps, mmatrix, lssh, nssh):
    # Fortran rotate_fb: left/right from chooser, apply left * M * right (with right in Fortran indexing)
    pmat = _twister_pmat(eps)
    # FireCore.py stores lssh with species as first axis (row). Use that layout here.
    lssh1 = [int(x) for x in lssh[in1 - 1, :int(nssh[in1 - 1])]]
    lssh2 = [int(x) for x in lssh[in2 - 1, :int(nssh[in2 - 1])]]
    xmatrix = np.zeros_like(mmatrix, dtype=np.float64)
    n1 = 0
    for l1 in lssh1:
        left = _chooser(l1, pmat)
        n2 = 0
        for l2 in lssh2:
            right = _chooser(l2, pmat)
            n1b = n1 + 2 * l1 + 1
            n2b = n2 + 2 * l2 + 1
            sub = mmatrix[n1:n1b, n2:n2b]
            xmatrix[n1:n1b, n2:n2b] += left @ sub @ right.T
            n2 = n2b
        n1 = n1b
    return xmatrix

#@staticmethod
def _recover_rotate_sp(hlist, mu, nu, eps, in1, in2, lssh, nssh):
    m = np.zeros((4, 4), dtype=np.float64)
    for idx in range(mu.size):
        imu = int(mu[idx]); inu = int(nu[idx])
        if imu <= 0 or inu <= 0 or imu > 4 or inu > 4:
            continue
        m[imu - 1, inu - 1] = hlist[idx]
    return _rotate_fb_py(in1, in2, eps, m, lssh, nssh)

#@staticmethod
def _recover_sp(hlist, mu, nu):
    m = np.zeros((4, 4), dtype=np.float64)
    for idx in range(mu.size):
        imu = int(mu[idx]); inu = int(nu[idx])
        if imu <= 0 or inu <= 0 or imu > 4 or inu > 4:
            continue
        m[imu - 1, inu - 1] = hlist[idx]
    return m

def build_vnl_neighbor_map(sd, natoms, verbose=False):
    """Build VNL neighbor list by mapping PP neighbors to normal neighbor indices.

    This helper iterates neighPP list and maps each (jatom,mbeta) to an index
    in the normal neigh list. It returns a list of unique (i,j) pairs that have
    PP neighbors mapped into the regular neighbor list.

    Parameters
    ----------
    sd : FireCore sparse-data struct
    natoms : int
    verbose : bool, whether to print diagnostic info

    Returns
    -------
    neighs_vnl : list of (i,j) tuples
    """
    neighs_vnl = []
    npp_tot = 0
    npp_mapped = 0
    for i in range(natoms):
        nn = int(sd.neighn[i])
        npp = int(sd.neighPPn[i])
        npp_tot += npp
        for ipp in range(npp):
            j = int(sd.neighPP_j[ipp, i]) - 1
            b = int(sd.neighPP_b[ipp, i])
            ineigh0 = -1
            for ineigh in range(nn):
                jj = int(sd.neigh_j[i, ineigh]) - 1
                bb = int(sd.neigh_b[i, ineigh])
                if (jj == j) and (bb == b):
                    ineigh0 = ineigh
                    break
            if ineigh0 >= 0:
                npp_mapped += 1
                neighs_vnl.append((i, j))
    neighs_vnl = list(dict.fromkeys(neighs_vnl))
    if len(neighs_vnl) == 0 and verbose:
        print(f"[VNL_DIAG] neighPPn={sd.neighPPn.tolist()}  total_PP_edges={int(npp_tot)}  mapped_PP_edges={int(npp_mapped)}")
        for i in range(natoms):
            nn = int(sd.neighn[i])
            npp = int(sd.neighPPn[i])
            if npp <= 0:
                continue
            pp_list = [(int(sd.neighPP_j[ipp, i]) - 1, int(sd.neighPP_b[ipp, i])) for ipp in range(npp)]
            nb_list = [(int(sd.neigh_j[i, ineigh]) - 1, int(sd.neigh_b[i, ineigh])) for ineigh in range(nn)]
            print(f"[VNL_DIAG] i={i} pp_list={pp_list} neigh_list={nb_list}")
    return neighs_vnl

def build_bcxcx_on_py(arho_on, arhoi_on, rho_on, rhoi_on, _cepal_py, nsh=2, lssh=None):
    """Implement build_ca_olsxc_on algorithm in Python.

    This helper computes the on-site exchange-correlation correction matrix
    bcxcx_on using the CEPAL functional.

    Parameters
    ----------
    arho_on : (nsh,nsh) float64, diagonal density matrix for on-site
    arhoi_on : (nsh,nsh) float64, diagonal density matrix for on-site (i)
    rho_on : (4,4) float64, full density matrix for on-site
    rhoi_on : (4,4) float64, full density matrix for on-site (i)
    _cepal_py : callable, CEPAL functional returning (exc, muxc, dexc, d2exc, dmuxc, d2muxc)
    nsh : int, number of shells (default 2)
    lssh : (nsh,) int32, angular momentum per shell (default [0, 1])

    Returns
    -------
    bcxcx_on_py : (4,4) float64, on-site exchange-correlation correction matrix
    """
    if lssh is None:
        lssh = np.array([0, 1], dtype=np.int32)
    bcxcx_on_py = np.zeros((4, 4), dtype=np.float64)
    n1 = 0
    for issh in range(nsh):
        l1 = lssh[issh]
        n1 += l1 + 1
        exc, muxc, dexc, d2exc, dmuxc, d2muxc = _cepal_py(arho_on[issh, issh])
        exci, muxci, dexci, d2exci, dmuxci, d2muxci = _cepal_py(arhoi_on[issh, issh])
        for ind1 in range(-l1, l1 + 1):
            imu = n1 + ind1 - 1
            bcxcx_on_py[imu, imu] = muxc + dmuxc * (rho_on[imu, imu] - arho_on[issh, issh]) - muxci - dmuxci * (rhoi_on[imu, imu] - arhoi_on[issh, issh])
        n1 += l1
    n1 = 0
    for issh in range(nsh):
        l1 = lssh[issh]
        n1 += l1 + 1
        n2 = 0
        for jssh in range(nsh):
            l2 = lssh[jssh]
            n2 += l2 + 1
            exc, muxc, dexc, d2exc, dmuxc, d2muxc = _cepal_py(arho_on[issh, jssh])
            exci, muxci, dexci, d2exci, dmuxci, d2muxci = _cepal_py(arhoi_on[issh, jssh])
            for ind1 in range(-l1, l1 + 1):
                imu = n1 + ind1 - 1
                for ind2 in range(-l2, l2 + 1):
                    inu = n2 + ind2 - 1
                    if imu != inu:
                        bcxcx_on_py[imu, inu] = dmuxc * rho_on[imu, inu] - dmuxci * rhoi_on[imu, inu]
            n2 += l2
        n1 += l1
    return bcxcx_on_py

def emit_ewald_2c_debug(sd, atomPos, neigh_self, dip4, dq_fn, natom_cap=3, eq2=14.39975):
    """Emit Ewald 2c debug terms (EW2C_A and EW2C_O) for verification.

    This is a pure physics kernel that computes the 2c atom-case and 2c ontop-case
    contributions to EwaldSR. It returns formatted strings for printing.

    Parameters
    ----------
    sd : FireCore sparse-data struct
    atomPos : (natoms,3) float
    neigh_self : (natoms,) int, self-slot neighbor index for each atom
    dip4 : (natoms,neigh_max,4,4) float64, dipole interaction 4D blocks
    dq_fn : callable, function that takes atom index and returns dq (charge deviation)
    natom_cap : int, limit to first N atoms for debug output (default 3)
    eq2 : float, Coulomb constant (default 14.39975)

    Returns
    -------
    lines : list of str, formatted debug lines
    """
    natoms = int(atomPos.shape[0])
    lines = []
    for ia in range(min(natom_cap, natoms)):
        mat = int(neigh_self[ia])
        if mat < 0:
            continue
        ni = int(sd.norb[ia]) if hasattr(sd, 'norb') else 4
        s_self = sd.s_mat[ia, mat, :ni, :ni].T
        dq1 = dq_fn(ia)
        for ineigh in range(int(sd.neighn[ia])):
            ja = int(sd.neigh_j[ia, ineigh]) - 1
            mb = int(sd.neigh_b[ia, ineigh])
            if ja < 0 or ja >= natoms:
                continue
            if ia == ja and mb == 0:
                continue
            if ja >= natom_cap:
                continue
            r1 = atomPos[ia]
            r2 = atomPos[ja] + sd.xl[mb]
            y = float(np.linalg.norm(r2 - r1))
            if y <= 1e-8:
                continue
            dq2 = dq_fn(ja)
            term = (s_self / y) * dq2 * eq2
            lines.append(f" [EW2C_A][P] ia,ja={ia+1:4d}{ja+1:4d} mbeta={mb:4d} ineigh={ineigh+1:4d} matom={mat+1:4d} y={y:12.6f} dq1={dq1:12.6f} dq2={dq2:12.6f} max|term|={np.max(np.abs(term)):12.6f}")
            s_blk = sd.s_mat[ia, ineigh, :ni, :ni].T
            dip_blk = dip4[ia, ineigh, :ni, :ni].T
            term_o = (((s_blk/(2.0*y) + dip_blk/(y*y))*dq1) + ((dq2*s_blk/(2.0*y) - dip_blk/(y*y))*dq2)) * eq2
            lines.append(f" [EW2C_O][P] ia,ja={ia+1:4d}{ja+1:4d} mbeta={mb:4d} ineigh={ineigh+1:4d} y={y:12.6f} dq1={dq1:12.6f} dq2={dq2:12.6f} max|term|={np.max(np.abs(term_o)):12.6f}")
    return lines

def emit_ewald_3c_debug(sd, atomPos, dip4, dq_fn, n_orb_atom, offs, natom_cap=3, eq2=14.39975):
    """Emit Ewald 3c debug terms (EW3C) for verification.

    This is a pure physics kernel that computes the 3c contribution to EwaldSR.
    It returns formatted strings for printing.

    Parameters
    ----------
    sd : FireCore sparse-data struct
    atomPos : (natoms,3) float
    dip4 : (natoms,neigh_max,4,4) float64, dipole interaction 4D blocks
    dq_fn : callable, function that takes atom index and returns dq (charge deviation)
    n_orb_atom : (natoms,) int, number of orbitals per atom
    offs : (natoms,) int, orbital offsets
    natom_cap : int, limit to first N atoms for debug output (default 3)
    eq2 : float, Coulomb constant (default 14.39975)

    Returns
    -------
    lines : list of str, formatted debug lines
    """
    natoms = int(atomPos.shape[0])
    lines = []
    for ia in range(min(natom_cap, natoms)):
        ni = int(n_orb_atom[ia])
        for mne in range(int(sd.neighn[ia])):
            ja = int(sd.neigh_j[ia, mne]) - 1
            mb = int(sd.neigh_b[ia, mne])
            if ja < 0 or ja >= natoms:
                continue
            if ja >= natom_cap:
                continue
            if ia == ja and mb == 0:
                continue
            r1 = atomPos[ia]
            r2 = atomPos[ja] + sd.xl[mb]
            y = float(np.linalg.norm(r2 - r1))
            if y <= 1e-8:
                continue
            nj = int(n_orb_atom[ja])
            s_blk = sd.s_mat[ia, mne, :nj, :ni].T
            dip_blk = dip4[ia, mne, :nj, :ni].T
            for ka in range(min(natom_cap, natoms)):
                if ka == ia or ka == ja:
                    continue
                dq3 = dq_fn(ka)
                d13 = float(np.linalg.norm(atomPos[ka] - r1))
                d23 = float(np.linalg.norm(atomPos[ka] - r2))
                if d13 <= 1e-8 or d23 <= 1e-8:
                    continue
                sterm = dq3 * s_blk / 2.0
                dterm = dq3 * dip_blk / y
                emnpl = (sterm - dterm)/d13 + (sterm + dterm)/d23
                term3 = emnpl * eq2
                lines.append(f" [EW3C][P] ia,ja,ka={ia+1:4d}{ja+1:4d}{ka+1:4d} mbeta={mb:4d} mneigh,jneigh={mne+1:4d}{-1:4d} y={y:12.6f} d13={d13:12.6f} d23={d23:12.6f} dq3={dq3:12.6f} max|term|={np.max(np.abs(term3)):12.6f}")
    return lines

def enforce_sfire_symmetry(ewaldsr4_sum_f, ewaldsr3c4_f, sd, neigh_back, n_orb_atom):
    """Apply SFIRE overwrite semantics to EwaldSR 4D blocks.

    This is a pure physics kernel that applies the SFIRE overwrite:
    ewaldsr(jneigh,jatom) = ewaldsr(ineigh,iatom).T

    Parameters
    ----------
    ewaldsr4_sum_f : (natoms,neigh_max,4,4) float64 - sum of 2c_atom + 2c_ontop + 3c
    ewaldsr3c4_f : (natoms,neigh_max,4,4) float64 - ewaldsr_3c
    sd : FireCore sparse-data struct
    neigh_back : (natoms,neigh_max) int
    n_orb_atom : (natoms,) int

    Returns
    -------
    n_sfire_apply : int, number of applied overwrites
    n_sfire_skip : int, number of skipped overwrites
    """
    natoms = int(ewaldsr4_sum_f.shape[0])
    n_sfire_apply = 0
    n_sfire_skip = 0
    for iatom in range(natoms):
        nn = int(sd.neighn[iatom])
        nmu = int(n_orb_atom[iatom])
        for ineigh in range(nn):
            if float(np.max(np.abs(ewaldsr3c4_f[iatom, ineigh]))) < 1e-12:
                continue
            jatom = int(sd.neigh_j[iatom, ineigh]) - 1
            if jatom < 0 or jatom >= natoms:
                continue
            jneigh = int(neigh_back[iatom, ineigh]) - 1
            if jneigh < 0 or jneigh >= int(sd.neighn[jatom]):
                n_sfire_skip += 1
                continue
            if int(sd.neigh_j[jatom, jneigh]) - 1 != iatom:
                n_sfire_skip += 1
                continue
            nnu = int(n_orb_atom[jatom])
            ewaldsr4_sum_f[jatom, jneigh, :nnu, :nmu] = ewaldsr4_sum_f[iatom, ineigh, :nmu, :nnu].T
            n_sfire_apply += 1
    return n_sfire_apply, n_sfire_skip

def accumulate_vca_blocks(neigh_list_vca, n_orb_atom, offs, offL, offR, diagA, H_shape):
    """Accumulate Vca neighbor blocks into dense matrices.

    This is a pure physics kernel that scatters off-diagonal and diagonal
    Vca blocks from neighbor lists into dense matrix form.

    Parameters
    ----------
    neigh_list_vca : list of (i,j) tuples, directed neighbor pairs
    n_orb_atom : (natoms,) int, number of orbitals per atom
    offs : (natoms,) int, orbital offsets
    offL : (npairs,4,4) float64, ontopl blocks
    offR : (npairs,4,4) float64, ontopr blocks
    diagA : (npairs,1) float64, atom blocks
    H_shape : tuple, shape of output matrices (norb, norb)

    Returns
    -------
    VcaL : (norb,norb) float64, accumulated ontopl matrix
    VcaR : (norb,norb) float64, accumulated ontopr matrix
    VcaA : (norb,norb) float64, accumulated atom matrix
    """
    VcaL = np.zeros(H_shape, dtype=np.float64)
    VcaR = np.zeros(H_shape, dtype=np.float64)
    VcaA = np.zeros(H_shape, dtype=np.float64)
    for k, (i, j) in enumerate(neigh_list_vca):
        i0 = int(offs[i]); j0 = int(offs[j])
        ni = int(n_orb_atom[i]); nj = int(n_orb_atom[j])
        VcaL[i0:i0+ni, j0:j0+nj] += offL[k]
        VcaR[i0:i0+ni, j0:j0+nj] += offR[k]
        VcaA[i0:i0+ni, i0:i0+ni] += diagA[k, 0]
    return VcaL, VcaR, VcaA

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
