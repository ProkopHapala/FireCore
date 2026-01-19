import numpy as np
import os

class FdataParser:
    """
    Fdata helper:
    - 2c files (overlap, kinetic, vna, vnl, …) store columns of mu-nu matrix elements.
      For vnl (non-local PP) each column is a projector-channel pair; count equals
      num_nonzero_pp computed from lsshPP of both species. Plotters take a channels
      argument only to avoid drawing every column when many exist.
    - basis/*.pp are one-body radial pseudopotentials (no channels); use read_wf/read_wf
      to load wavefunctions and parse_info to get rc_PP cutoffs.
    """
    def __init__(self, fdata_dir):
        self.fdata_dir = fdata_dir

    def parse_info(self):
        """Parses info.dat to get basis information for species."""
        path = os.path.join(self.fdata_dir, "info.dat")
        if not os.path.exists(path):
            raise FileNotFoundError(f"info.dat not found in {self.fdata_dir}")
        
        with open(path, 'r') as f:
            lines = f.readlines()
        
        # Line 2: Number of species
        num_species = int(lines[1].split()[0])
        
        self.species_info = {}
        idx = 3
        for _ in range(num_species):
            # Skip until "Information for this species"
            while "Information for this species" not in lines[idx]: idx += 1
            idx += 2 # Skip element name
            nz = int(lines[idx].split()[0])
            idx += 2 # Skip mass
            nssh = int(lines[idx].split()[0])
            idx += 1
            lssh = [int(x) for x in lines[idx].split()]
            idx += 1
            # Pseudopotential shells
            nsshPP = int(lines[idx].split()[0])
            idx += 1
            # DEBUG keep loud if missing PP shells
            if lines[idx].strip():
                lsshPP = [int(x) for x in lines[idx].split()]
            else:
                print(f"[DEBUG] parse_info: empty lsshPP for nz={nz}, nsshPP={nsshPP}")
                lsshPP = []
            idx += 1
            rc_PP = float(lines[idx].split()[0])
            idx += 3  # skip Qneutral and rcuts
            # Store
            self.species_info[nz] = {
                'nssh': nssh, 'lssh': lssh,
                'nsshPP': nsshPP, 'lsshPP': lsshPP,
                'rc_PP': rc_PP,
            }
            idx += 1
    def get_num_nonzero(self, nz1, nz2):
        """Calculates num_nonzero for a pair of species based on make_munu.f90."""
        if not hasattr(self, 'species_info'):
            self.parse_info()
        
        info1 = self.species_info[nz1]
        info2 = self.species_info[nz2]
        
        index = 0
        for l1 in info1['lssh']:
            for l2 in info2['lssh']:
                index += 2 * min(l1, l2) + 1
        return index

    def get_num_nonzero_pp(self, nz1, nz2):
        """Num_nonzero for PP (index_maxPP) using lsshPP."""
        if not hasattr(self, 'species_info'):
            self.parse_info()
        info1 = self.species_info[nz1]
        info2 = self.species_info[nz2]
        index = 0
        for l1 in info1['lsshPP']:
            for l2 in info2['lsshPP']:
                index += 2 * min(l1, l2) + 1
        return index

    def read_2c(self, fname):
        """
        Reads 2-center integral data from a .dat file.
        Uses get_num_nonzero to read correct number of values per point.
        """
        with open(fname, 'r') as f:
            lines = f.readlines()
        
        # Skip 9 header lines
        idx = 9
        
        # Line 10: nucz1, rc1
        while not lines[idx].strip(): idx += 1
        parts1 = lines[idx].split()
        nucz1 = int(parts1[0])
        rc1 = float(parts1[1])
        idx += 1
        
        # Line 11: nucz2, rc2
        while not lines[idx].strip(): idx += 1
        parts2 = lines[idx].split()
        nucz2 = int(parts2[0])
        rc2 = float(parts2[1])
        idx += 1
        
        # Detect interaction type from filename prefix
        base = os.path.basename(fname)
        is_pp = base.startswith("vnl")
        npseudo = None
        cl_pseudo = None

        # Fortran readheader_2c has extra 2 header lines for interaction=5:
        #   read npseudo
        #   read cl_pseudo(1:npseudo)
        if is_pp:
            # npseudo
            while idx < len(lines) and not lines[idx].strip():
                idx += 1
            if idx >= len(lines):
                raise ValueError(f"read_2c: unexpected EOF while reading npseudo in {fname}")
            npseudo = int(float(lines[idx].split()[0].replace('D','E')))
            idx += 1

            # cl_pseudo values (can be on one line)
            cl_vals = []
            while idx < len(lines) and len(cl_vals) < npseudo:
                if lines[idx].strip():
                    for t in lines[idx].split():
                        try:
                            cl_vals.append(float(t.replace('D','E')))
                        except ValueError:
                            pass
                idx += 1
            if len(cl_vals) < npseudo:
                raise ValueError(f"read_2c: failed to read cl_pseudo(1:{npseudo}) in {fname}, got {len(cl_vals)}")
            cl_pseudo = np.array(cl_vals[:npseudo], dtype=np.float64)

        # zmax, numz (Fortran: read (iounit,*) zmax, numz)
        while idx < len(lines) and not lines[idx].strip():
            idx += 1
        if idx >= len(lines):
            raise ValueError(f"read_2c: unexpected EOF while reading zmax/numz in {fname}")
        parts3 = lines[idx].split()
        if len(parts3) < 2:
            raise ValueError(f"read_2c: cannot parse zmax/numz in {fname}, line='{lines[idx].rstrip()}'")
        zmax = float(parts3[0].replace('D','E'))
        numz = int(float(parts3[1].replace('D','E')))
        idx += 1
        
        num_nonzero = self.get_num_nonzero_pp(nucz1, nucz2) if is_pp else self.get_num_nonzero(nucz1, nucz2)
        
        # Data: Read num_nonzero values for each of numz points
        # One point's data can be spread across multiple lines.
        data = []
        all_values = []
        # Flatten all remaining lines into one list of values
        for line in lines[idx:]:
            all_values.extend([float(x.replace('D', 'E')) for x in line.split()])
        
        total = len(all_values)
        expected = numz * num_nonzero
        if total < expected:
            # DEBUG: fall back to inferred num_nonzero from file length
            if total % numz == 0:
                inferred = total // numz
                print(f"[WARN] read_2c: values short for {fname}, expected {expected} got {total}; "
                      f"using inferred num_nonzero={inferred}")
                num_nonzero = inferred
            else:
                raise ValueError(f"read_2c: cannot reshape {total} values into ({numz}, {num_nonzero}) for {fname}")
        data = np.array(all_values[:numz * num_nonzero]).reshape(numz, num_nonzero)
        
        return {
            'nucz1': nucz1, 'rc1': rc1,
            'nucz2': nucz2, 'rc2': rc2,
            'zmax': zmax, 'numz': numz,
            'num_nonzero': num_nonzero,
            'data': data,
            'npseudo': npseudo,
            'cl_pseudo': cl_pseudo,
        }

    def build_spline_1d(self, y, xmax):
        """
        Replicates buildspline_1d.f90 logic.
        Constructs clamped cubic spline coefficients.
        y: array of values at points i=0..N-1
        xmax: maximum distance
        Returns: [4, N] array where [0,:]=a, [1,:]=b, [2,:]=c, [3,:]=d
        """
        n = len(y)
        h = xmax / (n - 1)
        
        # a is the input values
        a = y.copy()
        norder = n - 1
        
        # Get estimates of derivatives (fpo, fpn)
        fpo = (a[1] - a[0]) / h
        fpn = (a[norder] - a[norder-1]) / h
        
        alpha = np.zeros(n)
        alpha[0] = 3.0*(a[1] - a[0])/h - 3.0*fpo
        for i in range(1, norder):
            alpha[i] = 3.0*(a[i+1] - 2.0*a[i] + a[i-1])/h
        alpha[norder] = 3.0*fpn - 3.0*(a[norder] - a[norder-1])/h
        
        L = np.zeros(n)
        mu = np.zeros(n)
        Z = np.zeros(n)
        
        L[0] = 2.0*h
        mu[0] = 0.5
        Z[0] = alpha[0]/L[0]
        
        for i in range(1, norder):
            L[i] = (4.0 - mu[i-1])*h
            mu[i] = h/L[i]
            Z[i] = (alpha[i] - h*Z[i-1])/L[i]
            
        L[norder] = (2.0 - mu[norder-1])*h
        mu[norder] = 0.0
        Z[norder] = (alpha[norder] - h*Z[norder-1])/L[norder]
        
        c = np.zeros(n)
        b = np.zeros(n)
        d = np.zeros(n)
        
        c[norder] = Z[norder]
        for i in range(norder - 1, -1, -1):
            c[i] = Z[i] - mu[i]*c[i+1]
            b[i] = (a[i+1] - a[i])/h - h*(c[i+1] + 2.0*c[i])/3.0
            d[i] = (c[i+1] - c[i])/(3.0*h)
            
        # Match splineint_2c structure: (4, numz)
        spline = np.zeros((4, n))
        spline[0, :] = a
        spline[1, :] = b
        spline[2, :] = c
        spline[3, :] = d
        return spline

    def get_num_nonzero_3c(self, nz1, nz2):
        """Calculates num_nonzero_3c for a pair of species based on make_munu.f90."""
        if not hasattr(self, 'species_info'):
            self.parse_info()
        
        info1 = self.species_info[nz1]
        info2 = self.species_info[nz2]
        
        index = 0
        # 2-center interactions (sigma, pi, delta...)
        for l1 in info1['lssh']:
            for l2 in info2['lssh']:
                index += 2 * min(l1, l2) + 1
        
        # Extra 3-center interactions from getIsorpRange etc.
        # FireCore uses complex logic. Let's replicate make_munu.f90 lines 244-344.
        # This part depends on m_i and m_j mixing.
        for l1 in info1['lssh']:
            for l2 in info2['lssh']:
                # M=1 cases
                if l1 == 0 and l2 != 0: index += 1
                if l1 == 1:
                    index += 1
                    if l2 != 0: index += 1
                    if l2 == 2: index += 2
                if l1 == 2:
                    index += 1
                    if l2 != 0: index += 3
                    if l2 == 2: index += 2
                # M=2 cases
                if l1 == 2: index += 1
                if l2 == 2: index += 1
                
        return index

    def _try_float(self, s):
        try:
            return float(s.replace('D', 'E'))
        except ValueError:
            return None

    def read_3c(self, fname):
        """
        Reads 3-center integral data from a .dat file.
        Matches readheader_3c and readdata_3c.
        """
        # Parse nz from filename to determine num_nonzero
        # fname: ...root_it_is.nz1.nz2.nz3.dat
        basename = os.path.basename(fname)
        parts = basename.replace('.dat', '').split('.')
        nz1, nz2, nz3 = int(parts[1]), int(parts[2]), int(parts[3])
        num_nonzero = self.get_num_nonzero_3c(nz1, nz2)

        with open(fname, 'r') as f:
            lines = f.readlines()
            
        # Parse header exactly like readheader_3c.f90
        # - skip 10 lines
        # - read: nphi2 nr ntheta
        # - read: ymax numy
        # - read: xmax numx
        # - skip 1 line
        # - read 3 lines: nucZ rc
        # - skip 1 line
        idx = 10
        parts_n = lines[idx].split(); nphi2 = int(parts_n[0]); nr = int(parts_n[1]); ntheta_in = int(parts_n[2]); idx += 1
        parts_y = lines[idx].split(); ymax = self._try_float(parts_y[0]); numy = int(parts_y[1]); idx += 1
        parts_x = lines[idx].split(); xmax = self._try_float(parts_x[0]); numx = int(parts_x[1]); idx += 1
        # skip decorative line
        idx += 1
        # read 3 charge lines (we don't use values here)
        idx += 3
        # skip decorative line
        idx += 1
        
        all_values = []
        for line in lines[idx:]:
            for x in line.split():
                val = self._try_float(x)
                if val is not None:
                    all_values.append(val)
        
        total_points = numy * numx
        if len(all_values) < total_points * num_nonzero:
            # Handle potential mismatch or truncated files
            pass
            
        # Fortran file read loops (readdata_3c.f90):
        #   do jpoint=1,numy
        #     do ipoint=1,numx
        #       read ... (integral=1,num_nonzero)
        # So the natural C-order reshape is [y,x,integral].
        data = np.array(all_values[:total_points * num_nonzero]).reshape(numy, numx, num_nonzero)
        
        return {
            'ymax': ymax, 'numy': numy,
            'xmax': xmax, 'numx': numx,
            'num_nonzero': num_nonzero,
            'data': data
        }

    def find_2c(self, root, nz1, nz2):
        """Finds 2-center data file."""
        fname = f"{root}.{nz1:02d}.{nz2:02d}.dat"
        return os.path.join(self.fdata_dir, fname)

    def find_3c(self, root, itheta, isorp, nz1, nz2, nz3):
        """Finds 3-center data file."""
        fname = f"{root}_{itheta:02d}_{isorp:02d}.{nz1:02d}.{nz2:02d}.{nz3:02d}.dat"
        return os.path.join(self.fdata_dir, fname)

    def read_wf(self, fname):
        """
        Reads radial wavefunction from a .wf file.
        Matches read_wf.f90.
        """
        with open(fname, 'r') as f:
            lines = f.readlines()
        
        # Header (trash filename, nzx, mesh, rcutoff, rcmax, xnocc, l)
        nzxwf = int(lines[1].split()[0])
        mesh = int(lines[2].split()[0])
        parts_r = lines[3].split()
        rcutoffwf = float(parts_r[0].replace('D','E'))
        rcmax = float(parts_r[1].replace('D','E'))
        xnoccwf = float(parts_r[2].replace('D','E'))
        lqnwf = int(lines[4].split()[0])
        
        all_values = []
        for line in lines[5:]:
            all_values.extend([float(x.replace('D', 'E')) for x in line.split()])
            
        data = np.array(all_values[:mesh])
        return {
            'nzx': nzxwf, 'mesh': mesh, 
            'rcutoff': rcutoffwf, 'rcmax': rcmax, 'xnocc': xnoccwf, 'l': lqnwf,
            'data': data
        }

    def find_wf(self, species_nz):
        """Finds all .wf-style files for a given species (supports .wf, .wf1, etc.)."""
        import glob
        patterns = [
            os.path.join(self.fdata_dir, f"*.{species_nz:02d}.wf"),
            os.path.join(self.fdata_dir, f"**/*{species_nz:03d}*.wf*"),
        ]
        found = []
        for pat in patterns:
            found.extend(glob.glob(pat, recursive=True))
        return sorted(set(found))

    def parse_info(self):
        """Parses info.dat to populate species_info (nssh, lssh, element, Z)."""
        info_path = os.path.join(self.fdata_dir, "info.dat")
        if not os.path.exists(info_path):
            # try parent dir if we are inside basis/
            alt = os.path.join(os.path.dirname(self.fdata_dir), "info.dat")
            if os.path.exists(alt):
                info_path = alt
        if not os.path.exists(info_path):
            raise FileNotFoundError(f"info.dat not found under {self.fdata_dir} or parent")

        species_info = {}
        with open(info_path, 'r') as f:
            lines = [ln.strip() for ln in f.readlines()]

        i = 0
        while i < len(lines):
            line = lines[i]
            if "- Information for this species" in line:
                # Expect next lines: element, Z, mass, shells, L list, PP shells, rc_PP
                if i + 8 >= len(lines): break
                element = lines[i+1].split()[0]
                z = int(float(lines[i+2].split()[0]))
                # skip mass line (i+3)
                nssh = int(lines[i+4].split()[0])
                lssh_line = lines[i+5].strip()
                lssh = [int(x) for x in lssh_line.split()] if lssh_line else []
                nsshPP = int(lines[i+6].split()[0])
                lsshPP_line = lines[i+7].strip()
                if lsshPP_line:
                    lsshPP = [int(x) for x in lsshPP_line.split()]
                else:
                    print(f"[DEBUG] parse_info (trimmed): empty lsshPP for z={z}, nsshPP={nsshPP}")
                    lsshPP = []
                rc_PP = None
                try:
                    rc_PP = float(lines[i+8].split()[0])
                except Exception as e:
                    print(f"[DEBUG] parse_info (trimmed): failed rc_PP for z={z} line='{lines[i+8]}' err={e}")
                species_info[z] = {
                    'element': element,
                    'nssh': nssh, 'lssh': lssh,
                    'nsshPP': nsshPP, 'lsshPP': lsshPP,
                    'rc_PP': rc_PP,
                }
                # Advance to next species block (skip until separator or move past parsed lines)
                i = i + 9
                while i < len(lines) and lines[i].strip() and "====" not in lines[i]:
                    i += 1
            else:
                i += 1
        self.species_info = species_info

    def load_species_data(self, species_nz):
        """
        Loads all relevant Fdata for a list of species nuclear charges.
        species_nz: list of unique nuclear charges in the system.
        """
        if not hasattr(self, 'species_info'):
            self.parse_info()
        
        data_2c = {}
        for nz1 in species_nz:
            for nz2 in species_nz:
                for root in ['overlap', 'kinetic', 'vna', 'vnl', 'vxc', 'vna_atom_00', 'vna_ontopl_00', 'vna_ontopr_00', 'dipole_z', 'dipole_y', 'dipole_x']:
                    path = self.find_2c(root, nz1, nz2)
                    if os.path.exists(path):
                        # NOTE: Fortran assemble_2c uses interaction=2 (ontopl) and interaction=3 (ontopr).
                        # Our scans show the *files* vna_ontopl_* and vna_ontopr_* are swapped relative to
                        # these interaction numbers, so we swap the semantic keys here.
                        store_root = root
                        if root.startswith('vna_ontopl_'):
                            store_root = root.replace('vna_ontopl_', 'vna_ontopr_', 1)
                        elif root.startswith('vna_ontopr_'):
                            store_root = root.replace('vna_ontopr_', 'vna_ontopl_', 1)
                        data_2c[(store_root, nz1, nz2)] = self.read_2c(path)

                # For vna we also need all shell-resolved files like vna_atom_XX, vna_ontopl_XX, vna_ontopr_XX
                nssh1 = self.species_info.get(nz1, {}).get('nssh', 0)
                nssh2 = self.species_info.get(nz2, {}).get('nssh', 0)
                # atom case uses shells of atom 2 (matches read_2c logic: interaction 4 uses nssh(in2))
                for isorp in range(nssh2 + 1):
                    root = f"vna_atom_{isorp:02d}"
                    path = self.find_2c(root, nz1, nz2)
                    if os.path.exists(path):
                        data_2c[(root, nz1, nz2)] = self.read_2c(path)
                # ontop left uses shells of atom 1 (interaction 2 uses nssh(in1))
                for isorp in range(nssh1 + 1):
                    root = f"vna_ontopl_{isorp:02d}"
                    path = self.find_2c(root, nz1, nz2)
                    if os.path.exists(path):
                        store_root = f"vna_ontopr_{isorp:02d}"
                        data_2c[(store_root, nz1, nz2)] = self.read_2c(path)
                # ontop right uses shells of atom 2 (interaction 3 uses nssh(in2))
                for isorp in range(nssh2 + 1):
                    root = f"vna_ontopr_{isorp:02d}"
                    path = self.find_2c(root, nz1, nz2)
                    if os.path.exists(path):
                        store_root = f"vna_ontopl_{isorp:02d}"
                        data_2c[(store_root, nz1, nz2)] = self.read_2c(path)

                # For average_rho off-site 2c density seed we need shell-resolved density tables:
                #   interaction=15: den_ontopl <i|n_i|j> weighted by Qneutral(isorp,in1), isorp=1..nssh(in1)
                #   interaction=16: den_ontopr <i|n_j|j> weighted by Qneutral(isorp,in2), isorp=1..nssh(in2)
                # File naming in fdata uses explicit isorp suffix (01,02,...) for these.
                for isorp in range(1, nssh1 + 1):
                    root = f"den_ontopl_{isorp:02d}"
                    path = self.find_2c(root, nz1, nz2)
                    if os.path.exists(path):
                        data_2c[(root, nz1, nz2)] = self.read_2c(path)
                for isorp in range(1, nssh2 + 1):
                    root = f"den_ontopr_{isorp:02d}"
                    path = self.find_2c(root, nz1, nz2)
                    if os.path.exists(path):
                        data_2c[(root, nz1, nz2)] = self.read_2c(path)
        
        # 3-center data (example for den3)
        data_3c = {}
        # We need itheta and isorp ranges. 
        # Typically itheta=1..5, isorp=0..isorpmax.
        # FireCore uses getIsorpRange to determine isorp range.
        # For now, let's just search for what exists.
        import glob
        for root in ['den3', 'bcna']:
            pattern = os.path.join(self.fdata_dir, f"{root}_*.*.*.*.dat")
            for path in glob.glob(pattern):
                fname = os.path.basename(path)
                # Parse: root_itheta_isorp.nz1.nz2.nz3.dat
                parts = fname.replace('.dat', '').split('.')
                p_root_th_is = parts[0].split('_')
                if len(p_root_th_is) < 3: continue
                root_f = p_root_th_is[0]
                it = int(p_root_th_is[1])
                isorp = int(p_root_th_is[2])
                nz1 = int(parts[1])
                nz2 = int(parts[2])
                nz3 = int(parts[3])
                
                if nz1 in species_nz and nz2 in species_nz and nz3 in species_nz:
                    data_3c[(root_f, it, isorp, nz1, nz2, nz3)] = self.read_3c(path)
                    
        return data_2c, data_3c




def read_pp(path: str):
    """
    Minimal parser for radial pseudopotential basis/*.pp.
    Returns r (Angstrom) and V(r) arrays using lines that contain two floats.
    """
    if not os.path.exists(path):
        raise RuntimeError(f"PP file not found: {path}")
    r_list, v_list = [], []
    with open(path, "r") as f:
        for line in f:
            parts = line.split()
            if len(parts) != 2:
                continue
            try:
                r_val = float(parts[0].replace('D','E'))
                v_val = float(parts[1].replace('D','E'))
            except ValueError:
                continue
            r_list.append(r_val)
            v_list.append(v_val)
    if not r_list:
        raise RuntimeError(f"Failed to parse any (r,V) pairs from {path}")
    return np.array(r_list), np.array(v_list)

    def read_xc1c_tables(self):
        """
        Reads Vxc_1c tables from xc1c_dqi.dat and nuxc1crho.dat.
        Returns dict with exc1c0, nuxc1c, dexc1c, dnuxc1c arrays.
        Follows read_1c.f90 logic for itheory_xc=2/4.
        """
        if not hasattr(self, 'species_info'):
            self.parse_info()
        
        # Get nspecies from species_info
        nspecies = len(self.species_info)
        nsh_max = max(info['nssh'] for info in self.species_info.values())
        
        # Allocate arrays
        exc1c0 = np.zeros((nspecies, nsh_max, nsh_max), dtype=np.float64)
        nuxc1c = np.zeros((nspecies, nsh_max, nsh_max), dtype=np.float64)
        dexc1c = np.zeros((nspecies, nsh_max, nsh_max, nsh_max), dtype=np.float64)
        dnuxc1c = np.zeros((nspecies, nsh_max, nsh_max, nsh_max), dtype=np.float64)
        
        # Read xc1c_dqi.dat
        path_dqi = os.path.join(self.fdata_dir, 'xc1c_dqi.dat')
        if not os.path.exists(path_dqi):
            raise FileNotFoundError(f"xc1c_dqi.dat not found in {self.fdata_dir}")
        
        with open(path_dqi, 'r') as f:
            lines = f.readlines()
        
        # Skip 4 header lines
        idx = 4
        
        # Skip species header lines
        for _ in range(nspecies + 1):
            while idx < len(lines) and not lines[idx].strip():
                idx += 1
            if idx < len(lines):
                idx += 1
        
        # Read data
        in2 = 1
        for in1 in range(1, nspecies + 1):
            if in1 > nspecies:
                break
            while idx < len(lines) and not lines[idx].strip():
                idx += 1
            if idx >= len(lines):
                break
            
            parts = lines[idx].split()
            itype = int(parts[0])
            numsh = int(parts[1])
            idx += 1
            
            if numsh != self.species_info[in2]['nssh']:
                raise ValueError(f"numsh mismatch: expected {self.species_info[in2]['nssh']}, got {numsh}")
            
            # Read exc1c0
            for issh in range(numsh):
                while idx < len(lines) and not lines[idx].strip():
                    idx += 1
                if idx >= len(lines):
                    break
                values = [float(x.replace('D', 'E')) for x in lines[idx].split()]
                for jssh in range(numsh):
                    exc1c0[in2-1, issh, jssh] = values[jssh]
                idx += 1
            
            # Skip blank line
            while idx < len(lines) and not lines[idx].strip():
                idx += 1
            
            # Read nuxc1c
            for issh in range(numsh):
                while idx < len(lines) and not lines[idx].strip():
                    idx += 1
                if idx >= len(lines):
                    break
                values = [float(x.replace('D', 'E')) for x in lines[idx].split()]
                for jssh in range(numsh):
                    nuxc1c[in2-1, issh, jssh] = values[jssh]
                idx += 1
            
            in2 += 1
            if in2 > nspecies:
                break
        
        # Read nuxc1crho.dat
        path_rho = os.path.join(self.fdata_dir, 'nuxc1crho.dat')
        if not os.path.exists(path_rho):
            raise FileNotFoundError(f"nuxc1crho.dat not found in {self.fdata_dir}")
        
        with open(path_rho, 'r') as f:
            lines = f.readlines()
        
        # Skip 4 header lines
        idx = 4
        
        # Skip species header lines
        for _ in range(nspecies + 1):
            while idx < len(lines) and not lines[idx].strip():
                idx += 1
            if idx < len(lines):
                idx += 1
        
        # Read data
        in2 = 1
        for in1 in range(1, nspecies + 1):
            if in1 > nspecies:
                break
            while idx < len(lines) and not lines[idx].strip():
                idx += 1
            if idx >= len(lines):
                break
            
            parts = lines[idx].split()
            itype = int(parts[0])
            numsh = int(parts[1])
            kkssh = int(parts[2])
            idx += 1
            
            if numsh != self.species_info[in2]['nssh']:
                raise ValueError(f"numsh mismatch: expected {self.species_info[in2]['nssh']}, got {numsh}")
            
            # Read dnuxc1c for each kssh
            for kssh in range(numsh):
                while idx < len(lines) and not lines[idx].strip():
                    idx += 1
                if idx >= len(lines):
                    break
                parts = lines[idx].split()
                itype = int(parts[0])
                numsh_k = int(parts[1])
                kkssh_k = int(parts[2])
                idx += 1
                
                for issh in range(numsh):
                    while idx < len(lines) and not lines[idx].strip():
                        idx += 1
                    if idx >= len(lines):
                        break
                    values = [float(x.replace('D', 'E')) for x in lines[idx].split()]
                    for jssh in range(numsh):
                        dnuxc1c[in2-1, issh, jssh, kssh] = values[jssh]
                    idx += 1
            
            in2 += 1
            if in2 > nspecies:
                break
        
        return {
            'exc1c0': exc1c0,
            'nuxc1c': nuxc1c,
            'dexc1c': dexc1c,
            'dnuxc1c': dnuxc1c,
        }


if __name__ == "__main__":
    parser = FdataParser("/home/prokophapala/git/FireCore/tests/pyFireball/Fdata")
    
    # Test for H and C
    species = [1, 6]
    d2c, d3c = parser.load_species_data(species)
    
    print("Loaded 2C interactions:", len(d2c))
    for k in list(d2c.keys())[:5]:
        print(f"  {k}: {d2c[k]['data'].shape}")
        
    print("Loaded 3C interactions:", len(d3c))
    for k in list(d3c.keys())[:5]:
        print(f"  {k}: {d3c[k]['data'].shape}")
