import numpy as np
import os

class FdataParser:
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
            self.species_info[nz] = {'nssh': nssh, 'lssh': lssh}
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
        
        # Line 12: zmax, numz
        while not lines[idx].strip(): idx += 1
        parts3 = lines[idx].split()
        zmax = float(parts3[0].replace('D','E'))
        numz = int(parts3[1])
        idx += 1
        
        num_nonzero = self.get_num_nonzero(nucz1, nucz2)
        
        # Data: Read num_nonzero values for each of numz points
        # One point's data can be spread across multiple lines.
        data = []
        all_values = []
        # Flatten all remaining lines into one list of values
        for line in lines[idx:]:
            all_values.extend([float(x.replace('D', 'E')) for x in line.split()])
        
        # Reshape into [numz, num_nonzero]
        if len(all_values) < numz * num_nonzero:
            # Maybe it's a special interaction (like interaction 12 coulomb)?
            # In read_2c.f90, some have different num_nonzero.
            # But overlapping, kinetic, vna usually follow index_max2c.
            pass
            
        data = np.array(all_values[:numz * num_nonzero]).reshape(numz, num_nonzero)
        
        return {
            'nucz1': nucz1, 'rc1': rc1,
            'nucz2': nucz2, 'rc2': rc2,
            'zmax': zmax, 'numz': numz,
            'num_nonzero': num_nonzero,
            'data': data
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
            
        # 10 header lines
        idx = 10
        while idx < len(lines) and not lines[idx].strip(): idx += 1
        parts_n = lines[idx].split()
        nphi2 = int(parts_n[0])
        nr = int(parts_n[1])
        ntheta_in = int(parts_n[2])
        idx += 1
        
        while idx < len(lines) and not lines[idx].strip(): idx += 1
        parts_y = lines[idx].split()
        ymax = self._try_float(parts_y[0])
        numy = int(parts_y[1])
        idx += 1
        
        while idx < len(lines) and not lines[idx].strip(): idx += 1
        parts_x = lines[idx].split()
        xmax = self._try_float(parts_x[0])
        numx = int(parts_x[1])
        idx += 1
        
        # Skip decorative and charge lines
        while idx < len(lines) and "===" not in lines[idx]: idx += 1
        idx += 1 # skip first ===
        while idx < len(lines) and "===" not in lines[idx]: idx += 1
        idx += 1 # skip second ===
        
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

    def load_species_data(self, species_nz):
        """
        Loads all relevant Fdata for a list of species nuclear charges.
        species_nz: list of unique nuclear charges in the system.
        """
        data_2c = {}
        for nz1 in species_nz:
            for nz2 in species_nz:
                for root in ['overlap', 'kinetic', 'vna', 'vxc', 'vna_atom_00', 'vna_ontopl_00', 'vna_ontopr_00']:
                    path = self.find_2c(root, nz1, nz2)
                    if os.path.exists(path):
                        data_2c[(root, nz1, nz2)] = self.read_2c(path)
                # Fallback: if vna exists as vna_atom_00, map it to 'vna' for simplicity if needed
                if ('vna_atom_00', nz1, nz2) in data_2c and ('vna', nz1, nz2) not in data_2c:
                    data_2c[('vna', nz1, nz2)] = data_2c[('vna_atom_00', nz1, nz2)]
        
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
