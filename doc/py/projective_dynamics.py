import numpy as np

def build_grid_2d(nx, ny, m=1.0, m_end=1000.0, l=1.0, k=1.0, k_diag=-1.0,  l_rnd=0.0, k_rnd=0.0, m_rnd=0.0 ):
    """Build a 2D grid of points connected by springs"""
    np_total = (nx + 1) * (ny + 1)
    masses = np.ones(np_total) * m
    masses[0] = m_end
    masses[nx] = m_end
    
    # Create points
    points = np.zeros((np_total, 3))
    for iy in range(ny + 1):
        for ix in range(nx + 1):
            i = iy * (nx + 1) + ix
            points[i, 0] = ix * l
            points[i, 1] = -iy * l
    
    # Create bonds
    bonds = []
    ks    = []
    fixed = [0, nx]  # Fixed points (0-based indexing)
    
    for iy in range(ny + 1):
        for ix in range(nx + 1):
            i = iy * (nx + 1) + ix
            # Horizontal bonds
            if ix < nx:
                bonds.append((i, i + 1))
                ks.append(k)
            # Vertical bonds
            if iy < ny:
                bonds.append((i, i + nx + 1))
                ks.append(k)
            # Diagonal bonds
            if k_diag > 0:
                if ix < nx and iy < ny:
                    bonds.append((i, i + nx + 2))
                    ks.append(k_diag)
                if ix > 0 and iy < ny:
                    bonds.append((i, i + nx))
                    ks.append(k_diag)
    
    return np.array(bonds), points, masses, np.array(ks), fixed

def build_neighbor_list(bonds, n_points):
    """Build list of neighboring bonds for each point"""
    neighbs = [[] for _ in range(n_points)]
    for i, (i_, j_) in enumerate(bonds):
        neighbs[i_].append(i)
        neighbs[j_].append(i)
    return neighbs

def make_pd_matrix(neighbs, bonds, masses, dt, ks):
    """Create the system matrix for projective dynamics"""
    np_total = len(masses)
    A = np.zeros((np_total, np_total))
    idt2 = 1.0 / (dt * dt)
    
    for i in range(np_total):
        Aii = masses[i] * idt2
        for ib in neighbs[i]:
            k = ks[ib]
            i_, j_ = bonds[ib]
            Aii += k
            if j_ > i:
                A[i, j_] = -k
                A[j_, i] = -k
            elif i_ > i:
                A[i, i_] = -k
                A[i_, i] = -k
        A[i, i] = Aii
    
    Mt = masses * idt2
    return A, Mt

def make_pd_rhs(neighbs, bonds, masses, dt, ks, points, l0s, pnew):
    """Build the right-hand side of the system following the Julia implementation"""
    np_total = len(masses)
    b = np.zeros((np_total, 3))
    idt2 = 1.0 / (dt * dt)
    
    for i in range(np_total):
        # Mass term (inertial prediction)
        bi = pnew[i] * (masses[i] * idt2)
        
        # Spring terms
        for ib in neighbs[i]:
            k = ks[ib]
            i_, j_ = bonds[ib]
            j = j_ if i_ == i else i_
            
            # Using predicted positions for better propagation
            d = pnew[i] - pnew[j]
            d_norm = np.linalg.norm(d)
            if d_norm > 1e-10:  # Avoid division by zero
                d *= k * l0s[ib] / d_norm
                bi += d
        
        b[i] = bi
    
    return b

def update_velocity(ps_cor, points, velocity, dt):
    """Update velocity based on position change"""
    return (ps_cor - points) / dt

def solve_pd(points, velocity, bonds, masses, ks, dt=0.1, n_iter=100, gravity=np.array([0, -9.81, 0]), fixed_points=None, call_back=None, damping=0.01 ):
    """Solve the system using projective dynamics following the Julia implementation"""
    n_points = len(points)
    neighbs = build_neighbor_list(bonds, n_points)
    A, Mt = make_pd_matrix(neighbs, bonds, masses, dt, ks)
    
    # Initialize
    pos = points.copy()
    pos_cor = points.copy()
    pos_pred = points.copy()
    l0s = np.array([np.linalg.norm(points[j] - points[i]) for i, j in bonds])
    
    # Main simulation loop
    for itr in range(n_iter):
        
        velocity*=(1-damping)
        velocity[:,:] += gravity[None,:] * dt
        pos_pred[:] = pos + velocity * dt
        
        # Apply fixed point constraints
        if fixed_points is not None:
            pos_pred[fixed_points] = points[fixed_points]
        
        # Build right-hand side
        b = make_pd_rhs(neighbs, bonds, masses, dt, ks, pos, l0s, pos_pred)
        
        # Solve system for each coordinate
        for i in range(3):
            pos_cor[:, i] = np.linalg.solve(A, b[:, i])
        
        # Apply fixed point constraints
        if fixed_points is not None:
            pos_cor[fixed_points] = points[fixed_points]
        
        # Update velocity and position
        velocity = update_velocity(pos_cor, pos, velocity, dt)
        pos[:] = pos_cor[:]
        
        if call_back is not None:
            call_back(pos)
        
        # Print debug info
        v_norm = np.linalg.norm(velocity)
        pos_diff = np.linalg.norm(pos - pos_pred)
        print(f"iter:{itr} dt={dt:.3e} |v|={v_norm:.3e} dp={pos_diff:.3e}")
        
        # Optional: check bond lengths for debugging
        if False:  # Set to True to enable
            bond_errors = []
            for (i, j), l0 in zip(bonds, l0s):
                current_length = np.linalg.norm(pos[j] - pos[i])
                bond_errors.append(abs(current_length - l0))
            max_bond_error = max(bond_errors)
            print(f"    max bond length error: {max_bond_error:.3e}")
    
    return pos, velocity

# Example usage
if __name__ == "__main__":
    # Create a simple 5x5 grid
    nx, ny = 5, 5
    bonds, points, masses, ks, fixed = build_grid_2d(nx, ny)
    velocity = np.zeros_like(points)
    
    # Run simulation
    new_points, new_velocity = solve_pd(points, velocity, bonds, masses, ks, fixed_points=fixed)
