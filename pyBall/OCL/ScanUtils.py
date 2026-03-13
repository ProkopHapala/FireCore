"""
ScanUtils.py - Trajectory generators for molecule-substrate interaction scanning.
Generates arrays of rigid-body transforms (rotation + translation) for various scan types.
No GUI or OpenCL dependency - pure numpy.
"""

import numpy as np

# ======== Quaternion utilities ========

def quat_from_axis_angle(axis, angle):
    """Quaternion (w,x,y,z) from axis-angle."""
    axis = np.asarray(axis, dtype=np.float64)
    axis = axis / np.linalg.norm(axis)
    s = np.sin(angle * 0.5)
    return np.array([np.cos(angle * 0.5), axis[0]*s, axis[1]*s, axis[2]*s])

def quat_multiply(q1, q2):
    """Hamilton product of two quaternions (w,x,y,z)."""
    w1,x1,y1,z1 = q1
    w2,x2,y2,z2 = q2
    return np.array([
        w1*w2 - x1*x2 - y1*y2 - z1*z2,
        w1*x2 + x1*w2 + y1*z2 - z1*y2,
        w1*y2 - x1*z2 + y1*w2 + z1*x2,
        w1*z2 + x1*y2 - y1*x2 + z1*w2,
    ])

def quat_to_rotmat(q):
    """Convert quaternion (w,x,y,z) to 3x3 rotation matrix."""
    w,x,y,z = q
    return np.array([
        [1-2*(y*y+z*z),   2*(x*y-z*w),   2*(x*z+y*w)],
        [  2*(x*y+z*w), 1-2*(x*x+z*z),   2*(y*z-x*w)],
        [  2*(x*z-y*w),   2*(y*z+x*w), 1-2*(x*x+y*y)],
    ], dtype=np.float64)

def quat_slerp(q0, q1, t):
    """Spherical linear interpolation between quaternions q0,q1 at parameter t in [0,1]."""
    q0 = np.asarray(q0, dtype=np.float64)
    q1 = np.asarray(q1, dtype=np.float64)
    dot = np.dot(q0, q1)
    if dot < 0:  # take shortest path
        q1 = -q1
        dot = -dot
    dot = min(dot, 1.0)
    if dot > 0.9995:  # very close, use linear interpolation
        result = q0 + t * (q1 - q0)
        return result / np.linalg.norm(result)
    theta = np.arccos(dot)
    sin_theta = np.sin(theta)
    return (np.sin((1-t)*theta) * q0 + np.sin(t*theta) * q1) / sin_theta

def euler_to_rotmat(yaw, pitch, roll):
    """ZYX Euler angles (in radians) to 3x3 rotation matrix."""
    cy, sy = np.cos(yaw),   np.sin(yaw)
    cp, sp = np.cos(pitch), np.sin(pitch)
    cr, sr = np.cos(roll),  np.sin(roll)
    Rz = np.array([[cy,-sy,0],[sy,cy,0],[0,0,1]])
    Ry = np.array([[cp,0,sp],[0,1,0],[-sp,0,cp]])
    Rx = np.array([[1,0,0],[0,cr,-sr],[0,sr,cr]])
    return Rz @ Ry @ Rx

# ======== Transform packing ========

def pack_transform(R, T):
    """Pack 3x3 rotation R and (3,) translation T into 3x float4 = (12,) float32 array.
    Layout matches the OpenCL Transform struct: row0=(R[0,:], Tx), row1=(R[1,:], Ty), row2=(R[2,:], Tz)."""
    out = np.zeros(12, dtype=np.float32)
    out[0:3] = R[0,:]; out[3]  = T[0]
    out[4:7] = R[1,:]; out[7]  = T[1]
    out[8:11]= R[2,:]; out[11] = T[2]
    return out

def pack_transforms(rotmats, translations):
    """Pack arrays of rotation matrices and translations into flat transform buffer.
    rotmats:      (N,3,3) or list of (3,3)
    translations: (N,3) or list of (3,)
    Returns: (N,12) float32 array."""
    N = len(rotmats)
    buf = np.zeros((N, 12), dtype=np.float32)
    for i in range(N):
        buf[i] = pack_transform(rotmats[i], translations[i])
    return buf

# ======== Scan generators ========
# Each returns (transforms, scan_info) where:
#   transforms: (N, 12) float32 array ready for GPU
#   scan_info:  dict with metadata (coordinates, labels, shape, etc.)

def scan_z_approach(pos_xy, z_range, R=None, nz=50):
    """1D scan: vary z at fixed (x,y) and rotation.
    pos_xy:  (x, y) lateral position
    z_range: (z_min, z_max) 
    R:       3x3 rotation matrix (default: identity)
    nz:      number of z points"""
    if R is None: R = np.eye(3)
    zs = np.linspace(z_range[0], z_range[1], nz)
    rots = [R] * nz
    trans = [np.array([pos_xy[0], pos_xy[1], z]) for z in zs]
    return pack_transforms(rots, trans), {'type': '1D_z', 'z': zs, 'shape': (nz,)}

def scan_lateral_2d(z, x_range, y_range, R=None, nx=50, ny=50):
    """2D scan: vary (x,y) at fixed z and rotation.
    Returns transforms for nx*ny poses, row-major (y changes fastest)."""
    if R is None: R = np.eye(3)
    xs = np.linspace(x_range[0], x_range[1], nx)
    ys = np.linspace(y_range[0], y_range[1], ny)
    rots, trans = [], []
    for x in xs:
        for y in ys:
            rots.append(R)
            trans.append(np.array([x, y, z]))
    return pack_transforms(rots, trans), {'type': '2D_xy', 'x': xs, 'y': ys, 'shape': (nx, ny)}

def scan_xz_slice(y, x_range, z_range, R=None, nx=50, nz=50):
    """2D scan: vary (x,z) at fixed y — vertical slice."""
    if R is None: R = np.eye(3)
    xs = np.linspace(x_range[0], x_range[1], nx)
    zs = np.linspace(z_range[0], z_range[1], nz)
    rots, trans = [], []
    for ix in range(nx):
        for iz in range(nz):
            rots.append(R)
            trans.append(np.array([xs[ix], y, zs[iz]]))
    return pack_transforms(rots, trans), {'type': '2D_xz', 'x': xs, 'z': zs, 'shape': (nx, nz)}

def scan_rotation_1d(pos, axis, angle_range, nrot=36):
    """1D scan: rotate around axis at fixed position.
    axis:        (3,) rotation axis
    angle_range: (a_min, a_max) in radians"""
    angles = np.linspace(angle_range[0], angle_range[1], nrot)
    rots = [quat_to_rotmat(quat_from_axis_angle(axis, a)) for a in angles]
    trans = [np.array(pos)] * nrot
    return pack_transforms(rots, trans), {'type': '1D_rot', 'angles': np.degrees(angles), 'shape': (nrot,)}

def scan_rotation_z_2d(pos_xy, z_range, axis, angle_range, nz=30, nrot=36):
    """2D scan: rotation angle vs z-distance."""
    angles = np.linspace(angle_range[0], angle_range[1], nrot)
    zs = np.linspace(z_range[0], z_range[1], nz)
    rots, trans = [], []
    for a in angles:
        R = quat_to_rotmat(quat_from_axis_angle(axis, a))
        for z in zs:
            rots.append(R)
            trans.append(np.array([pos_xy[0], pos_xy[1], z]))
    return pack_transforms(rots, trans), {'type': '2D_rot_z', 'angles': np.degrees(angles), 'z': zs, 'shape': (nrot, nz)}

def scan_slerp_path(q0, q1, t0, t1, npts=50):
    """1D scan along quaternion SLERP path with linear position interpolation.
    q0, q1: start/end quaternions (w,x,y,z)
    t0, t1: start/end translation (3,)"""
    ts = np.linspace(0, 1, npts)
    rots, trans = [], []
    t0, t1 = np.asarray(t0), np.asarray(t1)
    for t in ts:
        q = quat_slerp(q0, q1, t)
        rots.append(quat_to_rotmat(q))
        trans.append(t0 + t * (t1 - t0))
    return pack_transforms(rots, trans), {'type': '1D_slerp', 't': ts, 'shape': (npts,)}

def scan_monte_carlo(pos_center, pos_spread, nsamples=1000, seed=42):
    """Random sampling of position and rotation around a center point.
    pos_center: (3,) center position
    pos_spread: (3,) half-width for uniform random shift in x,y,z"""
    rng = np.random.RandomState(seed)
    rots, trans = [], []
    pos_center = np.asarray(pos_center)
    pos_spread = np.asarray(pos_spread)
    for _ in range(nsamples):
        # Random rotation via random quaternion
        q = rng.randn(4)
        q /= np.linalg.norm(q)
        rots.append(quat_to_rotmat(q))
        trans.append(pos_center + pos_spread * (2*rng.rand(3) - 1))
    return pack_transforms(rots, trans), {'type': 'MC', 'shape': (nsamples,)}

def scan_multi_dof(dof_specs):
    """General N-dimensional scan over arbitrary DOF combinations.
    dof_specs: list of dicts, each with:
        'type': 'x'|'y'|'z'|'yaw'|'pitch'|'roll'
        'values': array of values for this DOF
    Builds Cartesian product of all DOF values.
    Returns transforms and scan_info with shape = product of all DOF sizes."""
    from itertools import product as cartprod
    dof_names  = [d['type'] for d in dof_specs]
    dof_values = [d['values'] for d in dof_specs]
    shape = tuple(len(v) for v in dof_values)
    defaults = {'x': 0.0, 'y': 0.0, 'z': 0.0, 'yaw': 0.0, 'pitch': 0.0, 'roll': 0.0}

    rots, trans = [], []
    for combo in cartprod(*dof_values):
        d = dict(defaults)
        for name, val in zip(dof_names, combo):
            d[name] = val
        R = euler_to_rotmat(d['yaw'], d['pitch'], d['roll'])
        T = np.array([d['x'], d['y'], d['z']])
        rots.append(R)
        trans.append(T)
    info = {'type': 'multi_dof', 'dof_names': dof_names, 'dof_values': dof_values, 'shape': shape}
    return pack_transforms(rots, trans), info
