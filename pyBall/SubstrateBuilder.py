#!/usr/bin/env python3
"""
Substrate builder for ionic crystals (NaCl-type simple cubic with alternating ions).
Generates flat slabs and step edges, returning AtomicSystem instances.

The NaCl structure here is simple cubic with alternating Na/Cl ions,
nearest-neighbor distance a = 5.6413/2 = 2.82 A along cartesian axes.
This is NOT the conventional FCC rock-salt unit cell (a=5.64A).
"""

import numpy as np

def gen_nacl_slab(a=2.82065, nx=10, ny=10, nz=3, Q0=0.7, bCharges=True):
    """Generate flat NaCl slab (simple cubic, alternating Na/Cl).
    
    Args:
        a:  nearest-neighbor distance (Na-Cl bond length) in Angstrom
        nx, ny, nz: number of atoms along x, y, z
        Q0: magnitude of ionic charge (sign alternates)
        bCharges: if True, assign +-Q0 charges
    Returns:
        AtomicSystem with apos, enames, atypes, qs, lvec set
    """
    assert nx > 0 and ny > 0 and nz > 0, f"gen_nacl_slab: invalid nx,ny,nz = {nx},{ny},{nz}"
    natoms = nx * ny * nz
    apos   = np.empty((natoms, 3))
    enames = []
    qs     = np.empty(natoms)
    atypes = np.empty(natoms, dtype=np.int32)
    idx = 0
    for iz in range(nz):
        for ix in range(nx):
            for iy in range(ny):
                x, y, z = ix * a, iy * a, iz * a
                i = ix + iy + iz
                if i % 2 == 0:
                    enames.append('Na'); atypes[idx] = 11; qs[idx] = +Q0
                else:
                    enames.append('Cl'); atypes[idx] = 17; qs[idx] = -Q0
                apos[idx] = (x, y, z)
                idx += 1
    if not bCharges: qs = None
    lvec = np.array([[nx*a, 0, 0], [0, ny*a, 0], [0, 0, nz*a]])
    from . import AtomicSystem as AS
    return AS.AtomicSystem(apos=apos, atypes=atypes, enames=enames, lvec=lvec, qs=qs, bPreinit=False)


def gen_nacl_step(a=2.82065, nx=15, ny=8, nz=3, Q0=0.7, bCharges=True,
                  nsteps=1, step_dir=0):
    """Generate NaCl slab with step edge(s).
    
    The step is created by a linear skew of z vs x (tilt), plus a discrete
    height offset (+a) for atoms past the midline. Multiple steps split the
    slab into nsteps+1 terraces.

    Args:
        a:  nearest-neighbor distance
        nx, ny, nz: grid size
        Q0: charge magnitude
        bCharges: assign charges
        nsteps: number of step edges (terraces = nsteps+1)
        step_dir: 0 = step along x, 1 = step along y
    Returns:
        AtomicSystem
    """
    assert nx > 0 and ny > 0 and nz > 0, f"gen_nacl_step: invalid nx,ny,nz = {nx},{ny},{nz}"
    assert nsteps >= 1, f"gen_nacl_step: nsteps must be >= 1, got {nsteps}"
    natoms = nx * ny * nz
    apos   = np.empty((natoms, 3))
    enames = []
    qs     = np.empty(natoms)
    atypes = np.empty(natoms, dtype=np.int32)

    if step_dir == 0:
        nL = nx
    else:
        nL = ny
    Ldir    = a * nL
    terrace = Ldir / (nsteps + 1)   # width of each terrace along step_dir
    ax      = -1.0 / nL             # tilt slope

    idx = 0
    for iz in range(nz):
        for ix in range(nx):
            for iy in range(ny):
                x, y, z = ix * a, iy * a, iz * a
                i = ix + iy + iz

                if step_dir == 0:
                    coord = x
                else:
                    coord = y

                # apply linear skew (tilt)
                z_ = z + ax * coord
                if step_dir == 0:
                    x -= ax * z; z = z_
                else:
                    y -= ax * z; z = z_

                # discrete step offsets
                step_idx = int(coord / terrace)
                if step_idx > nsteps: step_idx = nsteps
                if step_idx > 0:
                    z += step_idx * a
                    i += step_idx

                if i % 2 == 0:
                    enames.append('Na'); atypes[idx] = 11; qs[idx] = +Q0
                else:
                    enames.append('Cl'); atypes[idx] = 17; qs[idx] = -Q0
                apos[idx] = (x, y, z)
                idx += 1

    if not bCharges: qs = None
    lvec = np.array([[nx*a, 0, 0], [0, ny*a, 0], [0, 0, nz*a]])
    from . import AtomicSystem as AS
    return AS.AtomicSystem(apos=apos, atypes=atypes, enames=enames, lvec=lvec, qs=qs, bPreinit=False)
