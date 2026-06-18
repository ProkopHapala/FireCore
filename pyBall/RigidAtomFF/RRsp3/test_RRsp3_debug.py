import os, sys, re, subprocess
import numpy as np

_THIS_DIR = os.path.abspath(os.path.dirname(__file__))
_SHARED_DIR = os.path.abspath(os.path.join(_THIS_DIR, '..', 'shared'))
for _p in (_THIS_DIR, _SHARED_DIR):
    if _p not in sys.path:
        sys.path.insert(0, _p)

from XPTB_utils import pack_molecules_contiguous, make_h2o_geometry, masses_from_elems
from RRsp3 import build_neighs_bk_from_bonds, make_bk_slots_clustered, make_exclusions_1st_2nd
from RRsp3_utils import make_ports_from_neighs

# Defines expected interaction outcomes for synthetic verification.

# Format: (a, b, action)
# Actions: "EXCLUDE_BOND", "EXCLUDE_ANGLE", "COLLIDE", "IGNORE_FAR"

INTERACTION_TRUTH_H2O2 = [
    # Molecule A (cluster 0): O(0), H(1), H(2)
    (0, 1, "EXCLUDE_BOND"),
    (0, 2, "EXCLUDE_BOND"),
    (1, 2, "EXCLUDE_ANGLE"),

    # Molecule B (cluster 1): O(64), H(65), H(66) in packed layout
    (64, 65, "EXCLUDE_BOND"),
    (64, 66, "EXCLUDE_BOND"),
    (65, 66, "EXCLUDE_ANGLE"),

    # Inter-molecular (to be configured by test geometry)
    # (0, 64, "COLLIDE"),
]

def expected_actions_for_two_h2o(group_size=64):
    # Using packed layout: group0 O,H,H -> 0,1,2 ; group1 -> 64,65,66
    truth = {}
    def add(a,b,act):
        truth[(min(a,b), max(a,b))] = act

    add(0,1,'EXCLUDE_BOND'); add(0,2,'EXCLUDE_BOND'); add(1,2,'EXCLUDE_ANGLE')
    add(group_size+0, group_size+1,'EXCLUDE_BOND'); add(group_size+0, group_size+2,'EXCLUDE_BOND'); add(group_size+1, group_size+2,'EXCLUDE_ANGLE')

    return truth


def validate_exclusion_mapping(neighs, excl1, excl2, truth):
    # For every expected EXCLUDE_* pair, ensure it appears in excl1 or excl2 for both endpoints
    neighs = np.asarray(neighs, dtype=np.int32)
    excl1 = np.asarray(excl1, dtype=np.int32)
    excl2 = np.asarray(excl2, dtype=np.int32)
    for (a,b), act in truth.items():
        if not act.startswith('EXCLUDE'):
            continue
        ok_ab = (b in excl1[a]) or (b in excl2[a])
        ok_ba = (a in excl1[b]) or (a in excl2[b])
        if not (ok_ab and ok_ba):
            raise RuntimeError(f"exclusion mapping missing for pair {a}-{b} act={act}: a_has={ok_ab} b_has={ok_ba}")


def main():
    group_size = 64

    pos_h2o, _ = make_h2o_geometry(add_angle=False)
    elems = ['O','H','H']
    bonds = [(0,1),(0,2)]
    nnode = 1

    # place molecule B close to A to force collisions and ghost discovery
    shift = np.array([1.5, 0.0, 0.0], dtype=np.float32)

    mols = [
        {'elems': elems, 'pos': pos_h2o, 'bonds': bonds, 'nnode': nnode},
        {'elems': elems, 'pos': pos_h2o + shift, 'bonds': bonds, 'nnode': nnode},
    ]

    packed = pack_molecules_contiguous(mols, group_size=group_size, nodes_first=True, pad_to_group=True)
    natoms = int(packed['natoms_total'])
    pos = packed['pos']
    elems_all = packed['elems']
    is_pad = packed['is_padding']
    real = ~is_pad

    # neighs/exclusions
    neighs, _ = build_neighs_bk_from_bonds(natoms, packed['bonds'], max_deg=4)
    excl1, excl2 = make_exclusions_1st_2nd(neighs)

    truth = expected_actions_for_two_h2o(group_size=group_size)
    validate_exclusion_mapping(neighs, excl1, excl2, truth)

    # ports/stiffness
    port_local_atoms, K_atoms = make_ports_from_neighs(pos, neighs, K=200.0)

    # masses / inv masses (padding atoms: invm=0)
    elems_real = [e for e in elems_all if e != 'X']
    m_real = masses_from_elems(elems_real)
    m = np.ones((natoms,), dtype=np.float32)
    m[real] = m_real
    invm = np.zeros_like(m)
    invm[real] = 1.0 / m[real]

    # radius
    rad = np.zeros((natoms,), dtype=np.float32)
    rad[real] = 1.0

    nnode_per_group = int(packed['nnode_group'][0])
    bkSlots = make_bk_slots_clustered(neighs, group_size=group_size, nnode_per_group=nnode_per_group, natoms=natoms)

    # Build with debug prints and collision-only
    build_opts = [
        "-DENABLE_DEBUG_PRINTS",
        "-DDEBUG_VERBOSITY=3",
        "-DDEBUG_GID_START=0",
        "-DDEBUG_GID_END=8",
        "-DENABLE_COLL=1",
        "-DENABLE_PORT=0",
    ]
    dbg_start = 0
    dbg_end = 8

    # Run a child python that instantiates RRsp3 and runs one step, capturing stdout
    runner = os.path.join(os.path.dirname(__file__), "test_RRsp3_debug_runner.py")
    cmd = [sys.executable, runner]
    env = os.environ.copy()
    env['RRSP3_BUILD_OPTS'] = " ".join(build_opts)

    # Pass data via npy temp files to avoid inline scripting
    tmpdir = os.path.join(os.path.dirname(__file__), "_tmp_rrsp3")
    os.makedirs(tmpdir, exist_ok=True)
    np.save(os.path.join(tmpdir, "pos.npy"), pos)
    np.save(os.path.join(tmpdir, "invm.npy"), invm)
    np.save(os.path.join(tmpdir, "quat.npy"), np.array([[0,0,0,1]], dtype=np.float32).repeat(natoms, axis=0))
    np.save(os.path.join(tmpdir, "rad.npy"), rad)
    np.save(os.path.join(tmpdir, "neighs.npy"), neighs)
    np.save(os.path.join(tmpdir, "excl1.npy"), excl1)
    np.save(os.path.join(tmpdir, "excl2.npy"), excl2)
    np.save(os.path.join(tmpdir, "port.npy"), port_local_atoms)
    np.save(os.path.join(tmpdir, "K.npy"), K_atoms)
    np.save(os.path.join(tmpdir, "bkSlots.npy"), bkSlots)
    env['RRSP3_TMPDIR'] = tmpdir
    env['RRSP3_NNODE_PER_GROUP'] = str(nnode_per_group)

    out = subprocess.run(cmd, cwd=os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..")), capture_output=True, text=True, env=env)
    print(out.stdout)
    if out.returncode != 0:
        print(out.stderr)
        raise SystemExit(out.returncode)

    # Parse COLL lines
    rgx = re.compile(r"^COLL: MeG=(\d+) OtherL=(\d+) OtherG=(\d+) Dist=([0-9eE+\.-]+) Action=(\w+)")
    seen = {}
    for line in out.stdout.splitlines():
        m = rgx.match(line.strip())
        if not m:
            continue
        me = int(m.group(1)); og = int(m.group(3)); act = m.group(5)
        k = (min(me,og), max(me,og))
        # keep strongest action: COLLIDE > SKIP_EXCL > TOO_FAR
        prev = seen.get(k, None)
        if prev is None:
            seen[k] = act
        else:
            rank = {'COLLIDE':3,'SKIP_EXCL':2,'TOO_FAR':1}
            if rank.get(act,0) > rank.get(prev,0):
                seen[k] = act

    # Validate key intramolecular exclusions were logged as SKIP_EXCL at least once.
    # Note: kernel printf is gated by DEBUG_GID_START/END, so we only require logs
    # for pairs where at least one endpoint is inside the printed GID window.
    for (a,b), exp in truth.items():
        if exp.startswith('EXCLUDE'):
            if not ((dbg_start <= a < dbg_end) or (dbg_start <= b < dbg_end)):
                continue
            got = seen.get((a,b), None)
            if got is None:
                raise RuntimeError(f"missing any log for expected excluded pair {a}-{b}")
            if got != 'SKIP_EXCL':
                raise RuntimeError(f"excluded pair {a}-{b} expected SKIP_EXCL got {got}")

    # At least one inter-molecular collision should happen due to close shift
    inter = []
    for (a,b), act in seen.items():
        if (a < group_size and b >= group_size) and act == 'COLLIDE':
            inter.append((a,b))
    if len(inter) == 0:
        raise RuntimeError("expected at least one inter-molecular COLLIDE, got none")
    print("[PASS] debug collision/exclusion logging looks consistent")
    print("inter COLLIDE examples:", inter[:8])


if __name__ == '__main__':
    main()
