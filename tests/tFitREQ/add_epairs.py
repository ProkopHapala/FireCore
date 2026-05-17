#!/usr/bin/env python3
"""
Add electron pairs and/or sigma holes to multi-frame XYZ dimer files.

Processes each frame separately, splitting at n0 (first n0 atoms = molecule A,
remaining = molecule B). Updates n0 in the header after adding atoms.

AtomTypes.dat is read to determine which atoms get epairs (nepair column)
and which have sigma-hole epair types. The epair distance (Lepair) and
sigma-hole distance are configurable.

Usage examples:
  # Process one file:
  python tests/tFitREQ/add_epairs.py -i scan.xyz -o scan_ep.xyz

  # Process all .xyz files in a directory:
  python tests/tFitREQ/add_epairs.py -i input_dir/ -o output_dir/

  # Custom distances and sigma holes:
  python tests/tFitREQ/add_epairs.py -i input_dir/ -o out_dir/ \\
      --mode both --lepair 1.2 --sigma-dist 0.6

  # Use a specific AtomTypes.dat:
  python tests/tFitREQ/add_epairs.py -i scan.xyz -o out.xyz \\
      --atypes tests/tFitREQ_PN/data/AtomTypes.dat

  # Show help:
  python tests/tFitREQ/add_epairs.py -h
"""
import argparse
import os
import re
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))
from pyBall import atomicUtils as au
from pyBall import elements
from pyBall.AtomicSystem import AtomicSystem


# ── Data structures for AtomTypes.dat / ElementTypes.dat ──

class AtomType:
    def __init__(self, name, parent_name="*", element_name="", epair_name="",
                 valence=0, nepair=0, npi=0, sym=0,
                 Ruff=0.0, RvdW=0.0, EvdW=0.0, Qbase=0.0, Hb=0.0):
        self.name = name
        self.parent_name = parent_name
        self.element_name = element_name
        self.epair_name = epair_name
        self.valence = valence
        self.nepair = nepair
        self.npi = npi
        self.sym = sym
        self.Ruff = Ruff
        self.RvdW = RvdW
        self.EvdW = EvdW
        self.Qbase = Qbase
        self.Hb = Hb


def read_atom_types(filepath):
    atom_types = {}
    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('#') or not line:
                continue
            parts = line.split()
            if len(parts) < 5:
                continue
            name = parts[0]
            at = AtomType(name=name, parent_name=parts[1], element_name=parts[2],
                          epair_name=parts[3])
            try:
                at.valence = int(parts[4])
                at.nepair = int(parts[5])
                at.npi = int(parts[6])
                at.sym = int(parts[7])
                at.Ruff = float(parts[8])
                at.RvdW = float(parts[9])
                at.EvdW = float(parts[10])
                at.Qbase = float(parts[11])
                at.Hb = float(parts[12])
            except (ValueError, IndexError):
                continue
            atom_types[name] = at
    return atom_types


# ── Core logic ──

def parse_n0_from_comment(comment):
    m = re.search(r"\bn0\s+(\d+)", comment)
    return int(m.group(1)) if m else None


def update_n0_in_comment(comment, new_n0):
    return re.sub(r"\bn0\s+\d+", f"n0 {new_n0}", comment)


def element_from_typename(tname):
    return tname.split("_")[0] if "_" in tname else tname


def count_epairs_from_type(tname, atom_types):
    """Return how many epairs this atom type should have, or 0 if unknown."""
    at = atom_types.get(tname)
    if at is None:
        return 0
    return max(0, at.nepair)


def has_sigma_hole_type(tname, atom_types):
    """Return True if this atom type is a sigma-hole donor (has epair_name, nepair==0)."""
    at = atom_types.get(tname)
    if at is None:
        return False
    return at.nepair == 0 and at.epair_name not in ("*", "", "E")


def _z_to_qs(z):
    try:
        return elements.ELEMENTS[z - 1][9]
    except (IndexError, KeyError):
        return 0.0


def _z_to_rs(z):
    try:
        return elements.ELEMENTS[z - 1][7]
    except (IndexError, KeyError):
        return 1.0


def build_fragment(apos, enames, qs, atom_types, lepair, sigma_dist, do_epairs, do_sigma):
    """
    Build an AtomicSystem for one fragment, add epairs and/or sigma holes.
    Returns (sys, n_added).
    """
    elem = []
    atypes_list = []
    valid = []
    for i, e in enumerate(enames):
        en = element_from_typename(e)
        if en in elements.ELEMENT_DICT:
            elem.append(en)
            atypes_list.append(elements.ELEMENT_DICT[en][0])
            valid.append(i)
        else:
            # Unknown type (e.g., 'E' dummy from previous run) — keep but use Z=0
            elem.append(en)
            atypes_list.append(0)
            valid.append(i)

    atypes = np.array(atypes_list, dtype=np.int32)
    elem_list = list(elem)

    sys = AtomicSystem(
        apos=apos.copy(),
        atypes=atypes,
        enames=elem_list,
        qs=qs.copy() if qs is not None else None,
        bPreinit=False,
    )
    if sys.qs is None:
        sys.qs = np.array([_z_to_qs(z) for z in sys.atypes])
    if sys.Rs is None:
        sys.Rs = np.array([_z_to_rs(z) for z in sys.atypes])
    sys.neighs()

    n_orig = len(sys.apos)
    n_added = 0

    # ── Electron pairs ──
    if do_epairs:
        # Override VALENCE_DICT in the AtomicSystem instance so it knows
        # which atoms get epairs from AtomTypes.dat (not just O/N).
        # VALENCE_DICT is a module-level dict; we temporarily replace it.
        from pyBall.AtomicSystem import VALENCE_DICT as _orig_vd

        valence_map = {}
        for i, tname in enumerate(enames):
            at = atom_types.get(tname)
            if at is not None and at.nepair > 0:
                elem_name = elem[i]
                nb = at.valence
                nsigma = len(sys.ngs[i]) if (sys.ngs is not None and i < len(sys.ngs)) else 0
                npi = nb - nsigma
                valence_map[elem_name] = (nb, at.nepair, npi)

        # Patch VALENCE_DICT so add_electron_pairs uses our data
        import pyBall.AtomicSystem as _asmod
        _backup = dict(_asmod.VALENCE_DICT)
        for ename, (nb, nep, _) in valence_map.items():
            _asmod.VALENCE_DICT[ename] = (nb, nep)

        # Need to also set npi per atom. add_electron_pairs/dd_epair uses npi from
        # difference between valence and sigma neighbors, which we already did above.
        # The VALENCE_DICT only gives (nb, nep), npi is derived.
        # So we just need the correct (nb, nep) entries.
        sys.add_electron_pairs(distance=lepair)
        _asmod.VALENCE_DICT.clear()
        _asmod.VALENCE_DICT.update(_backup)

        n_added += len(sys.apos) - n_orig

    # ── Sigma holes ──
    if do_sigma:
        for i in range(n_orig):
            if not has_sigma_hole_type(enames[i], atom_types):
                continue
            neighs = (list(sys.ngs[i].keys())
                      if (sys.ngs is not None and i < len(sys.ngs)) else [])
            if len(neighs) != 1:
                continue
            j = neighs[0]
            direction = au.normalize(sys.apos[i] - sys.apos[j])
            sys.place_electron_pair(
                i, direction, distance=sigma_dist,
                ename="E", atype=200, qs=0.0, Rs=1.0,
            )
            n_added += 1

    # Restore original type names
    for i in range(len(enames)):
        sys.enames[i] = enames[i]

    return sys, n_added


def process_frame(es, apos, qs, rs, comment, atom_types, lepair, sigma_dist, do_epairs, do_sigma):
    comment = comment.strip()
    n0 = parse_n0_from_comment(comment)
    if n0 is None:
        return None, None
    natoms = len(es)
    if n0 > natoms:
        return None, None

    esA, esB = es[:n0], es[n0:]
    aposA, aposB = apos[:n0].copy(), apos[n0:].copy()
    qsA = qs[:n0].copy() if qs is not None else None
    qsB = qs[n0:].copy() if qs is not None else None

    sysA, nA = build_fragment(aposA, esA, qsA, atom_types, lepair, sigma_dist, do_epairs, do_sigma)
    sysB, nB = build_fragment(aposB, esB, qsB, atom_types, lepair, sigma_dist, do_epairs, do_sigma)

    new_es = list(sysA.enames) + list(sysB.enames)
    new_apos = np.vstack([sysA.apos, sysB.apos])
    new_qs = (np.concatenate([sysA.qs, sysB.qs])
              if (sysA.qs is not None and sysB.qs is not None) else None)
    new_comment = update_n0_in_comment(comment, n0 + nA)
    return (new_es, new_apos, new_qs, None, new_comment), (nA, nB)


def process_one_file(inpath, outpath, atom_types, lepair, sigma_dist, do_epairs, do_sigma, verbose, simple_names):
    trj = au.load_xyz_movie(inpath)
    if not trj:
        return 0, 0, 0
    if verbose:
        print(f"    {len(trj)} frames")

    fout = open(outpath, "w")
    total_frames = 0
    total_A = 0
    total_B = 0
    for es, apos, qs, rs, comment in trj:
        res = process_frame(es, apos, qs, rs, comment, atom_types, lepair, sigma_dist, do_epairs, do_sigma)
        if res is None:
            continue
        (new_es, new_apos, new_qs, _, new_comment), (nA, nB) = res
        out_es = [element_from_typename(e) for e in new_es] if simple_names else new_es
        au.writeToXYZ(fout, out_es, new_apos, qs=new_qs, comment=new_comment, bHeader=True)
        total_frames += 1
        total_A += nA
        total_B += nB
        if verbose:
            print(f"      frame -> {len(new_es)} atoms, n0→{parse_n0_from_comment(new_comment)}, +{nA}/+{nB}")
    fout.close()
    return total_frames, total_A, total_B


def find_data_path(script_dir):
    candidates = [
        os.path.join(script_dir, "..", "..", "tests", "tFitREQ_PN", "data"),
        os.path.join(script_dir, "..", "..", "tests", "tFitREQ", "data"),
        os.path.join(script_dir, "data"),
    ]
    for d in candidates:
        d = os.path.abspath(d)
        if os.path.isfile(os.path.join(d, "AtomTypes.dat")):
            return d
    return None


# ── CLI ──

def main():
    parser = argparse.ArgumentParser( description="Add electron pairs and sigma holes to multi-frame XYZ dimer files", formatter_class=argparse.RawTextHelpFormatter, )
    parser.add_argument("-i", "--input", required=True, help="Input .xyz file or directory of .xyz files")
    parser.add_argument("-o", "--output", default="output",  help="Output file or directory (default: output)")
    parser.add_argument("--mode", choices=["epairs", "sigma", "both"], default="epairs",  help="What to add (default: epairs)")
    parser.add_argument("--lepair", type=float, default=1.0,  help="Epair distance from host in Angstrom (default: 1.0)")
    parser.add_argument("--sigma-dist", type=float, default=0.5,  help="Sigma-hole distance from host in Angstrom (default: 0.5)")
    parser.add_argument("--atypes", default=None,  help="Path to AtomTypes.dat (auto-searched if not given)")
    parser.add_argument("--simple-names", action="store_true",   help="Output simple element names (N, O, H) instead of full type names (N_3, O_2, H_O)")
    parser.add_argument("--verbose", "-v", action="store_true")
    args = parser.parse_args()

    do_epairs = args.mode in ("epairs", "both")
    do_sigma = args.mode in ("sigma", "both")

    # ── Find AtomTypes.dat ──
    script_dir = os.path.dirname(os.path.abspath(__file__))
    if args.atypes:
        atypes_path = args.atypes
    else:
        data_dir = find_data_path(script_dir)
        if data_dir is None:
            print("Error: cannot find AtomTypes.dat. Use --atypes to specify path.")
            sys.exit(1)
        atypes_path = os.path.join(data_dir, "AtomTypes.dat")

    if not os.path.isfile(atypes_path):
        print(f"Error: {atypes_path} not found")
        sys.exit(1)

    atom_types = read_atom_types(atypes_path)
    print(f"Loaded {len(atom_types)} atom types from {atypes_path}")

    # ── Collect input files ──
    inpath = args.input
    outpath = args.output

    if os.path.isfile(inpath):
        files = [(inpath, outpath)]
    elif os.path.isdir(inpath):
        os.makedirs(outpath, exist_ok=True)
        files = []
        for fn in sorted(os.listdir(inpath)):
            if not fn.endswith(".xyz"):
                continue
            if "_bak" in fn or "_Epairs" in fn:
                continue
            files.append((os.path.join(inpath, fn), os.path.join(outpath, fn)))
        if not files:
            print(f"No .xyz files found in {inpath}")
            sys.exit(1)
        print(f"Found {len(files)} .xyz files in {inpath}")
    else:
        print(f"Error: {inpath} is not a file or directory")
        sys.exit(1)

    # ── Process ──
    total_frames = 0
    total_A = 0
    total_B = 0

    for src, dst in files:
        print(f"Processing: {os.path.basename(src)}")
        nf, nA, nB = process_one_file(
            src, dst, atom_types, args.lepair, args.sigma_dist,
            do_epairs, do_sigma, args.verbose, args.simple_names,
        )
        if nf > 0:
            print(f"  → {nf} frames, +{nA}/+{nB} (molA/molB) → {dst}")
        total_frames += nf
        total_A += nA
        total_B += nB

    print(f"\nTotal: {total_frames} frames across all files")
    print(f"  Added atoms: molA={total_A}, molB={total_B}")
    print(f"  Mode: {args.mode}, Lepair={args.lepair}, sigma-dist={args.sigma_dist}")


if __name__ == "__main__":
    main()
