charges={
   "HF":  { "H":0.41821209, "F" :-0.41821209},
   "HCl": { "H":0.21775824, "Cl":-0.21775824},
   "HBr": { "H":0.16290409, "Br":-0.16290409},
}

import os
import sys
from typing import List, Optional, Tuple
import re
import argparse


def append_charges_to_xyz_lines(
    lines: List[str],
    system_name: str,
    angle_atom: Optional[int] = 3,
    angle_pair: Optional[Tuple[int, int]] = None,
) -> List[str]:
    """Append RESP charges (from `charges`) to each atom line of every XYZ frame.

    Assumes frames of the form (repeated):
      natoms\n
      comment (e.g., energy)\n
      natoms lines: "Sym x y z"\n

    We preserve natoms lines EXACTLY. We REWRITE the comment to the format:
      '# n0 <n0> Etot <Etot> x0 <r> z <ang> <system_name>'
    where n0 = natoms//2, Etot is parsed from the old comment (numeric),
    r is the COM distance between halves, and ang is the angle [deg] between the
    axis of the second half and the COM separation vector. For atom lines, we append
    a fourth numeric field: the RESP partial charge. For diatomics HX (X in F/Cl/Br),
    we infer the H charge from the most recent halogen line in the frame.
    """
    out: List[str] = []
    i = 0
    n = len(lines)
    while i < n:
        line = lines[i]
        s = line.strip()
        if s == "":
            out.append(line)
            i += 1
            continue

        # Try to parse natoms; if not an int, pass-through and continue
        try:
            nat = int(s)
        except ValueError:
            out.append(line)
            i += 1
            continue

        # Frame header: natoms and comment
        out.append(line)  # natoms (preserve exactly)
        if i + 1 >= n:
            raise ValueError("Malformed XYZ: missing comment after natoms at line %d" % (i + 1))
        old_comment = lines[i + 1]

        # Parse atom block first to compute r/angle and build modified atom lines with charges
        coords: List[List[float]] = []
        types_block: List[str] = []
        mod_atom_lines: List[str] = []
        current_mol = None  # 'HF' | 'HCl' | 'HBr'
        for j in range(nat):
            k = i + 2 + j
            if k >= n:
                raise ValueError("Malformed XYZ: incomplete atom block, expected %d atoms" % nat)
            raw = lines[k].rstrip("\n")
            parts = raw.split()
            if len(parts) < 4:
                raise ValueError(f"Malformed atom line (need symbol + 3 coords): '{raw}'")
            sym = parts[0]
            try:
                x, y, z = map(float, parts[1:4])
            except Exception:
                raise ValueError(f"Malformed coordinates in line: '{raw}'")
            types_block.append(sym)
            coords.append([x, y, z])

            # Determine molecule and charge
            if sym in ("F", "Cl", "Br"):
                current_mol = "H" + sym  # HF/HCl/HBr
                q = charges[current_mol][sym]
            elif sym == "H":
                if current_mol is None:
                    raise ValueError("Encountered H before halogen; cannot infer molecule for charge")
                q = charges[current_mol]["H"]
            else:
                raise ValueError(f"Unsupported atom symbol '{sym}' — expected H, F, Cl, or Br")

            # Rebuild line from symbol and first 3 coordinate tokens to avoid duplicating existing charge
            mod_atom_lines.append(f"{sym}  {parts[1]}  {parts[2]}  {parts[3]}  {q:.8f}\n")

        # Compute r and angle from coords
        # Radial separation r := distance between COMs of first and second halves
        h = nat // 2
        import math
        cA = [sum(p[d] for p in coords[:h]) / h for d in range(3)]
        cB = [sum(p[d] for p in coords[h:]) / (nat - h) for d in range(3)]
        R = [cB[d] - cA[d] for d in range(3)]
        r = math.sqrt(R[0] * R[0] + R[1] * R[1] + R[2] * R[2])
        # Angle definition
        # Priority: angle_pair (j - i) in xz-plane -> atan2(dx, dz); else angle_atom -> atan2(x, z)
        ang = float('nan')
        try:
            if angle_pair is not None:
                i0, i1 = angle_pair
                if not (0 <= i0 < nat and 0 <= i1 < nat):  raise IndexError("angle_pair indices out of range")
                vx = coords[i1][0] - coords[i0][0]
                vz = coords[i1][2] - coords[i0][2]
                ang = math.degrees(math.atan2(vx, vz))
            else:
                ia = angle_atom if angle_atom is not None else 3
                if not (0 <= ia < nat):
                    # Fallback to last atom if default is invalid
                    ia = nat - 1
                x, z = coords[ia][0], coords[ia][2]
                ang = math.degrees(math.atan2(x, z))
        except Exception:
            ang = float('nan')
        n0 = h
        # Extract Etot numeric from old comment if present (pattern: '# E = <number> [units]')
        m = re.search(r"E\s*=\s*([+-]?(?:\d*\.\d+|\d+)(?:[Ee][+-]?\d+)?)", old_comment)
        Etot_str = m.group(1) if m else None
        if Etot_str is None:
            m2 = re.search(r"\bEtot\s+([^\s]+)", old_comment)
            if m2:
                v = m2.group(1)
                if v.lower() != "nan":
                    Etot_str = v
        if Etot_str is None:
            Etot_str = "nan"
        new_comment = f"# n0 {n0} Etot {Etot_str} x0 {r:5.2f} z {ang:5.2f} {system_name}\n"
        out.append(new_comment)

        # Append modified atom lines
        out.extend(mod_atom_lines)

        i += 2 + nat

    return out


def process_xyz_file(
    in_path: str,
    out_path: str,
    angle_atom: Optional[int] = 3,
    angle_pair: Optional[Tuple[int, int]] = None,
) -> None:
    with open(in_path, "r") as f:
        lines = f.readlines()
    system_name = os.path.splitext(os.path.basename(in_path))[0]
    new_lines = append_charges_to_xyz_lines(lines, system_name, angle_atom=angle_atom, angle_pair=angle_pair)
    with open(out_path, "w") as f:
        f.writelines(new_lines)


def process_dir(
    in_dir: str,
    out_dir: Optional[str] = None,
    angle_atom: Optional[int] = 3,
    angle_pair: Optional[Tuple[int, int]] = None,
) -> None:
    if out_dir is None:
        out_dir = os.path.join(in_dir, "porcessed")  # keep historical name
    os.makedirs(out_dir, exist_ok=True)
    files = [
        os.path.join(in_dir, fn)
        for fn in sorted(os.listdir(in_dir))
        if fn.lower().endswith(".xyz")
    ]
    if not files:
        print(f"No .xyz files found in {in_dir}")
        return
    for in_fp in files:
        out_fp = os.path.join(out_dir, os.path.basename(in_fp))
        process_xyz_file(in_fp, out_fp, angle_atom=angle_atom, angle_pair=angle_pair)
        print(f"Wrote: {out_fp}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Append RESP charges to XYZ files and add r/angle metadata.")
    parser.add_argument("input_dir", nargs="?", help="Directory with input .xyz files", default=os.path.abspath(os.path.join(os.path.dirname(__file__), "HHalogens")))
    parser.add_argument("--out", "--output-dir", dest="out_dir", default=None, help="Output directory (default: <input_dir>/porcessed)")
    parser.add_argument("--angle-atom", dest="angle_atom", type=int, default=3, help="Atom index (0-based) to define angle as atan2(x,z); default=3")
    parser.add_argument("--angle-pair", dest="angle_pair", type=int, nargs=2, default=None,  help="Two atom indices i j (0-based) to define angle from vector r_j - r_i in xz-plane")
    args = parser.parse_args()

    in_dir = args.input_dir
    if not os.path.isabs(in_dir):
        in_dir = os.path.abspath(in_dir)
    process_dir(in_dir, out_dir=args.out_dir, angle_atom=args.angle_atom, angle_pair=tuple(args.angle_pair) if args.angle_pair else None)