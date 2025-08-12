charges={
   "HF":  { "H":0.41821209, "F" :-0.41821209},
   "HCl": { "H":0.21775824, "Cl":-0.21775824},
   "HBr": { "H":0.16290409, "Br":-0.16290409},
}

import os
import sys
from typing import List
import re


def append_charges_to_xyz_lines(lines: List[str], system_name: str) -> List[str]:
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
        # Using first half vs second half as fragments
        h = nat // 2
        import math
        cA = [sum(p[d] for p in coords[:h]) / h for d in range(3)]
        cB = [sum(p[d] for p in coords[h:]) / (nat - h) for d in range(3)]
        R = [cB[d] - cA[d] for d in range(3)]
        r = math.sqrt(R[0] * R[0] + R[1] * R[1] + R[2] * R[2])
        vB = [coords[nat - 1][d] - coords[h][d] for d in range(3)]  # axis of B: last - first in second half
        # angle between vB and R
        def dot(a, b):
            return a[0]*b[0] + a[1]*b[1] + a[2]*b[2]
        def norm(a):
            return math.sqrt(dot(a, a))
        nb = norm(vB)
        if r < 1e-12 or nb < 1e-12:
            ang = float('nan')
        else:
            # Signed angle: atan2(|vB×R|, vB·R) with sign from (vB×R)_z
            dotvr = dot(vB, R)
            cx = vB[1]*R[2] - vB[2]*R[1]
            cy = vB[2]*R[0] - vB[0]*R[2]
            cz = vB[0]*R[1] - vB[1]*R[0]
            sn = math.sqrt(cx*cx + cy*cy + cz*cz)
            a0 = math.degrees(math.atan2(sn, dotvr))  # [0,180]
            ang = math.copysign(a0, cz)               # [-180,180]
        # Format fields
        r_str = f"{r:05.2f}"
        ang_str = str(int(round(ang))) if math.isfinite(ang) else "nan"
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
        new_comment = f"# n0 {n0} Etot {Etot_str} x0 {r_str} z {ang_str} {system_name}\n"
        out.append(new_comment)

        # Append modified atom lines
        out.extend(mod_atom_lines)

        i += 2 + nat

    return out


def process_xyz_file(path: str) -> None:
    with open(path, "r") as f:
        lines = f.readlines()
    system_name = os.path.splitext(os.path.basename(path))[0]
    new_lines = append_charges_to_xyz_lines(lines, system_name)
    with open(path, "w") as f:
        f.writelines(new_lines)


def process_dir(dirpath: str) -> None:
    files = [
        os.path.join(dirpath, fn)
        for fn in sorted(os.listdir(dirpath))
        if fn.lower().endswith(".xyz")
    ]
    if not files:
        print(f"No .xyz files found in {dirpath}")
        return
    for fp in files:
        process_xyz_file(fp)
        print(f"Updated: {fp}")


if __name__ == "__main__":
    if len(sys.argv) > 1:
        target_dir = sys.argv[1]
        if not os.path.isabs(target_dir):
            target_dir = os.path.abspath(target_dir)
    else:
        # Default to HHalogens/ next to this script
        target_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), "HHalogens"))
    process_dir(target_dir)