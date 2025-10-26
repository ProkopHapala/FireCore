import argparse
import sys
from pathlib import Path

import numpy as np

from ..AtomicSystem import AtomicSystem
from .. import atomicUtils as au
from .MMFFL import MMFFL


def _default_input() -> Path:
    base = Path(__file__).resolve().parents[2]
    return base / "tests" / "tDFT" / "data" / "mol" / "formaldehyde.mol2"


def _derive_output(input_path: Path, suffix: str = "_linearized") -> Path:
    stem = input_path.stem + suffix
    return input_path.with_name(stem + input_path.suffix)


def _collect_bonds(mol, mmffl) -> list[tuple[int, int]]:
    bond_set = set()
    if mol.bonds is not None:
        for a, b in mol.bonds:
            bond_set.add(tuple(sorted((int(a), int(b)))))
    for (a, b, _l0, _k, _tag) in mmffl.linear_bonds:
        bond_set.add(tuple(sorted((int(a), int(b)))))
    bonds = sorted(bond_set)
    return bonds


def _build_linearized(mmffl: MMFFL, mol: AtomicSystem):
    apos = np.asarray(mmffl.apos[:mmffl.natoms, :3], dtype=float)
    orig_natoms = int(mol.apos.shape[0])
    extra = int(mmffl.natoms - orig_natoms)

    base_names = [str(e) for e in mol.enames]
    if len(base_names) < orig_natoms:
        base_names = list(base_names) + ["X"] * (orig_natoms - len(base_names))
    for _ in range(extra):
        base_names.append("Du")

    bonds = _collect_bonds(mol, mmffl)
    return apos, base_names, bonds, extra

'''
#python -m pyBall.OCL.mmffl_export --input tests/tDFT/data/mol/formaldehyde.mol2 --output formaldehyde_linearized.mol2 --two-pi-dummies --L-pi 1.2
python -m pyBall.OCL.mmffl_export 
'''
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate linearized MMFF topology with pi dummies.")
    parser.add_argument("--input",     type=Path,  default=_default_input(), help="Input molecule file (mol2/mol/xyz).")
    parser.add_argument("--output",    type=Path,  default=None,  help="Output MOL2 path for linearized topology.")
    parser.add_argument("--Lpi",       type=float, default=1.0,   help="Distance for pi dummy atoms (default 1.0).")
    parser.add_argument("--twopi",     type=int,   default=0,       help="Place pi dummies on both sides (default off).")
    parser.add_argument("--Kang",      type=float, default=0.0,    help="Stiffness for angle-replacement bonds.")
    parser.add_argument("--KpiA",      type=float, default=0.0,    help="Stiffness for host-pi bonds.")
    parser.add_argument("--KpiB",      type=float, default=0.0,    help="Stiffness for pi orthogonality bonds.")
    parser.add_argument("--useuff",    type=int,   default=0,     help="Use UFF bond parameters when available.")
    parser.add_argument("--comment",   type=str,   default="",    help="Header comment stored in the output mol2.")
    parser.add_argument("--verbosity", type=int,   default=0,   help="Verbosity forwarded to MMFFL constructor.")
    args = parser.parse_args()

    input_path = args.input.resolve()
    if not input_path.exists():
        print(f"ERROR: input file '{input_path}' not found", file=sys.stderr)
        sys.exit(1)

    output_path = args.output
    if output_path is None:
        output_path = _derive_output(input_path)
    else:
        output_path = output_path.resolve()

    mol = AtomicSystem(fname=str(input_path))
    if not hasattr(mol, "natoms") or mol.natoms is None:
        mol.natoms = int(mol.apos.shape[0])
    if getattr(mol, "enames", None) is not None and not isinstance(mol.enames, list):
        mol.enames = list(mol.enames)

    mmffl = MMFFL(
        L_pi=float(args.Lpi),
        two_pi_dummies=bool(args.twopi),
        Kang=float(args.Kang),
        Kpi_host=float(args.KpiA),
        Kpi_orth=float(args.KpiB),
        verbosity=args.verbosity,
    )

    use_uff = bool(args.useuff)
    mmffl.build_linearized(mol, bUFF=use_uff)

    apos, enames, bonds, n_dummies = _build_linearized(mmffl, mol)

    au.save_mol2(str(output_path), enames, apos, bonds, qs=None, comment=args.comment)

    print(f"Input atoms: {mol.apos.shape[0]}")
    print(f"Dummy atoms added: {n_dummies}")
    print(f"Total atoms exported: {len(enames)}")
    print(f"Total bonds exported: {len(bonds)}")
    print(f"Linearized topology written to: {output_path}")
