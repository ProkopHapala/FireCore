#!/usr/bin/env python3
"""Utility CLI for loading molecules into MoleculeEditor2D.

Run with a backbone or standalone molecule, optionally overlay endgroups.

Examples
--------
python tests/tAttach/run_editor.py --molecule backbones/porphirin.mol2
python tests/tAttach/run_editor.py --backbone backbones/porphirin.mol2 --endgroup endgroups/guanine-SeCl.mol2
"""

import argparse
import sys
from pathlib import Path

sys.path.append(str(Path(__file__).resolve().parents[2]))

from pyBall.GUI.MoleculeEditor2D import MoleculeDocument, launch_editor
from pyBall.AtomicSystem import AtomicSystem


def parse_args(argv=None):
    p = argparse.ArgumentParser(description="Launch MoleculeEditor2D with preloaded systems.")
    src_group = p.add_mutually_exclusive_group(required=True)
    src_group.add_argument("--molecule", type=Path, help="Single molecule file (.xyz/.mol/.mol2)")
    src_group.add_argument("--backbone", type=Path, help="Backbone molecule file")
    p.add_argument("--endgroup", type=Path, action="append", default=[], help="Endgroup molecule file (repeatable)")
    p.add_argument("--out", type=Path, help="Optional path to save combined result before editing")
    return p.parse_args(argv)


def load_system(path: Path) -> AtomicSystem:
    sys = AtomicSystem(fname=str(path))
    sys.preinitialize_atomic_properties()
    return sys


def make_document(args) -> MoleculeDocument:
    if args.molecule:
        doc = MoleculeDocument()
        doc.load(args.molecule)
        return doc

    backbone = load_system(args.backbone)
    for eg_path in args.endgroup:
        endgroup = load_system(eg_path)
        backbone.attach_group_by_marker(endgroup, markerX="Se", markerY="Cl", _0=1)
    if args.out:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        suffix = args.out.suffix.lower()
        if suffix == ".mol2":
            backbone.save_mol2(str(args.out))
        elif suffix == ".xyz":
            backbone.saveXYZ(str(args.out))
        else:
            raise ValueError(f"Unsupported --out extension: {suffix}")
    return MoleculeDocument(system=backbone)


def main(argv=None):
    args = parse_args(argv)
    doc = make_document(args)
    if args.molecule:
        sys.exit(launch_editor(str(args.molecule)))
    if args.out:
        doc.save(str(args.out))
        sys.exit(launch_editor(str(args.out)))
    sys.exit(launch_editor())


if __name__ == "__main__":
    main()
