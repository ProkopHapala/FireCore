#!/usr/bin/env python3

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import sys
import time

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.append(str(ROOT))
sys.path.append(str(ROOT / "tests" / "tGridFF"))

import all_scan
from pyBall import MMFF as mmff


def _parse_vec3(text: str) -> tuple[float, float, float]:
    values = tuple(float(part) for part in text.split(","))
    if len(values) != 3:
        raise ValueError(f"Expected 3 comma-separated values, got '{text}'")
    return values


def parse_args():
    parser = argparse.ArgumentParser(description="Run a headless GridFF scan and save raw energies/forces")
    parser.add_argument("--molecule", required=True, help="Molecule path without or with .xyz")
    parser.add_argument("--substrate", required=True, help="Surface path without or with .xyz")
    parser.add_argument("--output-dir", required=True, help="Directory for scan outputs")
    parser.add_argument("--scan-type", choices=["total", "morse", "coulomb", "pauli", "london"], default="total")
    parser.add_argument("--nscan", type=int, default=200)
    parser.add_argument("--span-min", type=float, default=2.5)
    parser.add_argument("--span-max", type=float, default=15.0)
    parser.add_argument("--scan-dir", default="0,0,1")
    parser.add_argument("--scan-origin", default="0,0,0")
    return parser.parse_args()


def _strip_xyz_suffix(path_text: str) -> str:
    path = Path(path_text)
    if path.suffix == ".xyz":
        return str(path.with_suffix(""))
    return str(path)


def main():
    args = parse_args()
    os.environ["GRIDFF_SKIP_PLOTS"] = "1"
    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    molecule = _strip_xyz_suffix(args.molecule)
    substrate = _strip_xyz_suffix(args.substrate)
    direction = _parse_vec3(args.scan_dir)
    origin = _parse_vec3(args.scan_origin)

    t0 = time.perf_counter()
    mmff.init(xyz_name=molecule, surf_name=substrate, bUFF=True, bSimple=True)
    mmff.getBuffs()
    init_time = time.perf_counter() - t0

    scan_prefix = all_scan.SCAN_TYPES[args.scan_type]()
    mol_name = Path(molecule).name

    ts = np.linspace(args.span_min, args.span_max, args.nscan, endpoint=False)
    poss = np.zeros((args.nscan, 3), dtype=np.float64)
    poss[:, 0] = origin[0] + ts * direction[0]
    poss[:, 1] = origin[1] + ts * direction[1]
    poss[:, 2] = origin[2] + ts * direction[2]

    t1 = time.perf_counter()
    energies, forces, positions = mmff.scan_rigid_uff(poss, bF=True, bP=True)
    scan_time = time.perf_counter() - t1

    energies = np.asarray(energies, dtype=np.float64)
    forces = np.asarray(forces, dtype=np.float64)
    positions = np.asarray(positions, dtype=np.float64)
    energies_norm = energies - energies[-1]

    data_path = out_dir / f"{mol_name}_{scan_prefix}.dat"
    np.savetxt(
        data_path,
        np.column_stack((ts, energies_norm)),
        header="ts\tEnergy(eV)",
        comments="# ",
    )
    np.savez(
        out_dir / f"{mol_name}_{scan_prefix}.npz",
        ts=ts,
        energies=energies,
        energies_norm=energies_norm,
        forces=forces,
        positions=positions,
        scan_positions=poss,
    )

    payload = {
        "molecule": molecule,
        "substrate": substrate,
        "scan_type": args.scan_type,
        "scan_prefix": scan_prefix,
        "nscan": int(args.nscan),
        "span": [float(args.span_min), float(args.span_max)],
        "direction": list(direction),
        "origin": list(origin),
        "timing_s": {
            "init": float(init_time),
            "scan": float(scan_time),
            "total": float(init_time + scan_time),
        },
    }
    with (out_dir / f"{mol_name}_{scan_prefix}_metrics.json").open("w", encoding="utf8") as fout:
        json.dump(payload, fout, indent=2)

    sys.stdout.flush()
    os._exit(0)


if __name__ == "__main__":
    main()
