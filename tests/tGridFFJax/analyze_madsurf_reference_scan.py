#!/usr/bin/env python3
"""Analyze real DFT-labeled MAD-SURF Ag interface data against the local model.

This script does three things for one Ag+molecule composition:

1. Finds a near-1D interface subset inside the real DFT-labeled interface data by
   keeping lateral position and molecular pose nearly fixed and sorting by height.
2. Compares REF energy / REF forces against MAD-SURF model predictions on those
   exact DFT-labeled structures.
3. Builds an approximate interaction decomposition
      E_int ~= E_interface - E_slab - E_molecule
   using matching slab-only and pure-molecule reference data from the same bundle.

The decomposition is approximate unless the isolated-molecule geometry matches the
interface geometry exactly. The script reports that geometry mismatch explicitly.
"""

from __future__ import annotations

import argparse
from collections import Counter
import csv
from dataclasses import dataclass
from pathlib import Path
import sys

import matplotlib

matplotlib.use("Agg")

from ase.io import iread
from ase.data import atomic_numbers, covalent_radii
import matplotlib.pyplot as plt
from mace.calculators import MACECalculator
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.append(str(ROOT))

from pyBall.gridff_jax.utils import cartesian_to_fractional, cluster_z_layers, ensure_dir, save_json


@dataclass
class FrameRecord:
    index: int
    config_type: str
    composition: str
    atoms: object
    ref_energy: float
    ref_forces: np.ndarray
    slab_indices: np.ndarray
    mol_indices: np.ndarray
    mol_com: np.ndarray
    mol_centered: np.ndarray
    mol_pairwise: np.ndarray
    graph_signature: str
    height: float
    frac_xy: np.ndarray


def parse_args():
    base = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser(description="Compare real Ag DFT-labeled MAD-SURF frames against the local model")
    parser.add_argument(
        "--dataset-path",
        default=str(base / "mad-surf_data" / "full_train_test_std_config_types.extxyz"),
    )
    parser.add_argument(
        "--model-path",
        default=str(base / "mad-surf_data" / "models" / "full_dataset_config_weights" / "MACE_model.model"),
    )
    parser.add_argument("--device", default="cpu")
    parser.add_argument("--metal", default="Ag")
    parser.add_argument("--config-type", default="MD_interface_Ag")
    parser.add_argument("--slab-config-type", default="Slab_Ag")
    parser.add_argument("--composition", default="C5H5N1")
    parser.add_argument("--xy-cutoff", type=float, default=0.60, help="Max lateral COM drift in Angstrom for one extracted scan")
    parser.add_argument("--pose-cutoff", type=float, default=0.30, help="Max direct centered-coordinate RMSD in Angstrom for one extracted scan")
    parser.add_argument("--min-points", type=int, default=16)
    parser.add_argument("--max-points", type=int, default=120)
    parser.add_argument(
        "--pure-config-types",
        default="NMS,AISS_clusters",
        help="Comma-separated pure-molecule config types allowed for isolated-molecule matching",
    )
    parser.add_argument("--out-dir", default=str(base / "madsurf_reference_scan"))
    return parser.parse_args()


def composition_string(symbols):
    counts = Counter(symbols)
    return "".join(f"{key}{counts[key]}" for key in sorted(counts))


def split_indices(symbols, metal):
    symbols = np.asarray(symbols)
    slab_indices = np.where(symbols == metal)[0]
    mol_indices = np.where(symbols != metal)[0]
    return slab_indices.astype(int), mol_indices.astype(int)


def pairwise_descriptor(positions):
    positions = np.asarray(positions, dtype=float)
    if len(positions) < 2:
        return np.zeros(0, dtype=float)
    diffs = positions[:, None, :] - positions[None, :, :]
    dists = np.linalg.norm(diffs, axis=-1)
    iu = np.triu_indices(len(positions), 1)
    return np.sort(dists[iu])


def bond_adjacency(symbols, positions, scale=1.20):
    positions = np.asarray(positions, dtype=float)
    nat = len(symbols)
    adjacency = np.zeros((nat, nat), dtype=bool)
    if nat < 2:
        return adjacency
    radii = np.array([covalent_radii[atomic_numbers[sym]] for sym in symbols], dtype=float)
    diffs = positions[:, None, :] - positions[None, :, :]
    dists = np.linalg.norm(diffs, axis=-1)
    for i in range(nat):
        for j in range(i + 1, nat):
            cutoff = scale * (radii[i] + radii[j])
            if dists[i, j] <= cutoff:
                adjacency[i, j] = True
                adjacency[j, i] = True
    return adjacency


def graph_signature(symbols, positions, iterations=4):
    symbols = list(symbols)
    if len(symbols) == 0:
        return "empty"
    adjacency = bond_adjacency(symbols, positions)
    labels = [str(sym) for sym in symbols]
    for _ in range(iterations):
        new_labels = []
        for i, label in enumerate(labels):
            neighbors = sorted(labels[j] for j in np.where(adjacency[i])[0])
            new_labels.append(label + "|" + ",".join(neighbors))
        labels = new_labels
    return "||".join(sorted(labels))


def direct_centered_rmsd(a, b):
    if a.shape != b.shape:
        return 1.0e9
    return float(np.sqrt(np.mean((a - b) ** 2)))


def kabsch_rmsd(a, b):
    if a.shape != b.shape:
        return 1.0e9
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    a0 = a - a.mean(axis=0, keepdims=True)
    b0 = b - b.mean(axis=0, keepdims=True)
    cov = a0.T @ b0
    u, _, vt = np.linalg.svd(cov)
    rot = vt.T @ u.T
    if np.linalg.det(rot) < 0.0:
        vt[-1, :] *= -1.0
        rot = vt.T @ u.T
    a_rot = a0 @ rot.T
    return float(np.sqrt(np.mean((a_rot - b0) ** 2)))


def nearest_surface_height(metal_positions, mol_com):
    if len(metal_positions) == 0:
        return float("nan")
    zs = np.asarray(metal_positions[:, 2], dtype=float)
    layers = cluster_z_layers(zs, tol=0.50)
    bottom = float(np.mean(zs[layers[0]]))
    top = float(np.mean(zs[layers[-1]]))
    return float(min(abs(mol_com[2] - bottom), abs(mol_com[2] - top)))


def fractional_xy(position, cell):
    cell = np.asarray(cell, dtype=float)
    if abs(np.linalg.det(cell)) < 1.0e-12:
        return np.zeros(2, dtype=float)
    frac = cartesian_to_fractional(np.asarray(position, dtype=float), cell)
    return frac[:2]


def periodic_xy_distance(frac_xy_a, frac_xy_b, cell):
    duv = np.asarray(frac_xy_a, dtype=float) - np.asarray(frac_xy_b, dtype=float)
    duv = duv - np.round(duv)
    lateral = duv[0] * np.asarray(cell[0], dtype=float) + duv[1] * np.asarray(cell[1], dtype=float)
    return float(np.linalg.norm(lateral[:2]))


def frame_record(index, atoms, config_type, composition, metal):
    symbols = np.asarray(atoms.get_chemical_symbols())
    slab_indices, mol_indices = split_indices(symbols, metal)
    positions = np.asarray(atoms.get_positions(), dtype=float)
    mol_positions = positions[mol_indices]
    slab_positions = positions[slab_indices]
    if len(mol_positions) > 0:
        mol_com = np.mean(mol_positions, axis=0)
        mol_centered = mol_positions - mol_com[None, :]
        mol_pairwise = pairwise_descriptor(mol_positions)
        graph_code = graph_signature(symbols[mol_indices], mol_positions)
        frac_xy = fractional_xy(mol_com, atoms.cell.array)
    else:
        mol_com = np.zeros(3, dtype=float)
        mol_centered = np.zeros((0, 3), dtype=float)
        mol_pairwise = np.zeros(0, dtype=float)
        graph_code = "empty"
        frac_xy = np.zeros(2, dtype=float)
    ref_forces = np.asarray(atoms.arrays["REF_forces"], dtype=float)
    return FrameRecord(
        index=int(index),
        config_type=str(config_type),
        composition=str(composition),
        atoms=atoms.copy(),
        ref_energy=float(atoms.info["REF_energy"]),
        ref_forces=ref_forces,
        slab_indices=slab_indices,
        mol_indices=mol_indices,
        mol_com=mol_com,
        mol_centered=mol_centered,
        mol_pairwise=mol_pairwise,
        graph_signature=graph_code,
        height=nearest_surface_height(slab_positions, mol_com),
        frac_xy=frac_xy,
    )


def collect_matching_frames(dataset_path, composition, config_type, slab_config_type, pure_config_types, metal):
    interface_comp = f"{metal}75{composition}"
    interface_frames = []
    pure_frames = []
    slab_frames = []
    pure_config_types = set(pure_config_types)
    for index, atoms in enumerate(iread(str(dataset_path), index=":")):
        info = atoms.info
        if "REF_energy" not in info or "REF_forces" not in atoms.arrays:
            continue
        ct = info.get("config_type", "")
        comp = composition_string(atoms.get_chemical_symbols())
        if ct == config_type and comp == interface_comp:
            interface_frames.append(frame_record(index, atoms, ct, comp, metal))
        elif ct == slab_config_type and comp == f"{metal}75":
            slab_frames.append(frame_record(index, atoms, ct, comp, metal))
        elif ct in pure_config_types and comp == composition:
            pure_frames.append(frame_record(index, atoms, ct, comp, metal))
    return interface_frames, slab_frames, pure_frames


def choose_scan_subset(interface_frames, xy_cutoff, pose_cutoff, min_points, max_points):
    if not interface_frames:
        raise RuntimeError("No interface frames found for the requested composition")
    best_seed = None
    best_indices = None
    best_score = None
    for seed_idx, seed in enumerate(interface_frames):
        selected = []
        for idx, frame in enumerate(interface_frames):
            lateral = periodic_xy_distance(seed.frac_xy, frame.frac_xy, seed.atoms.cell.array)
            pose = direct_centered_rmsd(seed.mol_centered, frame.mol_centered)
            if lateral <= xy_cutoff and pose <= pose_cutoff:
                selected.append(idx)
        if len(selected) < min_points:
            continue
        heights = np.array([interface_frames[idx].height for idx in selected], dtype=float)
        score = (float(np.max(heights) - np.min(heights)), len(selected))
        if best_score is None or score > best_score:
            best_score = score
            best_seed = seed_idx
            best_indices = selected
    if best_indices is None:
        raise RuntimeError("Could not find a near-1D interface subset with the current thresholds")
    selected_frames = [interface_frames[idx] for idx in best_indices]
    selected_frames.sort(key=lambda frame: frame.height)
    if len(selected_frames) > max_points:
        sample_indices = np.linspace(0, len(selected_frames) - 1, max_points).round().astype(int)
        selected_frames = [selected_frames[idx] for idx in sample_indices]
    heights = np.array([frame.height for frame in selected_frames], dtype=float)
    return selected_frames, {
        "seed_frame_index": int(interface_frames[best_seed].index),
        "n_points": int(len(selected_frames)),
        "height_min": float(np.min(heights)),
        "height_max": float(np.max(heights)),
        "height_span": float(np.max(heights) - np.min(heights)),
        "xy_cutoff": float(xy_cutoff),
        "pose_cutoff": float(pose_cutoff),
    }


def choose_slab_reference(slab_frames, target_cell):
    if not slab_frames:
        raise RuntimeError("No slab-only Ag75 frames found in the dataset")
    target_cell = np.asarray(target_cell, dtype=float)
    best = None
    best_delta = None
    for frame in slab_frames:
        delta = float(np.linalg.norm(np.asarray(frame.atoms.cell.array, dtype=float) - target_cell))
        if best is None or delta < best_delta:
            best = frame
            best_delta = delta
    return best, {"cell_distance": float(best_delta)}


def match_pure_frames(scan_frames, pure_frames):
    if not pure_frames:
        raise RuntimeError("No pure-molecule frames found for the requested composition")
    pure_desc = np.stack([frame.mol_pairwise for frame in pure_frames], axis=0)
    matched = []
    for frame in scan_frames:
        delta = pure_desc - frame.mol_pairwise[None, :]
        dist2 = np.sum(delta * delta, axis=1)
        pure_idx = int(np.argmin(dist2))
        pure_frame = pure_frames[pure_idx]
        matched.append(
            {
                "frame": pure_frame,
                "pairwise_distance": float(np.sqrt(np.min(dist2))),
                "kabsch_rmsd": kabsch_rmsd(frame.mol_centered, pure_frame.mol_centered),
            }
        )
    return matched


def evaluate_atoms_list(calc, atoms_list):
    energies = []
    forces = []
    for atoms in atoms_list:
        test_atoms = atoms.copy()
        test_atoms.calc = calc
        energies.append(float(test_atoms.get_potential_energy()))
        forces.append(np.asarray(test_atoms.get_forces(), dtype=float))
    return np.asarray(energies, dtype=float), forces


def error_metrics(values):
    values = np.asarray(values, dtype=float)
    return {
        "rmse": float(np.sqrt(np.mean(values * values))),
        "mae": float(np.mean(np.abs(values))),
        "max_abs": float(np.max(np.abs(values))),
        "mean_signed": float(np.mean(values)),
    }


def plot_two_series(x, y_ref, y_model, xlabel, ylabel, title, out_path, legend_ref="REF (DFT)", legend_model="MAD-SURF"):
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(x, y_ref, "-o", ms=3.5, lw=1.6, label=legend_ref, color="#1f77b4")
    ax.plot(x, y_model, "-s", ms=3.2, lw=1.4, label=legend_model, color="#d62728", alpha=0.85)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.grid(True, alpha=0.25)
    ax.legend()
    fig.tight_layout()
    fig.savefig(out_path, dpi=180)
    plt.close(fig)


def plot_residuals(x, residual, xlabel, ylabel, title, out_path):
    fig, axes = plt.subplots(2, 1, figsize=(8, 7), sharex=True)
    axes[0].plot(x, residual, "-o", ms=3.5, lw=1.4, color="#9467bd")
    axes[0].axhline(0.0, color="black", lw=0.8, alpha=0.7)
    axes[0].set_ylabel(ylabel)
    axes[0].set_title(title)
    axes[0].grid(True, alpha=0.25)
    axes[1].semilogy(x, np.maximum(np.abs(residual), 1.0e-12), "-o", ms=3.5, lw=1.4, color="#2ca02c")
    axes[1].set_xlabel(xlabel)
    axes[1].set_ylabel(f"|{ylabel}|")
    axes[1].grid(True, alpha=0.25)
    fig.tight_layout()
    fig.savefig(out_path, dpi=180)
    plt.close(fig)


def plot_reference_components(x, interface_ref, slab_ref, mol_ref, interface_model, slab_model, mol_model, out_path):
    fig, ax = plt.subplots(figsize=(9, 5))
    ax.plot(x, interface_ref, "-o", ms=3.5, lw=1.5, label="Interface REF", color="#1f77b4")
    ax.plot(x, interface_model, "-s", ms=3.2, lw=1.4, label="Interface MAD-SURF", color="#d62728")
    ax.plot(x, slab_ref * np.ones_like(x), "--", lw=1.3, label="Slab REF", color="#17becf")
    ax.plot(x, slab_model * np.ones_like(x), "--", lw=1.3, label="Slab MAD-SURF", color="#ff9896")
    ax.plot(x, mol_ref, "-^", ms=3.0, lw=1.2, label="Matched molecule REF", color="#2ca02c")
    ax.plot(x, mol_model, "-v", ms=3.0, lw=1.2, label="Matched molecule MAD-SURF", color="#8c564b")
    ax.set_xlabel("Adsorbate COM height above nearest Ag surface (Å)")
    ax.set_ylabel("Absolute energy (eV)")
    ax.set_title("Three-set decomposition inputs on the same extracted Ag scan")
    ax.grid(True, alpha=0.25)
    ax.legend(ncol=2, fontsize=9)
    fig.tight_layout()
    fig.savefig(out_path, dpi=180)
    plt.close(fig)


def write_scan_table(path, rows):
    path = Path(path)
    ensure_dir(path.parent)
    with path.open("w", newline="", encoding="utf8") as fout:
        writer = csv.DictWriter(fout, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def signature_histogram(frames):
    counts = Counter(frame.graph_signature for frame in frames)
    return counts


def main():
    args = parse_args()
    out_dir = ensure_dir(args.out_dir)
    plots_dir = ensure_dir(out_dir / "plots")
    pure_config_types = [item.strip() for item in args.pure_config_types.split(",") if item.strip()]

    interface_frames, slab_frames, pure_frames = collect_matching_frames(
        dataset_path=args.dataset_path,
        composition=args.composition,
        config_type=args.config_type,
        slab_config_type=args.slab_config_type,
        pure_config_types=pure_config_types,
        metal=args.metal,
    )
    raw_interface_count = len(interface_frames)
    raw_pure_count = len(pure_frames)
    interface_signature_counts = signature_histogram(interface_frames)
    pure_signature_counts = signature_histogram(pure_frames)
    common_signatures = [sig for sig in interface_signature_counts if sig in pure_signature_counts]
    if not common_signatures:
        raise RuntimeError("No overlapping graph signature exists between the interface and pure-molecule reference sets")
    selected_signature = max(
        common_signatures,
        key=lambda sig: (min(interface_signature_counts[sig], pure_signature_counts[sig]), interface_signature_counts[sig]),
    )
    interface_frames = [frame for frame in interface_frames if frame.graph_signature == selected_signature]
    pure_frames = [frame for frame in pure_frames if frame.graph_signature == selected_signature]
    scan_frames, scan_meta = choose_scan_subset(
        interface_frames=interface_frames,
        xy_cutoff=args.xy_cutoff,
        pose_cutoff=args.pose_cutoff,
        min_points=args.min_points,
        max_points=args.max_points,
    )
    slab_frame, slab_meta = choose_slab_reference(slab_frames, scan_frames[0].atoms.cell.array)
    matched_pure = match_pure_frames(scan_frames, pure_frames)

    calc = MACECalculator(model_paths=args.model_path, device=args.device, default_dtype="float64")
    interface_atoms = [frame.atoms for frame in scan_frames]
    interface_model_energy, interface_model_forces = evaluate_atoms_list(calc, interface_atoms)
    slab_model_energy, slab_model_forces = evaluate_atoms_list(calc, [slab_frame.atoms])

    unique_pure_indices = sorted({entry["frame"].index for entry in matched_pure})
    pure_lookup = {frame.index: frame for frame in pure_frames if frame.index in unique_pure_indices}
    pure_atoms = [pure_lookup[idx].atoms for idx in unique_pure_indices]
    pure_model_energy, pure_model_forces = evaluate_atoms_list(calc, pure_atoms)
    pure_model_map = {idx: (pure_model_energy[i], pure_model_forces[i]) for i, idx in enumerate(unique_pure_indices)}

    heights = np.array([frame.height for frame in scan_frames], dtype=float)
    interface_ref_energy = np.array([frame.ref_energy for frame in scan_frames], dtype=float)
    interface_ref_fz = np.array([np.sum(frame.ref_forces[frame.mol_indices, 2]) for frame in scan_frames], dtype=float)
    interface_model_fz = np.array([np.sum(forces[frame.mol_indices, 2]) for forces, frame in zip(interface_model_forces, scan_frames)], dtype=float)

    slab_ref_energy = float(slab_frame.ref_energy)
    slab_model_energy = float(slab_model_energy[0])

    pure_ref_energy = np.array([entry["frame"].ref_energy for entry in matched_pure], dtype=float)
    pure_ref_fz = np.array([np.sum(entry["frame"].ref_forces[entry["frame"].mol_indices, 2]) for entry in matched_pure], dtype=float)
    pure_model_energy = np.array([pure_model_map[entry["frame"].index][0] for entry in matched_pure], dtype=float)
    pure_model_fz = np.array([np.sum(pure_model_map[entry["frame"].index][1][:, 2]) for entry in matched_pure], dtype=float)

    interaction_ref_energy = interface_ref_energy - slab_ref_energy - pure_ref_energy
    interaction_model_energy = interface_model_energy - slab_model_energy - pure_model_energy
    interaction_ref_fz = interface_ref_fz - pure_ref_fz
    interaction_model_fz = interface_model_fz - pure_model_fz
    interaction_energy_residual = interaction_model_energy - interaction_ref_energy
    interaction_energy_residual_centered = interaction_energy_residual - np.mean(interaction_energy_residual)

    interface_energy_rel_ref = interface_ref_energy - np.min(interface_ref_energy)
    interface_energy_rel_model = interface_model_energy - np.min(interface_ref_energy)

    plot_two_series(
        heights,
        interface_energy_rel_ref,
        interface_energy_rel_model,
        xlabel="Adsorbate COM height above nearest Ag surface (Å)",
        ylabel="Relative total interface energy (eV)",
        title=f"Real DFT vs MAD-SURF on extracted {args.composition}/Ag scan",
        out_path=plots_dir / "interface_energy_trace.png",
    )
    plot_residuals(
        heights,
        interface_model_energy - interface_ref_energy,
        xlabel="Adsorbate COM height above nearest Ag surface (Å)",
        ylabel="ΔE_model-ref (eV)",
        title="Total interface energy residual",
        out_path=plots_dir / "interface_energy_residual.png",
    )
    plot_two_series(
        heights,
        interface_ref_fz,
        interface_model_fz,
        xlabel="Adsorbate COM height above nearest Ag surface (Å)",
        ylabel="Total adsorbate Fz on interface (eV/Å)",
        title=f"Real DFT vs MAD-SURF adsorbate Fz on extracted {args.composition}/Ag scan",
        out_path=plots_dir / "interface_force_trace.png",
    )
    plot_residuals(
        heights,
        interface_model_fz - interface_ref_fz,
        xlabel="Adsorbate COM height above nearest Ag surface (Å)",
        ylabel="ΔFz_model-ref (eV/Å)",
        title="Total adsorbate Fz residual on the interface scan",
        out_path=plots_dir / "interface_force_residual.png",
    )
    plot_reference_components(
        heights,
        interface_ref_energy,
        slab_ref_energy,
        pure_ref_energy,
        interface_model_energy,
        slab_model_energy,
        pure_model_energy,
        plots_dir / "three_set_energy_inputs.png",
    )
    plot_two_series(
        heights,
        interaction_ref_energy,
        interaction_model_energy,
        xlabel="Adsorbate COM height above nearest Ag surface (Å)",
        ylabel="Approximate interaction energy (eV)",
        title="Approximate interaction energy from interface - slab - molecule",
        out_path=plots_dir / "interaction_energy_trace.png",
    )
    plot_residuals(
        heights,
        interaction_energy_residual,
        xlabel="Adsorbate COM height above nearest Ag surface (Å)",
        ylabel="ΔEint_model-ref (eV)",
        title="Approximate interaction-energy residual",
        out_path=plots_dir / "interaction_energy_residual.png",
    )
    plot_residuals(
        heights,
        interaction_energy_residual_centered,
        xlabel="Adsorbate COM height above nearest Ag surface (Å)",
        ylabel="ΔEint-centered (eV)",
        title="Approximate interaction-energy residual after removing the mean offset",
        out_path=plots_dir / "interaction_energy_residual_centered.png",
    )
    plot_two_series(
        heights,
        interaction_ref_fz,
        interaction_model_fz,
        xlabel="Adsorbate COM height above nearest Ag surface (Å)",
        ylabel="Approximate interaction Fz (eV/Å)",
        title="Approximate interaction Fz from interface - isolated molecule",
        out_path=plots_dir / "interaction_force_trace.png",
    )
    plot_residuals(
        heights,
        interaction_model_fz - interaction_ref_fz,
        xlabel="Adsorbate COM height above nearest Ag surface (Å)",
        ylabel="ΔFz_int_model-ref (eV/Å)",
        title="Approximate interaction-force residual",
        out_path=plots_dir / "interaction_force_residual.png",
    )
    plot_two_series(
        heights,
        pure_ref_fz,
        pure_model_fz,
        xlabel="Adsorbate COM height above nearest Ag surface (Å)",
        ylabel="Matched isolated-molecule total Fz (eV/Å)",
        title="Matched isolated-molecule net force check",
        out_path=plots_dir / "matched_molecule_force_trace.png",
        legend_ref="Matched molecule REF",
        legend_model="Matched molecule MAD-SURF",
    )
    plot_residuals(
        heights,
        np.array([entry["kabsch_rmsd"] for entry in matched_pure], dtype=float),
        xlabel="Adsorbate COM height above nearest Ag surface (Å)",
        ylabel="Kabsch RMSD to matched pure molecule (Å)",
        title="Geometry mismatch in the isolated-molecule matching step",
        out_path=plots_dir / "matched_molecule_geometry_mismatch.png",
    )

    scan_rows = []
    for i, frame in enumerate(scan_frames):
        pure_match = matched_pure[i]
        scan_rows.append(
            {
                "scan_order": i,
                "interface_frame_index": frame.index,
                "height_A": float(heights[i]),
                "interface_ref_energy_eV": float(interface_ref_energy[i]),
                "interface_model_energy_eV": float(interface_model_energy[i]),
                "interface_energy_residual_eV": float(interface_model_energy[i] - interface_ref_energy[i]),
                "interface_ref_total_fz_eVA": float(interface_ref_fz[i]),
                "interface_model_total_fz_eVA": float(interface_model_fz[i]),
                "interface_force_residual_eVA": float(interface_model_fz[i] - interface_ref_fz[i]),
                "slab_frame_index": slab_frame.index,
                "slab_ref_energy_eV": float(slab_ref_energy),
                "slab_model_energy_eV": float(slab_model_energy),
                "pure_frame_index": int(pure_match["frame"].index),
                "pure_config_type": pure_match["frame"].config_type,
                "pure_pairwise_distance": float(pure_match["pairwise_distance"]),
                "pure_kabsch_rmsd_A": float(pure_match["kabsch_rmsd"]),
                "pure_ref_energy_eV": float(pure_ref_energy[i]),
                "pure_model_energy_eV": float(pure_model_energy[i]),
                "pure_ref_total_fz_eVA": float(pure_ref_fz[i]),
                "pure_model_total_fz_eVA": float(pure_model_fz[i]),
                "interaction_ref_energy_eV": float(interaction_ref_energy[i]),
                "interaction_model_energy_eV": float(interaction_model_energy[i]),
                "interaction_energy_residual_eV": float(interaction_model_energy[i] - interaction_ref_energy[i]),
                "interaction_ref_total_fz_eVA": float(interaction_ref_fz[i]),
                "interaction_model_total_fz_eVA": float(interaction_model_fz[i]),
                "interaction_force_residual_eVA": float(interaction_model_fz[i] - interaction_ref_fz[i]),
            }
        )
    write_scan_table(out_dir / "scan_table.csv", scan_rows)

    summary = {
        "dataset_path": str(Path(args.dataset_path).resolve()),
        "model_path": str(Path(args.model_path).resolve()),
        "device": args.device,
        "composition": args.composition,
        "config_type": args.config_type,
        "slab_config_type": args.slab_config_type,
        "pure_config_types": pure_config_types,
        "training_source_note": {
            "full_dataset_file": str((Path(__file__).resolve().parent / "mad-surf_data" / "full_train_test_std_config_types.extxyz").resolve()),
            "small_train_file": str((Path(__file__).resolve().parent / "mad-surf_data" / "dataset" / "train_small_dataset_std.extxyz").resolve()),
            "small_test_file": str((Path(__file__).resolve().parent / "mad-surf_data" / "dataset" / "test_small_dataset_std.extxyz").resolve()),
            "full_training_script": str((Path(__file__).resolve().parent / "MAD-SURF-main" / "madsurf" / "mace_train" / "full_dataset" / "generate_training_configs-full.py").resolve()),
            "small_training_script": str((Path(__file__).resolve().parent / "MAD-SURF-main" / "madsurf" / "mace_train" / "small_dataset" / "generate_training_configs-small.py").resolve()),
            "dft_generation_script": str((Path(__file__).resolve().parent / "MAD-SURF-main" / "madsurf" / "boss" / "optimize" / "DFT_cal_aims.py").resolve()),
            "dft_backend": "FHI-aims with PBE + vdW_ts (tssurf), as indicated by the local generation scripts",
        },
        "counts": {
            "interface_frames_found_raw": int(raw_interface_count),
            "interface_frames_found": len(interface_frames),
            "slab_frames_found": len(slab_frames),
            "pure_frames_found_raw": int(raw_pure_count),
            "pure_frames_found": len(pure_frames),
        },
        "graph_signatures": {
            "selected_signature_interface_count": int(interface_signature_counts[selected_signature]),
            "selected_signature_pure_count": int(pure_signature_counts[selected_signature]),
            "n_common_signatures": int(len(common_signatures)),
            "n_interface_signatures": int(len(interface_signature_counts)),
            "n_pure_signatures": int(len(pure_signature_counts)),
            "pure_frames_after_signature_filter": int(len(pure_frames)),
        },
        "scan_selection": scan_meta,
        "slab_reference": {
            "frame_index": int(slab_frame.index),
            **slab_meta,
        },
        "pure_matching": {
            "mean_pairwise_distance": float(np.mean([entry["pairwise_distance"] for entry in matched_pure])),
            "max_pairwise_distance": float(np.max([entry["pairwise_distance"] for entry in matched_pure])),
            "mean_kabsch_rmsd_A": float(np.mean([entry["kabsch_rmsd"] for entry in matched_pure])),
            "max_kabsch_rmsd_A": float(np.max([entry["kabsch_rmsd"] for entry in matched_pure])),
        },
        "metrics": {
            "slab_energy_residual_eV": float(slab_model_energy - slab_ref_energy),
            "pure_energy": error_metrics(pure_model_energy - pure_ref_energy),
            "interface_total_energy": error_metrics(interface_model_energy - interface_ref_energy),
            "interface_total_fz": error_metrics(interface_model_fz - interface_ref_fz),
            "interaction_energy": error_metrics(interaction_energy_residual),
            "interaction_energy_centered": error_metrics(interaction_energy_residual_centered),
            "interaction_fz": error_metrics(interaction_model_fz - interaction_ref_fz),
            "matched_pure_total_fz": error_metrics(pure_model_fz - pure_ref_fz),
        },
    }
    save_json(summary, out_dir / "summary.json")
    print(f"analysis -> {out_dir}")


if __name__ == "__main__":
    main()
