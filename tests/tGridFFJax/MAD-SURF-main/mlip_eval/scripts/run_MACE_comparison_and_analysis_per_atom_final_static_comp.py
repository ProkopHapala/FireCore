import os
import sys
import glob
import numpy as np
from collections import Counter
import pandas as pd
from ase.io import read, write, Trajectory
from ase.optimize import BFGS
from ase.constraints import FixAtoms
from mace.calculators import MACECalculator

# -------------------------
# User settings
# -------------------------
input_file = sys.argv[1]        # input structures
reference_file = sys.argv[2]    # reference structures
output_dir = "relaxed_MACE_results_final_comparison_static"
os.makedirs(output_dir, exist_ok=True)

results_file = os.path.join(output_dir, "comparison_results_static.csv")

# Define your models
models = {
    "MACE_MPA-0_foundational_model": "/scratch/project_2008059/GEN_MLIP_EVAL/models/descriptor_filtering/full_set/mace-mpa-0-medium.model",
    "MACE_MPA-0_foundational_on_test": "/scratch/project_2008059/GEN_MLIP_EVAL/models/finetuned/foundational_model/MACE_finetuned.model",
    "MACE_MPA-0_foundational_on_train": "/scratch/project_2008059/GEN_MLIP_EVAL/models/finetuned/foundational_model/on_training_subset/MACE_finetuned.model",
    #"small_only_energy_stageone": "/scratch/project_2008059/GEN_MLIP_EVAL/models/small_set_lambda/only_energy/MACE_model.model",
    "small_lambda1_stageone": "/scratch/project_2008059/GEN_MLIP_EVAL/models/small_set_lambda/lambda_1/MACE_model.model",
    "small_lambda10_stageone": "/scratch/project_2008059/GEN_MLIP_EVAL/models/small_set_lambda/lambda_10/MACE_model.model",
    #"small_lambda100_stageone": "/scratch/project_2008059/GEN_MLIP_EVAL/models/small_set_lambda/lambda_100/MACE_model.model",
    #"small_lambda1000_stageone": "/scratch/project_2008059/GEN_MLIP_EVAL/models/small_set_lambda/lambda_1000/MACE_model.model",
    #"small_only_forces_stageone": "/scratch/project_2008059/GEN_MLIP_EVAL/models/small_set_lambda/only_forces/MACE_model.model",
    #"full_only_energy_stageone": "/scratch/project_2008059/GEN_MLIP_EVAL/models/full_set_lambda/only_energy/MACE_model.model",
    "full_lambda1_stageone": "/scratch/project_2008059/GEN_MLIP_EVAL/models/full_set_lambda/lambda_1/MACE_model.model",
    #"full_lambda10_stageone": "/scratch/project_2008059/GEN_MLIP_EVAL/models/full_set_lambda/lambda_10/MACE_model.model",
    "full_lambda10_stageone_scaleshift": "/scratch/project_2008059/GEN_MLIP_EVAL/models/full_set_lambda/lambda_10_ScaleShift/MACE_model.model",
    #"full_lambda10_stagetwo_scaleshift": "/scratch/project_2008059/GEN_MLIP_EVAL/models/full_set_lambda/lambda_10_ScaleShift/MACE_model_stagetwo.model",
    #"full_lambda100_stageone": "/scratch/project_2008059/GEN_MLIP_EVAL/models/full_set_lambda/lambda_100/MACE_model.model",
    #"full_lambda1000_stageone": "/scratch/project_2008059/GEN_MLIP_EVAL/models/full_set_lambda/lambda_1000/MACE_model.model",
    #"full_only_forces_stageone": "/scratch/project_2008059/GEN_MLIP_EVAL/models/full_set_lambda/only_forces/MACE_model.model",
    "df_0.01_l10": "/scratch/project_2008059/GEN_MLIP_EVAL/models/descriptor_filtering/full_set/threshold_0.01/MACE_model.model",
    "df_0.08_l10": "/scratch/project_2008059/GEN_MLIP_EVAL/models/descriptor_filtering/full_set/threshold_0.08/MACE_model.model",    
    #"df_0.08_l10_stageone_on_test_allrel": "/scratch/project_2008059/GEN_MLIP_EVAL/models/finetuned/descfilter_0.08_stageone/all_rel_data_test_mols/MACE_finetuned.model",
    #"df_0.08_l10_stagetwo_on_test_allopt": "/scratch/project_2008059/GEN_MLIP_EVAL/models/finetuned/descfilter_0.08/all_opt_test_mols/MACE_finetuned.model",
    #"df_0.08_l10_stagetwo_on_test_allrel": "/scratch/project_2008059/GEN_MLIP_EVAL/models/finetuned/descfilter_0.08/all_rel_data_test_mols/MACE_finetuned.model",    
    #"df_0.08_l10_stagetwo_on_test_filtered": "/scratch/project_2008059/GEN_MLIP_EVAL/models/finetuned/descfilter_0.08/filtered_rel_data_test_mols/MACE_finetuned.model",
}

relax_steps = 2000
FixIndex = 72  # Fix bottom N atoms of the substrate (set None if no constraint)

# -------------------------
# Helper functions
# -------------------------
def rmsd_per_element(a, b):
    """
    RMSD where each element type contributes equally,
    independent of how many atoms of that element appear.
    """
    # squared distances
    diff2 = (a.get_positions() - b.get_positions())**2
    dist2 = diff2.sum(axis=1)

    # count atoms per element
    elems = a.get_chemical_symbols()
    counts = Counter(elems)

    # weight per atom = 1 / (# atoms of that atom’s element)
    weights = np.array([1.0 / counts[e] for e in elems], dtype=float)

    # normalize so total weight = number of elements (optional but consistent)
    weights /= weights.sum()

    return np.sqrt((weights * dist2).sum())

def rmsd(a, b):
    """Simple RMSD between two ASE Atoms objects (no alignment)."""
    return np.sqrt(((a.get_positions() - b.get_positions())**2).sum(axis=1).mean())

def force_mae(forces_pred, forces_ref):
    """Mean absolute error per atom between predicted and reference forces."""
    return np.mean(np.linalg.norm(forces_pred - forces_ref, axis=1))

def adsorption_height(atoms):
    """
    Compute adsorption height = avg(z of carbons) - avg(z of topmost metal layer).
    Metals considered: Au, Cu, Ag.
    """
    pos = atoms.get_positions()
    syms = atoms.get_chemical_symbols()

    carbon_z = [pos[i, 2] for i, s in enumerate(syms) if s == "C"]
    metal_indices = [i for i, s in enumerate(syms) if s in ["Au", "Cu", "Ag"]]

    if len(carbon_z) == 0 or len(metal_indices) == 0:
        return None

    metal_z = pos[metal_indices, 2]

    # Define "top layer" as metals within 0.5 Å of the max z
    z_max = metal_z.max()
    top_layer = metal_z[metal_z > z_max - 0.5]

    if len(top_layer) == 0:
        return None

    return np.mean(carbon_z) - np.mean(top_layer)

def recompute_metrics_from_extxyz(extxyz_path, model_name, refs):
    frames = read(extxyz_path, ":")
    assert len(frames) == len(refs), f"Frame mismatch in {model_name}"
    metrics = []
    for i, atoms_copy in enumerate(frames):
        ref = refs[i]
        energy_diff = atoms_copy.get_potential_energy() - ref.get_potential_energy()
        energy_diff_per_atom = energy_diff/len(atoms_copy)
        pos_rmsd = rmsd(atoms_copy, ref)
        f_mae = force_mae(atoms_copy.get_forces(), ref.get_forces())
        ads_height_relaxed = adsorption_height(atoms_copy)
        ads_height_ref = adsorption_height(ref)
        metrics.append({
            "model": model_name,
            "index": i,
            "energy_diff": energy_diff,
            "energy_diff_per_atom": energy_diff_per_atom,
            "rmsd": pos_rmsd,
            "force_mae": f_mae,
            "force_mae_static" : f_mae_static, 
            "ads_height_relaxed": ads_height_relaxed,
            "ads_height_reference": ads_height_ref,
            "ads_height_diff": None if (ads_height_relaxed is None or ads_height_ref is None)
                               else ads_height_relaxed - ads_height_ref,
        })
    return metrics

# -------------------------
# Load data
# -------------------------
inputs = read(input_file, ":")
refs = read(reference_file, ":")
assert len(inputs) == len(refs), "Input and reference extxyz must have the same number of frames."

# -------------------------
# Recover previous progress
# -------------------------
results = []
completed_models = set()

if os.path.exists(results_file):
    print(f"Recovering results from {results_file}")
    results = pd.read_csv(results_file).to_dict(orient="records")
    completed_models = {r["model"] for r in results}
else:
    for f in glob.glob(os.path.join(output_dir, "relaxed_*.extxyz")):
        filename = os.path.basename(f)  # FIX: only take filename
        model_name = filename.replace("relaxed_", "").replace(".extxyz", "")
        print(f"Recomputing metrics for {model_name} from existing extxyz")
        metrics = recompute_metrics_from_extxyz(f, model_name, refs)
        results.extend(metrics)
        completed_models.add(model_name)

print(f"Already completed models: {completed_models}")

# -------------------------
# Run relaxation for unfinished models
# -------------------------
for model_name, model_path in models.items():
    if model_name in completed_models:
        print(f"Skipping {model_name} (already done)")
        continue

    print(f" Running relaxation with model: {model_name}")
    
    relaxed_structures = []
    for i, atoms in enumerate(inputs):
        atoms_copy = atoms.copy()

        # Attach calculator
        macemp = MACECalculator(
            model_path=model_path,
            device="cuda",
            enable_cueq=True,
            default_dtype="float64"
        )
        atoms_copy.calc = macemp

        # Apply constraint (fix lower part of substrate)
        if FixIndex is not None:
            fixed_indices = list(range(FixIndex))
            fix_constraint = FixAtoms(indices=fixed_indices)
            atoms_copy.set_constraint(fix_constraint)

        # Optimizer
        save_file = os.path.join(output_dir, f"{model_name}_idx{i}")
        #traj = Trajectory(f"{save_file}_rel_BFGS.traj", 'w', atoms_copy)
        pred_static_forces = atoms_copy.get_forces()
        pred_static_energy = atoms_copy.get_potential_energy()
        qn = BFGS(atoms_copy, logfile=None)
        #qn.attach(traj)
        qn.run(fmax=0.01, steps=relax_steps)

        relaxed_structures.append(atoms_copy)

        # Compare with reference
        ref = refs[i]
        energy_diff = atoms_copy.get_potential_energy() - ref.get_potential_energy()
        energy_diff_per_atom = energy_diff/len(atoms_copy)
        energy_diff_static = pred_static_energy - ref.get_potential_energy()
        energy_diff_per_atom_static = energy_diff_static/len(atoms_copy)
        pos_rmsd = rmsd(atoms_copy, ref)
        pos_rmsd_per_elem = rmsd_per_element(atoms_copy, ref)

        # Force comparison
        ref_forces = ref.get_forces()
        pred_forces = atoms_copy.get_forces()
        f_mae = force_mae(pred_forces, ref_forces)
        f_mae_static = force_mae(pred_static_forces, ref_forces)

        # Adsorption height
        ads_height_relaxed = adsorption_height(atoms_copy)
        ads_height_ref = adsorption_height(ref)

        results.append({
            "model": model_name,
            "index": i,
            "energy_diff": energy_diff,
            "energy_diff_per_atom": energy_diff_per_atom,
            "energy_diff_static": energy_diff_static,
            "energy_diff_per_atom_static": energy_diff_per_atom_static,            
            "rmsd": pos_rmsd,
            "rmsd_per_element": pos_rmsd_per_elem,
            "force_mae": f_mae,
            "force_mae_static" : f_mae_static, 
            "ads_height_relaxed": ads_height_relaxed,
            "ads_height_reference": ads_height_ref,
            "ads_height_diff": None if (ads_height_relaxed is None or ads_height_ref is None)
                               else ads_height_relaxed - ads_height_ref,
        })

    # Save relaxed structures
    out_file = os.path.join(output_dir, f"relaxed_{model_name}.extxyz")
    write(out_file, relaxed_structures)

# -------------------------
# Save comparison results
# -------------------------
df = pd.DataFrame(results)
df.to_csv(results_file, index=False)
print("Done. Relaxations + comparisons saved.")

