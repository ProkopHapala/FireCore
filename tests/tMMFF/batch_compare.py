import argparse
import subprocess
from pathlib import Path

# --- Configuration ---
# The comparison script is located in the doc/py directory, relative to the project root.
COMPARISON_SCRIPT = Path(__file__).parent.parent.parent / "doc/py/energy_comparison_1d.py"
FIRECORE_BASE = Path(__file__).parent
LAMMPS_BASE = Path("/home/indranil/Documents/Project_1/Lammps")

# --- Mappings ---
# Maps FireCore directory names to their LAMMPS counterparts.

DIRECTION_MAP = {
    "dir_1.0_1.0_0.0": "angle45",
    "dir_2.0_1.0_0.0": "angle26.565051177078",
}

CONSTRAINT_MAP = {
    "cons_24": "fixedatom25",
    "cons_26": "fixedatom27",
}

GRID_MAP = {
    "20x20": "nx20",
}

def find_and_run_comparisons(dry_run=False):
    """
    Finds FireCore data files and runs the comparison script against
    the corresponding LAMMPS data files.
    """
    print("Starting batch comparison...")
    firecore_files = list(FIRECORE_BASE.glob("relax_*/dir_*/cons_*/PTCDA_20x20_*_total.dat"))
    print(f"Found {len(firecore_files)} FireCore data files to process.")

    for fc_file in firecore_files:
        try:
            parts = fc_file.relative_to(FIRECORE_BASE).parts
            scan_type_fc, direction_fc, constraint_fc = parts[0], parts[1], parts[2]
            grid_fc = "20x20"  # Based on the glob pattern

            # Get LAMMPS components from mappings
            lammps_angle_dir = DIRECTION_MAP[direction_fc]
            lammps_fixed_atom_dir = CONSTRAINT_MAP[constraint_fc]
            lammps_grid_dir = GRID_MAP[grid_fc]

            # Determine base path and construct final path based on scan type
            if "perfect" in scan_type_fc:
                lammps_base_scan_path = "4-relaxed_linescan"
                lammps_path = LAMMPS_BASE / lammps_base_scan_path / lammps_angle_dir / lammps_fixed_atom_dir / lammps_grid_dir / "total.dat"
            elif "defect" in scan_type_fc:
                lammps_base_scan_path = "5-relaxed_linescan_defect"
                if "aligned" in scan_type_fc:
                    defect_type_dir = "defect_aligned"
                elif "perpendicular" in scan_type_fc:
                    defect_type_dir = "defect_perpendicular"
                else:
                    print(f"\n[SKIP] Unknown defect type in {scan_type_fc}")
                    continue
                lammps_path = LAMMPS_BASE / lammps_base_scan_path / lammps_angle_dir / defect_type_dir / lammps_fixed_atom_dir / lammps_grid_dir / "total.dat"
            else:
                print(f"\n[SKIP] Unknown scan type {scan_type_fc}")
                continue

            if not lammps_path.exists():
                print(f"\n[SKIP] LAMMPS file not found for {fc_file.relative_to(FIRECORE_BASE)}")
                print(f"       Searched for: {lammps_path}")
                continue

            output_png = fc_file.with_name(f"comparison_{fc_file.stem}.png")
            title = f"{scan_type_fc.replace('_', ' ')}\n{direction_fc} | {constraint_fc}"

            cmd = [
                "python", str(COMPARISON_SCRIPT.resolve()),
                "--lammps", str(lammps_path),
                "--firecore", str(fc_file),
                "--title", title,
                "--out", str(output_png)
            ]

            if dry_run:
                print("\n[DRY-RUN] " + " ".join(cmd))
            else:
                print(f"\nProcessing: {fc_file.relative_to(FIRECORE_BASE)}")
                print(f"-> Plotting against: {lammps_path}")
                print(f"-> Output: {output_png.name}")
                result = subprocess.run(cmd, check=True, capture_output=True, text=True)
                if result.stdout:
                    print(result.stdout.strip())

        except KeyError as e:
            print(f"\n[SKIP] Mapping not found for a component in path: {fc_file.relative_to(FIRECORE_BASE)}")
            print(f"       (Missing key: {e})")
        except subprocess.CalledProcessError as e:
            print(f"[ERROR] Subprocess failed for {fc_file.name}")
            print("STDOUT:", e.stdout)
            print("STDERR:", e.stderr)
        except Exception as e:
            print(f"[ERROR] An unexpected error occurred while processing {fc_file.name}: {e}")

def main():
    parser = argparse.ArgumentParser(description="Batch run energy comparison between FireCore and LAMMPS.")
    parser.add_argument("--dry-run", action="store_true", help="Print commands without executing them.")
    args = parser.parse_args()
    find_and_run_comparisons(dry_run=args.dry_run)

if __name__ == "__main__":
    main()
