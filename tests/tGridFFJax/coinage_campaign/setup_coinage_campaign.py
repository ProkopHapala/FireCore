#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[3]
sys.path.append(str(Path(__file__).resolve().parent.parent))

from coinage_campaign.config import build_default_campaign_config
from coinage_campaign.workflow import (
    create_clean_slab_final_phase,
    create_local_pilot_campaign,
    create_production_adsorption_phases,
    save_json,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Set up the local-first coinage-metal DFT campaign tree.")
    parser.add_argument(
        "--campaign-root",
        default="/home/niel/git/coinage_gridff_dft",
        help="External campaign root. Use /tmp/... for sandbox validation.",
    )
    parser.add_argument(
        "--pilot-only",
        action="store_true",
        help="Reserved switch for future use. Current implementation sets up the local pilot tree.",
    )
    parser.add_argument(
        "--selected-slab",
        choices=("3x3x4",),
        help="Optional promotion target. If provided, also materialize 02_clean_slab_final for the selected slab.",
    )
    return parser.parse_args()


def build_root_readme(campaign_root: str, selected_slab: str | None) -> str:
    selected_text = selected_slab if selected_slab is not None else "not promoted yet"
    return "\n".join(
        [
            "# Coinage GridFF DFT Campaign",
            "",
            f"- Campaign root: `{campaign_root}`",
            "- Primary protocol: `PBE + vdWsurf` for slabs/adsorbates",
            "- Bulk protocol: no `vdWsurf` tags",
            "- Fixed slab: `3x3x4`",
            "- HPC default: `1` node, `96` CPUs, `NCORE = 16`",
            f"- Selected slab status: `{selected_text}`",
            "",
            "## What Is Generated Initially",
            "- `00_bulk`: Ag/Cu/Au bulk multi-step optimization",
            "- `01_universal_slab_selection`: Ag `3x3x4` clean slab and representative adsorption pilots",
            "- `03_gas_refs`: isolated gas-phase molecule references",
            "",
            "## Stage Pattern",
            "- Bulk: `relax_stage1 -> relax_stage2 -> final_static`",
            "- Slab/adsorbates: `relax_stage1_nodipole -> relax_stage2_cycle1_dipole -> relax_stage2_cycle2_dipole -> relax_stage2_cycle3_dipole -> final_static -> workfunction -> bader`",
            "",
            "## Representative Adsorbates",
            "- `H`, `CO`, `H2O`, `HCONH2`, `pyridine`",
            "",
            "## Production Adsorbates",
            "- `H`, `CO`, `H2O`, `NH3`, `methanol`, `HCONH2`, `pyridine`",
            "",
            "## How To Run On HPC",
            "- Set `VASP_BIN_ON_HPC` if `vasp_std` is not already on `PATH`.",
            "- Submit each job with `submit_hpc.sh` or `qsub run.pbs` inside the job directory.",
            "- Recommended first wave: `00_bulk/*`, `03_gas_refs/*`, then `01_universal_slab_selection/Ag/3x3x4/*`.",
            "",
            "## After Bulk Finishes",
            "- Promote the fixed production slab with:",
            "```bash",
            f"python3 tests/tGridFFJax/coinage_campaign/setup_coinage_campaign.py --campaign-root {campaign_root} --selected-slab 3x3x4",
            "```",
            "- That materializes `02_clean_slab_final`, `04_ads_seed_library`, and `05_ads_relax` from real bulk outputs.",
            "",
            "## Scan Generation",
            "- Rigid and relaxed scans are generated later from converged minima using `setup_scans_from_minima.py`.",
            "- Interaction reference jobs are generated later from converged scans using `setup_interaction_references.py`.",
            "",
        ]
    )


def main() -> None:
    args = parse_args()
    config = build_default_campaign_config(ROOT, Path(args.campaign_root))
    manifest = create_local_pilot_campaign(config, Path(args.campaign_root))
    promoted = None
    production = None
    if args.selected_slab:
        promoted = create_clean_slab_final_phase(config, Path(args.campaign_root), args.selected_slab)
        production = create_production_adsorption_phases(config, Path(args.campaign_root), args.selected_slab)
    save_json(
        {
            "campaign_root": args.campaign_root,
            "job_count": len(manifest["jobs"]),
            "selected_slab": args.selected_slab,
            "promoted_job_count": 0 if promoted is None else len(promoted["jobs"]),
            "production_seed_count": 0 if production is None else len(production["seed_entries"]),
            "production_relax_job_count": 0 if production is None else len(production["relax_jobs"]),
            "notes": manifest["notes"],
        },
        Path(args.campaign_root) / "setup_summary.json",
    )
    (Path(args.campaign_root) / "README.md").write_text(
        build_root_readme(args.campaign_root, args.selected_slab)
    )
    print(f"campaign root: {args.campaign_root}")
    print(f"jobs planned : {len(manifest['jobs'])}")


if __name__ == "__main__":
    main()
