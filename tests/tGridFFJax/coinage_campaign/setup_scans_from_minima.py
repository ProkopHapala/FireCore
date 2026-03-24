#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[3]
sys.path.append(str(Path(__file__).resolve().parent.parent))

from coinage_campaign.config import build_default_campaign_config
from coinage_campaign.workflow import create_scan_jobs_from_minimum


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Create rigid, molecule-relaxed, and/or slab+molecule-relaxed scan jobs from converged adsorption minima.")
    parser.add_argument("--minimum-dir", help="Path to one converged minimum directory containing CONTCAR and job_manifest.json.")
    parser.add_argument("--minima-root", help="Root containing many converged minimum directories (final_static stage directories).")
    parser.add_argument("--out-root", required=True, help="Output root for scan jobs.")
    parser.add_argument(
        "--family",
        choices=["rigid", "relaxed", "relaxed_slab", "both", "all"],
        default="all",
        help="Scan family to generate. 'both' keeps the legacy two-family behavior; 'all' generates all three families.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if bool(args.minimum_dir) == bool(args.minima_root):
        raise SystemExit("Provide exactly one of --minimum-dir or --minima-root.")
    config = build_default_campaign_config(ROOT)
    if args.family == "both":
        families = ("rigid", "relaxed")
    elif args.family == "all":
        families = ("rigid", "relaxed", "relaxed_slab")
    else:
        families = (args.family,)
    minima_dirs: list[Path]
    if args.minimum_dir:
        minima_dirs = [Path(args.minimum_dir)]
    else:
        minima_root = Path(args.minima_root)
        minima_dirs = sorted(path for path in minima_root.rglob("final_static") if (path / "job_manifest.json").exists())
    all_summaries = []
    for minimum_dir in minima_dirs:
        for family in families:
            all_summaries.append(create_scan_jobs_from_minimum(minimum_dir, Path(args.out_root), config, family))
    phase_jobs = []
    for summary in all_summaries:
        phase_jobs.extend(summary["jobs"])
    payload = {
        "minimum_count": len(minima_dirs),
        "family_mode": args.family,
        "jobs": phase_jobs,
        "summaries": all_summaries,
    }
    (Path(args.out_root) / "scan_suite_manifest.json").write_text(json.dumps(payload, indent=2) + "\n")
    (Path(args.out_root) / "phase_manifest.json").write_text(json.dumps({"jobs": phase_jobs}, indent=2) + "\n")
    print(f"scan family : {args.family}")
    print(f"minima count : {len(minima_dirs)}")
    print(f"scan groups  : {len(all_summaries)}")
    print(f"scan root    : {args.out_root}")


if __name__ == "__main__":
    main()
