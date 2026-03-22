#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[3]
sys.path.append(str(Path(__file__).resolve().parent.parent))

from coinage_campaign.config import build_default_campaign_config
from coinage_campaign.workflow import create_interaction_reference_jobs


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Create slab-only and molecule-only reference jobs from scan geometries.")
    parser.add_argument("--scan-root", required=True, help="Root containing scan point job directories.")
    parser.add_argument("--out-root", required=True, help="Output root for reference jobs.")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    config = build_default_campaign_config(ROOT)
    summary = create_interaction_reference_jobs(Path(args.scan_root), Path(args.out_root), config)
    print(f"reference jobs: {len(summary['jobs'])}")
    print(f"reference root: {args.out_root}")


if __name__ == "__main__":
    main()
