#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import os
import shutil
import subprocess
import sys
import time
from dataclasses import dataclass
from pathlib import Path


SUCCESS_MARKERS = (
    "General timing and accounting informations for this job:",
    "reached required accuracy - stopping structural energy minimisation",
)

BUSY_PATTERNS = (
    "run_pipeline_local.sh",
    "run_local.sh",
    "/vasp_std",
    " vasp_std",
)

PHASE_PRIORITY = {
    "00_bulk": 0,
    "03_gas_refs": 1,
    "01_universal_slab_selection": 2,
    "02_clean_slab_final": 3,
}


@dataclass(frozen=True)
class QueueItem:
    kind: str
    path: Path
    phase: str
    description: str


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run coinage campaign jobs one-by-one on the local workstation.")
    parser.add_argument("--campaign-root", required=True, help="Generated campaign root.")
    parser.add_argument("--pilot-only", action="store_true", help="Restrict execution to the local-pilot phases.")
    parser.add_argument(
        "--phases",
        default="00_bulk,01_universal_slab_selection,02_clean_slab_final,03_gas_refs,04_ads_seed_library,05_ads_relax,06_scan_rigid,07_scan_relaxed,08_volumetrics,09_analysis",
        help="Comma-separated phases to consider.",
    )
    parser.add_argument("--phase", action="append", default=[], help="Additional phase filter. Can be passed multiple times.")
    parser.add_argument("--max-items", type=int, default=0, help="Optional limit on launched queue items.")
    parser.add_argument("--status-only", action="store_true", help="Only print the queued work; do not launch anything.")
    parser.add_argument("--dry-run", action="store_true", help="Alias for --status-only.")
    parser.add_argument("--wait-if-busy", action="store_true", help="Wait until no other VASP/pipeline job is active before launching.")
    parser.add_argument("--poll-seconds", type=int, default=30, help="Polling interval for --wait-if-busy.")
    return parser.parse_args()


def load_jobs(campaign_root: Path) -> list[dict]:
    jobs: list[dict] = []
    campaign_manifest = campaign_root / "campaign_manifest.json"
    if campaign_manifest.exists():
        jobs.extend(json.loads(campaign_manifest.read_text())["jobs"])
    for phase_manifest in sorted(campaign_root.glob("*/phase_manifest.json")):
        payload = json.loads(phase_manifest.read_text())
        jobs.extend(payload.get("jobs", []))
    return jobs


def phase_filter(raw: str) -> set[str]:
    return {part.strip() for part in raw.split(",") if part.strip()}


def outcar_complete(job_dir: Path) -> bool:
    outcar = job_dir / "OUTCAR"
    if not outcar.exists():
        return False
    text = outcar.read_text(errors="ignore")
    return any(marker in text for marker in SUCCESS_MARKERS)


def parse_pipeline_stage_order(pipeline_root: Path) -> list[str]:
    script = pipeline_root / "run_pipeline_local.sh"
    lines = script.read_text().splitlines()
    stages: list[str] = []
    for line in lines:
        stripped = line.strip()
        if stripped.startswith('echo "Running ') and stripped.endswith('"'):
            stages.append(stripped[len('echo "Running '):-1])
    return stages


def pipeline_stage_dirs(pipeline_root: Path) -> list[Path]:
    return [pipeline_root / stage for stage in parse_pipeline_stage_order(pipeline_root)]


def pipeline_complete(pipeline_root: Path) -> bool:
    stage_dirs = pipeline_stage_dirs(pipeline_root)
    if not stage_dirs:
        return False
    return outcar_complete(stage_dirs[-1])


def pipeline_started(pipeline_root: Path) -> bool:
    return any((stage_dir / "OUTCAR").exists() for stage_dir in pipeline_stage_dirs(pipeline_root))


def copy_restart_files(source_dir: Path, target_dir: Path) -> None:
    copies = (
        ("CONTCAR", "POSCAR"),
        ("WAVECAR", "WAVECAR"),
        ("CHGCAR", "CHGCAR"),
        ("AECCAR0", "AECCAR0"),
        ("AECCAR2", "AECCAR2"),
    )
    for source_name, target_name in copies:
        source_path = source_dir / source_name
        target_path = target_dir / target_name
        if source_path.exists():
            shutil.copy2(source_path, target_path)


def ps_lines() -> list[str]:
    proc = subprocess.run(
        ["ps", "-eo", "pid=,args="],
        check=True,
        text=True,
        capture_output=True,
    )
    return [line.rstrip() for line in proc.stdout.splitlines() if line.strip()]


def system_busy() -> bool:
    this_pid = os.getpid()
    for line in ps_lines():
        pid_text, _, args = line.strip().partition(" ")
        if not pid_text.isdigit():
            continue
        if int(pid_text) == this_pid:
            continue
        if any(pattern in args for pattern in BUSY_PATTERNS):
            return True
    return False


def wait_until_idle(poll_seconds: int) -> None:
    while system_busy():
        print(f"[queue] another VASP/pipeline job is active; sleeping {poll_seconds}s", flush=True)
        time.sleep(poll_seconds)


def relative_sort_key(campaign_root: Path, path: Path) -> tuple:
    rel = path.relative_to(campaign_root)
    phase = rel.parts[0]
    priority = PHASE_PRIORITY.get(phase, 99)
    clean_slab_priority = 0 if "clean_slab" in rel.parts else 1
    representative_priority = 1 if "representative_adsorbates" in rel.parts else 0
    return (priority, clean_slab_priority, representative_priority, rel.as_posix())


def enqueue_items(campaign_root: Path, phases: set[str]) -> list[QueueItem]:
    jobs = load_jobs(campaign_root)
    items: list[QueueItem] = []
    pipeline_roots: set[Path] = set()
    stage_paths: set[Path] = set()
    for job in jobs:
        if job["phase"] not in phases:
            continue
        stage_path = campaign_root / job["path"]
        stage_paths.add(stage_path)
        pipeline_root = stage_path.parent
        if (pipeline_root / "run_pipeline_local.sh").exists():
            pipeline_roots.add(pipeline_root)

    queued_pipelines: set[Path] = set()
    for pipeline_root in sorted(pipeline_roots, key=lambda path: relative_sort_key(campaign_root, path)):
        if pipeline_complete(pipeline_root):
            continue
        if not pipeline_started(pipeline_root):
            phase = pipeline_root.relative_to(campaign_root).parts[0]
            items.append(
                QueueItem(
                    kind="pipeline",
                    path=pipeline_root,
                    phase=phase,
                    description=str(pipeline_root.relative_to(campaign_root)),
                )
            )
            queued_pipelines.add(pipeline_root)

    grouped_stage_paths: dict[Path, list[Path]] = {}
    for stage_path in stage_paths:
        grouped_stage_paths.setdefault(stage_path.parent, []).append(stage_path)

    for parent_root in sorted(grouped_stage_paths, key=lambda path: relative_sort_key(campaign_root, path)):
        if parent_root in queued_pipelines:
            continue
        stage_group = grouped_stage_paths[parent_root]
        if (parent_root / "run_pipeline_local.sh").exists():
            ordered_paths = [parent_root / stage for stage in parse_pipeline_stage_order(parent_root)]
        else:
            ordered_paths = sorted(stage_group)
        for stage_path in ordered_paths:
            if stage_path not in stage_group:
                continue
            if outcar_complete(stage_path):
                continue
            manifest = json.loads((stage_path / "job_manifest.json").read_text())
            items.append(
                QueueItem(
                    kind="stage",
                    path=stage_path,
                    phase=manifest["phase"],
                    description=str(stage_path.relative_to(campaign_root)),
                )
            )
    return items


def campaign_root_for(stage_dir: Path, phase: str) -> Path:
    current = stage_dir
    while current.name != phase:
        current = current.parent
    return current.parent


def prepare_stage_inputs(stage_dir: Path) -> None:
    manifest = json.loads((stage_dir / "job_manifest.json").read_text())
    restart_from = manifest.get("restart_files_from")
    if restart_from:
        campaign_root = campaign_root_for(stage_dir, manifest["phase"])
        source_dir = campaign_root / restart_from
        copy_restart_files(source_dir, stage_dir)


def run_queue_item(item: QueueItem) -> int:
    if item.kind == "pipeline":
        cmd = ["./run_pipeline_local.sh"]
    else:
        prepare_stage_inputs(item.path)
        cmd = ["./run_local.sh"]
    proc = subprocess.run(cmd, cwd=item.path, check=False)
    return proc.returncode


def main() -> None:
    args = parse_args()
    campaign_root = Path(args.campaign_root).resolve()
    phases = phase_filter(args.phases)
    if args.pilot_only:
        phases = {"00_bulk", "01_universal_slab_selection", "02_clean_slab_final", "03_gas_refs"}
    phases.update(args.phase)
    items = enqueue_items(campaign_root, phases)
    if not items:
        print("queue is empty")
        return

    print(f"queue items: {len(items)}")
    for idx, item in enumerate(items, start=1):
        print(f"{idx:04d}  {item.kind:<8}  {item.description}")
    if args.status_only or args.dry_run:
        return

    launched = 0
    skipped = 0
    for item in items:
        if args.max_items and launched >= args.max_items:
            break
        if item.kind == "pipeline" and pipeline_complete(item.path):
            skipped += 1
            continue
        if item.kind == "stage" and outcar_complete(item.path):
            skipped += 1
            continue
        if args.wait_if_busy:
            wait_until_idle(args.poll_seconds)
        print(f"[queue] starting {item.kind}: {item.description}", flush=True)
        returncode = run_queue_item(item)
        print(f"[queue] return code {returncode}: {item.description}", flush=True)
        if returncode != 0:
            sys.exit(returncode)
        if item.kind == "stage" and not outcar_complete(item.path):
            print(f"[queue] stage did not finish cleanly: {item.description}", flush=True)
            sys.exit(2)
        if item.kind == "pipeline" and not pipeline_complete(item.path):
            print(f"[queue] pipeline did not finish cleanly: {item.description}", flush=True)
            sys.exit(3)
        launched += 1
    print(f"[queue] launched={launched} skipped={skipped}", flush=True)


if __name__ == "__main__":
    main()
