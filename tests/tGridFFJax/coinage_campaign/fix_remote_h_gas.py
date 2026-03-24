#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import shutil
import subprocess
from datetime import datetime
from pathlib import Path


def patch_incar(path: Path) -> None:
    lines = path.read_text().splitlines()
    out: list[str] = []
    seen = {"NCORE": False, "ISPIN": False, "MAGMOM": False, "NUPDOWN": False}
    for line in lines:
        key = line.split("=")[0].strip() if "=" in line else ""
        if key == "NCORE":
            out.append("NCORE     = 1")
            seen["NCORE"] = True
        elif key == "ISPIN":
            out.append("ISPIN     = 2")
            seen["ISPIN"] = True
        elif key == "MAGMOM":
            out.append("MAGMOM    = 1")
            seen["MAGMOM"] = True
        elif key == "NUPDOWN":
            out.append("NUPDOWN   = 1")
            seen["NUPDOWN"] = True
        else:
            out.append(line)
    if not seen["NCORE"]:
        out.append("NCORE     = 1")
    if not seen["ISPIN"]:
        out.append("ISPIN     = 2")
    if not seen["MAGMOM"]:
        out.append("MAGMOM    = 1")
    if not seen["NUPDOWN"]:
        out.append("NUPDOWN   = 1")
    path.write_text("\n".join(out) + "\n")


def patch_pbs(path: Path) -> None:
    lines = path.read_text().splitlines()
    out: list[str] = []
    for line in lines:
        if line.startswith("#PBS -l select="):
            out.append("#PBS -l select=1:ncpus=4:mem=8gb")
        else:
            out.append(line)
    path.write_text("\n".join(out) + "\n")


def main() -> None:
    parser = argparse.ArgumentParser(description="Fix atomic-H gas reference on remote HPC campaign.")
    parser.add_argument("--campaign-root", required=True)
    parser.add_argument("--state-path", required=True)
    parser.add_argument("--stdout-path", required=True)
    parser.add_argument("--manager-path", required=True)
    parser.add_argument("--vasp-bin", required=True)
    args = parser.parse_args()

    root = Path(args.campaign_root).resolve() / "H"
    state_path = Path(args.state_path).resolve()
    stdout_path = Path(args.stdout_path).resolve()
    manager_path = Path(args.manager_path).resolve()

    state = json.loads(state_path.read_text())
    old_jobid = state["H"].get("jobid")
    if isinstance(old_jobid, str) and old_jobid:
        subprocess.run(["qdel", old_jobid], check=False)

    for stage in ("relax_stage1", "final_static"):
        patch_incar(root / stage / "INCAR")
        patch_pbs(root / stage / "run.pbs")

    stage1 = root / "relax_stage1"
    backup_dir = stage1 / datetime.now().strftime("failed_H_spinfix_%Y%m%d_%H%M%S")
    backup_dir.mkdir(exist_ok=True)
    for name in (
        "CHG",
        "CHGCAR",
        "CONTCAR",
        "DOSCAR",
        "EIGENVAL",
        "IBZKPT",
        "OSZICAR",
        "OUTCAR",
        "PCDAT",
        "REPORT",
        "WAVECAR",
        "XDATCAR",
        "run_luna.jobout",
        "run_luna.joberr",
        "vasp.out",
    ):
        path = stage1 / name
        if path.exists():
            shutil.move(str(path), str(backup_dir / name))

    jobid = subprocess.check_output(
        ["qsub", "-v", f"VASP_BIN_ON_HPC={args.vasp_bin}", "run.pbs"],
        cwd=stage1,
        text=True,
    ).strip()
    state["H"] = {"stage": "relax_stage1", "jobid": jobid, "cycle": 1, "done": False, "failed": False}
    state_path.write_text(json.dumps(state, indent=2) + "\n")

    subprocess.run(["pkill", "-f", "coin_gas_manager.py"], check=False)
    if stdout_path.exists():
        stdout_path.unlink()
    with stdout_path.open("w") as handle:
        subprocess.Popen(
            ["nohup", "python3", str(manager_path)],
            stdout=handle,
            stderr=subprocess.STDOUT,
            cwd=str(root.parent.parent),
            start_new_session=True,
        )

    print(jobid)


if __name__ == "__main__":
    main()
