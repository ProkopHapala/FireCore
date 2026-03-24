#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import os
import shutil
import subprocess
from datetime import datetime
from pathlib import Path


MOLECULES = ("CO", "H", "H2O", "HCONH2", "NH3", "methanol", "pyridine")


def patch_incar(incar_path: Path, ncore: int, spin_polarized: bool = False) -> None:
    lines = incar_path.read_text().splitlines()
    patched: list[str] = []
    seen_ncore = False
    seen_ispin = False
    seen_magmom = False
    seen_nupdown = False
    for line in lines:
        stripped = line.strip()
        if stripped.startswith("NCORE"):
            patched.append(f"NCORE     = {ncore}")
            seen_ncore = True
        elif stripped.startswith("ISPIN"):
            patched.append(f"ISPIN     = {2 if spin_polarized else 1}")
            seen_ispin = True
        elif stripped.startswith("MAGMOM"):
            if spin_polarized:
                patched.append("MAGMOM    = 1")
            seen_magmom = True
        elif stripped.startswith("NUPDOWN"):
            if spin_polarized:
                patched.append("NUPDOWN   = 1")
            seen_nupdown = True
        else:
            patched.append(line)
    if not seen_ncore:
        patched.append(f"NCORE     = {ncore}")
    if not seen_ispin:
        patched.append(f"ISPIN     = {2 if spin_polarized else 1}")
    if spin_polarized and not seen_magmom:
        patched.append("MAGMOM    = 1")
    if spin_polarized and not seen_nupdown:
        patched.append("NUPDOWN   = 1")
    incar_path.write_text("\n".join(patched) + "\n")


def patch_pbs_resources(run_pbs_path: Path, ncpus: int, mem: str) -> None:
    lines = run_pbs_path.read_text().splitlines()
    patched: list[str] = []
    for line in lines:
        if line.startswith("#PBS -l select="):
            patched.append(f"#PBS -l select=1:ncpus={ncpus}:mem={mem}")
        else:
            patched.append(line)
    run_pbs_path.write_text("\n".join(patched) + "\n")


def patch_pbs_binary(run_pbs_path: Path, vasp_bin: str) -> None:
    lines = run_pbs_path.read_text().splitlines()
    patched: list[str] = []
    for line in lines:
        if line.startswith('VASP_BIN="${'):
            patched.append(f'VASP_BIN="${{VASP_BIN_ON_HPC:-{vasp_bin}}}"')
        else:
            patched.append(line)
    run_pbs_path.write_text("\n".join(patched) + "\n")


def patch_gas_manager(manager_path: Path) -> None:
    old = '    return ("reached required accuracy" in text) and (stage_dir / "CONTCAR").exists()\n'
    new = (
        '    contcar = stage_dir / "CONTCAR"\n'
        '    return ((("reached required accuracy" in text) or ("aborting loop because EDIFF is reached" in text)) '
        'and contcar.exists() and contcar.stat().st_size > 0)\n'
    )
    text = manager_path.read_text()
    if old in text:
        manager_path.write_text(text.replace(old, new))


def backup_stage_outputs(stage_dir: Path, backup_tag: str) -> None:
    backup_dir = stage_dir / backup_tag
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
        path = stage_dir / name
        if path.exists():
            shutil.move(str(path), str(backup_dir / name))


def submit(stage_dir: Path, vasp_bin: str) -> str:
    return subprocess.check_output(
        ["qsub", "-v", f"VASP_BIN_ON_HPC={vasp_bin}", "run.pbs"],
        cwd=stage_dir,
        text=True,
    ).strip()


def stop_existing_manager(manager_path: Path) -> None:
    this_pid = os.getpid()
    proc = subprocess.run(
        ["ps", "-eo", "pid=,args="],
        check=True,
        text=True,
        capture_output=True,
    )
    for line in proc.stdout.splitlines():
        stripped = line.strip()
        if not stripped:
            continue
        pid_text, _, args = stripped.partition(" ")
        if not pid_text.isdigit():
            continue
        pid = int(pid_text)
        if pid == this_pid:
            continue
        if "python3" in args and str(manager_path) in args:
            subprocess.run(["kill", str(pid)], check=False)


def main() -> None:
    parser = argparse.ArgumentParser(description="Recover and resubmit remote gas-reference jobs.")
    parser.add_argument("--campaign-root", required=True)
    parser.add_argument("--state-path", required=True)
    parser.add_argument("--manager-path", required=True)
    parser.add_argument("--stdout-path", required=True)
    parser.add_argument("--vasp-bin", required=True)
    parser.add_argument("--ncore", type=int, default=1)
    parser.add_argument("--ncpus", type=int, default=4)
    parser.add_argument("--mem", default="8gb")
    args = parser.parse_args()

    campaign_root = Path(args.campaign_root).resolve()
    state_path = Path(args.state_path).resolve()
    manager_path = Path(args.manager_path).resolve()
    stdout_path = Path(args.stdout_path).resolve()

    patch_gas_manager(manager_path)
    stop_existing_manager(manager_path)

    if state_path.exists():
        old_state = json.loads(state_path.read_text())
        for payload in old_state.values():
            jobid = payload.get("jobid")
            if isinstance(jobid, str) and jobid:
                subprocess.run(["qdel", jobid], check=False)

    backup_tag = datetime.now().strftime("failed_gas_%Y%m%d_%H%M%S")
    state: dict[str, dict[str, object]] = {}
    for mol in MOLECULES:
        mol_root = campaign_root / mol
        spin_polarized = mol == "H"
        for stage in ("relax_stage1", "final_static"):
            patch_incar(mol_root / stage / "INCAR", args.ncore, spin_polarized=spin_polarized)
            patch_pbs_resources(mol_root / stage / "run.pbs", args.ncpus, args.mem)
            patch_pbs_binary(mol_root / stage / "run.pbs", args.vasp_bin)
        stage1 = mol_root / "relax_stage1"
        backup_stage_outputs(stage1, backup_tag)
        jobid = submit(stage1, args.vasp_bin)
        state[mol] = {
            "stage": "relax_stage1",
            "jobid": jobid,
            "cycle": 1,
            "done": False,
            "failed": False,
        }

    state_path.write_text(json.dumps(state, indent=2) + "\n")
    if stdout_path.exists():
        stdout_path.unlink()

    with stdout_path.open("w") as handle:
        subprocess.Popen(
            ["nohup", "python3", str(manager_path)],
            stdout=handle,
            stderr=subprocess.STDOUT,
            cwd=str(campaign_root.parent.parent),
            start_new_session=True,
        )

    print(json.dumps(state, indent=2))


if __name__ == "__main__":
    main()
