#!/usr/bin/env python3
from __future__ import annotations

import argparse
import getpass
import json
import math
import os
import subprocess
import time
import traceback
from datetime import datetime
from pathlib import Path


STAGE_SEQUENCE = (
    "relax_stage1_nodipole",
    "relax_stage2_cycle1_dipole",
    "relax_stage2_cycle2_dipole",
    "relax_stage2_cycle3_dipole",
    "final_static",
    "workfunction",
    "bader",
)
RELAX_STAGES = set(STAGE_SEQUENCE[:4])
OUTPUT_REQUIREMENTS = {
    "workfunction": ("LOCPOT",),
    "bader": ("CHGCAR", "AECCAR0", "AECCAR2"),
}
ADSORBATE_PRIORITY = {
    "H": 0,
    "CO": 1,
    "H2O": 2,
    "HCONH2": 3,
    "pyridine": 4,
}


def parse_args() -> argparse.Namespace:
    return argparse.ArgumentParser(description="Manage the Ag universal-slab-selection pilot on Luna.").parse_args([])


def default_phase_root() -> Path:
    root_env = os.environ.get("COINAGE_CAMPAIGN_ROOT")
    if root_env:
        campaign_root = Path(root_env).expanduser().resolve()
    else:
        home = Path.home()
        campaign_root = home / "work" / "Vasp_work" / "coinage_gridff_dft_hpc"
    return campaign_root / "01_universal_slab_selection" / "Ag" / "3x3x4"


def build_parser() -> argparse.ArgumentParser:
    default_root = default_phase_root()
    manager_dir = default_root.parents[2] / "launchers"
    parser = argparse.ArgumentParser(description="Manage the Ag universal-slab-selection pilot on Luna.")
    parser.add_argument("--root", type=Path, default=default_root)
    parser.add_argument("--state-path", type=Path, default=manager_dir / "coin_ag_selection_manager_state.json")
    parser.add_argument("--log-path", type=Path, default=manager_dir / "coin_ag_selection_manager.log")
    parser.add_argument("--vasp-bin", default="/storage/praha1/home/indranil/src/vasp_std")
    parser.add_argument("--poll-seconds", type=int, default=300)
    parser.add_argument("--queue-grace-seconds", type=int, default=7200)
    parser.add_argument("--disp-threshold", type=float, default=5.0e-7)
    parser.add_argument("--max-active", type=int, default=24)
    parser.add_argument("--max-relax-cycles", type=int, default=8)
    parser.add_argument(
        "--adsorbates",
        default="",
        help="Comma-separated representative adsorbates to manage under representative_adsorbates/. Empty means all.",
    )
    return parser


def log(msg: str, log_path: Path) -> None:
    line = "[%s] %s" % (datetime.now().strftime("%Y-%m-%d %H:%M:%S"), msg)
    print(line, flush=True)
    log_path.parent.mkdir(parents=True, exist_ok=True)
    with log_path.open("a") as handle:
        handle.write(line + "\n")


def qstat_live_jobids(user: str | None = None) -> set[str] | None:
    target_user = user or os.environ.get("USER") or getpass.getuser()
    proc = subprocess.run(
        ["qstat", "-u", target_user],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )
    if proc.returncode != 0:
        return None
    live: set[str] = set()
    for line in (proc.stdout or "").splitlines():
        parts = line.split()
        if len(parts) >= 1 and ".pbs-" in parts[0]:
            live.add(parts[0].split(".", 1)[0])
    return live


def read_text(path: Path) -> str:
    return path.read_text(errors="ignore") if path.exists() else ""


def parse_poscar(path: Path):
    lines = [line.rstrip() for line in path.read_text().splitlines() if line.strip()]
    scale = float(lines[1].split()[0])
    cell = []
    for i in range(3):
        cell.append([float(x) * scale for x in lines[2 + i].split()[:3]])
    counts_line = 5
    try:
        int(lines[5].split()[0])
    except ValueError:
        counts_line = 6
    counts = [int(x) for x in lines[counts_line].split()]
    natoms = sum(counts)
    idx = counts_line + 1
    if lines[idx].lower().startswith("s"):
        idx += 1
    coord_mode = lines[idx].strip().lower()[0]
    idx += 1
    positions = []
    for i in range(natoms):
        positions.append([float(x) for x in lines[idx + i].split()[:3]])
    return cell, coord_mode, positions


def is_valid_structure(path: Path) -> bool:
    if not path.exists() or path.stat().st_size == 0:
        return False
    try:
        parse_poscar(path)
    except Exception:
        return False
    return True


def det3(m):
    return (
        m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
        - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
        + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0])
    )


def inv3(m):
    d = det3(m)
    if abs(d) < 1e-16:
        raise ValueError("Singular cell matrix")
    return [
        [(m[1][1] * m[2][2] - m[1][2] * m[2][1]) / d, -(m[0][1] * m[2][2] - m[0][2] * m[2][1]) / d, (m[0][1] * m[1][2] - m[0][2] * m[1][1]) / d],
        [-(m[1][0] * m[2][2] - m[1][2] * m[2][0]) / d, (m[0][0] * m[2][2] - m[0][2] * m[2][0]) / d, -(m[0][0] * m[1][2] - m[0][2] * m[1][0]) / d],
        [(m[1][0] * m[2][1] - m[1][1] * m[2][0]) / d, -(m[0][0] * m[2][1] - m[0][1] * m[2][0]) / d, (m[0][0] * m[1][1] - m[0][1] * m[1][0]) / d],
    ]


def matvec(m, v):
    return [sum(m[i][j] * v[j] for j in range(3)) for i in range(3)]


def to_frac(cell, mode, positions):
    if mode == "d":
        return positions
    inv = inv3(cell)
    return [matvec(inv, p) for p in positions]


def max_displacement(poscar: Path, contcar: Path) -> float:
    cell1, mode1, pos1 = parse_poscar(poscar)
    cell2, mode2, pos2 = parse_poscar(contcar)
    frac1 = to_frac(cell1, mode1, pos1)
    frac2 = to_frac(cell2, mode2, pos2)
    max_d = 0.0
    for a, b in zip(frac1, frac2):
        df = [b[i] - a[i] for i in range(3)]
        df = [x - round(x) for x in df]
        dc = matvec(cell1, df)
        dist = math.sqrt(sum(x * x for x in dc))
        if dist > max_d:
            max_d = dist
    return max_d


def stage_success(stage_dir: Path, stage: str) -> bool:
    outcar = stage_dir / "OUTCAR"
    if not outcar.exists():
        return False
    text = read_text(outcar)
    if stage in RELAX_STAGES:
        return (
            ("reached required accuracy" in text) or ("aborting loop because EDIFF is reached" in text)
        ) and is_valid_structure(stage_dir / "CONTCAR")
    if not (("aborting loop because EDIFF is reached" in text) or ("General timing and accounting informations for this job" in text)):
        return False
    required = OUTPUT_REQUIREMENTS.get(stage, ())
    return all((stage_dir / name).exists() and (stage_dir / name).stat().st_size > 0 for name in required)


def stage_restartable(stage_dir: Path, stage: str) -> bool:
    if stage not in RELAX_STAGES:
        return False
    outcar_text = read_text(stage_dir / "OUTCAR")
    has_contcar = is_valid_structure(stage_dir / "CONTCAR")
    return has_contcar and ("ZBRENT" in outcar_text or "reached required accuracy" in outcar_text or "aborting loop because EDIFF is reached" in outcar_text)


def copy_atomic(src: Path, dst: Path) -> None:
    tmp = dst.with_name(dst.name + ".tmp")
    tmp.write_bytes(src.read_bytes())
    tmp.replace(dst)


def snapshot_outputs(stage_dir: Path, tag: str) -> None:
    backup_dir = stage_dir / tag
    backup_dir.mkdir(exist_ok=True)
    for name in (
        "CONTCAR",
        "OSZICAR",
        "OUTCAR",
        "INCAR",
        "KPOINTS",
        "POSCAR",
        "run_luna.jobout",
        "run_luna.joberr",
        "vasp.out",
    ):
        path = stage_dir / name
        if path.exists() and path.stat().st_size > 0:
            copy_atomic(path, backup_dir / name)


def copy_forward(root: Path, source_stage: str, target_stage: str) -> None:
    src = root / source_stage
    dst = root / target_stage
    contcar = src / "CONTCAR"
    if not is_valid_structure(contcar):
        raise FileNotFoundError(f"Invalid CONTCAR for restart: {contcar}")
    copy_atomic(contcar, dst / "POSCAR")
    if src == dst:
        return
    wavecar = src / "WAVECAR"
    if wavecar.exists() and wavecar.stat().st_size > 0:
        target_wavecar = dst / "WAVECAR"
        if target_wavecar.exists():
            target_wavecar.unlink()
        wavecar.replace(target_wavecar)


def submit(stage_dir: Path, vasp_bin: str) -> str:
    attempts = 3
    for attempt in range(1, attempts + 1):
        try:
            proc = subprocess.run(
                ["qsub", "-v", f"VASP_BIN_ON_HPC={vasp_bin}", "run.pbs"],
                cwd=stage_dir,
                text=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                timeout=120,
            )
        except subprocess.TimeoutExpired:
            if attempt == attempts:
                raise RuntimeError(f"qsub timed out in {stage_dir} after {attempts} attempts")
            time.sleep(20)
            continue
        if proc.returncode == 0:
            return proc.stdout.strip()
        if attempt == attempts:
            stderr = (proc.stderr or "").strip()
            stdout = (proc.stdout or "").strip()
            raise RuntimeError(
                f"qsub failed in {stage_dir} after {attempts} attempts; stdout={stdout!r} stderr={stderr!r}"
            )
        time.sleep(30)
    raise RuntimeError(f"qsub failed in {stage_dir}")


def stage_manifest(root: Path, stage: str) -> dict:
    return json.loads((root / stage / "job_manifest.json").read_text())


def next_stage(stage: str) -> str | None:
    idx = STAGE_SEQUENCE.index(stage)
    if idx + 1 >= len(STAGE_SEQUENCE):
        return None
    return STAGE_SEQUENCE[idx + 1]


def discover_roots(root: Path, allowed_adsorbates: set[str] | None = None) -> list[Path]:
    roots = [root / "clean_slab"]
    rep_root = root / "representative_adsorbates"
    if rep_root.exists():
        for ads_root in sorted(p for p in rep_root.iterdir() if p.is_dir()):
            if allowed_adsorbates and ads_root.name not in allowed_adsorbates:
                continue
            roots.extend(sorted(p for p in ads_root.glob("*") if p.is_dir()))
    return roots


def root_key(root: Path, phase_root: Path) -> str:
    return str(root.relative_to(phase_root))


def root_priority(key: str) -> tuple[int, int, str]:
    if key == "clean_slab":
        return (0, -1, key)
    parts = Path(key).parts
    adsorbate = parts[1] if len(parts) > 1 else ""
    return (1, ADSORBATE_PRIORITY.get(adsorbate, 99), key)


def initial_state(phase_root: Path, allowed_adsorbates: set[str] | None = None) -> dict[str, dict[str, object]]:
    state: dict[str, dict[str, object]] = {}
    for root in discover_roots(phase_root, allowed_adsorbates):
        key = root_key(root, phase_root)
        state[key] = {
            "stage": STAGE_SEQUENCE[0],
            "jobid": None,
            "cycle": 0,
            "done": False,
            "failed": False,
        }
    return state


def load_state(state_path: Path, phase_root: Path, allowed_adsorbates: set[str] | None = None) -> dict[str, dict[str, object]]:
    if state_path.exists():
        state = json.loads(state_path.read_text())
    else:
        state = initial_state(phase_root, allowed_adsorbates)
    current_keys = {root_key(root, phase_root) for root in discover_roots(phase_root, allowed_adsorbates)}
    state = {key: value for key, value in state.items() if key in current_keys}
    for root in discover_roots(phase_root, allowed_adsorbates):
        key = root_key(root, phase_root)
        state.setdefault(
            key,
            {"stage": STAGE_SEQUENCE[0], "jobid": None, "cycle": 0, "done": False, "failed": False, "submitted_at": None},
        )
    now = time.time()
    for entry in state.values():
        entry.setdefault("submitted_at", None)
        if entry.get("jobid") and entry.get("submitted_at") is None:
            entry["submitted_at"] = now
    return state


def save_state(state_path: Path, state: dict[str, dict[str, object]]) -> None:
    state_path.write_text(json.dumps(state, indent=2) + "\n")


def maybe_resubmit_relax(root: Path, key: str, entry: dict[str, object], stage: str, stage_dir: Path, vasp_bin: str, max_relax_cycles: int, log_path: Path) -> bool:
    cycle = int(entry.get("cycle", 0))
    if cycle >= max_relax_cycles:
        entry["failed"] = True
        log(f"{key}: exceeded max relax cycles in {stage}", log_path)
        return True
    snapshot_outputs(stage_dir, f"retry_{stage}_{cycle + 1:02d}")
    copy_forward(root, stage, stage)
    jobid = submit(stage_dir, vasp_bin)
    entry.update({"jobid": jobid, "cycle": cycle + 1, "submitted_at": time.time()})
    log(f"{key}: resubmitted {stage} as {jobid}", log_path)
    return True


def maybe_resubmit_relax_limited(
    root: Path,
    key: str,
    entry: dict[str, object],
    stage: str,
    stage_dir: Path,
    vasp_bin: str,
    max_relax_cycles: int,
    log_path: Path,
    active: list[int],
    max_active: int,
) -> bool:
    if active[0] >= max_active:
        log(f"{key}: deferred resubmission of {stage}; active cap {max_active} reached", log_path)
        return True
    handled = maybe_resubmit_relax(root, key, entry, stage, stage_dir, vasp_bin, max_relax_cycles, log_path)
    if handled:
        active[0] += 1
    return handled


def advance_entry(
    phase_root: Path,
    key: str,
    entry: dict[str, object],
    live_jobids: set[str],
    active: list[int],
    max_active: int,
    vasp_bin: str,
    disp_threshold: float,
    max_relax_cycles: int,
    queue_grace_seconds: int,
    log_path: Path,
) -> None:
    if entry.get("done") or entry.get("failed"):
        return
    root = phase_root / key
    stage = str(entry["stage"])
    jobid = entry.get("jobid")
    live = isinstance(jobid, str) and jobid.split(".", 1)[0] in live_jobids
    if live:
        return

    stage_dir = root / stage
    outcar_exists = (stage_dir / "OUTCAR").exists()
    if jobid is not None and not live and not outcar_exists:
        submitted_at = float(entry.get("submitted_at") or time.time())
        if time.time() - submitted_at < queue_grace_seconds:
            return
        entry["jobid"] = None
        entry["submitted_at"] = None
        log(f"{key}: cleared stale jobid for {stage} after queue grace timeout", log_path)
        return
    if jobid is None and not outcar_exists:
        return

    if stage_success(stage_dir, stage):
        if stage in RELAX_STAGES:
            disp = max_displacement(stage_dir / "POSCAR", stage_dir / "CONTCAR")
            log(f"{key}: {stage} cycle {entry.get('cycle', 0)} max displacement {disp:.6e} A", log_path)
            if disp > disp_threshold:
                maybe_resubmit_relax_limited(
                    root,
                    key,
                    entry,
                    stage,
                    stage_dir,
                    vasp_bin,
                    max_relax_cycles,
                    log_path,
                    active,
                    max_active,
                )
                return
        nxt = next_stage(stage)
        if nxt is None:
            entry.update({"done": True, "jobid": None})
            entry["submitted_at"] = None
            log(f"{key}: workflow complete", log_path)
            return
        manifest = stage_manifest(root, nxt)
        source_stage = manifest.get("source_stage", stage)
        if active[0] >= max_active:
            log(f"{key}: deferred submission of {nxt}; active cap {max_active} reached", log_path)
            return
        copy_forward(root, str(source_stage), nxt)
        try:
            new_job = submit(root / nxt, vasp_bin)
        except Exception as exc:
            log(f"{key}: submit failed for {nxt}: {exc}", log_path)
            return
        entry.update({"stage": nxt, "jobid": new_job, "cycle": 0, "submitted_at": time.time()})
        active[0] += 1
        log(f"{key}: submitted {nxt} as {new_job}", log_path)
        return

    if stage_restartable(stage_dir, stage):
        maybe_resubmit_relax_limited(
            root,
            key,
            entry,
            stage,
            stage_dir,
            vasp_bin,
            max_relax_cycles,
            log_path,
            active,
            max_active,
        )
        return

    entry["failed"] = True
    log(f"{key}: {stage} finished without success marker", log_path)


def count_active(state: dict[str, dict[str, object]], live_jobids: set[str]) -> int:
    total = 0
    for entry in state.values():
        jobid = entry.get("jobid")
        if isinstance(jobid, str) and jobid.split(".", 1)[0] in live_jobids:
            total += 1
    return total


def submit_new_work(
    phase_root: Path,
    state: dict[str, dict[str, object]],
    active: list[int],
    vasp_bin: str,
    max_active: int,
    log_path: Path,
) -> None:
    for key in sorted(state, key=root_priority):
        if active[0] >= max_active:
            break
        entry = state[key]
        if entry.get("done") or entry.get("failed"):
            continue
        if entry.get("jobid"):
            continue
        stage = str(entry["stage"])
        stage_dir = phase_root / key / stage
        if (stage_dir / "OUTCAR").exists():
            continue
        try:
            jobid = submit(stage_dir, vasp_bin)
        except Exception as exc:
            log(f"{key}: submit failed for {stage}: {exc}", log_path)
            continue
        if stage in RELAX_STAGES and int(entry.get("cycle", 0)) == 0:
            entry["cycle"] = 1
        entry["jobid"] = jobid
        entry["submitted_at"] = time.time()
        active[0] += 1
        log(f"{key}: submitted {stage} as {jobid}", log_path)


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()
    phase_root = args.root.resolve()
    state_path = args.state_path.resolve()
    log_path = args.log_path.resolve()
    allowed_adsorbates = {item.strip() for item in args.adsorbates.split(",") if item.strip()}
    state = load_state(state_path, phase_root, allowed_adsorbates or None)
    log("Ag selection manager started", log_path)
    while True:
        try:
            live_jobids = qstat_live_jobids()
            if live_jobids is None:
                log("qstat -u snapshot failed; skipping this polling cycle", log_path)
                time.sleep(min(args.poll_seconds, 60))
                continue
            active = [count_active(state, live_jobids)]
            for key in sorted(state, key=root_priority):
                advance_entry(
                    phase_root,
                    key,
                    state[key],
                    live_jobids,
                    active,
                    args.max_active,
                    args.vasp_bin,
                    args.disp_threshold,
                    args.max_relax_cycles,
                    args.queue_grace_seconds,
                    log_path,
                )
            submit_new_work(phase_root, state, active, args.vasp_bin, args.max_active, log_path)
            save_state(state_path, state)
            if all(entry.get("done") or entry.get("failed") for entry in state.values()) and count_active(state, live_jobids) == 0:
                log("all managed Ag selection workflows finished", log_path)
                break
            time.sleep(args.poll_seconds)
        except KeyboardInterrupt:
            log("Ag selection manager interrupted", log_path)
            save_state(state_path, state)
            raise
        except Exception:
            log("manager exception:\n" + traceback.format_exc(), log_path)
            save_state(state_path, state)
            time.sleep(min(args.poll_seconds, 60))


if __name__ == "__main__":
    main()
