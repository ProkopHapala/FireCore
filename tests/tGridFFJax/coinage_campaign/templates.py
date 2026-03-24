from __future__ import annotations

from collections import OrderedDict
import hashlib
from pathlib import Path

from .config import HPCConfig, ProtocolConfig


def _format_value(value: str | int | float | bool) -> str:
    if isinstance(value, bool):
        return ".TRUE." if value else ".FALSE."
    if isinstance(value, float):
        return f"{value:.8g}"
    return str(value)


def _render_incar(tags: OrderedDict[str, str | int | float | bool]) -> str:
    lines = [f"{key:<10}= {_format_value(value)}" for key, value in tags.items()]
    return "\n".join(lines) + "\n"


def _protocol_tags(
    protocol: ProtocolConfig,
    surface_screened: bool | None = None,
    include_dispersion: bool = True,
) -> OrderedDict[str, str | int | float | bool]:
    tags: OrderedDict[str, str | int | float | bool] = OrderedDict()
    tags["GGA"] = protocol.gga
    if include_dispersion:
        tags["IVDW"] = protocol.ivdw
        use_surface_screened = protocol.ltssurf if surface_screened is None else surface_screened
        if use_surface_screened:
            tags["LTSSURF"] = True
    tags["LASPH"] = protocol.lasph
    tags["ENCUT"] = protocol.encut
    tags["PREC"] = "Accurate"
    tags["EDIFF"] = protocol.ediff
    return tags


def bulk_relax_incar(protocol: ProtocolConfig, system_label: str, stage: str, ncore: int | None = None) -> str:
    tags = _protocol_tags(protocol, surface_screened=False, include_dispersion=False)
    tags["SYSTEM"] = system_label
    tags["ISPIN"] = 1
    tags["ISMEAR"] = protocol.ismear_metal
    tags["SIGMA"] = protocol.sigma_metal
    tags["ALGO"] = "Normal"
    tags["ISYM"] = 0
    tags["IBRION"] = 2
    tags["ISIF"] = 3
    tags["LREAL"] = False
    tags["NSW"] = 50
    tags["EDIFFG"] = -0.01
    if ncore is not None:
        tags["NCORE"] = ncore
    tags["LWAVE"] = False
    tags["LCHARG"] = False
    return _render_incar(tags)


def slab_relax_incar(
    protocol: ProtocolConfig,
    system_label: str,
    stage: str,
    ncore: int | None = None,
    clean_slab: bool = False,
) -> str:
    tags = _protocol_tags(protocol, surface_screened=False)
    tags["SYSTEM"] = system_label
    tags["ISPIN"] = 2
    tags["ISMEAR"] = protocol.ismear_metal
    tags["SIGMA"] = protocol.sigma_metal
    tags["ISYM"] = 0
    tags["IBRION"] = 1 if clean_slab else 2
    tags["ISIF"] = 2
    tags["ALGO"] = "All"
    tags["TIME"] = 0.4
    tags["LREAL"] = False
    tags["NELM"] = 200
    tags["NELMIN"] = 8
    tags["IWAVPR"] = 11
    tags["NELMDL"] = -8
    tags["LORBIT"] = 11
    tags["NSW"] = 300
    tags["EDIFFG"] = protocol.ediffg_relax
    tags["POTIM"] = 0.3
    tags["MAXMIX"] = 40
    if ncore is not None:
        tags["NCORE"] = ncore
    tags["LWAVE"] = True
    tags["LCHARG"] = True
    if "stage1" in stage:
        tags["ISTART"] = 0
    else:
        tags["ISTART"] = 1
        tags["ICHARG"] = 1
        tags["LDIPOL"] = True
        tags["IDIPOL"] = 3
        tags["DIPOL"] = "0.5 0.5 0.5"
    return _render_incar(tags)


def final_static_incar(
    protocol: ProtocolConfig,
    system_label: str,
    metallic: bool = True,
    write_volumetrics: bool = True,
    restart_mode: str = "cold",
    apply_dipole: bool = True,
    surface_screened: bool = True,
    ncore: int | None = None,
    include_dispersion: bool = True,
    spin_polarized: bool = False,
    magmom: str | None = None,
    nupdown: int | None = None,
) -> str:
    tags = _protocol_tags(protocol, surface_screened=surface_screened, include_dispersion=include_dispersion)
    tags["SYSTEM"] = system_label
    tags["ISPIN"] = 2 if (metallic or spin_polarized) else 1
    tags["ISMEAR"] = protocol.ismear_metal if metallic else protocol.ismear_molecule
    tags["SIGMA"] = protocol.sigma_metal if metallic else protocol.sigma_molecule
    tags["ISYM"] = 0
    tags["ALGO"] = "All" if metallic else "Normal"
    tags["NELM"] = 500 if metallic else 200
    if metallic:
        tags["NELMIN"] = 8
    tags["IBRION"] = -1
    tags["NSW"] = 0
    tags["EDIFF"] = 1.0e-7
    tags["LREAL"] = False
    if metallic:
        tags["IWAVPR"] = 11
        tags["NELMDL"] = -8
        tags["LORBIT"] = 11
    if ncore is not None:
        tags["NCORE"] = ncore
    tags["LDIPOL"] = metallic and apply_dipole
    if metallic and apply_dipole:
        tags["IDIPOL"] = 3
        tags["DIPOL"] = "0.5 0.5 0.5"
    if restart_mode == "restart":
        tags["ISTART"] = 1
        tags["ICHARG"] = 1
    else:
        tags["ISTART"] = 0
        tags["ICHARG"] = 2
    tags["LWAVE"] = False
    tags["LCHARG"] = write_volumetrics
    tags["LVHAR"] = write_volumetrics
    tags["LAECHG"] = write_volumetrics
    if magmom is not None:
        tags["MAGMOM"] = magmom
    if nupdown is not None:
        tags["NUPDOWN"] = nupdown
    return _render_incar(tags)


def workfunction_incar(protocol: ProtocolConfig, system_label: str, ncore: int | None = None) -> str:
    tags = _protocol_tags(protocol, surface_screened=False)
    tags["SYSTEM"] = system_label
    tags["ISPIN"] = 2
    tags["ISMEAR"] = protocol.ismear_metal
    tags["SIGMA"] = 0.01
    tags["ALGO"] = "All"
    tags["NELM"] = 300
    tags["NELMIN"] = 8
    tags["IBRION"] = -1
    tags["NSW"] = 0
    tags["ISTART"] = 1
    tags["ICHARG"] = 11
    tags["LREAL"] = False
    tags["IWAVPR"] = 11
    tags["LDIPOL"] = True
    tags["IDIPOL"] = 3
    tags["DIPOL"] = "0.5 0.5 0.5"
    tags["ISYM"] = 0
    tags["LVHAR"] = True
    tags["LVTOT"] = True
    tags["LCHARG"] = True
    tags["LAECHG"] = True
    tags["LWAVE"] = False
    if ncore is not None:
        tags["NCORE"] = ncore
    return _render_incar(tags)


def bader_incar(protocol: ProtocolConfig, system_label: str, ncore: int | None = None) -> str:
    tags = _protocol_tags(protocol, surface_screened=False)
    tags["SYSTEM"] = system_label
    tags["ISPIN"] = 2
    tags["ISMEAR"] = protocol.ismear_metal
    tags["SIGMA"] = 0.01
    tags["ALGO"] = "All"
    tags["NELM"] = 300
    tags["NELMIN"] = 8
    tags["IBRION"] = -1
    tags["NSW"] = 0
    tags["ISTART"] = 1
    tags["ICHARG"] = 11
    tags["LREAL"] = False
    tags["IWAVPR"] = 11
    tags["LDIPOL"] = True
    tags["IDIPOL"] = 3
    tags["DIPOL"] = "0.5 0.5 0.5"
    tags["ISYM"] = 0
    tags["LCHARG"] = True
    tags["LAECHG"] = True
    tags["LWAVE"] = False
    if ncore is not None:
        tags["NCORE"] = ncore
    return _render_incar(tags)


def gas_phase_incar(
    protocol: ProtocolConfig,
    system_label: str,
    stage: str,
    ncore: int | None = None,
    spin_polarized: bool = False,
    magmom: str | None = None,
    nupdown: int | None = None,
) -> str:
    tags = _protocol_tags(protocol, surface_screened=False)
    tags["SYSTEM"] = system_label
    tags["ISPIN"] = 2 if spin_polarized else 1
    tags["ISMEAR"] = protocol.ismear_molecule
    tags["SIGMA"] = protocol.sigma_molecule
    tags["ISYM"] = 0
    tags["ALGO"] = "Normal"
    tags["LREAL"] = False
    tags["NELM"] = 200
    tags["IBRION"] = 2
    tags["ISIF"] = 2
    tags["NSW"] = 50
    tags["EDIFFG"] = -0.01
    if ncore is not None:
        tags["NCORE"] = ncore
    tags["LWAVE"] = False
    tags["LCHARG"] = False
    if magmom is not None:
        tags["MAGMOM"] = magmom
    if nupdown is not None:
        tags["NUPDOWN"] = nupdown
    return _render_incar(tags)


def rigid_scan_incar(protocol: ProtocolConfig, system_label: str, ncore: int | None = None) -> str:
    return final_static_incar(
        protocol,
        system_label,
        metallic=True,
        write_volumetrics=False,
        apply_dipole=True,
        surface_screened=True,
        ncore=ncore,
    )


def relaxed_scan_incar(protocol: ProtocolConfig, system_label: str, ncore: int | None = None) -> str:
    tags = _protocol_tags(protocol, surface_screened=False)
    tags["SYSTEM"] = system_label
    tags["ISPIN"] = 2
    tags["ISMEAR"] = protocol.ismear_metal
    tags["SIGMA"] = protocol.sigma_metal
    tags["ISYM"] = 0
    tags["ALGO"] = "All"
    tags["NELM"] = 200
    tags["NELMIN"] = 8
    tags["IBRION"] = 2
    tags["ISIF"] = 2
    tags["NSW"] = 120
    tags["EDIFFG"] = -0.01
    tags["LREAL"] = False
    tags["LDIPOL"] = True
    tags["IDIPOL"] = 3
    tags["DIPOL"] = "0.5 0.5 0.5"
    tags["LWAVE"] = False
    tags["LCHARG"] = False
    if ncore is not None:
        tags["NCORE"] = ncore
    return _render_incar(tags)


def kpoints_grid(mesh: tuple[int, int, int], gamma: bool = False) -> str:
    mode = "Gamma" if gamma else "Monkhorst-Pack"
    return "\n".join(
        [
            "Generated by coinage_campaign",
            "0",
            mode,
            f"{mesh[0]} {mesh[1]} {mesh[2]}",
            "0 0 0",
            "",
        ]
    )


def copy_forward_snippet(source_stage: str, target_stage: str) -> list[str]:
    return [
        f'cp "{source_stage}/CONTCAR" "{target_stage}/POSCAR"',
        f'if [ -f "{source_stage}/WAVECAR" ]; then cp "{source_stage}/WAVECAR" "{target_stage}/WAVECAR"; fi',
    ]


def local_run_script(vasp_bin: Path) -> str:
    gpu_root = vasp_bin.parent.parent
    gpu_wrapper = gpu_root / "tools" / "run_vasp_gpu.sh"
    variant = vasp_bin.name
    return "\n".join(
        [
            "#!/bin/bash",
            "set -euo pipefail",
            'export OMP_NUM_THREADS="${OMP_NUM_THREADS:-1}"',
            'export ACC_DEVICE_TYPE="${ACC_DEVICE_TYPE:-NVIDIA}"',
            'export ACC_SYNCHRONOUS="${ACC_SYNCHRONOUS:-1}"',
            'export NV_ACC_SYNCHRONOUS="${NV_ACC_SYNCHRONOUS:-1}"',
            f'VASP_GPU_WRAPPER="{gpu_wrapper}"',
            f'VASP_VARIANT="${{VASP_VARIANT:-{variant}}}"',
            'if [ -x "${VASP_GPU_WRAPPER}" ]; then',
            '  time -p "${VASP_GPU_WRAPPER}" "${VASP_VARIANT}" -np "${VASP_NP:-1}" --gpu "${CUDA_DEVICE:-0}" --omp "${OMP_NUM_THREADS}" > vasp.out 2>&1',
            "else",
            f'  VASP_CMD="${{VASP_CMD:-{vasp_bin}}}"',
            '  time -p "${VASP_CMD}" > vasp.out 2>&1',
            "fi",
            "",
        ]
    )


def run_custodian_template(
    hpc: HPCConfig,
    ncore: int,
    binary_env_name: str | None = None,
    binary_default: str | None = None,
) -> str:
    hours, minutes, seconds = (int(part) for part in hpc.walltime.split(":"))
    wall_seconds = max(1, hours * 3600 + minutes * 60 + seconds - 7200)
    env_name = binary_env_name or hpc.vasp_binary_env
    default_bin = binary_default or hpc.vasp_binary_default
    return "\n".join(
        [
            "#!/usr/bin/env python3",
            '"""Run VASP with PBS-aware ntask detection and optional custodian fallback."""',
            "import importlib.util",
            "import os",
            "import subprocess",
            "import sys",
            "",
            "def _detect_ntasks(default=24):",
            '    nodefile = os.environ.get("PBS_NODEFILE")',
            "    if nodefile and os.path.isfile(nodefile) and os.access(nodefile, os.R_OK):",
            "        with open(nodefile) as handle:",
            "            count = sum(1 for _ in handle)",
            "        if count > 0:",
            "            return count",
            '    for key in ("PBS_NTASKS", "PBS_NP", "PBS_NCPUS"):',
            "        value = os.environ.get(key)",
            "        if value and value.isdigit():",
            "            return int(value)",
            "    return default",
            "",
            "def _has_module(name):",
            "    return importlib.util.find_spec(name) is not None",
            "",
            f"ncore = {ncore}",
            f'default_bin = os.environ.get("{env_name}", "{default_bin}")',
            "ntasks = _detect_ntasks()",
            "ntasks = (ntasks // ncore) * ncore",
            "if ntasks < ncore:",
            "    ntasks = ncore",
            "vasp_cmd = [",
            '    "mpiexec",',
            f'    "-genv", "I_MPI_OFI_PROVIDER", "{hpc.mpi_provider}",',
            '    "-np", str(ntasks),',
            "    default_bin,",
            "]",
            'print("Running VASP with %d cores: %s" % (ntasks, " ".join(vasp_cmd)))',
            "",
            'if not (_has_module("custodian") and _has_module("custodian.vasp") and _has_module("custodian.vasp.handlers") and _has_module("custodian.vasp.jobs")):',
            '    print("WARNING: custodian not found; running VASP without custodian.", file=sys.stderr)',
            '    with open("vasp.out", "w") as handle:',
            "        proc = subprocess.run(vasp_cmd, stdout=handle, stderr=subprocess.STDOUT)",
            "    sys.exit(proc.returncode)",
            "",
            "from custodian import Custodian",
            "from custodian.vasp.handlers import (",
            "    FrozenJobErrorHandler,",
            "    NonConvergingErrorHandler,",
            "    PositiveEnergyErrorHandler,",
            "    StdErrHandler,",
            "    UnconvergedErrorHandler,",
            "    VaspErrorHandler,",
            "    WalltimeHandler,",
            ")",
            "from custodian.vasp.jobs import VaspJob",
            "",
            "handlers = [",
            "    VaspErrorHandler(),",
            "    UnconvergedErrorHandler(),",
            "    NonConvergingErrorHandler(),",
            "    PositiveEnergyErrorHandler(),",
            "    FrozenJobErrorHandler(),",
            "    StdErrHandler(),",
            f"    WalltimeHandler(wall_time={wall_seconds}),",
            "]",
            'jobs = [VaspJob(vasp_cmd=vasp_cmd, output_file="vasp.out")]',
            "Custodian(handlers, jobs, max_errors=5).run()",
            "",
        ]
    )


def hpc_pbs_template(
    job_name: str,
    hpc: HPCConfig,
    ncpus: int | None = None,
    mem: str | None = None,
    scratch_local: str | None = None,
    vasp_binary_default: str | None = None,
    use_custodian: bool = False,
) -> str:
    email_line = "#PBS -m bae" if hpc.email_notifications else ""
    module_lines = list(hpc.module_lines)
    digest = hashlib.sha1(job_name.encode("utf-8")).hexdigest()[:6]
    pbs_name = f"{job_name[:24]}_{digest}"
    requested_ncpus = ncpus or hpc.ncpus
    requested_mem = mem or hpc.mem
    requested_scratch = scratch_local or hpc.scratch_local
    binary_default = vasp_binary_default or hpc.vasp_binary_default
    lines = [
        "#!/bin/bash",
        f"#PBS -N {pbs_name}",
        f"#PBS -q {hpc.queue}",
        f"#PBS -l select={hpc.nodes}:ncpus={requested_ncpus}:mem={requested_mem}:scratch_local={requested_scratch}",
        f"#PBS -l walltime={hpc.walltime}",
        f"#PBS -o run_{hpc.queue}.jobout",
        f"#PBS -e run_{hpc.queue}.joberr",
        "#PBS -j oe",
    ]
    if email_line:
        lines.append(email_line)
    lines.extend(
        [
            "",
            "set -euo pipefail",
            "umask 002",
            'WORKDIR="${PBS_O_WORKDIR}"',
            'RUNDIR="${SCRATCHDIR:-${WORKDIR}/.scratch_${PBS_JOBID}}"',
            'mkdir -p "${RUNDIR}"',
            'export TMPDIR="${RUNDIR}"',
            "cleanup_jobdir() {",
            "  status=$?",
            '  shopt -s nullglob',
            '  for name in INCAR KPOINTS POSCAR CONTCAR OUTCAR OSZICAR vasp.out vasprun.xml CHGCAR WAVECAR DOSCAR EIGENVAL IBZKPT PCDAT REPORT XDATCAR AECCAR0 AECCAR2 LOCPOT PROCAR; do',
            '    if [ -e "${RUNDIR}/${name}" ]; then',
            '      cp -a "${RUNDIR}/${name}" "${WORKDIR}/"',
            "    fi",
            "  done",
            '  if [ -n "${SCRATCHDIR:-}" ] && command -v clean_scratch >/dev/null 2>&1; then',
            "    clean_scratch || true",
            '  elif [[ "${RUNDIR}" == "${WORKDIR}/.scratch_"* ]]; then',
            '    rm -rf "${RUNDIR}" || true',
            "  fi",
            "  exit ${status}",
            "}",
            "trap cleanup_jobdir EXIT",
            'for name in INCAR KPOINTS POSCAR POTCAR WAVECAR CHGCAR AECCAR0 AECCAR2 run_custodian.py job_manifest.json; do',
            '  if [ -e "${WORKDIR}/${name}" ]; then',
            '    cp -a "${WORKDIR}/${name}" "${RUNDIR}/"',
            "  fi",
            "done",
            'cd "${RUNDIR}"',
            "",
            *module_lines,
            "",
            f'VASP_BIN="${{{hpc.vasp_binary_env}:-{binary_default}}}"',
            "",
        ]
    )
    if use_custodian:
        lines.append("python3 run_custodian.py")
    else:
        lines.extend(
            [
                'NPROCS=$(wc -l < "${PBS_NODEFILE}")',
                'echo "Running on ${NPROCS} processors"',
                'time -p mpiexec -np "${NPROCS}" "${VASP_BIN}" > vasp.out',
                "",
            ]
        )
    return "\n".join(lines)


def hpc_submit_script() -> str:
    return "\n".join(
        [
            "#!/bin/bash",
            "set -euo pipefail",
            "qsub run.pbs",
            "",
        ]
    )


def relax_pipeline_script(stage_names: tuple[str, ...]) -> str:
    lines = ["#!/bin/bash", "set -euo pipefail"]
    for idx, stage in enumerate(stage_names):
        lines.append(f'echo "Running {stage}"')
        lines.append(f'( cd "{stage}" && ./run_local.sh )')
        if idx + 1 < len(stage_names):
            next_stage = stage_names[idx + 1]
            lines.extend(copy_forward_snippet(stage, next_stage))
    lines.append("")
    return "\n".join(lines)
