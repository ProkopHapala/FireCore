#!/usr/bin/env python3
"""One-time bootstrap for parallel vibration agent jobs. Writes gitignored fixtures under tests/tSiNCs/fixtures/vibration_parallel/."""
import hashlib, json, os, shutil, subprocess, sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
TEST_DIR = Path(__file__).resolve().parent
FIX = REPO / "tests/tSiNCs/fixtures/vibration_parallel"
STRUCT = FIX / "structures"
HESS = FIX / "hessian_mmff"
EXPECTED = FIX / "expected"


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def run(cmd, *, cwd=None):
    print(f"[bootstrap] {' '.join(map(str, cmd))}")
    subprocess.run(cmd, cwd=cwd or REPO, check=True)


def copy_ref(src: Path, dst: Path):
    dst.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(src, dst)
    print(f"[bootstrap] copied {src.relative_to(REPO)} -> {dst.relative_to(REPO)}")


def export_adamantane_hessian():
    import numpy as np
    tmmff = REPO / "tests/tMMFF"
    # MMFF resolves params via tests/tMMFF/data -> common_resources
    if not (tmmff / "data").exists() and not (tmmff / "data").is_symlink():
        run(["bash", "-c", "ln -sf ../../cpp/common_resources data"], cwd=tmmff)
    sys.path.insert(0, str(REPO))
    from pyBall import MMFF, FTIR

    mol2 = REPO / "cpp/common_resources/mol/adamantane.mol2"
    if not mol2.is_file():
        raise FileNotFoundError(mol2)
    orig = os.getcwd()
    try:
        os.chdir(tmmff)
        MMFF.setVerbosity(0, 0)
        MMFF.init(xyz_name=str(mol2.resolve()), bEpairs=False, bMMFF=True)
        MMFF.setSwitches(NonBonded=-1, MMFF=+1, Angles=+1, SurfAtoms=-1, GridFF=-1, PiSigma=-1, PiPiI=-1)
        MMFF.getBuffs()
        n = MMFF.natoms
        inds = np.arange(n, dtype=np.int32)
        H = MMFF.getHessian3Nx3N(inds, dx=1e-4)
        H = 0.5 * (H + H.T)
        M = FTIR.get_mass_matrix(MMFF, n)
        pos = MMFF.apos[:n].copy()
        Hp = FTIR.project_rigid_modes(H, M, pos, shift=1e6)
        m = np.diag(M).reshape(n, 3)[:, 0]
        m_sqrt = np.repeat(np.sqrt(m), 3)
        m_inv_sqrt = 1.0 / m_sqrt
        Hmw = (m_inv_sqrt[:, None] * Hp) * m_inv_sqrt[None, :]
        w = np.linalg.eigvalsh(Hmw)
        omegas = np.sqrt(np.maximum(w, 0.0))
        elements = []
        for i in range(n):
            tn = MMFF.getTypeName(i)
            elements.append(tn.split("_")[0] or "C")
    finally:
        os.chdir(orig)
    HESS.mkdir(parents=True, exist_ok=True)
    out = HESS / "adamantane_mmff_dense.npz"
    np.savez(out, H=H, H_projected=Hp, M=M, pos=pos, omegas_modes=omegas, elements=np.array(elements, dtype=object), natoms=n, source=str(mol2))
    print(f"[bootstrap] wrote {out.relative_to(REPO)} natoms={n} nmodes={len(omegas)}")
    EXPECTED.mkdir(parents=True, exist_ok=True)
    with open(EXPECTED / "adamantane_omegas_modes.txt", "w") as f:
        for i, om in enumerate(omegas):
            f.write(f"{i:4d} {om:.8f}\n")


def export_diamond_nc_hessian():
    tmmff = REPO / "tests/tMMFF"
    if not (tmmff / "data").exists() and not (tmmff / "data").is_symlink():
        run(["bash", "-c", "ln -sf ../../cpp/common_resources data"], cwd=tmmff)
    work = FIX / "work/diamond_nc"
    work.mkdir(parents=True, exist_ok=True)
    py = "/home/prokop/venvs/ML/bin/python3" if Path("/home/prokop/venvs/ML/bin/python3").is_file() else "python3"
    run(["bash", "run.sh", "test_nanocrystal_sparse_hessian.py",
         "--R", "6.0", "--nrep", "5", "--shells", "2", "--max-neigh", "48", "--dx", "1e-4",
         "--outdir", str(work.resolve())], cwd=tmmff)
    relaxed = work / "diamond_nc_R6.0_relaxed.xyz"
    blocks = work / "diamond_nc_R6.0_sparse_hessian.npz"
    init = work / "diamond_nc_R6.0_init.xyz"
    if not relaxed.is_file() or not blocks.is_file():
        raise FileNotFoundError(f"missing diamond_nc outputs under {work}")
    copy_ref(relaxed, STRUCT / "diamond_nc_R6_relaxed.xyz")
    copy_ref(init, STRUCT / "diamond_nc_R6_init.xyz")
    copy_ref(blocks, HESS / "diamond_nc_R6_sparse_blocks.npz")


def run_gen_nanocrystals(prefix: str, extra_args: list):
    out = FIX / "work" / prefix
    out.mkdir(parents=True, exist_ok=True)
    cmd = ["node", os.path.join(TEST_DIR, "gen_nanocrystals.mjs"),
           "--samples", "1", "--seed", "42", "--maxFiles", "1",
           "--caps", "H", "--insertProb", "0", "--collapseProb", "0",
           "--outDir", str(out), "--prefix", prefix] + extra_args
    run(cmd)
    mol2s = sorted(out.glob(f"{prefix}*.mol2"))
    if not mol2s:
        raise FileNotFoundError(f"no mol2 from gen_nanocrystals in {out}")
    src = mol2s[0]
    dst_mol2 = STRUCT / f"{prefix}.mol2"
    copy_ref(src, dst_mol2)
    dst_xyz = STRUCT / f"{prefix}.xyz"
    run(["node", "scripts/mol2_to_xyz.mjs", str(dst_mol2), str(dst_xyz)])


def write_manifest():
    FIX.mkdir(parents=True, exist_ok=True)
    files = {}
    for p in sorted(FIX.rglob("*")):
        if p.is_file() and "work/" not in str(p.relative_to(FIX)):
            rel = str(p.relative_to(REPO))
            files[rel] = {"sha256": sha256(p), "bytes": p.stat().st_size}
    manifest = {
        "repo_root": str(REPO),
        "fixture_root": str(FIX.relative_to(REPO)),
        "generated_by": "tests/tSiNCs/bootstrap_vibration_parallel_fixtures.py",
        "files": files,
    }
    man_path = FIX / "MANIFEST.json"
    with open(man_path, "w") as f:
        json.dump(manifest, f, indent=2)
    print(f"[bootstrap] wrote {man_path.relative_to(REPO)} ({len(files)} files)")


def main():
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--hessian-only", action="store_true", help="Skip structure generation; only MMFF Hessian exports + MANIFEST")
    args = ap.parse_args()
    for d in (STRUCT, HESS, EXPECTED, FIX / "work"):
        d.mkdir(parents=True, exist_ok=True)
    if not args.hessian_only:
        copy_ref(REPO / "cpp/common_resources/mol/adamantane.mol2", STRUCT / "adamantane.mol2")
        copy_ref(REPO / "cpp/common_resources/crystals/diamond_primitive.xyz", STRUCT / "diamond_primitive.xyz")
        run_gen_nanocrystals("si_G1_caps_only", [
            "--nx-range", "1,1", "--ny-range", "1,1", "--nz-range", "1,1", "--centered", "1",
            "--planeTemplates", "a111", "--planeSymC", "6", "--planeCScale", "0.50", "--planeCJitter", "0",
            "--minSiDegree", "2",
        ])
        run_gen_nanocrystals("si_G2_facet111_caps_only", [
            "--nx-range", "2,2", "--ny-range", "2,2", "--nz-range", "2,2", "--centered", "1",
            "--planeTemplates", "a111", "--planeSymC", "6", "--planeCScale", "0.40", "--planeCJitter", "0",
            "--minSiDegree", "2",
        ])
    export_adamantane_hessian()
    export_diamond_nc_hessian()
    write_manifest()
    print("[bootstrap] done.")


if __name__ == "__main__":
    main()
