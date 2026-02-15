
import sys
import os
import time
import numpy as np
import matplotlib.pyplot as plt
from pyBall import AtomicSystem
from pyBall.XPDB_AVBD.RRsp3 import RRsp3, build_neighs_bk_from_bonds, make_bk_slots_clustered, make_exclusions_1st_2nd


def _as_list_floats(s):
    if s is None:
        return []
    s = str(s).strip()
    if not s:
        return []
    return [float(x) for x in s.split(',')]


def _as_list_ints(s):
    if s is None:
        return []
    s = str(s).strip()
    if not s:
        return []
    return [int(x) for x in s.split(',')]


def _rotate_vec_by_quat(q, v):
    t = 2.0 * np.cross(q[:3], v)
    return v + q[3] * t + np.cross(q[:3], t)


def _measure_port_errors(p, q, port_local, neighs_pad, nnode):
    err_sum = 0.0
    err_max = 0.0
    count = 0
    for i in range(nnode):
        qi = q[i]
        for k in range(4):
            j = neighs_pad[i, k]
            if j == -1:
                continue
            vl = port_local[i, k, :3]
            rv = _rotate_vec_by_quat(qi, vl)
            target = p[i, :3] + rv
            diff = p[j, :3] - target
            err = float(np.linalg.norm(diff))
            err_sum += err
            if err > err_max:
                err_max = err
            count += 1
    if count <= 0:
        return 0.0, 0.0
    return err_sum / count, err_max


def _write_xyz_frame(mol_sys, perm, p, natoms, fname, comment, mode):
    mol_sys.apos[perm] = p[:natoms, :3]
    mol_sys.saveXYZ(fname, mode=mode, comment=comment)

def run_test(args):
    # 1. Load Molecule
    print(f"Loading molecule: {args.mol}")
    mol_sys = AtomicSystem.AtomicSystem(args.mol)
    
    # 2. Build Connectivity & Topology
    print("Building topology...")
    # Find bonds if not present (simple distance check or use existing)
    if len(mol_sys.bonds) == 0:
        print("Finding bonds...")
        mol_sys.findBonds(Rcut=1.6) # Standard covalent bond length
    
    # Pack into clusters (simple contiguous packing for now, assuming molecule fits in one or we just pack)
    # For this test, we might treat the whole molecule as one cluster or split it.
    # To be safe and simple, let's repack atoms to ensure nodes are first in groups if we were doing mixed groups.
    # But RRsp3 assumes [Node, Node, ..., Cap, Cap] per group. 
    # For now, let's assume the user provides a "Packable" molecule or we just treat all as nodes for simplicity if ports are not defined?
    # Actually, RRsp3 requires port definitions. 
    # Let's infer ports from geometry for this test.
    
    # For simplicity in this generic test, let's assume sp3 carbons/nitrogens/oxygens are Nodes and H are Caps.
    # But doing full port assignment here is complex. 
    # INSTEAD, let's use the provided 'common_resources/xyz/backbone.xyz' which is likely small.
    # AND, we need to generate valid ports.
    
    # ... Or we can reuse test_RRsp3_vispy's logic which does packing?
    # Let's simplify: Treat ALL heavy atoms as Nodes, H as Caps.
    # Assume we can just load standard ports.
    
    # ACTUALLY, for a robust CONVERGENCE test, we need a valid relaxed structure with valid ports.
    # If we infer ports from the CURRENT geometry, then the current geometry is by definition relaxed (0 strain).
    # This is perfect for a convergence test! 
    # 1. Load XYZ.
    # 2. Infer ports from neighbors (exact ports).
    # 3. Distort molecule.
    # 4. Relax back.
    
    # ... Implementation of "Infer Exact Ports" ...
    # For each atom i:
    #   For each neighbor j:
    #     v = pos[j] - pos[i]
    #     dist = length(v)
    #     dir = v / dist
    #     # We need to assign this bond to a port slot (0..3).
    #     # And we need an orientation Q[i]. 
    #     # For initial state, let's set Q[i] = Identity.
    #     # Then port_local[k] = v (rotated by Identity^-1 = v).
    #     # This guarantees 0 force at start.
    
    natoms = mol_sys.natoms
    pos = np.array(mol_sys.apos, dtype=np.float32)
    
    # Identify neighbors
    neighs = np.full((natoms, 4), -1, dtype=np.int32)
    counts = np.zeros(natoms, dtype=np.int32)
    bonds = []
    
    # Simple distance based bonds if empty
    if len(mol_sys.bonds) == 0:
        for i in range(natoms):
            for j in range(i+1, natoms):
                d2 = np.sum((pos[i]-pos[j])**2)
                if d2 < 1.7**2: # C-C is 1.54
                    mol_sys.bonds.append((i,j))
    
    for (i,j) in mol_sys.bonds:
        bonds.append((i,j))
        if counts[i]<4: neighs[i, counts[i]] = j; counts[i]+=1
        if counts[j]<4: neighs[j, counts[j]] = i; counts[j]+=1
    
    # Pack/Sort atoms: Nodes first?
    # To use RRsp3, we need to respect group_size=64.
    # And nnode_per_group. 
    # If we fit in one group, we can just sort nodes then caps.
    # Let's sort: Nodes (Heavy) then Caps (H).
    
    is_heavy = np.array([s != 'H' for s in mol_sys.enames], dtype=bool)
    ids_heavy = np.where(is_heavy)[0]
    ids_light = np.where(~is_heavy)[0]
    
    perm = np.concatenate([ids_heavy, ids_light])
    inv_perm = np.argsort(perm)
    
    pos = pos[perm]
    new_neighs = np.full((natoms, 4), -1, dtype=np.int32)
    
    # Remap neighbors
    for i in range(natoms):
        old_idx = perm[i]
        for k in range(4):
            old_j = neighs[old_idx, k]
            if old_j != -1:
                new_neighs[i, k] = inv_perm[old_j]
                
    nnode = len(ids_heavy)
    
    # Prepare System
    group_size = 64
    # Pad to group size
    n_padded = ((natoms + group_size - 1) // group_size) * group_size
    pos_pad = np.zeros((n_padded, 4), dtype=np.float32) # .w is invMass
    pos_pad[:natoms, :3] = pos
    
    # Masses (approximate)
    masses = np.ones(natoms, dtype=np.float32)
    masses[:nnode] = 12.0 # C
    masses[nnode:] = 1.0  # H
    pos_pad[:natoms, 3] = 1.0 / masses
    
    # Quats (Identity)
    quat_pad = np.zeros((n_padded, 4), dtype=np.float32)
    quat_pad[:, 3] = 1.0 # w=1
    
    # Ports (Exact)
    # For each node, calculate port_local such that:
    # pos[i] + port_local = pos[j]
    # port_local = pos[j] - pos[i] (since Q=Identity)
    
    # We need to map which neighbor maps to which port slot.
    # We already filled new_neighs.
    
    port_local = np.zeros((n_padded, 4, 4), dtype=np.float32)
    stiffness = np.zeros((n_padded, 4), dtype=np.float32)
    
    for i in range(nnode): # Only nodes have ports
        for k in range(4):
            j = new_neighs[i, k]
            if j != -1:
                v = pos[j] - pos[i]
                d = np.linalg.norm(v)
                port_local[i, k, :3] = v # Exact vector
                stiffness[i, k] = args.stiffness # K
    
    # Topology
    neighs_pad = np.full((n_padded, 4), -1, dtype=np.int32)
    neighs_pad[:natoms] = new_neighs
    
    bkSlots_pad = make_bk_slots_clustered(neighs_pad, group_size=group_size, nnode_per_group=nnode, natoms=n_padded)
    
    # Exclusions (simple 1-2)
    excl1, excl2 = make_exclusions_1st_2nd(neighs_pad)
    
    # --- Initialize Engine ---
    print("Initializing OpenCL...")
    rr = RRsp3(n_padded, group_size=group_size)
    rr.upload_state(pos_pad[:, :3], pos_pad[:, 3], quat=quat_pad, nan_padding=args.nan_padding)
    rr.upload_neighs_and_exclusions(neighs_pad, excl1, excl2)
    rr.upload_cluster_ports(port_local, stiffness, nnode_per_group=nnode) # Assuming 1 group for simplicity or nnode fits
    rr.upload_bkSlots(bkSlots_pad)
    
    # Upload Radius (for collisions)
    radii = np.zeros(n_padded, dtype=np.float32)
    radii[:nnode] = 1.5
    radii[nnode:natoms] = 1.0
    rr.upload_radius(radii)
    
    # --- Distort ---
    print(f"Applying distortion: {args.distort}...")
    pos_ref = pos_pad[:, :3].copy()
    pos_dist = pos_ref.copy()
    if args.distort == 'noise':
        pos_dist += np.random.normal(0.0, float(args.noise_sigma), pos_dist.shape)
    elif args.distort == 'stretch':
        pos_dist[:, 0] *= float(args.fstretch)
    elif args.distort == 'bend':
        x = pos_dist[:, 0]
        y = pos_dist[:, 1]
        z = pos_dist[:, 2]
        R = float(args.bend_R)
        alpha = x / R
        pos_dist[:, 0] = R * np.sin(alpha)
        pos_dist[:, 1] = y
        pos_dist[:, 2] = R * (1.0 - np.cos(alpha)) + z
    d = pos_dist[:natoms, :3] - pos_ref[:natoms, :3]
    dlen = np.linalg.norm(d, axis=1)
    print(f"Distortion displacement: mean={float(np.mean(dlen)):.6g} max={float(np.max(dlen)):.6g}")

    rr.upload_state(pos_dist, pos_pad[:, 3], quat=quat_pad, nan_padding=args.nan_padding)
    
    # --- Relax & Measure ---
    print("Running relaxation...")
    
    betas = [float(x) for x in args.betas.split(',')]
    results = {}
    
    # We need to reset for each beta, so we need to store distorted state
    pos_start = pos_dist.copy()
    quat_start = quat_pad.copy()
    
    for beta in betas:
        print(f" Testing beta={beta}")
        rr.upload_state(pos_start, pos_pad[:, 3], quat=quat_start, nan_padding=args.nan_padding)
        rr.reset_momentum()
        
        errors = []
        errors_max = []
        
        for it in range(args.niter):
            # Step
            rr.step_cluster(
                nnode_per_group=nnode,
                dt=float(args.dt),
                k_coll=float(args.k_coll),
                relaxation=float(args.relaxation),
                momentum_beta=beta
            )
            
            if it % 1 == 0: # Check every step
                # Download and measure constraints
                p, q = rr.download_pos_quat()
                
                avg_err, max_err = _measure_port_errors(p, q, port_local, neighs_pad, nnode)
                errors.append(avg_err)
                errors_max.append(max_err)
                
        results[beta] = (errors, errors_max)
        
        # Save trajectory for the last beta or specific one?
        if args.save_traj:
            # Download current state
            p, q = rr.download_pos_quat()
            # Map back to original order
            # pos_sim[i] corresponds to atom perm[i]
            # so sys.apos[perm[i]] = p[i]
            # We can use: sys.apos[perm] = p[:natoms]
            
            # Update system positions
            mol_sys.apos[perm] = p[:natoms, :3]
            
            fname = f"traj_beta_{beta}.xyz"
            # Append to file if it exists (trajectory) or overwrite? 
            # Usually trajectory files are one file per run.
            # But here we might want one file per beta.
            # Let's overwrite/start new file for each beta loop, but append frames within the loop?
            # Actually, we only save the FINAL state of the relaxation here (outside the loop).
            # If we want a trajectory of the relaxation, we should do it inside the loop.
            
            # Let's save the final relaxed state for now as 'relaxed_beta_X.xyz'
            mol_sys.saveXYZ(fname, mode="w", comment=f"Relaxed with beta={beta} iter={args.niter}")
            print(f"Saved relaxed state to {fname}")

    # --- Plotting ---
    plt.figure(figsize=(10, 6))
    for beta, (errs, errs_max) in results.items():
        plt.plot(errs, label=f'mean beta={beta}')
        plt.plot(errs_max, '--', label=f'max beta={beta}')
    
    plt.yscale('log')
    plt.xlabel('Iteration')
    plt.ylabel('Mean Constraint Violation (Angstrom)')
    plt.title(f'Convergence (distort={args.distort} dt={args.dt} relax={args.relaxation} k={args.stiffness} k_coll={args.k_coll})')
    plt.legend()
    plt.grid(True, which="both", ls="-")
    plt.savefig(args.output)
    print(f"Plot saved to {args.output}")


def run_scan(args):
    print(f"Loading molecule: {args.mol}")
    mol_sys = AtomicSystem.AtomicSystem(args.mol)
    if len(mol_sys.bonds) == 0:
        mol_sys.findBonds(Rcut=1.6)

    natoms = mol_sys.natoms
    pos0 = np.array(mol_sys.apos, dtype=np.float32)

    neighs = np.full((natoms, 4), -1, dtype=np.int32)
    counts = np.zeros(natoms, dtype=np.int32)
    if len(mol_sys.bonds) == 0:
        raise RuntimeError("No bonds found")
    for (i, j) in mol_sys.bonds:
        if counts[i] < 4:
            neighs[i, counts[i]] = j
            counts[i] += 1
        if counts[j] < 4:
            neighs[j, counts[j]] = i
            counts[j] += 1

    is_heavy = np.array([s != 'H' for s in mol_sys.enames], dtype=bool)
    ids_heavy = np.where(is_heavy)[0]
    ids_light = np.where(~is_heavy)[0]
    perm = np.concatenate([ids_heavy, ids_light])
    inv_perm = np.argsort(perm)
    pos0 = pos0[perm]

    new_neighs = np.full((natoms, 4), -1, dtype=np.int32)
    for i in range(natoms):
        old_idx = perm[i]
        for k in range(4):
            old_j = neighs[old_idx, k]
            if old_j != -1:
                new_neighs[i, k] = inv_perm[old_j]

    nnode = int(len(ids_heavy))
    group_size = 64
    n_padded = ((natoms + group_size - 1) // group_size) * group_size
    pos_pad = np.zeros((n_padded, 4), dtype=np.float32)
    pos_pad[:natoms, :3] = pos0

    masses = np.ones(natoms, dtype=np.float32)
    masses[:nnode] = 12.0
    masses[nnode:] = 1.0
    pos_pad[:natoms, 3] = 1.0 / masses

    quat_pad = np.zeros((n_padded, 4), dtype=np.float32)
    quat_pad[:, 3] = 1.0

    port_local = np.zeros((n_padded, 4, 4), dtype=np.float32)
    stiffness = np.zeros((n_padded, 4), dtype=np.float32)
    for i in range(nnode):
        for k in range(4):
            j = new_neighs[i, k]
            if j != -1:
                v = pos0[j] - pos0[i]
                port_local[i, k, :3] = v
                stiffness[i, k] = float(args.stiffness)

    neighs_pad = np.full((n_padded, 4), -1, dtype=np.int32)
    neighs_pad[:natoms] = new_neighs
    bkSlots_pad = make_bk_slots_clustered(neighs_pad, group_size=group_size, nnode_per_group=nnode, natoms=n_padded)
    excl1, excl2 = make_exclusions_1st_2nd(neighs_pad)

    rr = RRsp3(n_padded, group_size=group_size)
    rr.upload_state(pos_pad[:, :3], pos_pad[:, 3], quat=quat_pad, nan_padding=args.nan_padding)
    rr.upload_neighs_and_exclusions(neighs_pad, excl1, excl2)
    rr.upload_cluster_ports(port_local, stiffness, nnode_per_group=nnode)
    rr.upload_bkSlots(bkSlots_pad)
    radii = np.zeros(n_padded, dtype=np.float32)
    radii[:nnode] = 1.5
    radii[nnode:natoms] = 1.0
    rr.upload_radius(radii)

    pos_ref = pos_pad[:, :3].copy()
    pos_start = pos_ref.copy()
    if args.distort != 'stretch':
        raise ValueError("run_scan currently expects --distort stretch")
    pos_start[:, 0] *= float(args.fstretch)
    d = pos_start[:natoms, :3] - pos_ref[:natoms, :3]
    dlen = np.linalg.norm(d, axis=1)
    print(f"Distortion displacement: mean={float(np.mean(dlen)):.6g} max={float(np.max(dlen)):.6g}")

    out_root = os.path.abspath(args.out_dir)
    os.makedirs(out_root, exist_ok=True)

    betas = _as_list_floats(args.bmix_list)
    dts = _as_list_floats(args.dt_list)
    relaxs = _as_list_floats(args.relax_list)
    if not betas:
        betas = [0.0]
    if not dts:
        dts = [float(args.dt)]
    if not relaxs:
        relaxs = [float(args.relaxation)]

    save_every = int(args.save_every)
    if save_every < 1:
        save_every = 1

    thresh2 = float(args.thresh_1e2)
    thresh3 = float(args.thresh_1e3)
    blowup_err = float(args.blowup_err)

    summary_rows = []
    for beta in betas:
        for dt in dts:
            for relax in relaxs:
                tag = f"bmix{beta:g}_dt{dt:g}_rel{relax:g}_fs{float(args.fstretch):g}"
                run_dir = os.path.join(out_root, tag)
                os.makedirs(run_dir, exist_ok=True)
                traj_path = os.path.join(run_dir, "traj.xyz")
                csv_path = os.path.join(run_dir, "errors.csv")
                png_path = os.path.join(run_dir, "convergence.png")

                rr.upload_state(pos_start, pos_pad[:, 3], quat=quat_pad, nan_padding=args.nan_padding)
                rr.reset_momentum()

                it_1e2 = -1
                it_1e3 = -1
                status = "ok"
                errs_mean = []
                errs_max = []

                if args.save_traj:
                    if os.path.exists(traj_path):
                        os.remove(traj_path)
                    p, q = rr.download_pos_quat()
                    _write_xyz_frame(mol_sys, perm, p, natoms, traj_path, comment=f"start {tag}", mode="w")

                k_coll = float(args.k_coll)
                for it in range(int(args.niter)):
                    rr.step_cluster(
                        nnode_per_group=nnode,
                        dt=float(dt),
                        k_coll=float(k_coll),
                        relaxation=float(relax),
                        momentum_beta=float(beta),
                    )
                    p, q = rr.download_pos_quat()
                    e_mean, e_max = _measure_port_errors(p, q, port_local, neighs_pad, nnode)
                    errs_mean.append(e_mean)
                    errs_max.append(e_max)
                    if (not np.isfinite(e_mean)) or (not np.isfinite(e_max)) or (e_max > blowup_err):
                        status = "blowup"
                        break
                    if (it_1e2 < 0) and (e_max < thresh2):
                        it_1e2 = it
                    if (it_1e3 < 0) and (e_max < thresh3):
                        it_1e3 = it
                    if (it_1e3 >= 0) and args.stop_at_1e3:
                        break
                    if args.save_traj and ((it + 1) % save_every == 0):
                        _write_xyz_frame(mol_sys, perm, p, natoms, traj_path, comment=f"it={it+1} {tag} emax={e_max:.6g}", mode="a")

                with open(csv_path, 'w') as f:
                    f.write("it,mean_err,max_err\n")
                    for i, (em, ex) in enumerate(zip(errs_mean, errs_max)):
                        f.write(f"{i},{em:.9g},{ex:.9g}\n")

                plt.figure(figsize=(10, 6))
                plt.plot(errs_mean, label='mean')
                plt.plot(errs_max, '--', label='max')
                plt.yscale('log')
                plt.xlabel('Iteration')
                plt.ylabel('Constraint error [A]')
                plt.title(f"{tag} k={args.stiffness} k_coll={k_coll}")
                plt.legend()
                plt.grid(True, which="both", ls="-")
                plt.savefig(png_path)
                plt.close()

                nsteps = len(errs_max)
                last_err = float(errs_max[-1]) if nsteps else np.nan
                summary_rows.append((tag, beta, dt, relax, nsteps, it_1e2, it_1e3, last_err, status))

    sum_path = os.path.join(out_root, "summary.csv")
    with open(sum_path, 'w') as f:
        f.write("tag,bmix,dt,relax,nsteps,it_1e2,it_1e3,last_max_err,status\n")
        for row in summary_rows:
            f.write(",".join(str(x) for x in row) + "\n")

    print(f"Scan outputs saved to {out_root}")

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Headless Jacobi Convergence Test")
    parser.add_argument("--mol", type=str, required=True, help="Path to molecule .xyz/.mol2")
    parser.add_argument("--distort", type=str, default="noise", choices=['noise', 'stretch', 'bend'])
    parser.add_argument("--niter", type=int, default=100, help="Number of iterations")
    parser.add_argument("--betas", type=str, default="0.0,0.5,0.8,0.9", help="Comma separated betas")
    parser.add_argument("--stiffness", type=float, default=50.0, help="Bond stiffness")
    parser.add_argument("--output", type=str, default="convergence.png", help="Output plot filename")
    parser.add_argument("--save_traj", action="store_true", help="Save trajectory xyz")

    parser.add_argument("--nan_padding", action="store_true", help="Set padded atoms (inv_mass==0) to NaN (debug); default off (use finite padding)")

    parser.add_argument("--dt", type=float, default=0.1)
    parser.add_argument("--relaxation", type=float, default=1.0)
    parser.add_argument("--k_coll", type=float, default=0.0)
    parser.add_argument("--noise_sigma", type=float, default=0.2)
    parser.add_argument("--fstretch", type=float, default=1.2)
    parser.add_argument("--bend_R", type=float, default=10.0)

    parser.add_argument("--scan", action="store_true")
    parser.add_argument("--bmix_list", type=str, default="0.0,0.5,0.7,0.8,0.9,0.95")
    parser.add_argument("--dt_list", type=str, default="0.05,0.1,0.2,0.5,1.0")
    parser.add_argument("--relax_list", type=str, default="0.5,0.8,1.0,1.2")
    parser.add_argument("--save_every", type=int, default=1)
    parser.add_argument("--out_dir", type=str, default=None)
    parser.add_argument("--thresh_1e2", type=float, default=1e-2)
    parser.add_argument("--thresh_1e3", type=float, default=1e-3)
    parser.add_argument("--stop_at_1e3", action="store_true")
    parser.add_argument("--blowup_err", type=float, default=1e3)
    
    args = parser.parse_args()
    if args.out_dir is None:
        ts = time.strftime("scan_%Y%m%d_%H%M%S")
        args.out_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "scan_outputs", ts)
    if args.scan:
        run_scan(args)
    else:
        run_test(args)
