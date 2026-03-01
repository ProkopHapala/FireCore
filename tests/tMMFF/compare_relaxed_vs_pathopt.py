#!/usr/bin/env python3
import os, argparse
import numpy as np
import matplotlib.pyplot as plt

np.set_printoptions(precision=6, suppress=True)


def read_npz_seq(path):
    D = np.load(path)
    return D['A_scan'], D['tip_scan']


def read_xyz_frames(path):
    with open(path) as f:
        lines = f.read().splitlines()
    i = 0
    frames = []
    while i < len(lines):
        nat = int(lines[i]); i += 2  # skip comment
        pos = []
        for _ in range(nat):
            w = lines[i].split()
            pos.append([float(w[1]), float(w[2]), float(w[3])])
            i += 1
        frames.append(pos)
    return np.array(frames, dtype=float)


def plot_xz(seq, par, tip_seq, out_png):
    atoms = ['O','H1','H2'] if seq.shape[1] == 3 else [f'a{i}' for i in range(seq.shape[1])]
    plt.figure(figsize=(7,6))
    # tip
    plt.plot(tip_seq[:,0], tip_seq[:,2], 'k--', lw=2, label='tip (seq)')
    for ia in range(seq.shape[1]):
        plt.plot(seq[:,ia,0], seq[:,ia,2], '-', lw=1.5, label=f'seq {atoms[ia]}')
        plt.plot(par[:,ia,0], par[:,ia,2], '--', lw=1.2, label=f'par {atoms[ia]}')
    plt.xlabel('x (A)'); plt.ylabel('z (A)'); plt.axis('equal'); plt.grid(True)
    plt.legend(); plt.tight_layout(); plt.savefig(out_png); plt.close()


def plot_divergence(seq, par, out_png):
    n = min(len(seq), len(par))
    seq = seq[:n]; par = par[:n]
    atoms = ['O','H1','H2'] if seq.shape[1] == 3 else [f'a{i}' for i in range(seq.shape[1])]
    plt.figure(figsize=(7,4))
    for ia in range(seq.shape[1]):
        diff = par[:,ia,:] - seq[:,ia,:]
        dr = np.linalg.norm(diff, axis=1)
        plt.plot(range(n), dr, label=f'{atoms[ia]} |d|')
    plt.xlabel('scan step'); plt.ylabel('|par - seq| (A)'); plt.grid(True); plt.legend()
    plt.tight_layout(); plt.savefig(out_png); plt.close()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--seq_npz', default='out_relaxed_scan_long/scan_relaxed.npz')
    ap.add_argument('--par_xyz', default='out_pathopt_cmp3/traj_K_0.00.xyz')
    ap.add_argument('--out_dir', default='out_compare')
    args = ap.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)

    seq, tip_seq = read_npz_seq(args.seq_npz)
    par = read_xyz_frames(args.par_xyz)

    n = min(len(seq), len(par))
    seq = seq[:n]; par = par[:n]; tip_seq = tip_seq[:n]

    # basic stats
    diffs = par - seq
    rms = np.sqrt((diffs**2).sum(axis=2))  # shape (n, natoms)
    mean_rms = rms.mean(axis=0)
    print('Mean |par-seq| per atom:', mean_rms)
    print('Last step deltas per atom:', diffs[-1])

    plot_xz(seq, par, tip_seq, os.path.join(args.out_dir, 'xz_overlay.png'))
    plot_divergence(seq, par, os.path.join(args.out_dir, 'divergence.png'))

if __name__ == '__main__':
    main()
