#!/usr/bin/env python3

import os, argparse
import numpy as np
import matplotlib.pyplot as plt


def reconstruct_dense(neigh_idx, neigh_count, blocks, sym=True):
    natoms, max_neigh = neigh_idx.shape
    dim = natoms*3
    H = np.zeros((dim,dim), dtype=np.float64)
    for p in range(natoms):
        for j in range(int(neigh_count[p])):
            o = int(neigh_idx[p,j])
            if o < 0: continue
            H[o*3:(o+1)*3, p*3:(p+1)*3] = blocks[p,j]
    if sym: H = 0.5*(H + H.T)
    return H


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('-i','--input', required=True, help='Path to sparse Hessian .npz (contains neigh_idx, neigh_count, blocks)')
    ap.add_argument('--out', default=None, help='Output PNG path (default: <input>.png)')
    ap.add_argument('--dense', action='store_true', help='Also plot dense reconstructed Hessian (can be large)')
    ap.add_argument('--vmin', type=float, default=-12.0)
    ap.add_argument('--vmax', type=float, default=+2.0)
    args = ap.parse_args()

    d = np.load(args.input)
    neigh_idx = d['neigh_idx']
    neigh_count = d['neigh_count']
    blocks = d['blocks']
    natoms, max_neigh = neigh_idx.shape

    rect = np.zeros((natoms*3, max_neigh*3), dtype=np.float64)
    mask = np.zeros((natoms, max_neigh), dtype=np.int8)
    for p in range(natoms):
        n = int(neigh_count[p])
        if n > 0: mask[p,:n] = 1
        for j in range(n):
            rect[p*3:(p+1)*3, j*3:(j+1)*3] = blocks[p,j]

    nfig = 2 if args.dense else 1
    fig = plt.figure(figsize=(6*nfig, 6))

    ax = fig.add_subplot(1, nfig, 1)
    ax.set_title('Sparse Hessian blocks rectangle log10(|H|)')
    img = np.log10(np.abs(rect) + 1e-12)
    ax.imshow(img, cmap='viridis', origin='lower', interpolation='nearest', aspect='auto', vmin=args.vmin, vmax=args.vmax)
    ax.set_xlabel('neighbor slot (3x)')
    ax.set_ylabel('atom index (3x)')

    if args.dense:
        H = reconstruct_dense(neigh_idx, neigh_count, blocks, sym=True)
        ax2 = fig.add_subplot(1, nfig, 2)
        ax2.set_title('Dense reconstructed Hessian log10(|H|)')
        img2 = np.log10(np.abs(H) + 1e-12)
        ax2.imshow(img2, cmap='viridis', origin='lower', interpolation='nearest', vmin=args.vmin, vmax=args.vmax)
        ax2.set_xlabel('DOF')
        ax2.set_ylabel('DOF')

    fig.tight_layout()
    out = args.out
    if out is None:
        out = args.input + '.png'
    fig.savefig(out, dpi=200)
    print('saved', out)


if __name__ == '__main__':
    main()
