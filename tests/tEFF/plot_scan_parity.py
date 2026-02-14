#!/usr/bin/env python3
import argparse
import os
import numpy as np
import matplotlib.pyplot as plt

LABELS = ["Etot", "Ek", "Eee", "Eae", "Eaa"]


def load_arrays(cpu_path, gpu_path):
    cpu = np.load(cpu_path)
    gpu = np.load(gpu_path)
    n = min(len(cpu), len(gpu))
    if n == 0:
        raise ValueError("Empty arrays; nothing to plot")
    return cpu[:n], gpu[:n]


def plot_scan(cpu, gpu, out_path, title=None):
    diff = cpu - gpu
    n = cpu.shape[0]
    x = np.arange(n)

    fig, (ax_top, ax_bot) = plt.subplots(2, 1, figsize=(10, 7), sharex=True, gridspec_kw={"height_ratios": [2, 1]})

    # Top: Etot curves (same color, different linestyle/width)
    ax_top.plot(x, cpu[:, 0], label="CPU Etot", color="#1f77b4", ls=":", lw=1.5)
    ax_top.plot(x, gpu[:, 0], label="GPU Etot", color="#1f77b4", ls="-", lw=0.5)
    ax_top.set_ylabel("Energy [eV]")
    ax_top.legend(frameon=False)
    if title:
        ax_top.set_title(title)
    ax_top.grid(alpha=0.2)

    # Bottom: component diffs
    colors = ["#2ca02c", "#d62728", "#9467bd", "#8c564b", "#e377c2"]
    for i in range(5):
        ax_bot.plot(x, diff[:, i], label=f"d{LABELS[i]}", color=colors[i % len(colors)])
    ax_bot.axhline(0.0, color="k", lw=0.8, alpha=0.6)
    ax_bot.set_xlabel("Configuration index")
    ax_bot.set_ylabel("CPU - GPU [eV]")
    ax_bot.legend(frameon=False, ncol=3)
    ax_bot.grid(alpha=0.2)

    fig.tight_layout()
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    fig.savefig(out_path, dpi=180)
    plt.close(fig)

    diff_max = np.max(np.abs(diff), axis=0)
    return {
        "nconf": int(n),
        "Etot_max": float(diff_max[0]),
        "Ek_max": float(diff_max[1]),
        "Eee_max": float(diff_max[2]),
        "Eae_max": float(diff_max[3]),
        "Eaa_max": float(diff_max[4]),
    }


def main():
    ap = argparse.ArgumentParser(description="Plot CPU vs GPU scan parity (energies) and diffs.")
    ap.add_argument("--cpu", required=True, help="Path to Es5_cpu.npy")
    ap.add_argument("--gpu", required=True, help="Path to Es5_gpu.npy")
    ap.add_argument("--out", required=True, help="Output plot path (png/pdf)")
    ap.add_argument("--title", default=None, help="Optional title")
    args = ap.parse_args()

    cpu, gpu = load_arrays(args.cpu, args.gpu)
    stats = plot_scan(cpu, gpu, args.out, title=args.title)
    print("Saved", args.out)
    print("Stats:", stats)


if __name__ == "__main__":
    main()
