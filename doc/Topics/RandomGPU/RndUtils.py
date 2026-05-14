"""
RndUtils.py  –  Statistical analysis and plotting for normal-distribution RNG output.

Functions:
    moments(x)                -> dict
    chi_square(x, bins)       -> (chi2, p)
    autocorr(x, max_lag)      -> np.ndarray
    welch_psd(x, nperseg)     -> (freqs, psd)
    analyse(results)           -> prints a summary table
    plot_all(results, title, save_path)  -> combined figure
    plot_histograms(results, ...)
    plot_autocorr(results, ...)
    plot_psd(results, ...)
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy import stats as sp_stats
from scipy.signal import welch


# ---------------------------------------------------------------------------
# Statistical measures
# ---------------------------------------------------------------------------

def moments(x):
    """Return dict with mean, var, skew, kurtosis (excess)."""
    return {
        "mean": float(np.mean(x)),
        "var":  float(np.var(x)),
        "skew": float(sp_stats.skew(x)),
        "kurt": float(sp_stats.kurtosis(x)),
    }

def chi_square(x, bins=200, range_sigma=5.0):
    """
    Chi-square goodness-of-fit against N(0,1).
    Bins are placed symmetrically over ±range_sigma.
    Returns (chi2_stat, p_value).
    """
    lo, hi = -range_sigma, range_sigma
    hist, edges = np.histogram(x, bins=bins, range=(lo, hi), density=False)
    bin_probs   = np.diff(sp_stats.norm.cdf(edges))
    expected    = hist.sum() * bin_probs / bin_probs.sum()
    expected    = np.maximum(expected, 1e-10)
    chi2, p     = sp_stats.chisquare(hist, expected)
    return float(chi2), float(p)

def autocorr(x, max_lag=200):
    """Normalized autocorrelation at lags 0..max_lag."""
    x   = x[:max(10000, max_lag * 10)] - np.mean(x)
    ac  = np.correlate(x, x, mode='full')
    ac  = ac[len(ac) // 2:]
    return ac[:max_lag + 1] / ac[0]

def welch_psd(x, nperseg=1024):
    """Welch power spectral density estimate."""
    return welch(x[:min(len(x), 1_000_000)], nperseg=nperseg)

def lag1_autocorr(x):
    """Lag-1 autocorrelation (quick scalar)."""
    x0 = x[:-1] - np.mean(x)
    x1 = x[1:]  - np.mean(x)
    denom = np.sum(x0**2)
    return float(np.sum(x0 * x1) / denom) if denom > 0 else float('nan')


# ---------------------------------------------------------------------------
# Analysis table
# ---------------------------------------------------------------------------

def analyse(results, timings=None, N=None):
    """
    Print a summary table of statistical quality for each method.

    Parameters
    ----------
    results  : dict  name -> np.ndarray
    timings  : dict  name -> ms  (optional)
    N        : sample count used for throughput calculation
    """
    header = f"{'Method':<24} {'mean':>8} {'var':>8} {'skew':>8} {'kurt':>8} {'chi2-p':>10} {'lag1-ac':>9}"
    if timings:
        header += f"  {'ms':>7}  {'Msamp/s':>9}"
    print("\n" + "="*len(header))
    print(header)
    print("-"*len(header))

    for name, data in results.items():
        m    = moments(data)
        _, p = chi_square(data)
        ac1  = lag1_autocorr(data)
        row  = (f"{name:<24} {m['mean']:+8.4f} {m['var']:8.4f} "
                f"{m['skew']:+8.4f} {m['kurt']:+8.4f} {p:10.3e} {ac1:+9.5f}")
        if timings and name in timings:
            t = timings[name]
            n = N or len(data)
            row += f"  {t:7.2f}  {n/t*1e-6:9.1f}"
        print(row)
    print("="*len(header))


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------

def plot_histograms(results, bins=100, range_sigma=4.5, ncols=None,
                    title="Histograms vs N(0,1)", save_path=None, show=True):
    names  = list(results.keys())
    ncols  = ncols or len(names)
    fig, axes = plt.subplots(1, ncols, figsize=(4 * ncols, 4), squeeze=False)
    axes = axes[0]
    xref = np.linspace(-range_sigma, range_sigma, 400)
    pdf  = sp_stats.norm.pdf(xref)
    for ax, name in zip(axes, names):
        data = results[name]
        ax.hist(data, bins=bins, range=(-range_sigma, range_sigma),
                density=True, alpha=0.6, color='steelblue', label='generated')
        ax.plot(xref, pdf, 'r--', lw=1.5, label='N(0,1)')
        ax.set_title(name, fontsize=9)
        ax.set_xlim(-range_sigma, range_sigma)
        ax.legend(fontsize=7)
    fig.suptitle(title)
    fig.tight_layout()
    if save_path: plt.savefig(save_path, dpi=120)
    if show: plt.show()
    return fig


def plot_autocorr(results, max_lag=200, ncols=None,
                  title="Autocorrelation", save_path=None, show=True):
    names  = list(results.keys())
    ncols  = ncols or len(names)
    fig, axes = plt.subplots(1, ncols, figsize=(4 * ncols, 3.5), squeeze=False)
    axes = axes[0]
    lags = np.arange(max_lag + 1)
    for ax, name in zip(axes, names):
        ac = autocorr(results[name], max_lag)
        ax.plot(lags[1:], ac[1:], lw=0.8, color='steelblue')
        ax.axhline(0,  color='r', lw=0.8, ls='--')
        ax.axhline( 1.96 / np.sqrt(len(results[name])), color='gray', lw=0.7, ls=':')
        ax.axhline(-1.96 / np.sqrt(len(results[name])), color='gray', lw=0.7, ls=':')
        ax.set_title(name, fontsize=9)
        ax.set_xlabel("lag")
        ax.set_ylim(-0.05, 0.15)
    fig.suptitle(title)
    fig.tight_layout()
    if save_path: plt.savefig(save_path, dpi=120)
    if show: plt.show()
    return fig


def plot_psd(results, nperseg=1024, ncols=None,
             title="Power Spectral Density", save_path=None, show=True):
    names  = list(results.keys())
    ncols  = ncols or len(names)
    fig, axes = plt.subplots(1, ncols, figsize=(4 * ncols, 3.5), squeeze=False)
    axes = axes[0]
    for ax, name in zip(axes, names):
        freqs, psd = welch_psd(results[name], nperseg=nperseg)
        ax.semilogy(freqs, psd, lw=0.8, color='steelblue')
        ax.set_title(name, fontsize=9)
        ax.set_xlabel("normalized frequency")
        ax.set_ylabel("PSD")
        ax.set_xlim(0, 0.5)
    fig.suptitle(title)
    fig.tight_layout()
    if save_path: plt.savefig(save_path, dpi=120)
    if show: plt.show()
    return fig


def plot_all(results, timings=None, N=None, prefix="rng",
             save=True, show=True):
    """
    Produce all three plot types and print analysis table.
    If save=True, writes <prefix>_hist.png, <prefix>_acf.png, <prefix>_psd.png.
    """
    analyse(results, timings=timings, N=N)
    plot_histograms(results,
        save_path=f"{prefix}_hist.png" if save else None, show=show)
    plot_autocorr(results,
        save_path=f"{prefix}_acf.png"  if save else None, show=show)
    plot_psd(results,
        save_path=f"{prefix}_psd.png"  if save else None, show=show)


# ---------------------------------------------------------------------------
# CLI demo: load results from RndNormal and plot
# ---------------------------------------------------------------------------
if __name__ == "__main__":
    import sys
    sys.path.insert(0, ".")
    from RndNormal import RndNormal

    rng = RndNormal(N=2**20, buffer_size=2**20, lds_size=256)
    results, timings = rng.run_all(seed=42, n_warmup=3)
    plot_all(results, timings=timings, N=rng.N, prefix="rng", save=True, show=True)
