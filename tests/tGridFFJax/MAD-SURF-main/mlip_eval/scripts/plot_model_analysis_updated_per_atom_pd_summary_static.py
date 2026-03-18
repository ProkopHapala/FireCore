import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

def plot_comparisons(filepath, save_prefix="comparison",
                     x_limits=None, y_limits=(0, 2),
                     hist_RMSD_x_limits=(-0.1, 1.5), hist_RMSD_y_limits=None,
                     hist_fmae_x_limits=(-0.005, 0.005), hist_fmae_y_limits=None,
                     hist_x_limits=(-0.03, 0.03), hist_y_limits=None,
                     parity_limits=(2, 4)):

    # Load data
    df = pd.read_csv(filepath, sep="\t|,", engine="python")

    # --- Handle MACE_MPA-0 foundational model ---
    df.loc[df["model"] == "MACE_MPA-0_foundational_model", "energy_diff_per_atom"] = 0.0

    # --- Metrics for bar plot ---
    metrics_df = df.groupby("model").agg({
        "energy_diff_per_atom": lambda x: x.abs().mean(),
        "force_mae_static": lambda x: x.abs().mean(),
        "ads_height_diff": lambda x: x.abs().mean(),
        "rmsd_per_element": "mean"
    })

    # --- Consistent model ordering (by average RMSD, ascending) ---
    model_order = metrics_df.sort_values("rmsd_per_element", ascending=True).index.tolist()
    metrics_df = metrics_df.loc[model_order]

    # --- Consistent colors across plots ---
    colors = plt.cm.tab20(np.linspace(0, 1, len(model_order)))
    color_map = {model: colors[i] for i, model in enumerate(model_order)}

    # SUMMARY FIGURE 
    fig, axes = plt.subplots(1, 2, figsize=(24, 12))
    ax1, ax4 = axes.flat
    
    # --- Panel labels ---
    ax1.text(-0.05, 1.1, "a)", transform=ax1.transAxes, fontsize=28, fontweight='bold', va='top', ha='right')
    ax4.text(-0.05, 1.1, "b)", transform=ax4.transAxes, fontsize=28, fontweight='bold', va='top', ha='right')

    # --- Customize the axes boxes (spines) for all subplots ---
    for ax in axes.flat:
        for spine in ax.spines.values():
            spine.set_linewidth(3)       # thicker lines
    
    # --- Bar plot ---
    metrics_df.plot(kind="bar", ax=ax1)
    ax1.set_ylabel("Error value (eV/atom, eV/Å², Å, Å)", fontsize=20)
    if y_limits:
        ax1.set_ylim(y_limits)
    ax1.tick_params(axis="x", rotation=90)
    leg1 = ax1.legend(fontsize=20)
    leg1.get_frame().set_linewidth(2)
    ax1.tick_params(axis='x', labelsize=20)
    ax1.tick_params(axis='y', labelsize=16)

    # Color x-axis tick labels to match model colors
    for label in ax1.get_xticklabels():
        model_name = label.get_text()
        if model_name in color_map:
            label.set_color(color_map[model_name])
    
    ax1.xaxis.label.set_visible(False)

    # --- Parity plot (same model order and colors) ---
    for model in model_order:
        subset = df[df["model"] == model]
        x_vals = subset["ads_height_reference"].values
        y_vals = subset["ads_height_relaxed"].values

        mask = (x_vals >= 0) & (x_vals <= 5) & (y_vals >= 0) & (y_vals <= 5)
        x_fit = x_vals[mask]
        y_fit = y_vals[mask]

        if len(x_fit) > 1:
            slope, intercept = np.polyfit(x_fit, y_fit, 1)
            y_pred = slope * x_fit + intercept
            r2 = 1 - np.sum((y_fit - y_pred) ** 2) / np.sum((y_fit - np.mean(y_fit)) ** 2)
            ax4.scatter(x_vals, y_vals, s=150, alpha=0.5,
                        label=f"{model} (R²={r2:.2f})",
                        color=color_map[model])
            ax4.plot(x_fit, y_pred, color=color_map[model], linestyle="--", alpha=0.7)
        else:
            ax4.scatter(x_vals, y_vals, s=150, alpha=0.8,
                        label=model, color=color_map[model])

    lims = ax4.get_xlim()
    ax4.plot(lims, lims, "k--")
    if parity_limits:
        ax4.set_xlim(parity_limits)
        ax4.set_ylim(parity_limits)
    ax4.set_xlabel("Reference (DFT) Ads. Height (Å)", fontsize=20)
    ax4.set_ylabel("Relaxed (MACE) Ads. Height (Å)", fontsize=20)
    #ax4.legend(fontsize=16)
    leg2 = ax4.legend(
    fontsize=20,
    loc='upper center',         # position relative to the bbox
    bbox_to_anchor=(0.5, -0.15), # x=0.5 (center), y=-0.15 (below plot)
    ncol=1                      # number of columns in legend
)
    leg2.get_frame().set_linewidth(2)
    ax4.tick_params(axis='x', labelsize=16)
    ax4.tick_params(axis='y', labelsize=16)

    plt.tight_layout()
    plt.savefig(f"{save_prefix}_summary_only_bar_parity.png")
    plt.close()

if __name__ == "__main__":
    plot_comparisons("comparison_results_static.csv")

