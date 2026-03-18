#!/usr/bin/env python3
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def plot_training_strategies(filepath, save_name="figure2_training_strategies.pdf"):

    # ------------------------------------------------------------------
    # Load data
    # ------------------------------------------------------------------
    df = pd.read_csv(filepath, sep=r"\t|,", engine="python")

    # Set foundational model energy diff to zero
    df.loc[df["model"] == "MACE_MPA-0_foundational_model", "energy_diff_per_atom"] = 0.0

    # ------------------------------------------------------------------
    # Model ordering and labels
    # ------------------------------------------------------------------
    model_order = [
        
        # Full dataset, stage-one
        "full_lambda1_stageone",
        "full_lambda10_stageone_scaleshift",
        # Small subsets
        "small_lambda1_stageone",
        "small_lambda10_stageone",
        # Descriptor-filtered subsets
        "df_0.01_l10",
        "df_0.08_l10",
        # Pretrained + finetuning
        "MACE_MPA-0_foundational_model",
        "MACE_MPA-0_foundational_on_train",
        "MACE_MPA-0_foundational_on_test",
    ]

    label_map = {
        "MACE_MPA-0_foundational_model":      "MACE-MPA-0 (pretrained)",
        "MACE_MPA-0_foundational_on_train":   "MACE-MPA-0 (finetuned)",
        "MACE_MPA-0_foundational_on_test":    "MACE-MPA-0 (lifelong FT)",
        "full_lambda1_stageone":              "Full data λ = 1",
        "full_lambda10_stageone_scaleshift":  "Full data λ = 10 (scaleshift)",
        "small_lambda1_stageone":             "Small subset λ = 1",
        "small_lambda10_stageone":            "Small subset λ = 10",
        "df_0.01_l10":                        "Filtered 1% λ = 10",
        "df_0.08_l10":                        "Filtered 8% λ = 10",
    }

    # ------------------------------------------------------------------
    # Aggregate metrics
    # ------------------------------------------------------------------
    metrics = (
        df.groupby("model")
        .agg(
            energy_mae=("energy_diff_per_atom_static", lambda x: x.abs().mean()),
            force_mae=("force_mae_static", lambda x: x.abs().mean()),
            rmsd=("rmsd_per_element", "mean"),
        )
        .reindex(model_order)
    )

    # ------------------------------------------------------------------
    # Color palettes
    # ------------------------------------------------------------------
    soft_colors = {
        "energy": "#8CBED6",  # soft blue
        "force":  "#A8D5BA",  # soft green
        "rmsd":   "#F4C8A4",  # soft orange
    }

    # Select distinct colors for parity plots
    parity_color_map = {
        "full_lambda1_stageone":              "#1f77b4",  # blue
        "full_lambda10_stageone_scaleshift":  "#ff7f0e",  # orange
        "small_lambda1_stageone":             "#2ca02c",  # green
        "small_lambda10_stageone":            "#d62728",  # red
        "df_0.01_l10":                        "#9467bd",  # purple
        "df_0.08_l10":                        "#8c564b",  # brown
        "MACE_MPA-0_foundational_model":      "#e377c2",  # pink
        "MACE_MPA-0_foundational_on_train":   "#7f7f7f",  # gray
        "MACE_MPA-0_foundational_on_test":    "#bcbd22",  # olive
    }   

    # ------------------------------------------------------------------
    # Plot layout
    # ------------------------------------------------------------------
    plt.rcParams.update({
        "font.size": 17,
        "axes.labelsize": 17,
        "xtick.labelsize": 17,
        "ytick.labelsize": 17,
    })

    fig, axes = plt.subplots(2, 3, figsize=(22, 12))
    ax_energy, ax_force, ax_rmsd, ax_D, ax_E, ax_F = axes.flat

    y = np.arange(len(model_order))
    ylabels = [label_map[m] for m in model_order]

    # ========================= FIRST ROW ===============================

    # 1) Energy MAE per element
    ax_energy.barh(y, metrics["energy_mae"]*1000, color=soft_colors["energy"])
    ax_energy.set_yticks(y)
    ax_energy.set_yticklabels(ylabels)
    ax_energy.invert_yaxis()
    ax_energy.set_xlabel("Energy MAE per element (meV)")
    ax_energy.set_xlim(0, 30)
    ax_energy.tick_params(axis="y", labelsize=17)

    # 2) Force MAE
    ax_force.barh(y, metrics["force_mae"]*1000, color=soft_colors["force"])
    ax_force.set_yticks(y)
    ax_force.set_yticklabels([])
    ax_force.invert_yaxis()
    ax_force.set_xlabel("Force MAE (meV/Å)")
    ax_force.set_xlim(0, 200)

    # 3) RMSD
    ax_rmsd.barh(y, metrics["rmsd"], color=soft_colors["rmsd"])
    ax_rmsd.set_yticks(y)
    ax_rmsd.set_yticklabels([])
    ax_rmsd.invert_yaxis()
    ax_rmsd.set_xlabel("RMSD (Å)")
    ax_rmsd.set_xlim(0, 0.6)

    # ==================================================================
    # HELPER: parity plotting for groups
    # ==================================================================
    def plot_parity_group(ax, model_list):

        for model_key in model_list:
            sub = df[df["model"] == model_key]
            xvals = sub["ads_height_reference"].values
            yvals = sub["ads_height_relaxed"].values

            # Fit in 2–4 Å only
            mask = (
                (xvals >= 2.0) & (xvals <= 4.0) &
                (yvals >= 2.0) & (yvals <= 4.0)
            )
            x_fit = xvals[mask]
            y_fit = yvals[mask]

            if len(x_fit) >= 2:
                slope, intercept = np.polyfit(x_fit, y_fit, 1)
                y_pred = slope * x_fit + intercept
                r2 = 1.0 - np.sum((y_fit - y_pred)**2) / np.sum((y_fit - np.mean(y_fit))**2)

                # Regression line
                order = np.argsort(x_fit)
                ax.plot(
                    x_fit[order],
                    y_pred[order],
                    color=parity_color_map[model_key],
                    linewidth=2,
                )

                label = f"{label_map[model_key].replace(chr(10),' ')} (R²={r2:.2f})"

            else:
                label = label_map[model_key].replace("\n", " ")

            # Scatter points
            ax.scatter(xvals, yvals, s=60, alpha=0.7,
                       color=parity_color_map[model_key],
                       label=label)

        # Layout
        ax.plot([2, 4], [2, 4], "k--", linewidth=1.2)
        ax.set_xlim(2, 4)
        ax.set_ylim(2, 4)
        ax.set_xlabel("DFT ads. height (Å)", fontsize=17)
        ax.set_ylabel("MAD-SURF ads. height (Å)", fontsize=17)

    # ========================= SECOND ROW ===============================

    

    # Panel D — Full data
    full_data_models = [
        "full_lambda1_stageone",
        "full_lambda10_stageone_scaleshift",
    ]
    plot_parity_group(ax_D, full_data_models)
    ax_D.legend(loc="upper center", bbox_to_anchor=(0.5, -0.17), ncol=1)


    # Panel E — Small + filtered
    subset_models = [
        "small_lambda1_stageone",
        "small_lambda10_stageone",
        "df_0.01_l10",
        "df_0.08_l10",
    ]

    plot_parity_group(ax_E, subset_models)
    ax_E.legend(loc="upper center", bbox_to_anchor=(0.5, -0.17), ncol=1)

    # Panel F — Pretrained family
    pretrained_models = [
        "MACE_MPA-0_foundational_model",
        "MACE_MPA-0_foundational_on_train",
        "MACE_MPA-0_foundational_on_test",
    ]
    plot_parity_group(ax_F, pretrained_models)
    ax_F.legend(loc="upper center", bbox_to_anchor=(0.5, -0.17), ncol=1)

    

    # ------------------------------------------------------------------
    # Cosmetic: thicker spines
    # ------------------------------------------------------------------
    for ax in axes.flat:
        for spine in ax.spines.values():
            spine.set_linewidth(1.5)

    plt.tight_layout()

    outpath = os.path.join(".", save_name)
    plt.savefig(outpath, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved figure to: {outpath}")


if __name__ == "__main__":
    plot_training_strategies("relaxed_MACE_results_final_comparison/comparison_results_static.csv")
