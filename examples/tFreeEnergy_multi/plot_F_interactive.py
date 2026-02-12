#!/usr/bin/env python3
"""
Interactive Free Energy Plotter (TI and Jarzynski)

This script creates an interactive HTML plot using Plotly for exploring
Free Energy results (Thermodynamic Integration and Jarzynski Equation)
with zoom, pan, hover tooltips, and error analysis.
"""

import sys
import numpy as np
import argparse
from scipy import integrate
import re
import os

# Constants
kB = 8.617333262145e-5  # eV/K

try:
    import plotly.graph_objects as go
    from plotly.subplots import make_subplots
except ImportError:
    print("Error: plotly is not installed. Install it with: pip install plotly")
    sys.exit(1)

def read_ti_data(filename):
    """Read free energy data from file"""
    try:
        data = np.loadtxt(filename, comments='#')
        if data.ndim == 1:
            data = data.reshape(1, -1)
        return data
    except Exception as e:
        print(f"Error reading file {filename}: {e}")
        sys.exit(1)

def detect_columns(filename):
    """
    Parse header to determine column mapping.
    """
    col_map = {'lambda': 0, 'dE_dlambda': 1, 'sigma_dE': 2, 'cumulative_FE': 3, 'cumulative_err': 4, 'distance': 5}
    header_found = False
    try:
        with open(filename, 'r') as f:
            for line in f:
                line = line.strip()
                if line.startswith('#') and 'lambda' in line:
                    parts = line.replace('#','').split()
                    # Map standard names
                    temp_map = {}
                    for i, p in enumerate(parts):
                        if p == 'lambda': temp_map['lambda'] = i
                        elif p == 'dE/dlambda' or p == 'TI_dE/dlambda': temp_map['dE_dlambda'] = i
                        elif p == 'sigma_dE' or p == 'TI_sigma': temp_map['sigma_dE'] = i
                        elif p == 'cumulative_FE' or p == 'TI_F': temp_map['cumulative_FE'] = i
                        elif p == 'cumulative_err' or p == 'TI_err': temp_map['cumulative_err'] = i
                        elif p == 'distance': temp_map['distance'] = i
                        elif p == 'JE_F': temp_map['JE_F'] = i
                        elif p == 'JE_W_avg': temp_map['JE_W_avg'] = i
                        elif p == 'JE_W_sigma': temp_map['JE_W_sigma'] = i

                    if 'lambda' in temp_map:
                        col_map.update(temp_map)
                        header_found = True
                    break
    except Exception as e:
         print(f"Warning: Could not parse header, using default columns: {e}")

    return col_map, header_found

def plot_f_interactive(filename, output_prefix=None, N_segments=None, T=300.0, b=1.198):

    # Detect columns
    col_map, header_found = detect_columns(filename)
    if header_found:
        print(f"Detected columns: {col_map}")
    else:
        print("No header detected, using default column mapping.")

    # Read data
    data = read_ti_data(filename)
    if data.shape[0] > 1000:
        print(f"Downsampling data from {data.shape[0]} to 1000 points for performance.")
        indices = np.linspace(0, data.shape[0] - 1, 1000).astype(int)
        data = data[indices]
    ncols = data.shape[1]
    
    # Safe extraction helper
    def get_col(name):
        idx = col_map.get(name)
        if idx is not None and idx < ncols:
            return data[:, idx]
        return None

    lambda_vals = get_col('lambda')
    if lambda_vals is None:
        print("Error: Could not extract lambda column.")
        sys.exit(1)

    dE_dlambda = get_col('dE_dlambda')
    errors = get_col('sigma_dE')
    cumulative_F_ti = get_col('cumulative_FE')
    cumulative_err_ti = get_col('cumulative_err')
    distances = get_col('distance')
    je_F_col = get_col('JE_F')

    # Shift data to start at zero if requested
    # TI Data
    ti_plot = cumulative_F_ti
    if ti_plot is not None and not np.isnan(ti_plot).all():
        valid_ti = ~np.isnan(ti_plot)
        if np.any(valid_ti):
            first_val = ti_plot[valid_ti][0]
            ti_plot = ti_plot - first_val
            print(f"Shifted TI curve by {-first_val:.6f} eV to start at 0.")

    # JE Data
    je_plot = je_F_col
    if je_plot is not None and not np.isnan(je_plot).all():
        valid_je = ~np.isnan(je_plot)
        if np.any(valid_je):
            first_val = je_plot[valid_je][0]
            je_plot = je_plot - first_val
            print(f"Shifted JE curve by {-first_val:.6f} eV to start at 0.")

    # Prepare reference calculations if distance and N_segments are available
    k_spring = None
    dE_ref_dlambda = None
    F_ref = None
    FE_ref = None
    
    if distances is not None and N_segments is not None:
        k_spring = kB * T / (N_segments * b**2)
        F_ref = k_spring * distances
        FE_ref_abs = 0.5 * k_spring * distances**2
        FE_ref = FE_ref_abs - FE_ref_abs[0]
        
        # Estimate dE/dlambda for reference
        # dE/dlambda = d/dlambda (0.5 * k * d^2) = k * d * dd/dlambda
        if len(lambda_vals) > 1:
            dd_dlambda = np.gradient(distances, lambda_vals)
            dE_ref_dlambda = k_spring * distances * dd_dlambda

    # Try to load Jarzynski Delta F from work file as a backup/extra info
    je_dF_work = None
    je_filename = "jarzynski_work.dat"
    if os.path.exists(je_filename):
        try:
            work_data = np.loadtxt(je_filename)
            beta = 1.0 / (kB * T)
            avg_exp = np.mean(np.exp(-beta * work_data))
            je_dF_work = -np.log(avg_exp) / beta
            print(f"Loaded Jarzynski Work data from {je_filename}. Computed dF = {je_dF_work:.6f} eV")
        except Exception as e:
            print(f"Warning: Could not read {je_filename}: {e}")

    # Create subplots
    fig = make_subplots(
        rows=5, cols=1,
        subplot_titles=(
            'dE/dλ vs λ',
            'Force vs λ',
            'Cumulative Free Energy ΔF(λ)',
            'FE Difference (data vs ref)',
            'End-to-end distance vs λ'
        ),
        vertical_spacing=0.08,
        row_heights=[0.2, 0.2, 0.2, 0.2, 0.2],
        specs=[[{"secondary_y": False}], [{"secondary_y": False}], [{"secondary_y": False}], [{"secondary_y": True}], [{"secondary_y": False}]]
    )
    plot_rows = {'dE/dλ': 1, 'Force': 2, 'FE': 3, 'FE_diff': 4, 'Distance': 5}

    # --- Plot 1: dE/dlambda vs lambda ---
    # TI Data
    if dE_dlambda is not None and not np.isnan(dE_dlambda).all():
        errs = errors if errors is not None else np.zeros_like(dE_dlambda)
        hover_text_1 = [
            f"λ = {lam:.4f}<br>dE/dλ (TI) = {de:.4f} eV<br>σ = {err:.4f} eV"
            for lam, de, err in zip(lambda_vals, dE_dlambda, errs)
        ]
        fig.add_trace(
            go.Scatter(
                x=lambda_vals, y=dE_dlambda,
                mode='lines+markers', name='dE/dλ (TI)',
                line=dict(color='blue', width=2), marker=dict(size=8, symbol='circle'),
                error_y=dict(type='data', array=errors, visible=True, color='rgba(0,0,255,0.3)'),
                hovertext=hover_text_1, hoverinfo='text'
            ),
            row=plot_rows['dE/dλ'], col=1
        )
    
    # Reference
    if dE_ref_dlambda is not None:
        hover_text_ref_de = [
            f"λ = {lam:.4f}<br>dE/dλ (ref) = {de:.4f} eV<br>k = {k_spring:.4e} eV/Å²"
            for lam, de in zip(lambda_vals, dE_ref_dlambda)
        ]
        fig.add_trace(
            go.Scatter(
                x=lambda_vals, y=dE_ref_dlambda,
                mode='lines', name='dE/dλ (ref)',
                line=dict(color='purple', width=2, dash='dash'),
                opacity=0.7, hovertext=hover_text_ref_de, hoverinfo='text'
            ),
            row=plot_rows['dE/dλ'], col=1
        )

    fig.update_yaxes(title_text="dE/dλ [eV]", row=plot_rows['dE/dλ'], col=1)

    # --- Plot 2: Force vs lambda  ---
    # Force from TI
    if dE_dlambda is not None and distances is not None:
        dd_dlambda = np.gradient(distances, lambda_vals)
        safe_dd = np.copy(dd_dlambda)
        safe_dd[np.abs(safe_dd) < 1e-9] = np.nan
        force_ti = dE_dlambda / safe_dd
        
        hover_text_force_ti = [
            f"λ = {lam:.4f}<br>Force (TI) = {f:.4f} eV/Å<br>d = {d:.3f} Å"
            for lam, f, d in zip(lambda_vals, force_ti, distances)
        ]
        fig.add_trace(
            go.Scatter(
                x=lambda_vals, y=force_ti,
                mode='lines+markers', name='Force (TI)',
                line=dict(color='orange', width=2), marker=dict(size=8, symbol='hexagram'),
                hovertext=hover_text_force_ti, hoverinfo='text'
            ),
            row=plot_rows['Force'], col=1
        )

    # Reference Force
    if F_ref is not None:
        hover_text_force_ref = [
            f"λ = {lam:.4f}<br>Force (ref) = {f:.4f} eV/Å<br>d = {d:.3f} Å<br>k = {k_spring:.4e}"
            for lam, f, d in zip(lambda_vals, F_ref, distances)
        ]
        fig.add_trace(
            go.Scatter(
                x=lambda_vals, y=F_ref,
                mode='lines', name='Force (ref spring)',
                line=dict(color='purple', width=2, dash='dash'),
                opacity=0.7, hovertext=hover_text_force_ref, hoverinfo='text'
            ),
            row=plot_rows['Force'], col=1
        )

    fig.update_yaxes(title_text="Force [eV/Å]", row=plot_rows['Force'], col=1)
    
    # --- Plot 3: Cumulative Free Energy ---
    # TI Data
    if ti_plot is not None and not np.isnan(ti_plot).all():
        errs_ti = cumulative_err_ti if cumulative_err_ti is not None else np.zeros_like(ti_plot)
        hover_text_fe_ti = [
            f"λ = {lam:.4f}<br>ΔF (TI) = {cf:.4f} eV<br>σ = {ce:.4f} eV"
            for lam, cf, ce in zip(lambda_vals, ti_plot, errs_ti)
        ]
        fig.add_trace(
            go.Scatter(
                x=lambda_vals, y=ti_plot,
                mode='lines+markers', name='ΔF(λ) (TI)',
                line=dict(color='green', width=2), marker=dict(size=8, symbol='diamond'),
                error_y=dict(type='data', array=cumulative_err_ti, visible=True, color='rgba(0,255,0,0.3)'),
                hovertext=hover_text_fe_ti, hoverinfo='text',
            ),
            row=plot_rows['FE'], col=1
        )
        
    # JE Data
    if je_plot is not None and not np.isnan(je_plot).all():
        hover_text_fe_je = [
            f"λ = {lam:.4f}<br>ΔF (JE) = {cf:.4f} eV"
            for lam, cf in zip(lambda_vals, je_plot)
        ]
        fig.add_trace(
            go.Scatter(
                x=lambda_vals, y=je_plot,
                mode='lines+markers', name='ΔF(λ) (JE)',
                line=dict(color='red', width=2, dash='dash'), marker=dict(size=6, symbol='triangle-up'),
                hovertext=hover_text_fe_je, hoverinfo='text'
            ),
            row=plot_rows['FE'], col=1
        )

    # Reference FE
    if FE_ref is not None:
        # Since data is shifted to 0, reference is already relative (FE_ref = FE_ref_abs - FE_ref_abs[0])
        FE_ref_plot = FE_ref
        
        hover_text_fe_ref = [
            f"λ = {lam:.4f}<br>ΔF (ref) = {cf:.4f} eV<br>d = {d:.3f} Å"
            for lam, cf, d in zip(lambda_vals, FE_ref_plot, distances)
        ]
        fig.add_trace(
            go.Scatter(
                x=lambda_vals, y=FE_ref_plot,
                mode='lines', name='ΔF(λ) (ref spring)',
                line=dict(color='purple', width=2, dash='solid'),
                opacity=0.7, hovertext=hover_text_fe_ref, hoverinfo='text'
            ),
            row=plot_rows['FE'], col=1
        )

    # Annotation for Jarzynski work file
    if je_dF_work is not None:
        fig.add_annotation(
            x=0.5, y=0.05, xref="paper", yref="paper",
            xanchor='center', yanchor='bottom',
            text=f"Jarzynski (work file) ΔF: {je_dF_work:.6f} eV",
            showarrow=False,
            font=dict(size=12, color="red"),
            bgcolor="rgba(255, 255, 255, 0.8)",
            row=plot_rows['FE'], col=1
        )

    # Final Delta F results comparison
    final_text = []
    if ti_plot is not None and len(ti_plot) > 0: 
        final_text.append(f"TI: {ti_plot[-1]:.4f} eV")
    
    if je_plot is not None:
        valid_je_data = je_plot[~np.isnan(je_plot)]
        if len(valid_je_data) > 0:
            final_text.append(f"JE: {valid_je_data[-1]:.4f} eV")
    
    if FE_ref is not None and len(FE_ref) > 0: 
        final_text.append(f"Ref: {FE_ref[-1]:.4f} eV")
    
    if final_text:
        fig.add_annotation(
            x=0.98, y=0.02, xref="paper", yref="paper",
            xanchor='right', yanchor='bottom',
            text="Final ΔF: " + " | ".join(final_text),
            showarrow=False,
            font=dict(size=12, color="black"),
            bgcolor="rgba(255, 255, 255, 0.9)",
            row=plot_rows['FE'], col=1
        )

    fig.update_yaxes(title_text="ΔF(λ) [eV]", row=plot_rows['FE'], col=1)

    # --- Plot 4: FE Difference ---
    if FE_ref is not None:
        data_to_compare = None
        method_name = ""
        # Prefer TI for difference if available, else JE
        if ti_plot is not None and not np.isnan(ti_plot).all():
            data_to_compare = ti_plot
            method_name = "TI"
        elif je_plot is not None and not np.isnan(je_plot).all():
            data_to_compare = je_plot
            method_name = "JE"
        
        if data_to_compare is not None:
            # Data and ref are both zero-based now
            FE_ref_aligned = FE_ref 
            
            # Avoid division by zero
            safe_ref = np.copy(FE_ref_aligned)
            safe_ref[np.abs(safe_ref) < 1e-9] = np.nan # Use nan for percentage at zero
            
            fe_diff_abs = data_to_compare - FE_ref_aligned
            fe_diff_percent = (fe_diff_abs / safe_ref) * 100

            hover_text_diff_abs = [
                f"λ = {lam:.4f}<br>Abs Diff ({method_name}-ref) = {d:.4f} eV"
                for lam, d in zip(lambda_vals, fe_diff_abs)
            ]
            hover_text_diff_pct = [
                f"λ = {lam:.4f}<br>% Diff ({method_name}/ref) = {d:.2f}%"
                for lam, d in zip(lambda_vals, fe_diff_percent)
            ]

            # Percentage difference
            fig.add_trace(
                go.Scatter(
                    x=lambda_vals, y=fe_diff_percent,
                    mode='lines+markers', name=f'% Diff ({method_name}/ref)',
                    line=dict(color='cyan', width=2), marker=dict(size=8, symbol='cross'),
                    hovertext=hover_text_diff_pct, hoverinfo='text'
                ),
                row=plot_rows['FE_diff'], col=1,
                secondary_y=False
            )
            # Absolute difference
            fig.add_trace(
                go.Scatter(
                    x=lambda_vals, y=fe_diff_abs,
                    mode='lines+markers', name=f'Abs Diff ({method_name}-ref)',
                    line=dict(color='magenta', width=2, dash='dot'), marker=dict(size=8, symbol='x'),
                    hovertext=hover_text_diff_abs, hoverinfo='text'
                ),
                row=plot_rows['FE_diff'], col=1,
                secondary_y=True
            )

    fig.update_yaxes(title_text="% Difference", row=plot_rows['FE_diff'], col=1, secondary_y=False)
    fig.update_yaxes(title_text="Abs Difference [eV]", row=plot_rows['FE_diff'], col=1, secondary_y=True)

    # --- Plot 5: Distance ---
    if distances is not None:
        hover_text_dist = [
            f"λ = {lam:.4f}<br>Distance = {dist:.3f} Å"
            for lam, dist in zip(lambda_vals, distances)
        ]
        fig.add_trace(
            go.Scatter(
                x=lambda_vals, y=distances,
                mode='lines+markers', name='End-to-end Distance',
                line=dict(color='red', width=2), marker=dict(size=8, symbol='square'),
                hovertext=hover_text_dist, hoverinfo='text'
            ),
            row=plot_rows['Distance'], col=1
        )
        fig.update_yaxes(title_text="End-to-end Distance [Å]", row=plot_rows['Distance'], col=1)

    # Update all x-axes
    for i in range(1, 6):
        fig.update_xaxes(title_text="λ", row=i, col=1)

    # Layout
    total_height = 1600
    fig.update_layout(
        height=total_height,
        title=f"Free Energy Analysis: {os.path.basename(filename)}",
        showlegend=True,
        hovermode='closest',
        template='plotly_white'
    )

    # Output
    if output_prefix is None:
        output_prefix = filename.replace('.dat', '')
    output_file = f"{output_prefix}_F_interactive.html"
    fig.write_html(output_file)
    print(f"Interactive plot saved to: {output_file}")


def main():
    parser = argparse.ArgumentParser(description="Interactive Free Energy Plotter")
    parser.add_argument("--input", type=str, default="entropic_spring_20_TI.dat",
                       help="Input data file")
    parser.add_argument("--output", type=str, default=None,
                       help="Output prefix for plots")
    parser.add_argument("--N", type=int, default=None,
                       help="Number of segments (for reference). Guessed from filename if possible.")
    parser.add_argument("--T", type=float, default=300.0,
                       help="Temperature in Kelvin")
    parser.add_argument("--b", type=float, default=1.198,
                       help="Segment length in Angstroms")

    args = parser.parse_args()

    if args.N is None:
        match = re.search(r'_(\d+)_', args.input)
        if match:
             args.N = int(match.group(1))
             print(f"Guessed N={args.N} from filename.")

    # N_segments = N_atoms - 1 for a chain
    N_seg = args.N - 1 if args.N is not None else None
    
    plot_f_interactive(args.input, args.output, N_segments=N_seg, T=args.T, b=args.b)

if __name__ == "__main__":
    main()
