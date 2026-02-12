#!/usr/bin/env python3
"""
Interactive Thermodynamic Integration Plotter

This script creates an interactive HTML plot using Plotly for exploring
TI results with zoom, pan, hover tooltips, and error analysis.
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
    """Read thermodynamic integration data from file"""
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
    Returns a dictionary mapping 'lambda', 'dE_dlambda', 'sigma_dE', 'cumulative_FE', 'cumulative_err', 'distance' to indices.
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

def plot_ti_interactive(filename, output_prefix=None, N_segments=None, T=300.0, b=1.198):

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
    dE_dlambda = get_col('dE_dlambda')
    
    if lambda_vals is None:
        print("Error: Could not extract lambda columns.")
        sys.exit(1)
        
    errors = get_col('sigma_dE')
    cumulative_F_read = get_col('cumulative_FE')
    cumulative_err_read = get_col('cumulative_err')

    cumulative_F = cumulative_F_read
    cumulative_error = cumulative_err_read

    distances = get_col('distance')
    
    Force_from_data = None
    dE_ref_dlambda = None
    F_ref = None

    # JE columns
    je_F_col = get_col('JE_F')

    je_F_col_to_plot = None
    initial_je_F_val = 0.0
    je_final_dF = None # This will be the Delta F from the start of the curve

    if je_F_col is not None:
        valid_je_F_indices = np.where(~np.isnan(je_F_col))[0]
        if len(valid_je_F_indices) > 0:
            initial_je_F_val = je_F_col[valid_je_F_indices[0]]
            je_F_col_to_plot = je_F_col - initial_je_F_val
            je_final_dF = je_F_col_to_plot[valid_je_F_indices[-1]] # Final value of the shifted curve
            print(f"Loaded Jarzynski Free Energy from file. Final dF (curve shifted to zero) = {je_final_dF:.6f} eV")
        # If no valid data, je_F_col_to_plot remains None


    je_dF = je_final_dF # Use this for annotation if curve data is primary
    
    if je_dF is None: # Only try to load from work.dat if not already set from curve
        je_filename = "jarzynski_work.dat"
        if os.path.exists(je_filename):
            try:
                work_data = np.loadtxt(je_filename)
                beta = 1.0 / (kB * T)
                avg_exp = np.mean(np.exp(-beta * work_data))
                je_dF_from_work_file = -np.log(avg_exp) / beta
                je_dF = je_dF_from_work_file
                print(f"Loaded Jarzynski Work data from {je_filename}. Computed dF = {je_dF:.6f} eV")
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
    if dE_dlambda is not None and not np.isnan(dE_dlambda).all():
        hover_text_1 = [
            f"λ = {lam:.4f}<br>dE/dλ = {de:.4f} eV<br>σ(dE/dλ) = {err:.4f} eV"
            for lam, de, err in zip(lambda_vals, dE_dlambda, errors)
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
        if dE_ref_dlambda is not None:
            fig.add_trace(
                go.Scatter(
                    x=lambda_vals, y=dE_ref_dlambda,
                    mode='lines', name='dE/dλ (ref)',
                    line=dict(color='purple', width=2, dash='dash'),
                    opacity=0.7, hoverinfo='text'
                ),
                row=plot_rows['dE/dλ'], col=1
            )
        
        # Linear fit for dE/dlambda
        valid = ~np.isnan(dE_dlambda)
        if np.any(valid):
            p_dE = np.polyfit(lambda_vals[valid], dE_dlambda[valid], 1)
            fit_text_dE = f"Fit: {p_dE[0]:.4f}λ + {p_dE[1]:.4f}"
            fig.add_annotation(
                x=0.01, y=0.99, xref="paper", yref="paper",
                xanchor='left', yanchor='top',
                text=fit_text_dE, showarrow=False,
                font=dict(size=12, color="blue"),
                bgcolor="rgba(255, 255, 255, 0.8)",
                row=plot_rows['dE/dλ'], col=1
            )
    fig.update_yaxes(title_text="dE/dλ [eV]", row=plot_rows['dE/dλ'], col=1)

    # --- Plot 2: Force vs lambda  ---
    if 'Force' in plot_rows and Force_from_data is not None and not np.isnan(Force_from_data).all():
        hover_text_force = [
            f"λ = {lam:.4f}<br>Force = {f:.4f} eV/Å"
            for lam, f in zip(lambda_vals, Force_from_data)
        ]
        fig.add_trace(
            go.Scatter(
                x=lambda_vals, y=Force_from_data,
                mode='lines+markers', name='Force (from data)',
                line=dict(color='orange', width=2), marker=dict(size=8, symbol='hexagram'),
                hovertext=hover_text_force, hoverinfo='text'
            ),
            row=plot_rows['Force'], col=1
        )
        if F_ref is not None:
            fig.add_trace(
                go.Scatter(
                    x=lambda_vals, y=F_ref,
                    mode='lines', name='Force (ref spring)',
                    line=dict(color='purple', width=2, dash='dash'),
                    opacity=0.7, hoverinfo='text'
                ),
                row=plot_rows['Force'], col=1
            )
        
        # Linear fit for Force
        valid = ~np.isnan(Force_from_data)
        if np.any(valid):
            p_force = np.polyfit(lambda_vals[valid], Force_from_data[valid], 1)
            fit_text_force = f"Fit: {p_force[0]:.4f}λ + {p_force[1]:.4f}"
            fig.add_annotation(
                x=0.01, y=0.99, xref="paper", yref="paper",
                xanchor='left', yanchor='top',
                text=fit_text_force, showarrow=False,
                font=dict(size=12, color="orange"),
                bgcolor="rgba(255, 255, 255, 0.8)",
                row=plot_rows['Force'], col=1
            )
    fig.update_yaxes(title_text="Force [eV/Å]", row=plot_rows['Force'], col=1)
    
    # --- Plot 3: Cumulative Free Energy ---
    if cumulative_F is not None and not np.isnan(cumulative_F).all():
        hover_text_fe = [
            f"λ = {lam:.4f}<br>ΔF(λ) = {cf:.4f} eV<br>σ(ΔF) = {ce:.4f} eV"
            for lam, cf, ce in zip(lambda_vals, cumulative_F, cumulative_error)
        ]
        fig.add_trace(
            go.Scatter(
                x=lambda_vals, y=cumulative_F,
                mode='lines+markers', name='ΔF(λ) (TI)',
                line=dict(color='green', width=2), marker=dict(size=8, symbol='diamond'),
                error_y=dict(type='data', array=cumulative_error, visible=True, color='rgba(0,255,0,0.3)'),
                hovertext=hover_text_fe, hoverinfo='text',
            ),
            row=plot_rows['FE'], col=1
        )
        
    if je_F_col is not None and not np.isnan(je_F_col).all():
        fig.add_trace(
            go.Scatter(
                x=lambda_vals, y=je_F_col,
                mode='lines+markers', name='ΔF(λ) (JE)',
                line=dict(color='orange', width=2, dash='dash'), marker=dict(size=6, symbol='triangle-up'),
            ),
            row=plot_rows['FE'], col=1
        )

    FE_ref_shifted = None 
    FE_ref_smooth_plot_data = None  # Keep this initialized to None for the FE Difference plot
    if distances is not None and N_segments is not None and cumulative_F is not None:
        k_spring = kB * T / (N_segments * b**2)
        
        # Calculate FE_ref_shifted for the difference plot (using original lambda_vals points)
        FE_ref_original_points = 0.5 * k_spring * distances**2
        FE_ref_shifted = FE_ref_original_points - FE_ref_original_points[0] + cumulative_F[0]

        # Generate denser lambda values for plotting the smooth reference
        lambda_plot_vals = np.linspace(lambda_vals[0], lambda_vals[-1], 200)

        # Assume distance changes linearly with lambda for the reference parabola
        initial_distance = distances[0]
        final_distance = distances[-1]
        distances_smooth = initial_distance + (final_distance - initial_distance) * lambda_plot_vals

        # Calculate the smooth reference free energy for plotting
        FE_ref_smooth_plot = 0.5 * k_spring * distances_smooth**2
        FE_ref_smooth_plot_data = FE_ref_smooth_plot - FE_ref_smooth_plot[0] + cumulative_F[0]
        
        fig.add_trace(
            go.Scatter(
                x=lambda_plot_vals, y=FE_ref_smooth_plot_data, # Use smooth values for plotting
                mode='lines', name='ΔF(λ) (ref spring, smooth)', # Updated name
                line=dict(color='purple', width=2, dash='solid'), # Changed dash to solid for smooth
                opacity=0.7, hoverinfo='text'
            ),
            row=plot_rows['FE'], col=1
        )
        # Parabolic fit for data
        valid = ~np.isnan(cumulative_F)
        if np.any(valid):
            p_fe = np.polyfit(lambda_vals[valid], cumulative_F[valid], 2)
            fit_text_fe = f"Data Fit: {p_fe[0]:.4f}λ² + {p_fe[1]:.4f}λ + {p_fe[2]:.4f}"
            fig.add_annotation(
                x=0.01, y=0.99, xref="paper", yref="paper",
                xanchor='left', yanchor='top',
                text=fit_text_fe, showarrow=False,
                font=dict(size=12, color="green"),
                bgcolor="rgba(255, 255, 255, 0.8)",
                row=plot_rows['FE'], col=1
            )
        # Reference info
        ref_text_fe = f"Ref: 0.5*k*d², k={k_spring:.4f}"
        fig.add_annotation(
            x=0.01, y=0.85, xref="paper", yref="paper",
            xanchor='left', yanchor='top',
            text=ref_text_fe, showarrow=False,
            font=dict(size=12, color="purple"),
            bgcolor="rgba(255, 255, 255, 0.8)",
            row=plot_rows['FE'], col=1
        )

    # Add Jarzynski result if available
    if je_dF is not None:
        fig.add_annotation(
            x=0.5, y=0.1, xref="paper", yref="paper",
            xanchor='center', yanchor='bottom',
            text=f"Jarzynski ΔF: {je_dF:.4f} eV",
            showarrow=False,
            font=dict(size=14, color="red", weight="bold"),
            bgcolor="rgba(255, 255, 255, 0.9)",
            bordercolor="red",
            row=plot_rows['FE'], col=1
        )
        # Also add a horizontal line at the end
        fig.add_shape(
            type="line",
            x0=0, y0=je_dF, x1=1, y1=je_dF,
            line=dict(color="red", width=2, dash="dashdot"),
            row=plot_rows['FE'], col=1
        )

    fig.update_yaxes(title_text="ΔF(λ) [eV]", row=plot_rows['FE'], col=1)

    # --- Plot 4: FE Difference ---
    if FE_ref_shifted is not None and cumulative_F is not None:
        # Avoid division by zero, though FE_ref_shifted should be safe
        safe_ref = np.copy(FE_ref_shifted)
        safe_ref[np.abs(safe_ref) < 1e-9] = 1e-9
        fe_diff_percent = ((cumulative_F / safe_ref) - 1) * 100
        fe_diff_abs = cumulative_F - FE_ref_shifted

        # Percentage difference
        fig.add_trace(
            go.Scatter(
                x=lambda_vals, y=fe_diff_percent,
                mode='lines+markers', name='% Diff (data/ref)',
                line=dict(color='cyan', width=2), marker=dict(size=8, symbol='cross'),
            ),
            row=plot_rows['FE_diff'], col=1,
            secondary_y=False
        )
        # Absolute difference
        fig.add_trace(
            go.Scatter(
                x=lambda_vals, y=fe_diff_abs,
                mode='lines+markers', name='Abs Diff (data-ref)',
                line=dict(color='magenta', width=2, dash='dot'), marker=dict(size=8, symbol='x'),
            ),
            row=plot_rows['FE_diff'], col=1,
            secondary_y=True
        )

    fig.update_yaxes(title_text="% Difference", row=plot_rows['FE_diff'], col=1, secondary_y=False)
    fig.update_yaxes(title_text="Abs Difference [eV]", row=plot_rows['FE_diff'], col=1, secondary_y=True)

    # --- Plot 5: Distance ---
    if 'Distance' in plot_rows and distances is not None:
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
    for i in range(1, len(plot_rows) + 1):
        fig.update_xaxes(title_text="λ", row=i, col=1)

    # Add annotations with results
    annotation_text = (
        f"Number of λ windows: {len(lambda_vals)}"
    )

    fig.add_annotation(
        text=annotation_text,
        xref="paper", yref="paper",
        x=0.02, y=0.98,
        xanchor='left', yanchor='top',
        showarrow=False,
        bgcolor="rgba(255,255,255,0.8)",
        bordercolor="black",
        borderwidth=1,
        font=dict(size=11)
    )

    # Update layout
    total_height = 350 * len(plot_rows)
    fig.update_layout(
        height=total_height,
        showlegend=True,
        hovermode='closest',
        template='plotly_white'
    )

    # Determine output filename
    if output_prefix is None:
        output_prefix = filename.replace('.dat', '')

    output_file = f"{output_prefix}_interactive.html"

    # Save interactive HTML
    fig.write_html(output_file)
    print(f"\nInteractive plot saved to: {output_file}")
    print(f"Open this file in a web browser to explore the data interactively.")


def main():
    parser = argparse.ArgumentParser(description="Interactive Thermodynamic Integration Plotter")
    parser.add_argument("--input", type=str, default="entropic_spring_20_TI.dat",
                       help="Input TI data file")
    parser.add_argument("--output", type=str, default=None,
                       help="Output prefix for plots (default: derived from input)")
    parser.add_argument("--N", type=int, default=None,
                       help="Number of segments (for entropic spring reference). If None, tries to guess from filename.")
    parser.add_argument("--T", type=float, default=300.0,
                       help="Temperature in Kelvin (default: 300.0)")
    parser.add_argument("--b", type=float, default=1.198,
                       help="Segment length in Angstroms (default: 1.198)")

    args = parser.parse_args()

    # Try to guess N from filename if not provided
    if args.N is None:
        match = re.search(r'_(\d+)_', args.input)
        if match:
             args.N = int(match.group(1))
             print(f"Guessed N={args.N} from filename.")

    print(f"Reading TI data from: {args.input}")
    N_seg = args.N - 1 if args.N is not None else None
    plot_ti_interactive(args.input, args.output, N_segments=N_seg, T=args.T, b=args.b)

if __name__ == "__main__":
    main()