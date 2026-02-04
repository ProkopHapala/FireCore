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
                if line.startswith('#') and 'lambda' in line and 'dE' in line:
                    parts = line.replace('#','').split()
                    # Map standard names
                    temp_map = {}
                    for i, p in enumerate(parts):
                        if p == 'lambda': temp_map['lambda'] = i
                        elif p == 'dE/dlambda': temp_map['dE_dlambda'] = i
                        elif p == 'sigma_dE': temp_map['sigma_dE'] = i
                        elif p == 'cumulative_FE': temp_map['cumulative_FE'] = i
                        elif p == 'cumulative_err': temp_map['cumulative_err'] = i
                        elif p == 'distance': temp_map['distance'] = i

                    if 'lambda' in temp_map and 'dE_dlambda' in temp_map:
                        col_map.update(temp_map)
                        header_found = True
                    break
    except Exception as e:
         print(f"Warning: Could not parse header, using default columns: {e}")

    return col_map, header_found

def compute_cumulative_integral_and_error(lambda_vals, dE_dlambda, errors):
    """Compute cumulative integral and propagated error for visualization"""
    cumulative = np.zeros_like(lambda_vals)
    cumulative_error = np.zeros_like(lambda_vals)
    
    use_errors = errors is not None and np.any(errors > 0)
    
    for i in range(1, len(lambda_vals)):
        dx = lambda_vals[i] - lambda_vals[i-1]
        
        # Trapezoidal rule: area = 0.5 * (y1 + y2) * dx
        cumulative[i] = cumulative[i-1] + 0.5 * (dE_dlambda[i] + dE_dlambda[i-1]) * dx
        
        # Error propagation (assuming uncorrelated errors):
        if use_errors:
            var_contribution = (0.5 * dx)**2 * (errors[i]**2 + errors[i-1]**2)
            cumulative_error[i] = np.sqrt(cumulative_error[i-1]**2 + var_contribution)
        
    return cumulative, cumulative_error

def plot_ti_interactive(filename, output_prefix=None, N_segments=None, T=300.0, b=1.198):
    """Create interactive TI plots using Plotly"""
    
    # Detect columns
    col_map, header_found = detect_columns(filename)
    if header_found:
        print(f"Detected columns: {col_map}")
    else:
        print("No header detected, using default column mapping.")

    # Read data
    data = read_ti_data(filename)
    ncols = data.shape[1]
    
    # Safe extraction helper
    def get_col(name):
        idx = col_map.get(name)
        if idx is not None and idx < ncols:
            return data[:, idx]
        return None

    lambda_vals = get_col('lambda')
    dE_dlambda = get_col('dE_dlambda')
    
    if lambda_vals is None or dE_dlambda is None:
        print("Error: Could not extract lambda or dE/dlambda columns.")
        sys.exit(1)
        
    errors = get_col('sigma_dE')
    cumulative_F_read = get_col('cumulative_FE')
    cumulative_err_read = get_col('cumulative_err')
    distances = get_col('distance')

    # Fallback for old format if no header
    if not header_found:
        if ncols < 6:
             if ncols >= 3: errors = data[:, 2]
             else: errors = None

             if ncols >= 4: cumulative_F_read = data[:, 3]
             else: cumulative_F_read = None

             if ncols >= 5: distances = data[:, 4]
             else: distances = None

             cumulative_err_read = None

    # Handle missing errors array
    if errors is None:
        errors = np.zeros_like(lambda_vals)

    # Sort by lambda
    sort_idx = np.argsort(lambda_vals)
    lambda_vals = lambda_vals[sort_idx]
    dE_dlambda = dE_dlambda[sort_idx]
    if errors is not None:
        errors = errors[sort_idx]
    if distances is not None:
        distances = distances[sort_idx]
    if cumulative_F_read is not None:
        cumulative_F_read = cumulative_F_read[sort_idx]
    if cumulative_err_read is not None:
        cumulative_err_read = cumulative_err_read[sort_idx]

    # Compute cumulative integral and error
    cumulative_F, cumulative_error = compute_cumulative_integral_and_error(lambda_vals, dE_dlambda, errors)

    # Use read values if available (prioritize file content as requested)
    if cumulative_F_read is not None:
        cumulative_F = cumulative_F_read
    if cumulative_err_read is not None:
        cumulative_error = cumulative_err_read
    
    # Integrate to get total free energy (redundant with cumulative_F[-1] but good to keep)
    delta_F_trapz = np.trapz(dE_dlambda, lambda_vals)
    if len(lambda_vals) >= 3:
        delta_F_simps = integrate.simpson(y=dE_dlambda, x=lambda_vals)
    else:
        delta_F_simps = delta_F_trapz
    
    # Estimate total uncertainty from final cumulative error
    total_error = cumulative_error[-1]
    
    # Create subplots
    fig = make_subplots(
        rows=3, cols=1,
        subplot_titles=(
            'dE/dλ vs λ (Thermodynamic Integration)',
            'Cumulative Free Energy ΔF(λ)',
            'Si-Si Distance vs λ' if distances is not None else 'Error Analysis'
        ),
        vertical_spacing=0.12,
        row_heights=[0.35, 0.35, 0.30]
    )
    
    # Prepare hover text
    hover_text_1 = [
        f"λ = {lam:.4f}<br>dE/dλ = {de:.4f} eV<br>σ(F) = {err:.4f} eV" + 
        (f"<br>Distance = {dist:.3f} Å" if distances is not None else "")
        for lam, de, err, dist in zip(lambda_vals, dE_dlambda, errors, 
                                       distances if distances is not None else [0]*len(lambda_vals))
    ]
    
    # Plot 1: dE/dlambda vs lambda with error bars
    fig.add_trace(
        go.Scatter(
            x=lambda_vals, y=dE_dlambda,
            mode='lines+markers',
            name='dE/dλ',
            line=dict(color='blue', width=2),
            marker=dict(size=8, symbol='circle'),
            error_y=dict(
                type='data',
                array=errors,
                visible=True,
                color='rgba(0,0,255,0.3)'
            ) if np.any(errors > 0) else None,
            hovertext=hover_text_1,
            hoverinfo='text'
        ),
        row=1, col=1
    )

    # --- Linear Fit for dE/dlambda ---
    if len(lambda_vals) > 1:
        fit_coeffs = np.polyfit(lambda_vals, dE_dlambda, 1)
        fit_line = np.polyval(fit_coeffs, lambda_vals)
        fit_expression = f"dE/dλ = {fit_coeffs[0]:.4f}λ + {fit_coeffs[1]:.4f}"

        # Add fitted line to plot
        fig.add_trace(
            go.Scatter(
                x=lambda_vals, y=fit_line,
                mode='lines',
                name='Linear Fit',
                line=dict(color='cyan', width=2, dash='dot'),
                hoverinfo='skip'
            ),
            row=1, col=1
        )

        # Add annotation with the fit expression to the first subplot
        fig.add_annotation(
            text=fit_expression,
            xref="paper", yref="paper",
            x=0.05, y=0.95,
            xanchor='left', yanchor='top',
            showarrow=False,
            bgcolor="rgba(255, 255, 255, 0.8)",
            bordercolor="black",
            borderwidth=1,
            row=1, col=1
        )
    
    # Entropic Spring Reference (Force)
    if distances is not None and N_segments is not None:
        k_spring = kB * T / (N_segments * b**2)
        F_ref = k_spring * distances
        
        fig.add_trace(
            go.Scatter(
                x=lambda_vals, y=F_ref,
                mode='lines',
                name='Ref Force',
                line=dict(color='purple', width=2, dash='dash'),
                opacity=0.7
            ),
            row=1, col=1
        )
    
    # Plot 2: Cumulative integral
    hover_text_2 = [
        f"λ = {lam:.4f}<br>ΔF(λ) = {cf:.4f} eV<br>σ(ΔF) = {ce:.4f} eV"
        for lam, cf, ce in zip(lambda_vals, cumulative_F, cumulative_error)
    ]
    
    fig.add_trace(
        go.Scatter(
            x=lambda_vals, y=cumulative_F,
            mode='lines+markers',
            name='ΔF(λ)',
            line=dict(color='green', width=2),
            marker=dict(size=8, symbol='diamond'),
            error_y=dict(
                type='data',
                array=cumulative_error,
                visible=True,
                color='rgba(0,255,0,0.3)'
            ),
            hovertext=hover_text_2,
            hoverinfo='text',
        ),
        row=2, col=1
    )
    
    # Entropic Spring Reference (Free Energy)
    if distances is not None and N_segments is not None:
        k_spring = kB * T / (N_segments * b**2)
        FE_ref = 0.5 * k_spring * distances**2
        FE_ref_shifted = FE_ref - FE_ref[0] + cumulative_F[0]
        
        hover_text_ref_fe = [
            f"λ = {lam:.4f}<br>Ref FE = {fe:.4f} eV<br>Distance = {dist:.3f} Å"
            for lam, fe, dist in zip(lambda_vals, FE_ref_shifted, distances)
        ]
        
        fig.add_trace(
            go.Scatter(
                x=lambda_vals, y=FE_ref_shifted,
                mode='lines',
                name='Ref FE (spring)',
                line=dict(color='purple', width=2, dash='dash'),
                opacity=0.7,
                hovertext=hover_text_ref_fe,
                hoverinfo='text'
            ),
            row=2, col=1
        )
    
    # Plot 3: Distance or Error analysis
    if distances is not None:
        hover_text_3 = [
            f"λ = {lam:.4f}<br>Distance = {dist:.3f} Å"
            for lam, dist in zip(lambda_vals, distances)
        ]
        
        fig.add_trace(
            go.Scatter(
                x=lambda_vals, y=distances,
                mode='lines+markers',
                name='Si-Si Distance',
                line=dict(color='red', width=2),
                marker=dict(size=8, symbol='square'),
                hovertext=hover_text_3,
                hoverinfo='text'
            ),
            row=3, col=1
        )
    else:
        # Show error bars if no distance data
        hover_text_3 = [
            f"λ = {lam:.4f}<br>σ(F) = {err:.4f} eV"
            for lam, err in zip(lambda_vals, errors)
        ]
        
        fig.add_trace(
            go.Scatter(
                x=lambda_vals, y=errors,
                mode='lines+markers',
                name='σ(F)',
                line=dict(color='orange', width=2),
                marker=dict(size=8, symbol='triangle-up'),
                hovertext=hover_text_3,
                hoverinfo='text'
            ),
            row=3, col=1
        )

    # Update axes labels
    fig.update_xaxes(title_text="λ", row=1, col=1)
    fig.update_yaxes(title_text="dE/dλ (eV)", row=1, col=1)

    fig.update_xaxes(title_text="λ", row=2, col=1)
    fig.update_yaxes(title_text="ΔF(λ) (eV)", row=2, col=1)

    fig.update_xaxes(title_text="λ", row=3, col=1)
    if distances is not None:
        fig.update_yaxes(title_text="Si-Si Distance (Å)", row=3, col=1)
    else:
        fig.update_yaxes(title_text="σ(F) (eV)", row=3, col=1)

    # Add annotations with results
    annotation_text = (
        f"<b>Free Energy Results:</b><br>"
        f"ΔF (Trapezoidal) = {delta_F_trapz:.6f} eV<br>"
        f"ΔF (Simpson) = {delta_F_simps:.6f} eV<br>"
        f"Total σ(F) = {total_error:.6f} eV<br>"
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
    fig.update_layout(
        title=dict(
            text=f"Interactive Thermodynamic Integration Analysis<br><sub>{filename}</sub>",
            x=0.5,
            xanchor='center'
        ),
        height=1200,
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

    # Print summary
    print("\n" + "="*60)
    print("  Thermodynamic Integration Summary")
    print("="*60)
    print(f"  Number of λ windows: {len(lambda_vals)}")
    print(f"  λ range: {lambda_vals[0]:.4f} → {lambda_vals[-1]:.4f}")
    print(f"  Free energy change (trapezoidal): {delta_F_trapz:.6f} eV")
    print(f"  Free energy change (Simpson):     {delta_F_simps:.6f} eV")
    if total_error > 0:
        print(f"  Total uncertainty σ(F):            {total_error:.6f} eV")
    if distances is not None:
        print(f"  Distance range: {distances[0]:.3f} → {distances[-1]:.3f} Å")
    if N_segments is not None:
         print(f"  Entropic Spring Reference (N={N_segments}, T={T}K, b={b}A)")
    print("="*60)

    return delta_F_trapz, delta_F_simps

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
    plot_ti_interactive(args.input, args.output, N_segments=args.N-1, T=args.T, b=args.b)

if __name__ == "__main__":
    main()
