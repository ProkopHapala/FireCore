import pandas as pd
import os
import re
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import glob
import numpy as np

# Function to extract code from filename
def extract_code_from_filename(filename):
    match = re.search(r'minima__([0-9]+)_', filename)
    if match:
        return match.group(1)
    return None

# Code mapping to names and colors
code_mapping = {
    '11000': {'name': 'No GridFF', 'color': 'green'},
    '11110': {'name': 'GridFF', 'color': 'blue'}
}

# Format large numbers to compact form (e.g., 1M, 250M)
def millions_formatter(x, pos):
    if x >= 1e6:
        return f'{x*1e-6:.0f}M'
    elif x >= 1e3:
        return f'{x*1e-3:.0f}K'
    else:
        return f'{x:.0f}'

# Load data from files
def load_data_from_files(file_pattern):
    data_frames = []
    
    for file_path in glob.glob(file_pattern):
        filename = os.path.basename(file_path)
        code = extract_code_from_filename(filename)
        
        if code in code_mapping:
            # Column names based on data file header
            column_names = [
                "Nb Iteration", "Totally (conv+nonconv)", "Unique conv", "Time", 
                "Totally (conv)", "nbEvaluation", "nStepConvAvg", "nStepNonConvAvg",
                "nStepExploringAvg", "total evaluation"
            ]
            
            # Skip header rows and load data - filtering out lines that contain '[1]Nb Iteration'
            df = pd.read_csv(
                file_path, 
                sep='\s+', 
                comment='[', 
                names=column_names, 
                na_values=['-nan']  # Handle -nan values
            )
            
            # Add code information
            df['Code'] = code
            df['Type'] = code_mapping[code]['name']
            df['Color'] = code_mapping[code]['color']
            
            data_frames.append(df)
    
    return data_frames

# Create a plot with specified x and y scales (linear or log)
def create_plot(data_frames, x_scale='linear', y_scale='linear', y_column='Totally (conv+nonconv)', plot_unique=True):
    fig, ax = plt.subplots(figsize=(14, 10))
    
    # Line styles for distinguishing different data series
    line_styles = ['-', '--']
    
    # Track index for annotation positioning
    for i, df in enumerate(data_frames):
        # Total structures
        ax.plot(df['Time'], df[y_column], 
                color=df['Color'].iloc[0], 
                linestyle=line_styles[0],
                label=f"{df['Type'].iloc[0]}", 
                linewidth=4)
        
        # Add equation annotation next to line for evaluations graph
        if y_column == 'nbEvaluation':
            # Calculate linear fit coefficient (without drawing the line)
            mask = ~np.isnan(df['Time']) & ~np.isnan(df[y_column])  # Ignore NaN values
            coeffs = np.polyfit(df['Time'][mask], df[y_column][mask], 1)
            slope = coeffs[0]
            
            # Format slope for display
            if slope >= 1e3:
                slope_text = f"{slope/1e6:.1f}M"
            else:
                slope_text = f"{slope:.1f}"
            
            # Find a good position for the label (around 70% of x-axis)
            x_pos = df['Time'].max() * 0.7
            mask = (df['Time'] >= x_pos * 0.9) & (df['Time'] <= x_pos * 1.1)
            if mask.any():
                y_pos = df.loc[mask, y_column].mean()  # Y position near the line
                
                # Add annotation directly next to the line
                ax.annotate(f"y~{slope_text}t", 
                           xy=(x_pos, y_pos),  # Position next to the line
                           xytext=(10, 0),  # Small offset
                           textcoords='offset points',
                           fontsize=28,
                           color=df['Color'].iloc[0],
                           bbox=dict(boxstyle="round,pad=0.3", fc="white", alpha=0.8))
        
        # Unique structures (optional)
        if plot_unique and y_column == 'Totally (conv+nonconv)':
            ax.plot(df['Time'], df['Unique conv'], 
                    color=df['Color'].iloc[0], 
                    linestyle=line_styles[1],
                    label=f"{df['Type'].iloc[0]} - Unique structures", 
                    linewidth=4)
    
    # Set axis scales
    ax.set_xscale(x_scale)
    ax.set_yscale(y_scale)
    
    # Set x-axis range to 0-55
    ax.set_xlim(0.1, 55)
    
    # Format y-axis ticks
    if y_column == 'nbEvaluation':
        formatter = ticker.FuncFormatter(millions_formatter)
        ax.yaxis.set_major_formatter(formatter)
    
    # Set labels with larger font size
    ax.set_xlabel('Time [s]', fontsize=28)
    if y_column == 'nbEvaluation':
        ax.set_ylabel('Number of Evaluations', fontsize=28)
    else:
        ax.set_ylabel('Number of Structures', fontsize=28)
    
    # Increase tick label size
    ax.tick_params(axis='both', which='major', labelsize=28)
    
    # Completely reset the legend to match simple_plot.py style
    ax.legend().remove()  # Remove any existing legend
    
    # Create a new legend with explicit frameon parameter
    handles, labels = ax.get_legend_handles_labels()
    legend = ax.legend(handles, labels, loc='lower right', fontsize=20, 
                      frameon=True, framealpha=1.0, edgecolor='black')
    
    # Force update the legend frame
    legend.get_frame().set_linewidth(2.0)
    legend.get_frame().set_edgecolor('black')
    legend.get_frame().set_facecolor('white')
    
    # Add grid
    ax.grid(True, linestyle='--', alpha=0.7)
    
    # Set border styling to match simple_plot.py
    for spine in ax.spines.values():
        spine.set_linewidth(1.0)
        spine.set_color('black')
    
    plt.tight_layout()
    
    # Create descriptive filename
    x_scale_str = 'lin' if x_scale == 'linear' else 'log'
    y_scale_str = 'lin' if y_scale == 'linear' else 'log'
    
    if y_column == 'nbEvaluation':
        filename_prefix = 'evaluations'
    else:
        filename_prefix = 'structures'
    
    filename = f"{filename_prefix}_{x_scale_str}X_{y_scale_str}Y"
    
    # Save plot
    plt.savefig(f"{filename}.png", dpi=300, bbox_inches='tight')
    plt.savefig(f"{filename}.pdf", bbox_inches='tight')
    
    return fig

# Main function
def main():
    # Set style for plots
    plt.style.use('seaborn-v0_8-whitegrid')
    
    # Load data from all .dat files in current directory
    data_frames = load_data_from_files('*.dat')
    
    if not data_frames:
        print("No data files found!")
        return
    
    # Plot 1: Number of evaluations vs time (linear scales)
    fig1 = create_plot(data_frames, x_scale='linear', y_scale='linear', y_column='nbEvaluation', plot_unique=False)
    plt.close(fig1)  # Close the figure after saving
    
    # Plot 2-5: Number of structures vs time with various axis scales
    # Linear X, Linear Y
    fig2 = create_plot(data_frames, x_scale='linear', y_scale='linear')
    plt.close(fig2)  # Close the figure after saving
    
    # Linear X, Log Y
    fig3 = create_plot(data_frames, x_scale='linear', y_scale='log')
    plt.close(fig3)  # Close the figure after saving
    
    # Log X, Linear Y
    fig4 = create_plot(data_frames, x_scale='log', y_scale='linear')
    plt.close(fig4)  # Close the figure after saving
    
    # Log X, Log Y
    fig5 = create_plot(data_frames, x_scale='log', y_scale='log')
    plt.close(fig5)  # Close the figure after saving
    
    print("All plots saved as PNG and PDF files.")

if __name__ == "__main__":
    main()
