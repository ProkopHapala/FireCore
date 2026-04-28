import pandas as pd
import re
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import pandas as pd
import numpy as np
import argparse
import webbrowser
import os
import json

# Function to convert surface size to numeric value for sorting
def surface_size_to_numeric(size_str):
    """Convert a surface size string like '5x5' to a numeric value for proper sorting"""
    try:
        if 'x' in size_str:
            # Parse the dimensions from the string
            dimensions = size_str.split('x')
            if len(dimensions) == 2:
                x, y = map(int, dimensions)
                # For NxN grids, use N directly for more natural sorting
                if x == y:
                    return x
                # For non-square grids, use the first dimension as primary and second as secondary
                else:
                    return x + y/100.0  # This ensures 3x4 sorts after 3x3 but before 4x1
        # Handle case where it might just be a number
        return int(size_str)
    except (ValueError, TypeError):
        # If we can't convert, use a string comparison approach
        # This is a fallback for unexpected formats
        return 0
        
# Function to calculate the actual number of surface atoms
def calculate_surface_atoms(size_str, npbc_val):
    """Calculate the actual number of surface atoms given the surface size and nPBC value"""
    try:
        if 'x' in size_str:
            dimensions = size_str.split('x')
            if len(dimensions) == 2:
                N = int(dimensions[0])  # Using N from NxN
                
                # Parse nPBC to get the x component
                npbc_x = 0
                
                # Handle format like (1,1,0) or _1_1_0_ (escaped format)
                if '(' in npbc_val and ')' in npbc_val:
                    # Extract values within parentheses and split by comma
                    inner_values = npbc_val.strip('()').split(',')
                    if len(inner_values) >= 1:
                        npbc_x = int(inner_values[0])
                elif npbc_val.startswith('_') and npbc_val.endswith('_'):
                    # Handle escaped format like _1_1_0_
                    parts = npbc_val.strip('_').split('_')
                    if len(parts) >= 1 and parts[0].isdigit():
                        npbc_x = int(parts[0])
                # Handle format with 'x' like 1x1x0
                elif 'x' in npbc_val:
                    npbc_components = npbc_val.split('x')
                    npbc_x = int(npbc_components[0])
                # Handle single value format
                elif npbc_val != '--':
                    try:
                        npbc_x = int(npbc_val)  # Assume single value means all components are the same
                    except ValueError:
                        npbc_x = 0  # Default fallback
                    
                # Calculate using formula 6*N*N*(nPBC.x+2)^2
                return 6 * N * N * (2*npbc_x + 1)**2#(2*nPBC.x+1)*(2*nPBC.y+1)*nx*ny
        return 0  # Default return if format is unexpected
    except (ValueError, TypeError, IndexError):
        # Print debugging information
        print(f"Error calculating surface atoms for size={size_str}, npbc={npbc_val}")
        return 0  # Return 0 for any parsing errors

# Column names provided by the user
column_names = [
    "Code",
    "Surface Size",
    "nPBC",
    "Nloc_MMFF",
    "Nloc_move",
    "Nloc_NBFF",
    "Nloc_surf",
    "Nloc_gridFF",
    "Nloc_GridFF_bSpline",
    "Replica",
    "Perframe",
    "PerVF",
    "Nb Iteration",
    "Totally (conv+nonconv)",
    "Unique conv",
    "Time",
    "Totally (conv)",
    "nbEvaluation",
    "nStepConvAvg",
    "nStepNonConvAvg",
    "nStepExploringAvg",
    "total evaluation",
    "Evaluations per second",
]

# Read the data from testing_data.dat
with open("results/all/results.dat", "r") as f:
    lines = f.readlines()

# Prepare the data list
data = []
for line in lines:  # Process all lines including the first one
    values = line.strip().split()
    name = values[0]
    # Extract parameters from the name using regex
    # Updated regex to handle both old format with parentheses and new escaped format with double underscores
    match = re.match(r"minima__([0-9]+)_surf:_([^_]+)_([0-9]+x[0-9]+)_nPBC_(.+?)_nloc:_MMFF_([0-9-]+)_move_([0-9-]+)_NBFF_([^_]+)_surf_([^_]+)_gridFF_([^_]+)_gridFFbSpline_([^_]+)___replica:_([0-9]+)_perframe:_([0-9]+)_perVF:_([0-9]+)", name)
    if match:
        code = match.group(1)
        surface_type = match.group(2)  # Now captures the surface type (e.g., NaCl-Cl-hole)
        surface_size = match.group(3)  # Now specifically captures the NxN part
        npbc = match.group(4)
        nloc_MMFF = match.group(5)
        nloc_move = match.group(6)
        nloc_NBFF = match.group(7)
        nloc_surf = match.group(8)
        nloc_gridFF = match.group(9)
        nloc_GridFF_bSpline = match.group(10)
        replica = match.group(11)
        perframe = match.group(12)
        perVF = match.group(13)

        # Create a dictionary for the current row
        row = {
            "Code": code,
            "Surface Size": surface_size,
            "Surface Atoms": calculate_surface_atoms(surface_size, npbc),
            "nPBC": npbc,
            "Nloc_MMFF": nloc_MMFF,
            "Nloc_move": nloc_move,
            "Nloc_NBFF": nloc_NBFF,
            "Nloc_surf": nloc_surf,
            "Nloc_gridFF": nloc_gridFF,
            "Nloc_GridFF_bSpline": nloc_GridFF_bSpline,
            "Replica": replica,
            "Perframe": perframe,
            "PerVF": perVF,
        }

        # Add the measured data - ensure we don't go out of bounds
        for i, col_name in enumerate(column_names[12:]):
            if i+1 < len(values):  # +1 to account for the name at index 0
                row[col_name] = values[i+1]  # +1 to skip the name field
            else:
                row[col_name] = 0  # Default value if missing
                
        # Calculate and add 'Evaluations per second'
        try:
            if float(values[4]) != 0:
                row['Evaluations per second'] = float(values[6]) / float(values[4])
            else:
                row['Evaluations per second'] = 0
        except (IndexError, ValueError):
            row['Evaluations per second'] = 0
        data.append(row)
    else:
        print(f"No match for name: {name}")

# Create the Pandas DataFrame
df = pd.DataFrame(data, columns=column_names)

# Add a numeric column for surface size to enable proper sorting
df['Surface Size Numeric'] = df['Surface Size'].apply(surface_size_to_numeric)

# Add a column for the actual number of surface atoms
df['Surface Atoms'] = df.apply(lambda row: calculate_surface_atoms(row['Surface Size'], row['nPBC']), axis=1)

# Sort the DataFrame by the numeric surface size
df = df.sort_values('Surface Size Numeric')

# You can further process the DataFrame here, e.g., save it to a CSV file
df.to_csv("results_table.csv", index=False)

# Load the data from the CSV file
data = pd.read_csv("results_table.csv")

# Create a mapping between surface sizes and their numeric positions for proper ordering
size_order = {
    '0x0': 0,
    '1x1': 1,
    '2x2': 2,
    '3x3': 3,
    '4x4': 4,
    '5x5': 5,
    '6x6': 6,
    '7x7': 7,
    '8x8': 8,
    '9x9': 9,
    '10x10': 10,
    '11x11': 11,
    '12x12': 12,
    '13x13': 13,
    '14x14': 14,
    '15x15': 15,
    '16x16': 16
}

# Add numeric values for sorting and plotting
data['x_value'] = data['Surface Size'].map(size_order)

# Ensure data is sorted by these numeric values
data = data.sort_values('x_value')

# Calculate actual surface atoms count
data['Surface Atoms'] = data.apply(lambda row: calculate_surface_atoms(row['Surface Size'], row['nPBC']), axis=1)

# Create a list of unique surface sizes in correct order for axis ticks
ordered_sizes = []
size_to_pos = {}
pos = 0

for size in sorted([(k, v) for k, v in size_order.items() if k in data['Surface Size'].unique()], key=lambda x: x[1]):
    ordered_sizes.append(size[0])
    size_to_pos[size[0]] = pos
    pos += 1

# Map each surface size string to its position
data['position'] = data['Surface Size'].map(size_to_pos)

# Create a dictionary to map codes to labels
labels = {
    '10000': 'Without surface',
    '11000': 'No GridFF',
    '11110': 'GridFF with BSpline'
}

# Define specific colors for each category
color_map = {
    '10000': 'red',    # Without surface always red
    '11000': 'green',  # No GridFF always green
    '11110': 'blue'    # GridFF with BSpline always blue
}

# Find the largest surface size for each code group
largest_sizes = {}
for code in labels.keys():
    code_data = data[data['Code'] == int(code)]
    if not code_data.empty:
        largest_size = code_data['position'].max()
        largest_sizes[code] = largest_size

# Find the best performer for each code at its largest surface size
best_performers = {}
for code in labels.keys():
    code_data = data[data['Code'] == int(code)]
    if not code_data.empty:
        # Get data for this code at its largest surface size
        largest_size = largest_sizes[code]
        largest_size_data = code_data[code_data['position'] == largest_size]
        
        # Find the best performer
        if not largest_size_data.empty:
            best_idx = largest_size_data['Evaluations per second'].idxmax()
            best_performer = data.loc[best_idx]
            best_performers[code] = {
                'replica': best_performer['Replica'],
                'perframe': best_performer['Perframe'],
                'pervf': best_performer['PerVF'],
                'npbc': best_performer['nPBC'],
                'nloc_MMFF': best_performer['Nloc_MMFF'],
                'nloc_move': best_performer['Nloc_move'],
                'nloc_NBFF': best_performer['Nloc_NBFF'],
                'nloc_surf': best_performer['Nloc_surf'],
                'nloc_gridFF': best_performer['Nloc_gridFF'],
                'nloc_GridFF_bSpline': best_performer['Nloc_GridFF_bSpline'],
                'performance': best_performer['Evaluations per second'],
                'surface_size': best_performer['Surface Size']
            }

# # Print the best performers for each code group
# for code, info in best_performers.items():
#     print(f"Best performance for {labels[code]} at size {info['surface_size']}:")
#     print(f"  Replica: {info['replica']}, Perframe: {info['perframe']}, PerVF: {info['pervf']}")
#     print(f"  Evaluations per second: {info['performance']}")


# Create an interactive Plotly figure with light blue background and centered display
fig = go.Figure()
fig.update_layout(
    plot_bgcolor='rgba(240, 248, 255, 0.95)',  # Light blue background (aliceblue)
    paper_bgcolor='white',
    autosize=True, # Allow resizing while maintaining proportions
    margin=dict(l=80, r=80, b=100, t=100),  # Expanded margins to show all labels
    height=700,   # Fixed height to ensure enough vertical space
    xaxis=dict(
        showgrid=True,
        gridcolor='white',
        gridwidth=1,
        title_font=dict(size=16),  # Larger axis title
        tickfont=dict(size=14),    # Larger tick labels
    ),
    yaxis=dict(
        showgrid=True,
        gridcolor='white',
        gridwidth=1,
        title_font=dict(size=16),  # Larger axis title
        tickfont=dict(size=14),    # Larger tick labels
    ),
    legend=dict(
        font=dict(size=14),  # Larger legend text
        xanchor='center',     # Center the legend horizontally
        yanchor='top',        # Anchor to top
        y=-0.2                # Position below the plot
    ),
    font=dict(size=14)  # Larger global font size
)

# Parse command-line arguments
parser = argparse.ArgumentParser(description='Generate performance plot with customizable x-axis')
parser.add_argument('--nxn', action='store_true', help='Display NxN notation on x-axis instead of actual atom count')
args = parser.parse_args()

# Use NxN notation if specified, otherwise use actual atom count (default)
use_atom_count = not args.nxn

# Track lines added to the plot for each code
code_has_label = {code: False for code in labels.keys()}

# First group data by Code to ensure consistent appearance
for code in sorted(data['Code'].unique()):
    # Only process codes we have labels for
    if str(code) in labels:
        # Get all data for this code
        code_data = data[data['Code'] == code]
        
        # Create separate traces for each combination of parameters (including all Nloc parameters)
        for param_combo, group_data in code_data.groupby(['nPBC', 'Nloc_MMFF', 'Nloc_move', 'Nloc_NBFF', 
                                                    'Nloc_surf', 'Nloc_gridFF', 'Nloc_GridFF_bSpline',
                                                    'Replica', 'Perframe', 'PerVF']):
            # Sort by the numeric position
            group_data = group_data.sort_values(by='position')
            
            # Extract parameters from the parameter combination
            npbc, nloc_MMFF, nloc_move, nloc_NBFF, nloc_surf, nloc_gridFF, nloc_GridFF_bSpline, replica, perframe, pervf = param_combo
            
            # Use consistent color for this code
            color = color_map[str(code)]
            
            # Create hover text with the correct parameters for each point
            hover_text = []
            for _, row in group_data.iterrows():
                # Start with basic info
                hover_info = [f"<b>{labels[str(code)]}</b><br>",
                             f"Surface Size: {row['Surface Size']}<br>",
                             f"Surface Atoms: {row['Surface Atoms']}<br>",
                             f"nPBC: {row['nPBC']}<br>",
                             f"Nloc_MMFF: {row['Nloc_MMFF']}<br>",
                             f"Nloc_move: {row['Nloc_move']}<br>",
                             f"Nloc_NBFF: {row['Nloc_NBFF']}<br>",
                             f"Nloc_surf: {row['Nloc_surf']}<br>",
                             f"Nloc_gridFF: {row['Nloc_gridFF']}<br>",
                             f"Nloc_GridFF_bSpline: {row['Nloc_GridFF_bSpline']}<br>",
                             f"Performance: {row['Evaluations per second']:.2f} evals/sec<br>",
                             f"<b>Parameters:</b><br>",
                             f"Replica: {row['Replica']}<br>",
                             f"Perframe: {row['Perframe']}<br>",
                             f"PerVF: {row['PerVF']}"
                             ]
                
                # Join all parts into a single hover text
                hover_text.append(''.join(hover_info))
            
            # Check if this is the best performer for this code
            is_best_performer = False
            if str(code) in best_performers:
                best_info = best_performers[str(code)]
                if (str(replica) == str(best_info['replica']) and 
                    str(perframe) == str(best_info['perframe']) and 
                    str(pervf) == str(best_info['pervf']) and
                    # Also check if Nloc parameters match
                    str(npbc) == str(best_info['npbc']) and
                    str(nloc_MMFF) == str(best_info.get('nloc_MMFF', '--')) and
                    str(nloc_move) == str(best_info.get('nloc_move', '--')) and
                    str(nloc_NBFF) == str(best_info.get('nloc_NBFF', '--')) and
                    str(nloc_surf) == str(best_info.get('nloc_surf', '--')) and
                    str(nloc_gridFF) == str(best_info.get('nloc_gridFF', '--')) and
                    str(nloc_GridFF_bSpline) == str(best_info.get('nloc_GridFF_bSpline', '--'))):
                    is_best_performer = True
            
            # All lines start with the same appearance with low opacity
            line_width = 2
            opacity = 0.3  # Low opacity for all lines
            
            # Only show in legend if this is the first line for this code group
            show_legend = not code_has_label[str(code)]
            name = labels[str(code)] if show_legend else None
            
            # Add the trace to the figure using position values for correct ordering
            # But use either atom count or position for display depending on mode
            x_values = group_data["Surface Atoms"].tolist() if use_atom_count else group_data["position"].tolist()
            fig.add_trace(go.Scatter(
                x=x_values,
                y=group_data["Evaluations per second"].tolist(),
                mode="lines+markers",
                name=name,
                line=dict(color=color, width=line_width),
                marker=dict(
                    size=12 if is_best_performer else 10,  # Larger markers, even larger for best performers
                    symbol='circle',
                    line=dict(width=2, color=color)  # Add marker outline
                ),
                opacity=opacity,
                hoverinfo="text",
                hovertext=hover_text,
                showlegend=show_legend
            ))
            
            # Mark this code as having a label now
            if show_legend:
                code_has_label[str(code)] = True

# Enable full line highlighting on hover with custom JavaScript
fig.update_layout(
    title=dict(
        text="Evaluations per second vs. Surface Size",
        font=dict(size=24)  # Larger title font
    ),
    xaxis=dict(
        title=dict(
            text="Surface Atoms" if use_atom_count else "Surface Size",
            font=dict(size=18)  # Larger axis title font
        ),
        # Configure x-axis based on display mode
        **({
            'type': 'linear',
            'tickfont': dict(size=14),
            # For atom count mode, use auto-generated evenly spaced ticks
            'autorange': True
        } if use_atom_count else {
            'tickmode': 'array',
            'tickvals': list(range(len(ordered_sizes))),
            'ticktext': ordered_sizes,
            'type': 'category'
        })
    ),
    yaxis=dict(
        title=dict(
            text="Evaluations per second",
            font=dict(size=18)  # Larger axis title font
        ),
        tickfont=dict(size=14)  # Larger tick font
    ),
    hovermode="closest",
    legend=dict(
        orientation="h",
        yanchor="bottom",
        y=1.02,
        xanchor="right",
        x=1,
        font=dict(size=14)  # Larger legend font
    ),
    # Make the plot more interactive
    hoverdistance=20,  # Distance to highlight closest point
    # Add more interactive behavior with hover
    hoverlabel=dict(
        bgcolor="white",
        font_size=16,  # Larger hover text
        font_family="Arial"
    ),
    # Add custom JavaScript to highlight entire lines when hovering
    updatemenus=[
        dict(
            type="buttons",
            showactive=False,
            buttons=[
                dict(
                    label="",
                    method="relayout",
                    args=[{}, {}],
                    visible=False
                )
            ]
        )
    ]
)

# Add hover effects with custom behavior
fig.update_traces(
    hovertemplate='%{hovertext}',
    line=dict(width=2),  # Default line width
    hoverlabel=dict(
        bgcolor="white",
        font_size=12,
        font_family="Arial"
    )
)

# Determine which mode was used for the output filenames
mode_suffix = "_atoms" if use_atom_count else "_NxN"

# Create a fully self-contained interactive HTML plot with toggle controls
# Use full_html=True to include everything in one file
# Store both the NxN and atom count data for x-axis toggle functionality

# We need to create and store just the essential data for x-axis toggling
# Create a JSON-serializable object with only the data we need
toggle_data = {}

# Capture initial plot layout for background and styling
toggle_data['original_layout'] = {
    'plot_bgcolor': fig.layout.plot_bgcolor if hasattr(fig.layout, 'plot_bgcolor') else '#ffffff',
    'paper_bgcolor': fig.layout.paper_bgcolor if hasattr(fig.layout, 'paper_bgcolor') else '#ffffff',
    'font': {
        'family': fig.layout.font.family if hasattr(fig.layout, 'font') and hasattr(fig.layout.font, 'family') else 'Arial, sans-serif',
        'size': fig.layout.font.size if hasattr(fig.layout, 'font') and hasattr(fig.layout.font, 'size') else 12,
        'color': fig.layout.font.color if hasattr(fig.layout, 'font') and hasattr(fig.layout.font, 'color') else '#000000'
    },
    'margin': {
        'l': fig.layout.margin.l if hasattr(fig.layout, 'margin') and hasattr(fig.layout.margin, 'l') else 50,
        'r': fig.layout.margin.r if hasattr(fig.layout, 'margin') and hasattr(fig.layout.margin, 'r') else 50,
        't': fig.layout.margin.t if hasattr(fig.layout, 'margin') and hasattr(fig.layout.margin, 't') else 50,
        'b': fig.layout.margin.b if hasattr(fig.layout, 'margin') and hasattr(fig.layout.margin, 'b') else 50,
        'pad': fig.layout.margin.pad if hasattr(fig.layout, 'margin') and hasattr(fig.layout.margin, 'pad') else 4
    }
}

# Capture grid settings if available
if hasattr(fig.layout, 'xaxis'):
    toggle_data['original_layout']['xaxis'] = {
        'showgrid': fig.layout.xaxis.showgrid if hasattr(fig.layout.xaxis, 'showgrid') else True,
        'gridcolor': fig.layout.xaxis.gridcolor if hasattr(fig.layout.xaxis, 'gridcolor') else '#e0e0e0',
        'gridwidth': fig.layout.xaxis.gridwidth if hasattr(fig.layout.xaxis, 'gridwidth') else 1,
        'zeroline': fig.layout.xaxis.zeroline if hasattr(fig.layout.xaxis, 'zeroline') else False
    }

if hasattr(fig.layout, 'yaxis'):
    toggle_data['original_layout']['yaxis'] = {
        'showgrid': fig.layout.yaxis.showgrid if hasattr(fig.layout.yaxis, 'showgrid') else True,
        'gridcolor': fig.layout.yaxis.gridcolor if hasattr(fig.layout.yaxis, 'gridcolor') else '#e0e0e0',
        'gridwidth': fig.layout.yaxis.gridwidth if hasattr(fig.layout.yaxis, 'gridwidth') else 1,
        'zeroline': fig.layout.yaxis.zeroline if hasattr(fig.layout.yaxis, 'zeroline') else False
    }

# Store the x-axis tick settings for both modes
toggle_data['nxn_layout'] = {
    'xaxis.title.text': 'Surface Size',
    'xaxis.tickmode': 'array',
    'xaxis.tickvals': list(range(len(ordered_sizes))),
    'xaxis.ticktext': ordered_sizes,
    'xaxis.type': 'category'
}

# Settings for atom count mode
toggle_data['atom_layout'] = {
    'xaxis.title.text': 'Surface Atoms',
    'xaxis.tickmode': 'auto',
    'xaxis.type': 'linear',
    'xaxis.autorange': True
}

# Store trace data for x-axis toggle
toggle_data['traces'] = []

# Extract and store x-values and trace properties for both display modes for each trace
fig_traces = []

# First capture all the trace data from the current figure
for i, trace in enumerate(fig.data):
    # Get line properties
    line_width = trace.line.width if trace.line and hasattr(trace.line, 'width') else 2
    line_color = trace.line.color if trace.line and hasattr(trace.line, 'color') else None
    
    # Get marker properties
    marker_size = trace.marker.size if trace.marker and hasattr(trace.marker, 'size') else 6
    marker_color = trace.marker.color if trace.marker and hasattr(trace.marker, 'color') else line_color
    
    fig_traces.append({
        'index': i,
        'name': trace.name,
        'legendgroup': trace.legendgroup,
        'showlegend': trace.showlegend,
        'line': {
            'color': line_color,
            'width': line_width,
        },
        'marker': {
            'color': marker_color,
            'size': marker_size
        },
        'mode': trace.mode,
        'opacity': trace.opacity,
        'hovertext': trace.hovertext if hasattr(trace, 'hovertext') else None
    })

# Now extract the x-values for toggling from our grouped data
trace_index = 0
for code in sorted(data['Code'].unique()):
    if str(code) not in labels:
        continue
        
    # Get all data for this code
    code_data = data[data['Code'] == code]
    code_color = color_map[str(code)]
    
    # Create separate traces for each combination of parameters
    for param_combo, group_data in code_data.groupby(['nPBC', 'Nloc_MMFF', 'Nloc_move', 'Nloc_NBFF', 
                                             'Nloc_surf', 'Nloc_gridFF', 'Nloc_GridFF_bSpline',
                                             'Replica', 'Perframe', 'PerVF']):
        # Sort by the numeric position
        group_data = group_data.sort_values(by='position')
        
        # Create hover text with the correct parameters for each point
        hover_text = []
        npbc, nloc_MMFF, nloc_move, nloc_NBFF, nloc_surf, nloc_gridFF, nloc_GridFF_bSpline, replica, perframe, pervf = param_combo
        
        for _, row in group_data.iterrows():
            # Start with basic info
            hover_info = [f"<b>{labels[str(code)]}</b><br>",
                         f"Surface Size: {row['Surface Size']}<br>",
                         f"nPBC: {row['nPBC']}<br>",
                         f"Nloc_MMFF: {row['Nloc_MMFF']}<br>",
                         f"Nloc_move: {row['Nloc_move']}<br>",
                         f"Nloc_NBFF: {row['Nloc_NBFF']}<br>",
                         f"Nloc_surf: {row['Nloc_surf']}<br>",
                         f"Nloc_gridFF: {row['Nloc_gridFF']}<br>",
                         f"Nloc_GridFF_bSpline: {row['Nloc_GridFF_bSpline']}<br>",
                         f"Performance: {row['Evaluations per second']:.2f} evals/sec<br>",
                         f"<b>Parameters:</b><br>",
                         f"Replica: {row['Replica']}<br>",
                         f"Perframe: {row['Perframe']}<br>",
                         f"PerVF: {row['PerVF']}"
                        ]
            hover_text.append(''.join(hover_info))
        
        # Store the data for both x-axis modes with all relevant properties
        # Check if we're the first trace in this code group for legend visibility
        show_legend = not code in [t.get('code') for t in toggle_data['traces'] if t.get('showlegend', False)]
        
        toggle_data['traces'].append({
            'nxn_x': group_data["position"].tolist(),
            'atom_x': group_data["Surface Atoms"].tolist(),
            'y': group_data["Evaluations per second"].tolist(),
            'index': trace_index,
            'code': code,
            'name': labels[str(code)],
            'color': code_color,
            'showlegend': show_legend,
            'legendgroup': str(code),
            'hovertext': hover_text,
            'npbc': npbc,
            'replica': replica,
            'perframe': perframe,
            'pervf': pervf
        })
        
        trace_index += 1

# Helper function to convert NumPy types to standard Python types for JSON serialization
def convert_to_serializable(obj):
    if hasattr(obj, 'tolist'):
        return obj.tolist()
    elif hasattr(obj, 'item'):
        return obj.item()
    elif isinstance(obj, dict):
        return {k: convert_to_serializable(v) for k, v in obj.items()}
    elif isinstance(obj, list):
        return [convert_to_serializable(i) for i in obj]
    else:
        return obj

# Convert any NumPy types in our data for JSON serialization
toggle_data_serializable = convert_to_serializable(toggle_data)

# Convert the data to JSON for use in JavaScript
fig_data_json = json.dumps(toggle_data_serializable)

# Save the interactive HTML with controls and script injection
fig_html = fig.to_html(include_plotlyjs=True, full_html=True)

# Define the injection point for the control panel
injection_point = "<body>"

# Generate mode suffix based on x-axis display mode
mode_label = "NxN" if args.nxn else "AtomCount"

# Inject our custom control for y-axis toggle right after the <body> tag
control_panel_html = '''
<div class="plot-container">
    <!-- Control panel with centered content -->
    <div style="width: 100%; padding: 15px; background: #f5f5f5; border-radius: 5px; display: flex; align-items: center; justify-content: center; flex-wrap: wrap; margin-bottom: 20px;">
        <div style="margin: 0 30px 0 0; display: flex; align-items: center;">
            <span style="font-weight: bold; margin-right: 10px; font-size: 16px;">Y-Axis Scale:</span>
            <button id="toggleYScale" style="padding: 10px 15px; cursor: pointer; background: #4CAF50; color: white; border: none; border-radius: 4px; font-size: 16px;">
                Switch to Log Scale
            </button>
        </div>
        <div style="display: flex; align-items: center;">
            <span style="font-weight: bold; margin-right: 10px; font-size: 16px;">X-Axis Display:</span>
            <button id="toggleXAxis" style="padding: 10px 15px; cursor: pointer; background: #2196F3; color: white; border: none; border-radius: 4px; font-size: 16px;">
                Switch to ''' + ('NxN Notation' if use_atom_count else 'Atom Count') + '''
            </button>
        </div>
    </div>
</div>
'''

# Add the plot data as a hidden element and JavaScript before the </body> tag
script_injection = f'''
<script id="plotData" type="application/json">{fig_data_json}</script>
<script>'''

# Add the interactive JavaScript
script_injection += '''
    // Add y-axis scale toggle and hover highlighting functionality
    document.addEventListener('DOMContentLoaded', function() {
        // Setup for highlighting
        const plotDiv = document.getElementsByClassName('plotly-graph-div')[0];
        const plotDivId = plotDiv.id;
        let lockedTrace = null; // Track which trace is currently locked (if any)
        
        // Setup hover effect to highlight traces (only when no trace is locked)
        plotDiv.on('plotly_hover', function(data) {
            if (lockedTrace === null && data.points && data.points.length > 0) {
                const curveNumber = data.points[0].curveNumber;
                highlightTrace(curveNumber, false); // false = not locked
            }
        });
        
        // Reset highlighting when not hovering over any trace (only when no trace is locked)
        plotDiv.on('plotly_unhover', function() {
            if (lockedTrace === null) {
                resetOpacity();
            }
        });
        
        // Click event to lock/unlock a trace
        plotDiv.on('plotly_click', function(data) {
            if (data.points && data.points.length > 0) {
                const curveNumber = data.points[0].curveNumber;
                
                // If clicking on already locked trace, unlock it
                if (lockedTrace === curveNumber) {
                    lockedTrace = null;
                    resetOpacity();
                } else {
                    // Otherwise, lock the new trace
                    lockedTrace = curveNumber;
                    highlightTrace(curveNumber, true); // true = locked
                }
            }
        });
        
        // Function to highlight a specific trace
        function highlightTrace(curveNumber, isLocked) {
            const opacityUpdates = [];
            const traceIndices = [];
            const totalTraces = document.querySelectorAll('g.trace').length;
            
            // Set the selected trace to full opacity, all others to low
            for (let i = 0; i < totalTraces; i++) {
                opacityUpdates.push(i === curveNumber ? 1.0 : 0.15);
                traceIndices.push(i);
            }
            
            Plotly.restyle(plotDivId, 'opacity', opacityUpdates, traceIndices);
            
            // If this is a locked trace, also update line width to make it more prominent
            if (isLocked) {
                const lineWidthUpdates = [];
                for (let i = 0; i < totalTraces; i++) {
                    lineWidthUpdates.push(i === curveNumber ? 4 : 2);
                }
                Plotly.restyle(plotDivId, 'line.width', lineWidthUpdates, traceIndices);
            }
        }
        
        // Function to reset all traces to default opacity
        function resetOpacity() {
            const opacityUpdates = [];
            const lineWidthUpdates = [];
            const traceIndices = [];
            const totalTraces = document.querySelectorAll('g.trace').length;
            
            for (let i = 0; i < totalTraces; i++) {
                opacityUpdates.push(0.3);  // Default medium-low opacity
                lineWidthUpdates.push(2);   // Default line width
                traceIndices.push(i);
            }
            
            Plotly.restyle(plotDivId, {
                'opacity': opacityUpdates,
                'line.width': lineWidthUpdates
            }, traceIndices);
        }
        
        // SCALE TOGGLE FUNCTIONALITY
        let isLogScale = false;
        const toggleYButton = document.getElementById('toggleYScale');
        
        if (toggleYButton) {
            toggleYButton.addEventListener('click', function() {
                isLogScale = !isLogScale;
                
                const yAxisUpdate = {
                    'yaxis.type': isLogScale ? 'log' : 'linear'
                };
                
                Plotly.relayout(plotDivId, yAxisUpdate);
                
                // Update button text based on current state
                toggleYButton.textContent = isLogScale 
                    ? 'Switch to Linear Scale' 
                    : 'Switch to Log Scale';
            });
        } else {
            console.error('Y-axis toggle button not found');
        }
        
        // X-AXIS MODE TOGGLE FUNCTIONALITY
        // Get the plot data that was embedded in the HTML
        const plotData = JSON.parse(document.getElementById('plotData').textContent);
        let useAtomCount = false;
        const toggleXButton = document.getElementById('toggleXAxis');
        
        if (toggleXButton) {
            toggleXButton.addEventListener('click', function() {
                useAtomCount = !useAtomCount;
                
                // 1. First, reset any highlighting
                if (lockedTrace !== null) {
                    lockedTrace = null;
                    resetOpacity();
                }
                
                // 2. Create new data array for all traces
                const newData = [];
                
                for (let i = 0; i < plotData.traces.length; i++) {
                    const traceData = plotData.traces[i];
                    
                    // Create new traces with updated x values and enhance properties for atom mode
                    // Use larger line width and much larger marker size in both modes
                    const lineWidth = useAtomCount ? 3 : 3;      // Wider lines in atom count mode
                    const markerSize = useAtomCount ? 14 : 14;   // Large markers in both modes
                    
                    const updatedTrace = {
                        x: useAtomCount ? traceData.atom_x : traceData.nxn_x,
                        y: traceData.y,
                        type: 'scatter',
                        mode: 'lines+markers',
                        name: traceData.name,
                        legendgroup: traceData.legendgroup,
                        showlegend: traceData.showlegend,
                        line: {
                            color: traceData.color,
                            width: lineWidth  // Enhanced line width
                        },
                        marker: {
                            color: traceData.color,
                            size: markerSize,  // Enhanced marker size
                            line: {
                                width: 1,
                                color: 'white'  // Adding white outline for better visibility
                            }
                        },
                        hovertext: traceData.hovertext,
                        hoverinfo: 'text',
                        opacity: 0.3 // default opacity that can be changed by highlighting
                    };
                    
                    newData.push(updatedTrace);
                }
                
                // 3. Create a complete layout update that preserves original styling
                let fullLayout = {};
                
                // First set the background colors
                fullLayout.plot_bgcolor = plotData.original_layout.plot_bgcolor;
                fullLayout.paper_bgcolor = plotData.original_layout.paper_bgcolor;
                
                // Set font properties
                fullLayout.font = plotData.original_layout.font;
                
                // Set margins
                fullLayout.margin = plotData.original_layout.margin;
                
                // Setup xaxis and yaxis with preserved grid properties
                fullLayout.xaxis = {
                    showgrid: plotData.original_layout.xaxis.showgrid,
                    gridcolor: plotData.original_layout.xaxis.gridcolor,
                    gridwidth: plotData.original_layout.xaxis.gridwidth,
                    zeroline: plotData.original_layout.xaxis.zeroline
                };
                
                fullLayout.yaxis = {
                    showgrid: plotData.original_layout.yaxis.showgrid,
                    gridcolor: plotData.original_layout.yaxis.gridcolor,
                    gridwidth: plotData.original_layout.yaxis.gridwidth,
                    zeroline: plotData.original_layout.yaxis.zeroline
                };
                
                // Now apply the mode-specific settings
                if (useAtomCount) {
                    // For atom count mode - set to linear axis
                    fullLayout.xaxis.title = { 
                        text: 'Surface Atoms',
                        font: { size: 18 }  // Larger font for atom count mode
                    };
                    fullLayout.xaxis.type = 'linear';
                    fullLayout.xaxis.autorange = true;
                    fullLayout.xaxis.tickmode = 'auto';
                    fullLayout.xaxis.tickformat = ',d'; // Format as integers with commas
                    fullLayout.xaxis.tickfont = { size: 16 }; // Larger tick font for atom count
                    
                    // Make sure Y-axis title is also displayed with large font
                    fullLayout.yaxis.title = {
                        text: plotData.original_layout.yaxis.title?.text || 'Evaluations per second',
                        font: { size: 18 } // Larger font for y-axis title
                    };
                    fullLayout.yaxis.tickfont = { size: 16 }; // Larger tick font for y-axis
                } else {
                    // For NxN mode - set to categorical with proper ordering but use standard font sizes
                    fullLayout.xaxis.title = { 
                        text: 'Surface Size', 
                        font: { size: 16 }  // Standard size for NxN mode
                    };
                    fullLayout.xaxis.type = 'category';
                    fullLayout.xaxis.tickmode = 'array';
                    fullLayout.xaxis.tickvals = plotData.nxn_layout['xaxis.tickvals'];
                    fullLayout.xaxis.ticktext = plotData.nxn_layout['xaxis.ticktext'];
                    fullLayout.xaxis.categoryorder = 'array';
                    fullLayout.xaxis.categoryarray = plotData.nxn_layout['xaxis.ticktext'];
                    fullLayout.xaxis.tickfont = { size: 14 }; // Standard tick font for NxN mode
                    
                    // Make sure Y-axis title is consistently displayed too
                    fullLayout.yaxis.title = {
                        text: plotData.original_layout.yaxis.title?.text || 'Evaluations per second',
                        font: { size: 16 } // Standard font for y-axis title
                    };
                    fullLayout.yaxis.tickfont = { size: 14 }; // Standard tick font for y-axis
                }
                
                // 4. Completely redraw the plot with new data configuration
                Plotly.react(plotDivId, newData, fullLayout);
                
                // 5. Update button text
                toggleXButton.textContent = useAtomCount
                    ? 'Switch to NxN Notation'
                    : 'Switch to Atom Count';
                
                // 6. Re-register the hover and click events after redrawing
                const updatedPlotDiv = document.getElementById(plotDivId);
                
                // Re-register hover events
                updatedPlotDiv.on('plotly_hover', function(data) {
                    if (lockedTrace === null && data.points && data.points.length > 0) {
                        const curveNumber = data.points[0].curveNumber;
                        highlightTrace(curveNumber, false);
                    }
                });
                
                updatedPlotDiv.on('plotly_unhover', function() {
                    if (lockedTrace === null) {
                        resetOpacity();
                    }
                });
                
                // Re-register click events
                updatedPlotDiv.on('plotly_click', function(data) {
                    if (data.points && data.points.length > 0) {
                        const curveNumber = data.points[0].curveNumber;
                        
                        if (lockedTrace === curveNumber) {
                            lockedTrace = null;
                            resetOpacity();
                        } else {
                            lockedTrace = curveNumber;
                            highlightTrace(curveNumber, true);
                        }
                    }
                });
            });
        } else {
            console.error('X-axis toggle button not found');
        }
    });
</script>
'''

# Insert the control panel after <body> tag
fig_html_with_controls = fig_html.replace(injection_point, injection_point + control_panel_html)

# Add global styles for proper display and scrolling
global_styles = '''
<style>
    body, html {
        margin: 0;
        padding: 0;
        width: 100%;
        height: 100%;
        overflow-x: hidden;
    }
    .plot-container {
        width: 95%;
        max-width: 1200px;
        margin: 20px auto;
        padding: 15px;
    }
</style>
'''

# Insert global styles into head
head_end_tag = '</head>'
fig_html_with_controls = fig_html_with_controls.replace(head_end_tag, global_styles + head_end_tag)

# Wrap the plotly graph div in a container
# Find the plotly graph div which starts with <div id="plotly-"
plot_div_start = fig_html_with_controls.find('<div id="plotly-')
if plot_div_start != -1:
    # Find the end of the div by searching for a </div> after the plot div start
    plot_div_end = fig_html_with_controls.find('</div>', plot_div_start)
    if plot_div_end != -1:
        plot_div_end += 6  # Include the </div> tag itself
        # Extract the plot div
        plot_div = fig_html_with_controls[plot_div_start:plot_div_end]
        # Create the wrapper div with proper styling
        wrapper_start = '<div class="plot-container">'
        wrapper_end = '</div>'
        # Replace the plot div with the wrapped version
        fig_html_with_controls = fig_html_with_controls.replace(plot_div, wrapper_start + plot_div + wrapper_end)

# Insert the script before closing </body> tag
fig_html_with_controls = fig_html_with_controls.replace("</body>", script_injection + "</body>")

# Write the modified HTML to a file with appropriate name
mode_label = "NxN" if not use_atom_count else "AtomCount"
html_output_path = f"interactive_performance_plot_{mode_label}.html"
with open(html_output_path, "w") as f:
    f.write(fig_html_with_controls)

# Save as static image files in various formats
png_output_path = f"performance_plot_{mode_label}.png"
svg_output_path = f"performance_plot_{mode_label}.svg"
pdf_output_path = f"performance_plot_{mode_label}.pdf"
fig.write_image(png_output_path, width=1200, height=800, scale=2)  # High-resolution PNG
fig.write_image(svg_output_path)  # Vector SVG format
fig.write_image(pdf_output_path)  # PDF format

# Open the interactive HTML file directly in the browser
file_path = os.path.abspath(html_output_path)
webbrowser.open('file://' + file_path)

print(f"\nPlot generated with {mode_label} notation x-axis mode")

print("\nExported plot in the following formats:")
print(f"1. Interactive HTML: {html_output_path}")
print(f"2. High-resolution PNG: {png_output_path}")
print(f"3. Vector SVG: {svg_output_path}")
print(f"4. PDF: {pdf_output_path}")

print("\nTo switch between display modes, use:")
print("  python plot_data.py            # For actual surface atom count (default)")
print("  python plot_data.py --nxn      # For NxN notation")