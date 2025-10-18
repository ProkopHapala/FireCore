import pandas as pd
import re
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
import os

# Funkce pro formátování hodnot osy Y s příponami K, M
def y_fmt(x, pos):
    """Formátovací funkce pro osu Y"""
    if x >= 1e6:
        return f'{x*1e-6:.1f}M'
    elif x >= 1e3:
        return f'{x*1e-3:.1f}K'
    else:
        return f'{x:.0f}'

# Function to convert surface size to numeric value for sorting
# Mapování velikosti mřížky na počet atomů
size_to_atoms = {
    '1': 6,
    '2': 24,
    '3': 54,
    '4': 96,
    '5': 150,
    '6': 216,
    '7': 294,
    '8': 384,
    '9': 486,
    '10': 600,
    '11': 726,
    '12': 864,
    '13': 1014,
    '14': 1176,
    '15': 1350,
    '16': 1536
}

def surface_size_to_numeric(size_str):
    """Convert a surface size string like '5x5' to a numeric value for proper sorting"""
    try:
        if 'x' in size_str:
            dimensions = size_str.split('x')
            if len(dimensions) == 2:
                x, y = map(int, dimensions)
                if x == y:  # Square grid NxN
                    return x
        return int(size_str)  # For simple numeric size
    except (ValueError, TypeError):
        return 0
        
def get_npbc_x(npbc_str):
    """Extract the x component from nPBC string like '1,1,0'"""
    try:
        # Split by commas and get the first value
        return int(npbc_str.split(',')[0])
    except (IndexError, ValueError):
        return 0

def surface_size_to_atoms(size_str):
    """Convert a surface size string like '5x5' to number of atoms"""
    try:
        if 'x' in size_str:
            dimensions = size_str.split('x')
            if len(dimensions) == 2:
                x, y = map(int, dimensions)
                if x == y and str(x) in size_to_atoms:  # Square grid NxN
                    return size_to_atoms[str(x)]
        return 0  # Default for unknown sizes
    except (ValueError, TypeError):
        return 0

def calculate_scaled_atoms(size_str, npbc_str):
    """Calculate the number of atoms based on the formula 6*N**2*(2*npbc.x+1)*(2*npbc.y+1)-1"""
    try:
        # Extract N from NxN pattern
        if 'x' in size_str:
            N = int(size_str.split('x')[0])
        else:
            N = int(size_str)
            
        # Extract npbc.x and npbc.y
        npbc_parts = npbc_str.split(',')
        npbc_x = int(npbc_parts[0]) if len(npbc_parts) > 0 else 0
        npbc_y = int(npbc_parts[1]) if len(npbc_parts) > 1 else 0
        
        # Apply the formula
        return (6 * N**2-1) * (2*npbc_x+1) * (2*npbc_y+1)
    except (ValueError, IndexError, TypeError):
        return 0

# Funkce pro extrakci parametrů z prvního sloupce
def parse_first_column(first_column):
    # Extrakce kódu
    code_match = re.search(r'minima__(\d+)', first_column)
    if code_match:
        code = code_match.group(1)
    else:
        code = 'unknown'
    
    # Extrakce nPBC
    npbc_match = re.search(r'nPBC_\((\d+),(\d+),(\d+)\)', first_column)
    if npbc_match:
        npbc = tuple(map(int, npbc_match.groups()))
        npbc_str = f"{npbc[0]},{npbc[1]},{npbc[2]}"
    else:
        npbc = (0, 0, 0)
        npbc_str = "0,0,0"
    
    # Extrakce MMFF, move, NBFF, surf hodnot
    mmff_match = re.search(r'MMFF_(\d+)', first_column)
    mmff = mmff_match.group(1) if mmff_match else '--'
    
    move_match = re.search(r'move_(\d+)', first_column)
    move = move_match.group(1) if move_match else '--'
    
    nbff_match = re.search(r'NBFF_(\d+|--)', first_column)
    nbff = nbff_match.group(1) if nbff_match else '--'
    
    surf_match = re.search(r'surf_([0-9]+|--)', first_column)
    surf = surf_match.group(1) if surf_match else '--'
    
    # Extrakce velikosti povrchu
    surface_match = re.search(r'surf:_[^_]+_([0-9]+x[0-9]+)_nPBC', first_column)
    surface_size = surface_match.group(1) if surface_match else 'unknown'
    
    gridff_match = re.search(r'gridFF_(\d+|--)', first_column)
    gridff = gridff_match.group(1) if gridff_match else '--'
    
    gridffbspline_match = re.search(r'gridFFbSpline_(\d+|--)', first_column)
    gridffbspline = gridffbspline_match.group(1) if gridffbspline_match else '--'
    
    # Extrakce perframe a perVF
    perframe_match = re.search(r'perframe:_(\d+)', first_column)
    perframe = int(perframe_match.group(1)) if perframe_match else 0
    
    pervf_match = re.search(r'perVF:_(\d+)', first_column)
    pervf = int(pervf_match.group(1)) if pervf_match else 0
    
    # Extrakce replica
    replica_match = re.search(r'replica:_(\d+)', first_column)
    replica = int(replica_match.group(1)) if replica_match else 0
    
    # Extrakce velikosti povrchu
    surface_match = re.search(r'surf:_([^_]+)_([0-9]+x[0-9]+)_nPBC', first_column)
    if surface_match:
        surface_type = surface_match.group(1)  # e.g., NaCl-Cl-hole
        surface_size = surface_match.group(2)  # e.g., 16x16
    else:
        surface_size = 'unknown'
    
    # Vytvoříme parametrickou hash pro účely seskupování
    # Tato hash obsahuje všechny parametry kromě velikosti povrchu
    param_hash = f"{code}_{npbc_str}_{mmff}_{move}_{nbff}_{surf}_{gridff}_{gridffbspline}_{replica}_{perframe}_{pervf}"
    
    return {
        'Code': code,
        'nPBC': npbc_str,  # Použijeme celý nPBC string
        'MMFF': mmff,
        'Move': move,
        'NBFF': nbff,
        'Surf': surf,
        'GridFF': gridff,
        'GridFFbSpline': gridffbspline,
        'Perframe': perframe,
        'PerVF': pervf,
        'Replica': replica,
        'Surface Size': surface_size,
        'Parameter Hash': param_hash  # Hash pro seskupení dat podle všech parametrů
    }

# Read the data from testing_data.dat
with open("uff_data.dat", "r") as f:
    lines = f.readlines()

# Column indices for the test data
column_indices = {
    "Code": 0,  # Extracted from the name
    "Surface Size": 1,  # Extracted from the name
    "nPBC": 2,  # Extracted from the name
    "Nb Iteration": 1,  # First column after the name
    "Time": 4,  # Time column
    "Evaluations per second": 10  # Last column divided by the time
}

# Prepare the data list
data = []

print("Zpracování {0} řádků dat".format(len(lines)))

for line_idx, line in enumerate(lines):
    values = line.strip().split()
    if len(values) < 7:  # Potřebujeme min. 7 hodnot pro výpočet evaluací
        print(f"Příliš málo hodnot na řádku {line_idx+1}: {len(values)} hodnot")
        continue
        
    name = values[0]
    params = parse_first_column(name)
    
    # Create a dictionary for the current row
    row = {
        **params,
        "Nb Iteration": int(values[1]) if len(values) > 1 else 0,
        "Time": float(values[4]) if len(values) > 4 else 0,
    }
    
    # Calculate and add 'Evaluations per second'
    try:
        if float(values[4]) != 0:
            row['Evaluations per second'] = float(values[6]) / float(values[4])
        else:
            row['Evaluations per second'] = 0
    except (IndexError, ValueError):
        row['Evaluations per second'] = 0
    
    # Přidat řádek do seznamu dat - OPRAVENO ODSAZENÍ
    data.append(row)
    
    # Debugovací výpis pro první a poslední řádek
    if line_idx < 3 or line_idx >= len(lines) - 3:
        print(f"Zpracován řádek {line_idx+1}, parametry: kód={params['Code']}, povrch={params['Surface Size']}, hash={params['Parameter Hash']}")

print(f"Celkem zpracováno {len(data)} řádků dat")

# Create the Pandas DataFrame
df = pd.DataFrame(data)

# Debug výpis pro kontrolu existujících sloupců
print("Sloupce v DataFrame:")
print(df.columns.tolist())
print("První řádek dat:")
print(df.iloc[0] if not df.empty else "DataFrame je prázdný")

# Add numeric columns for surface size and atom count
df['Surface Size Numeric'] = df['Surface Size'].apply(surface_size_to_numeric)
df['Atoms Count'] = df['Surface Size'].apply(surface_size_to_atoms)

# Calculate Number of surface atoms for each row based on its nPBC value
df['Number of surface atoms'] = df.apply(lambda row: calculate_scaled_atoms(row['Surface Size'], str(row['nPBC'])), axis=1)

# Create a dict to map codes to labels
labels = {
    '10000': 'Without surface',
    '11000': 'No GridFF',
    '11110': 'GridFF with B-Spline',
}

# Create a dict to map codes to colors
color_map = {
    '10000': 'red',
    '11000': 'green',
    '11110': 'blue',
}

# Keep track of which codes have already been added to legend
legend_added = set()

# Create a matplotlib figure
plt.figure(figsize=(20, 20))

# Seskupit data podle všech parametrů kromě velikosti povrchu
for code in sorted(df['Code'].unique()):
    if str(code) in labels:
        # Získej data jen pro tento kód
        code_data = df[df['Code'] == code]
        
        # Seskupit podle parametrického hash - každý hash je unikátní kombinace všech parametrů
        param_groups = {}
        
        # Vytvoříme slovník skupin dat podle Parameter Hash
        for hash_value, group_data in code_data.groupby('Parameter Hash'):
            param_groups[hash_value] = group_data
        
        # Získáme seznam všech hash hodnot pro tento kód
        param_hashes = list(param_groups.keys())
        
        # Pro každou kombinaci parametrů zjistíme výkon na největší velikosti povrchu
        hash_performance = {}
        for hash_value in param_hashes:
            group_data = param_groups[hash_value]
            
            # Seřadit podle velikosti povrchu sestupně a vzít první řádek (největší povrch)
            if not group_data.empty:
                max_size_row = group_data.sort_values('Surface Size Numeric', ascending=False).iloc[0]
                hash_performance[hash_value] = max_size_row['Evaluations per second']
        
        # Seřadit kombinace parametrů podle výkonu (sestupně)
        sorted_hashes = sorted(hash_performance.items(), key=lambda x: x[1], reverse=True)
        
        # Vykreslit čáru pro každou kombinaci parametrů
        for i, (hash_value, performance) in enumerate(sorted_hashes):
            group_data = param_groups[hash_value]
            sorted_data = group_data.sort_values('Surface Size Numeric')
            
            # Nastavení stylu čáry - první je lepší (plná), ostatní jsou horší (čárkované)
            if i == 0:  # Nejlepší výkon
                line_style = '-'
                line_width = 4
                
                # Přidat popisek do legendy pouze pro nejlepší výkon každého kódu
                label = None
                if str(code) not in legend_added:
                    label = labels[str(code)]
                    legend_added.add(str(code))
            else:  # Horší výkon
                line_style = '--'
                line_width = 3
                label = None  # Bez popisku pro horší výkon
            
            # Vykreslit čáru
            plt.plot(
                sorted_data['Number of surface atoms'],
                sorted_data['Evaluations per second'],
                linestyle=line_style,
                linewidth=line_width,
                marker='o',
                markersize=10,
                color=color_map[str(code)],
                label=label
            )
            
            # Debug info - vypsat parametry pro tuto datovou řadu
            sample_row = group_data.iloc[0]
            print(f"Code: {code}, Performance: {performance:.2f}, Params: MMFF={sample_row['MMFF']}, Move={sample_row['Move']}, "
                  f"NBFF={sample_row['NBFF']}, Surf={sample_row['Surf']}, GridFF={sample_row['GridFF']}, "
                  f"GridFFbSpline={sample_row['GridFFbSpline']}, Replica={sample_row['Replica']}, "
                  f"Perframe={sample_row['Perframe']}, PerVF={sample_row['PerVF']}")

# Použijeme skutečně vykreslené hodnoty z dataframu pro popisky osy X
# Seřadím hodnoty podle škálovaného počtu atomů pro zajištění správného pořadí na ose x

# Získáme unikátní seřazené hodnoty škálovaného počtu atomů
sorted_atom_counts = sorted(df['Number of surface atoms'].unique())

# Funkce pro vykreslení grafu s různými typy os
def create_plot(log_x=False, log_y=False, uniform_ticks=False, nanoseconds=False):
    # Plot type - either "Evaluations per second" or "Nanoseconds per evaluation"
    # Vytvoříme nový graf
    plt.figure(figsize=(14, 10))
    
    # Determine which metric to use
    metric = "Nanoseconds per evaluation" if nanoseconds else "Evaluations per second"
    
    # Nastavíme větší velikost písma pro celý graf
    plt.rcParams.update({'font.size': 14})
    
    # Vykreslení dat pro každý kód
    legend_added = set()

    
    # Seskupit data podle všech parametrů kromě velikosti povrchu
    for code in sorted(df['Code'].unique()):
        if str(code) in labels:
            # Získej data jen pro tento kód
            code_data = df[df['Code'] == code]
            
            # Seskupit podle parametrického hash - každý hash je unikátní kombinace všech parametrů
            param_groups = {}
            
            # Vytvoříme slovník skupin dat podle Parameter Hash
            for hash_value, group_data in code_data.groupby('Parameter Hash'):
                param_groups[hash_value] = group_data
            
            # Získáme seznam všech hash hodnot pro tento kód
            param_hashes = list(param_groups.keys())
            
            # Pro každou kombinaci parametrů zjistíme výkon na největší velikosti povrchu
            hash_performance = {}
            for hash_value in param_hashes:
                group_data = param_groups[hash_value]
                
                # Seřadit podle velikosti povrchu sestupně a vzít první řádek (největší povrch)
                if not group_data.empty:
                    max_size_row = group_data.sort_values('Surface Size Numeric', ascending=False).iloc[0]
                    hash_performance[hash_value] = max_size_row['Evaluations per second']
            
            # Seřadit kombinace parametrů podle výkonu (sestupně)
            sorted_hashes = sorted(hash_performance.items(), key=lambda x: x[1], reverse=True)
            
            # Vykreslit čáru pro každou kombinaci parametrů
            for i, (hash_value, performance) in enumerate(sorted_hashes):
                group_data = param_groups[hash_value]
                sorted_data = group_data.sort_values('Surface Size Numeric')
                
                # Nastavení stylu čáry - první je lepší (plná), ostatní jsou horší (čárkované)
                if i == 0:  # Nejlepší výkon
                    line_style = '-'
                    line_width = 2.5
                    
                    # Přidat popisek do legendy pouze pro nejlepší výkon každého kódu
                    label = None
                    if str(code) not in legend_added:
                        label = labels[str(code)]
                        legend_added.add(str(code))
                else:  # Horší výkon
                    line_style = '--'
                    line_width = 1.5
                    label = None  # Bez popisku pro horší výkon
                
                if code == "10000":
                    # Získáme min a max hodnotu počtu atomů pro definici rozsahu vodorovných čar
                    sorted_atom_counts = sorted(df['Number of surface atoms'].unique())
                    xmin = sorted_atom_counts[1] if sorted_atom_counts else 0
                    xmax = max(sorted_atom_counts) if sorted_atom_counts else 10000
                    
                    # Vypočítat hodnotu y podle zvoleného metrika
                    if nanoseconds:
                        # Nanoseconds per evaluation = 10^9 / (Evaluations per second)
                        y_values = [1e9 / val for val in sorted_data['Evaluations per second']]
                    else:
                        y_values = sorted_data['Evaluations per second']
                    
                    # Vykreslit vodorovnou čáru
                    plt.hlines(
                        y=y_values,
                        xmin=xmin,
                        xmax=xmax,
                        color=color_map[str(code)],
                        linestyle=line_style,
                        linewidth=line_width,
                        label=label
                    )
                    
                    if nanoseconds:
                        print(f"Vykreslena vodorovná čára pro code 10000: Performance: {1e9/performance:.2f} ns/eval, Rozsah: {xmin}-{xmax}")
                    else:
                        print(f"Vykreslena vodorovná čára pro code 10000: Performance: {performance:.2f} eval/s, Rozsah: {xmin}-{xmax}")
                else:
                    # Vykreslit čáru - pro obě metriky
                    if nanoseconds:
                        # Nanoseconds per evaluation = 10^9 / (Evaluations per second)
                        y_values = [1e9 / val for val in sorted_data['Evaluations per second']]
                    else:
                        y_values = sorted_data['Evaluations per second']
                    
                    plt.plot(
                        sorted_data['Number of surface atoms'],
                        y_values,
                        linestyle=line_style,
                        linewidth=line_width,
                        marker='o',
                        markersize=10,
                        color=color_map[str(code)],
                        label=label
                    )
                    
                    if nanoseconds:
                        print(f"Vykreslena čára pro code {code}: (nanoseconds per evaluation)")
                    else:
                        print(f"Vykreslena čára pro code {code}: (evaluations per second)")
    
    # Získáme unikátní seřazené hodnoty škálovaného počtu atomů
    sorted_atom_counts = sorted(df['Number of surface atoms'].unique())
    
    if uniform_ticks:
        # Vytvoření rovnoměrných popisek na ose X začínajících od nuly
        max_val = max(sorted_atom_counts)
        # Vytvoříme 10 rovnoměrně rozdělených hodnot od 0 do max_val
        tick_positions = np.linspace(0, max_val, 10)
        plt.xticks(tick_positions, [str(int(val)) for val in tick_positions], rotation=45, fontsize=20)
        # Zvýrazníme mřížku
        plt.grid(True, linestyle='-', alpha=0.7, which='both')
    else:
        # Standardní popisky - každý datový bod
        plt.xticks(sorted_atom_counts, [str(int(val)) for val in sorted_atom_counts], rotation=45, fontsize=20)
    
    plt.yticks(fontsize=20)
    
    # Nastavíme popisky os s větší velikostí písma
    x_label = 'Number of surface atoms'
    y_label = metric  # Use the selected metric name
    
    plt.xlabel(x_label, fontsize=20)
    plt.ylabel(y_label, fontsize=20)
    
    # Použití formátovací funkce pro osu Y
    plt.gca().yaxis.set_major_formatter(ticker.FuncFormatter(y_fmt))
    
    # Přidat mřížku pouze pokud nebylo již nastaveno v uniform_ticks
    if not uniform_ticks:
        plt.grid(True, linestyle='--', alpha=0.7)
    
    # Nastavení logaritmických os podle parametrů
    if log_x:
        plt.xscale('log')
    if log_y:
        plt.yscale('log')
    
    # Upravit velikost a umístění legendy pro lepší čitelnost
    plt.legend(loc='upper right', fontsize=20, framealpha=1.0)
    
    # Create images directory if it doesn't exist
    images_dir = 'images'
    if not os.path.exists(images_dir):
        os.makedirs(images_dir)
    
    # Vytvoření vhodného názvu souboru
    x_type = "log_x" if log_x else "lin_x"
    y_type = "log_y" if log_y else "lin_y"
    tick_type = "uniform" if uniform_ticks else "data"
    metric_type = "nanoseconds" if nanoseconds else "evaluations"
    filename = f"performance_{metric_type}_{x_type}_{y_type}_{tick_type}.png"
    
    # Uložení grafu do souboru v images adresáři
    filepath = os.path.join(images_dir, filename)
    plt.savefig(filepath, dpi=300, bbox_inches='tight')
    
    return filename

# Export data to a file with performance metrics vs number of atoms
def export_performance_data():
    # Create a new DataFrame with just the columns we need
    export_df = df[['Code', 'Surface Size', 'Number of surface atoms', 'Evaluations per second',
                   'Parameter Hash', 'Replica', 'Perframe', 'PerVF']].copy()
    
    # Add nanoseconds per evaluation column
    export_df['Nanoseconds per evaluation'] = 1e9 / export_df['Evaluations per second']
    
    # Map code numbers to algorithm types (based on execution output and parameters)
    code_map = {
        '10000': 'No Surface',    # No NBFF, No Surf
        '11000': 'No GridFF',     # Has NBFF+Surf, No GridFF
        '11110': 'GridFF'         # Has GridFFbSpline
    }
    
    # Convert code column to integer to ensure matching
    # export_df['Code'] = export_df['Code'].astype(int)
    
    # Find best and worst performers for each algorithm type
    # First, get the baseline No Surface values (code 10000)
    no_surface_data = export_df[export_df['Code'] == '10000']
    no_surface_best = {
        'evals/sec': no_surface_data['Evaluations per second'].max() if not no_surface_data.empty else 0, 
        'ns/eval': (1e9 / no_surface_data['Evaluations per second'].max()) if not no_surface_data.empty else 0
    }
    no_surface_worst = {
        'evals/sec': no_surface_data['Evaluations per second'].min() if not no_surface_data.empty else 0, 
        'ns/eval': (1e9 / no_surface_data['Evaluations per second'].min()) if not no_surface_data.empty else 0
    }
    
    # Process data by surface size/atom count
    performance_data = {}
    
    # For each atom count, find the best and worst performance for each algorithm type
    for atom_count in sorted(export_df['Number of surface atoms'].unique()):
        atom_data = export_df[export_df['Number of surface atoms'] == atom_count]
        
        # Get the surface size for this atom count (should be the same for all rows)
        surface_size = atom_data['Surface Size'].iloc[0]
        
        # Initialize data structure
        performance_data[atom_count] = {
            'Surface Size': surface_size,
            'No Surface': {'best': no_surface_best, 'worst': no_surface_worst},
            'No GridFF': {'best': {'evals/sec': 0, 'ns/eval': 0}, 'worst': {'evals/sec': 0, 'ns/eval': 0}},
            'GridFF': {'best': {'evals/sec': 0, 'ns/eval': 0}, 'worst': {'evals/sec': 0, 'ns/eval': 0}}
        }
        
        # Find best and worst performers for the other code types
        for code in ['11000', '11110']:  # No GridFF and GridFF
            algo_type = code_map[code]
            code_data = atom_data[atom_data['Code'] == code]
            
            if not code_data.empty:
                # Find best performer (max evaluations per second)
                best_idx = code_data['Evaluations per second'].idxmax()
                best_row = code_data.loc[best_idx]
                performance_data[atom_count][algo_type]['best']['evals/sec'] = best_row['Evaluations per second']
                performance_data[atom_count][algo_type]['best']['ns/eval'] = best_row['Nanoseconds per evaluation']
                
                # Find worst performer (min evaluations per second)
                worst_idx = code_data['Evaluations per second'].idxmin()
                worst_row = code_data.loc[worst_idx]
                performance_data[atom_count][algo_type]['worst']['evals/sec'] = worst_row['Evaluations per second']
                performance_data[atom_count][algo_type]['worst']['ns/eval'] = worst_row['Nanoseconds per evaluation']

    
    # Write to file with header - matching the manual format
    with open('performance_vs_number_atoms.dat', 'w') as f:
        # Tab-separated header row matching manual format
        f.write('Number of atoms\tSize\t')
        f.write('Performance – No surface – best\tPerformance – No GridFF – best\tPerformance – GridFF – best\t')
        f.write('Performance – No surface – bad\tPerformance – No GridFF – bad\tPerformance – GridFF – bad\t\t')
        f.write('ns per eval – No surface – best\tns per eval – No GridFF – best\tns per eval – GridFF – best\t')
        f.write('ns per eval – No surface – bad\tns per eval – No GridFF – bad\tns per eval – GridFF – bad\n')
        
        # Data rows - sorted by atom count
        for atom_count in sorted(performance_data.keys()):
            # Skip the entry with atom_count = -1 (baseline) if it exists
            if atom_count < 0:
                continue
                
            data = performance_data[atom_count]
            
            # Format with tabs between columns
            f.write(f"{atom_count}\t{data['Surface Size']}\t")
            
            # Performance (evaluations per second) columns
            f.write(f"{data['No Surface']['best']['evals/sec']:.2f}\t")
            f.write(f"{data['No GridFF']['best']['evals/sec']:.2f}\t")
            f.write(f"{data['GridFF']['best']['evals/sec']:.2f}\t")
            
            f.write(f"{data['No Surface']['worst']['evals/sec']:.2f}\t")
            f.write(f"{data['No GridFF']['worst']['evals/sec']:.2f}\t")
            f.write(f"{data['GridFF']['worst']['evals/sec']:.2f}\t\t")
            
            # Nanoseconds per evaluation columns
            f.write(f"{data['No Surface']['best']['ns/eval']:.2f}\t")
            f.write(f"{data['No GridFF']['best']['ns/eval']:.2f}\t")
            f.write(f"{data['GridFF']['best']['ns/eval']:.2f}\t")
            
            f.write(f"{data['No Surface']['worst']['ns/eval']:.2f}\t")
            f.write(f"{data['No GridFF']['worst']['ns/eval']:.2f}\t")
            f.write(f"{data['GridFF']['worst']['ns/eval']:.2f}\n")
    
    print(f"\nExported performance data to: performance_vs_number_atoms.dat")
    return 'performance_vs_number_atoms.dat'

# Vytvoření všech šesti variant grafů
filenames = [
    # Original "Evaluations per second" plots
    create_plot(log_x=False, log_y=False, uniform_ticks=False, nanoseconds=False),
    create_plot(log_x=True, log_y=False, uniform_ticks=False, nanoseconds=False),
    create_plot(log_x=False, log_y=True, uniform_ticks=False, nanoseconds=False),
    create_plot(log_x=True, log_y=True, uniform_ticks=False, nanoseconds=False),
    create_plot(log_x=False, log_y=False, uniform_ticks=True, nanoseconds=False),
    create_plot(log_x=False, log_y=True, uniform_ticks=True, nanoseconds=False),
    
    # New "Nanoseconds per evaluation" plots
    create_plot(log_x=False, log_y=False, uniform_ticks=False, nanoseconds=True),
    create_plot(log_x=True, log_y=False, uniform_ticks=False, nanoseconds=True),
    create_plot(log_x=False, log_y=True, uniform_ticks=False, nanoseconds=True),
    create_plot(log_x=True, log_y=True, uniform_ticks=False, nanoseconds=True),
    create_plot(log_x=False, log_y=False, uniform_ticks=True, nanoseconds=True),
    create_plot(log_x=False, log_y=True, uniform_ticks=True, nanoseconds=True)
]

# Export performance data to file
export_performance_data()

# Výpis informace o vytvorených souborech
print("\nExported plots to current directory:")
for filename in filenames:
    print(f"- {filename}")

# plt.show()
