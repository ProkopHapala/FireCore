import numpy as np
import os
import subprocess
from pathlib import Path

# Default mappings as global variables
# Element to atom type mapping
DEFAULT_ATOM_TYPE_MAP = {
    '1': 2,   # Hydrogen
    '8': 1,   # Oxygen
    '92': -1,  # Spin-down electron
    '109': -3,  # Spin-up electron
}

# Atom type to mass mapping
DEFAULT_MASS_MAP = {
    1: 15.99903,  # Oxygen
    2: 1.000794,  # Hydrogen
    3: 1.0        # Electron
}

# Element number to symbol mapping for output parsing
DEFAULT_ELEMENT_MAPPING = {
    8: '8',    # Oxygen
    1: '1',    # Hydrogen
    0: {       # Electrons - map based on spin
        -1: '92',  # Spin down
        1: '109'   # Spin up
    }
}

# Default electron type
DEFAULT_ELECTRON_TYPE = 3

# Common LAMMPS utility functions

def transform_elements_for_eff(elements, atom_type_map=None, electron_type=None):
    """
    Transform elements for eFF simulation, handling electron spin states.
    
    Parameters:
    - elements: list of element symbols/numbers
    - atom_type_map: dictionary mapping element symbols to atom types
    - electron_type: atom type to use for electrons
    
    Returns a dictionary of transformed elements and spin states
    """
    if atom_type_map is None:
        atom_type_map = DEFAULT_ATOM_TYPE_MAP
        
    if electron_type is None:
        electron_type = DEFAULT_ELECTRON_TYPE
    
    transformed = []
    
    for element in elements:
        # Determine element number and initial atom_type
        element_str = element
        element_num = int(element_str) if element_str.isdigit() else 0
        atom_type = atom_type_map.get(element_str, 0)
        # Handle electrons: negative atom_type indicates spin state
        if atom_type < 0:
            spin = -atom_type  # spin encoded as positive integer
            atom_type = electron_type
            element_num = 0
        else:
            spin = 0
        
        transformed.append({
            'element': element_num,
            'atom_type': atom_type,
            'spin': spin
        })
    
    return transformed

def create_lammps_data_eff(elements, positions, charges, radii, box_size=500, mass_map=None, atom_type_map=None, electron_type=None):
    """
    Create LAMMPS data file content for eFF simulations.
    
    Parameters:
    - elements: list of element symbols/numbers
    - positions: list of [x,y,z] coordinates
    - charges: list of charges
    - radii: list of radii (for electrons)
    - box_size: simulation box size
    - mass_map: dictionary mapping atom types to masses
    - atom_type_map: dictionary mapping element symbols to atom types
    - electron_type: atom type to use for electrons
    
    Returns the content as a string
    """
    if mass_map is None:
        mass_map = DEFAULT_MASS_MAP
    
    if atom_type_map is None:
        atom_type_map = DEFAULT_ATOM_TYPE_MAP
        
    if electron_type is None:
        electron_type = DEFAULT_ELECTRON_TYPE
    
    # Transform elements and extract atom types
    atoms = []
    atom_types = set()
    
    for i, (element, pos, q, radius) in enumerate(zip(elements, positions, charges, radii)):
        atom_type = atom_type_map.get(element, 0)
        element_num = int(element) if element.isdigit() else 0
        
        if atom_type < 0:  # This is an electron
            element_num = 0
            spin = -atom_type
            atom_type = electron_type
        else:
            spin = 0
            
        atom_types.add(atom_type)
        atoms.append((atom_type, element_num, spin, radius, pos, q))
    
    # Sort atom types
    atom_types = sorted(atom_types)
    
    # Create data file content using a multi-line string template
    header = f"""LAMMPS data file generated from eFF output

{len(elements)} atoms
{len(atom_types)} atom types

{-box_size} {box_size} xlo xhi
{-box_size} {box_size} ylo yhi
{-box_size} {box_size} zlo zhi

Masses

"""
    
    # Add masses
    masses = ""
    for atype in atom_types:
        masses += f"{atype} {mass_map[atype]}\n"
    
    # Add atoms section
    atoms_section = "\nAtoms # electron\n\n"
    for atom_id, (atype, element, spin, radius, pos, q) in enumerate(atoms, 1):
        x, y, z = pos
        if np.isnan(radius):
            radius = 0.0
        atoms_section += f"{atom_id:5d} {atype:5d} {element:5d} {spin:5d} {radius:10.6f} {x:10.6f} {y:10.6f} {z:10.6f}\n"
    
    # Add velocities section
    velocities = "\nVelocities\n\n"
    for atom_id in range(1, len(elements)+1):
        velocities += f"{atom_id} 0 0 0 0\n"
    
    # Combine all sections
    return header + masses + atoms_section + velocities

def write_lammps_data(filename, content):
    """Write LAMMPS data file content to a file"""
    with open(filename, 'w') as f:
        f.write(content)

def read_lammps_data(filename, element_mapping=None):
    """
    Read geometry from LAMMPS data file
    
    Parameters:
    - filename: path to LAMMPS data file
    - element_mapping: dictionary to map element numbers back to symbols
    
    Returns:
    - elements: list of element symbols
    - positions: list of [x,y,z] coordinates
    - Optional: additional data extracted from file
    """
    if element_mapping is None:
        element_mapping = DEFAULT_ELEMENT_MAPPING
    
    elements = []
    positions = []
    
    with open(filename, 'r') as f:
        lines = f.readlines()
    
    # Find the Atoms section
    atom_section = False
    for i, line in enumerate(lines):
        if "Atoms" in line:
            atom_section = True
            # Skip the header line and any blank lines
            i += 1
            while i < len(lines) and not lines[i].strip():
                i += 1
            break
    
    if not atom_section:
        return None, None
    
    # Read atom data
    atom_data = []
    for j in range(i, len(lines)):
        line = lines[j].strip()
        if not line or "Velocities" in line:
            break
        
        parts = line.split()
        if len(parts) >= 8:
            atom_type = int(parts[1])
            element_num = int(parts[2])
            spin = int(parts[3])
            x, y, z = float(parts[5]), float(parts[6]), float(parts[7])
            
            # Convert atom type/element number back to element symbol
            # Handle electrons (element_num 0) first based on spin
            if element_num == 0 and spin in element_mapping.get(0, {}):
                element = element_mapping[0][spin]
            elif element_num in element_mapping and not isinstance(element_mapping[element_num], dict):
                element = element_mapping[element_num]
            else:
                # Default: just use the element number as string
                element = str(element_num)
                
            elements.append(element)
            positions.append([x, y, z])
            
    return elements, positions

def create_lammps_script(template, **kwargs):
    """
    Create a LAMMPS script from a template string with formatting parameters
    
    Parameters:
    - template: multiline string template with placeholders
    - kwargs: key-value pairs for template formatting
    
    Returns the formatted script as a string
    """
    return template.format(**kwargs)

def run_lammps(script_content, work_dir=".", lammps_executable="lmp_serial"):
    """
    Run LAMMPS as subprocess and capture output
    
    Parameters:
    - script_content: LAMMPS script content
    - work_dir: working directory
    - lammps_executable: path to LAMMPS executable
    
    Returns:
    - stdout: captured standard output
    - stderr: captured standard error
    - returncode: process return code
    """
    script_file = os.path.join(work_dir, "run.lammps")
    with open(script_file, "w") as f:
        f.write(script_content)
    
    # Run LAMMPS and capture output
    cmd = [lammps_executable, "-in", script_file]
    result = subprocess.run(cmd, cwd=work_dir, capture_output=True, text=True)
    
    return result.stdout, result.stderr, result.returncode

def process_xyz_with_lammps(xyz_frames, create_input_fn, process_output_fn, script_template, field_map=None, work_dir=".", box_size=500, max_frames=None, lammps_exec="lmp_serial", **kwargs):
    """
    Generic function to process XYZ file structures with LAMMPS
    
    Parameters:
    - xyz_frames: list of frames (each frame is a list of [elements, positions, ...])
    - create_input_fn: Function to create input data from XYZ frame
      Signature: fn(elements, positions, ...) -> data_content
    - process_output_fn: Function to process output data after simulation
      Signature: fn(output_file, stdout) -> processed_data
    - script_template: LAMMPS script template with placeholders
    - field_map: Dictionary mapping for extract_lammps_values
    - work_dir: Working directory
    - box_size: Simulation box size
    - max_frames: Maximum number of frames to process
    - lammps_exec: LAMMPS executable name/path
    - **kwargs: Additional arguments passed to create_input_fn
    
    Returns: List of result dictionaries for each frame
    """
    # Use provided list of frames
    imgs = xyz_frames
    # Limit number of frames if specified
    if max_frames is not None:
        imgs = imgs[:max_frames]
        print(f"Processing only first {max_frames} frames")
    
    results = []
    
    for i, img in enumerate(imgs):
        #rint(f"Processing frame {i}...")
        # Extract atom data from the frame - handle different frame formats
        if len(img) >= 4:  # Assume pyBall format
            elements = img[0]
            positions = img[1]
            charges = img[2] if len(img) > 2 else np.zeros(len(elements))
            radii = img[3] if len(img) > 3 else np.ones(len(elements))
        else:  # Simple XYZ format
            elements = img[0]
            positions = img[1]
            charges = np.zeros(len(elements))
            radii = np.ones(len(elements))
        
        # Create input and output file paths
        data_file = os.path.join(work_dir, f"frame_{i}.data")
        output_file = os.path.join(work_dir, f"frame_{i}_out.data")
        
        # Create LAMMPS data file content using the provided function
        data_content = create_input_fn(elements, positions, charges, radii, box_size=box_size, **kwargs)
        write_lammps_data(data_file, data_content)
        
        # Create LAMMPS script from template
        script = create_lammps_script(
            template=script_template,
            data_file=os.path.basename(data_file),
            output_file=os.path.basename(output_file)
        )
        
        # Run LAMMPS
        stdout, stderr, returncode = run_lammps(
            script_content=script,
            work_dir=work_dir,
            lammps_executable=lammps_exec
        )
        
        if returncode != 0:
            print(f"Error running LAMMPS for frame {i}:")
            print(stderr)
            continue
        
        # Process output using the provided function
        processed_data = process_output_fn(output_file, stdout)
        
        # Extract energy values from LAMMPS output if field_map provided
        if field_map:
            energies = extract_lammps_values(stdout, field_map)
        else:
            energies = {}
        
        # Store results
        results.append({
            'frame': i,
            'energies': energies,
            'processed_data': processed_data,
            'output': stdout
        })
        
        # Print basic results
        if energies and 'total' in energies:
            print(f"Frame {i}: Total energy = {energies['total']}")
        
    return results

def extract_lammps_values(lammps_output, field_map):
    """
    Extract values from LAMMPS output based on field mapping
    
    Parameters:
    - lammps_output: string output from LAMMPS
    - field_map: dictionary mapping output field names to info about how to extract them
      e.g. {'energy': {'pattern': 'TotEng', 'column': 1}}
    
    Returns: dictionary of extracted values
    """
    lines = lammps_output.split('\n')
    results = {}
    
    # Thermo output extraction - typically in columns with headers
    thermo_headers = None
    thermo_values = None
    
    for i, line in enumerate(lines):
        # Find thermo output headers (e.g. "Step TotEng etc")
        if any(pattern in line for pattern in ["Step", "TotEng", "Temp"]):
            thermo_headers = line.split()
            # Look for values in the next non-empty line
            for j in range(i+1, len(lines)):
                if lines[j].strip() and not lines[j].startswith("Loop"):
                    thermo_values = lines[j].split()
                    break
    
    # If we found thermo data, extract the requested fields
    if thermo_headers and thermo_values:
        for field, info in field_map.items():
            if 'pattern' in info and 'index' in info:
                # Extract by pattern and index
                pattern = info['pattern']
                idx = info['index']
                for h_idx, header in enumerate(thermo_headers):
                    if pattern in header and h_idx < len(thermo_values):
                        try:
                            results[field] = float(thermo_values[h_idx])
                        except (ValueError, IndexError):
                            results[field] = None
            elif 'column' in info and info['column'] < len(thermo_values):
                # Extract by column index directly
                try:
                    results[field] = float(thermo_values[info['column']])
                except (ValueError, IndexError):
                    results[field] = None
    
    # Custom extraction for patterns not in thermo output
    for field, info in field_map.items():
        if field not in results and 'regex' in info:
            import re
            pattern = info['regex']
            matches = re.findall(pattern, lammps_output)
            if matches:
                try:
                    results[field] = float(matches[-1])
                except (ValueError, IndexError):
                    results[field] = matches[-1]
    
    return results
