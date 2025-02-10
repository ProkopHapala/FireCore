import os
import traceback
import glob
import numpy as np
from .atomicUtils import AtomicSystem
from .plotUtils import render_POVray

def run_povray(pov_file, width=800, height=800, quality=9, antialias=0.3):
    """
    Run POVRay to render a .pov file to .png
    
    Args:
        pov_file: Path to .pov file
        width: Image width in pixels
        height: Image height in pixels
        quality: POVRay quality setting (0-9)
        antialias: Antialias amount (0.0-1.0)
    """
    import subprocess
    
    png_file = os.path.splitext(pov_file)[0] + '.png'
    
    cmd = [
        'povray',
        f'+I{pov_file}',      # Input file
        f'+O{png_file}',      # Output file
        f'+W{width}',         # Width
        f'+H{height}',        # Height
        f'+Q{quality}',       # Quality
        f'+A{antialias}',     # Antialias
        '+D',                 # Enable display
        '-V',                 # Verbose output
    ]
    
    try:
        subprocess.run(cmd, check=True, capture_output=True, text=True)
        print(f"Rendered {png_file}")
        return True
    except subprocess.CalledProcessError as e:
        print(f"Error running POVRay: {e.stderr}")
        return False

def render_molecule_files(directory, extension='.mol2', output_dir=None, cam_fw=(0.0,0.0,-1.0), cam_up=(0.0,1.0,0.0), cam_dist=100.0, zoom=20.0):
    """
    Find and render all molecule files of given extension in a directory
    
    Args:
        directory (str): Path to directory containing molecule files
        extension (str): File extension to look for ('.mol2' or '.xyz')
        output_dir (str): Directory for output POV files. If None, uses input directory
    """
    # Ensure directory path ends with /
    directory = os.path.abspath(directory)
    if not directory.endswith('/'):
        directory += '/'
        
    # Set up output directory
    if output_dir is None:
        output_dir = directory
    os.makedirs(output_dir, exist_ok=True)
    
    # Find all files with given extension
    pattern = directory + '*' + extension
    mol_files = glob.glob(pattern)
    
    if not mol_files:
        print(f"No {extension} files found in {directory}")
        return
        
    print(f"Found {len(mol_files)} {extension} files to render")

    cam_fw = np.array(cam_fw)
    cam_up = np.array(cam_up)

    cam_right  = np.cross(cam_fw, cam_up)
    cam_right /= np.linalg.norm(cam_right)*-1.
    cam_up     = np.cross(cam_fw,cam_right)
    light_pos  = ( cam_fw*-1.0 + cam_up ) * cam_dist

    print(f"cam_fw: {cam_fw}")
    print(f"cam_up: {cam_up}")
    
    # Process each file
    for mol_file in mol_files:
        basename = os.path.basename(mol_file)
        name = os.path.splitext(basename)[0]
        print(f"\render_molecule_files() {basename}...")
        
        # Load molecule
        try:
            mol = AtomicSystem(mol_file)
            
            # Set up camera based on molecule size
            pos      = mol.apos
            center   = np.mean(pos, axis=0)
            max_dist = np.max(np.linalg.norm(pos - center, axis=1))
            
            # Render with nice default settings
            pov_file = os.path.join(output_dir, f"{name}.pov")
            render_POVray(mol, pov_file,
                # Adjust view based on molecule size
                camera_pos  =center-cam_fw*cam_dist,
                camera_up   =cam_up,
                camera_right=cam_right,
                look_at=center,
                zoom=zoom,
                
                # Nice material settings
                atom_scale=0.5,
                bond_width=0.1,
                
                # Lighting for good visibility
                light_pos=light_pos,
                light_color=(1.0, 1.0, 1.0),
                #ambient_light=(0.5, 0.5, 0.5),
                
                # Enable shadows and color effects
                #shadows=True,
                #z_color_shift=True
                #viewAxis=True
            )
            # Run POVRay to generate PNG
            if run_povray(pov_file):
                print(f"render_molecule_files() Successfully rendered {pov_file} to PNG")
            else:
                print(f"render_molecule_files() Failed to render {pov_file} to PNG")
            #exit()
        except Exception as e:
            print(f"render_molecule_files() ERROR processing {basename}: {str(e)}")
            traceback.print_exc()
            #exit()
            continue