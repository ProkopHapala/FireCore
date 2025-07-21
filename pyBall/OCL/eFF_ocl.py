import os
import re
import numpy as np
import pyopencl as cl
from .OpenCLBase import OpenCLBase

class EFF_OCL(OpenCLBase):
    """
    PyOpenCL interface for running Electron Force Field (eFF) relaxations.
    
    This class handles parsing of multi-geometry XYZ files, setting up GPU buffers
    for batch processing, running the local memory eFF kernel, and retrieving
    the results. It separates buffer allocation from data upload for flexibility.
    """
    
    ATOM_PARAMS = {
    #                   Z_nuc , R_eff    , Zcore_eff , PA        , PB       , PC,        PD        ,  PE  
        'H': np.array([ 1.0   , 0.0      , 0.0       , 0.0       , 0.0      , 0.0      , 0.0       , 0.0      ], dtype=np.float32),
        'C': np.array([ 6.0   , 0.621427 , 2.0       , 22.721015 , 0.728733 , 1.103199 , 17.695345 , 6.693621 ], dtype=np.float32),
        'O': np.array([ 8.0   , 0.167813 , 2.0       , 25.080199 , 0.331574 , 1.276183 , 12.910142 , 3.189333 ], dtype=np.float32),
    }
    ELEMENT_Z = {name: params[0] for name, params in ATOM_PARAMS.items()}

    # def __init__(self, nloc=32, kernel_path="../../cpp/common_resources/cl/eFF.cl"):
    #     super().__init__(nloc=nloc)
    #     self.load_program(kernel_path=kernel_path, bPrint=True)
    #     self.systems = []

    def __init__(self, nloc: int = 32, device_index: int = 0):
        super().__init__(nloc=nloc, device_index=device_index)

        # Build the OpenCL program ------------------------------------------------
        base_path = os.path.dirname(os.path.abspath(__file__))
        rel_path  = "../../cpp/common_resources/cl/eFF.cl"
        if not self.load_program(rel_path=rel_path, base_path=base_path):
            raise RuntimeError("[EFF_OCL] Failed to load/compile eFF.cl")

    def load_xyzs(self, xyz_path, frozen_core_size=0.2):
        """
        Parses a multi-geometry XYZ file and prepares system definitions.
        Handles different core_modes and electron notations.
        """
        print(f"--- Loading systems from: {xyz_path}")
        self.systems = []
        with open(xyz_path, 'r') as f: lines = f.readlines()
        
        i = 0
        sys_count = 0
        while i < len(lines):
            try:
                ntot = int(lines[i].strip())
                comment_line = lines[i+1].strip()
                match = re.search(r"na,ne,core\s+(\d+)\s+(\d+)\s+(\w)", comment_line)
                if not match:
                    i += ntot + 2
                    continue
                _, _, core_mode = match.groups()
                ions, electrons = [], []
                for line in lines[i+2 : i+2+ntot]:
                    parts = line.split()
                    name, pos = parts[0], np.array(parts[1:4], dtype=np.float32)
                    if name in self.ELEMENT_Z: ions.append({'name': name, 'pos': pos})
                    elif name.startswith('e'):
                        size = float(parts[5]) if len(parts) > 5 else 1.0
                        if name == 'e2':
                            electrons.extend([{'pos': pos, 'size': size, 'spin': 1}, {'pos': pos, 'size': size, 'spin': -1}])
                        elif name == 'e+': electrons.append({'pos': pos, 'size': size, 'spin': 1})
                        elif name == 'e-': electrons.append({'pos': pos, 'size': size, 'spin': -1})
                if core_mode == 'f':
                    for ion in ions:
                        electrons.extend([{'pos': ion['pos'], 'size': frozen_core_size, 'spin': 1}, {'pos': ion['pos'], 'size': frozen_core_size, 'spin': -1}])
                system_def = {'na': len(ions), 'ne': len(electrons), 'ions': ions, 'electrons': electrons, 'core_mode': core_mode}
                if (system_def['na'] + system_def['ne']) > self.nloc:
                    print(f"Warning: System {sys_count} with {system_def['na']}+{system_def['ne']} particles exceeds workgroup size {self.nloc}. Skipping.")
                else:
                    self.systems.append(system_def)
                    sys_count += 1
                i += ntot + 2
            except (ValueError, IndexError) as e:
                print(f"Warning: Parsing error at line {i}. Stopping. Error: {e}"); 
                print("line: ", lines[i])
                break
        
        print(f"--- Successfully loaded and processed {len(self.systems)} systems.")

    def realloc_buffers(self):
        """
        (Re-)allocates GPU buffers based on the number of loaded systems.
        This separates allocation from data uploading.
        """
        if not self.systems:
            print("No systems loaded. Cannot allocate buffers."); return

        n_systems = len(self.systems)
        total_slots = n_systems * self.nloc
        
        print(f"Reallocating buffers for {n_systems} systems ({total_slots} total particle slots).")

        # Define buffers and their element sizes in bytes
        # (float4=16, float8=32, signed char=1)
        buffer_specs = {
            'g_pos_in': 16, 'g_vel_in': 16,
            'g_aparams': 32, 'g_spins': 1,
            'g_pos_out': 16, 'g_vel_out': 16,
        }
        
        for name, element_size in buffer_specs.items():
            self.try_make_buff(f"{name}_buff", total_slots * element_size)
            
        print("--- GPU buffers allocated.")

    def upload_data(self):
        """
        Populates host arrays with system data and uploads them to the GPU.
        Assumes realloc_buffers() has already been called.
        """
        if not self.systems:
            print("No systems loaded. Cannot upload data."); return
        
        n_systems = len(self.systems)
        total_slots = n_systems * self.nloc

        self.pos_h     = np.zeros((total_slots, 4), dtype=np.float32)
        self.vel_h     = np.zeros_like(self.pos_h)
        self.aparams_h = np.zeros((total_slots, 8), dtype=np.float32)
        self.spins_h   = np.zeros(total_slots, dtype=np.int8)
        
        for i, sys in enumerate(self.systems):
            offset = i * self.nloc
            for j, ion in enumerate(sys['ions']):
                idx = offset + j
                self.pos_h[idx, :3] = ion['pos']
                self.aparams_h[idx, :] = self.ATOM_PARAMS[ion['name']]
                if sys['core_mode'] == 'a': self.aparams_h[idx, 1:] = 0.0

            for j, ele in enumerate(sys['electrons']):
                idx = offset + sys['na'] + j
                self.pos_h[idx, :3] = ele['pos']
                self.pos_h[idx, 3]  = ele['size']
                self.spins_h[idx]   = ele['spin']

        self.toGPU_(self.g_pos_in_buff, self.pos_h)
        self.toGPU_(self.g_vel_in_buff, self.vel_h)
        self.toGPU_(self.g_aparams_buff, self.aparams_h)
        self.toGPU_(self.g_spins_buff, self.spins_h)
        self.queue.finish()
        print("--- Host data uploaded to GPU.")

    def relax_systems(self, n_steps=100, dt=0.1, damping=0.5):
        """
        Executes the eFF relaxation kernel on the GPU for all loaded systems.
        """
        na = np.int32(self.systems[0]['na'])
        ne = np.int32(self.systems[0]['ne'])
        
        kernel = self.prg.localMD
        krsrho = np.array([1.125, 0.9, -0.2, 1.0], dtype=np.float32)

        kernel.set_args(
            self.g_pos_in_buff,  
            self.g_vel_in_buff, 
            self.g_aparams_buff, 
            self.g_spins_buff,
            self.g_pos_out_buff, 
            self.g_vel_out_buff,
            na, ne, np.int32(n_steps), 
            np.float32(dt), 
            np.float32(damping), 
            krsrho
        )

        global_size = (len(self.systems) * self.nloc,)
        local_size = (self.nloc,)
        
        print(f"--- Running relaxation for {n_steps} steps...")
        cl.enqueue_nd_range_kernel(self.queue, kernel, global_size, local_size).wait()
        #print("--- Relaxation finished.")

        pos_out = np.empty_like(self.pos_h)
        self.fromGPU_(self.g_pos_out_buff, pos_out)
        self.queue.finish()
        return pos_out

# =================
#    MAIN USAGE
# =================
if __name__ == "__main__":
    xyz_content = """7
na,ne,core 3 4 f | A water molecule with frozen core oxygen
  O   -0.054841   0.000000  -0.117854 
  H    0.887105   0.000000   0.683742 
  H   -0.880789   0.000000   0.674609 
 e2     0.682333   0.000000   0.512430   2.0   1.0 
 e2    -0.670364   0.000000   0.512427   2.0   1.0 
 e2     0.009820   0.519605  -0.319574   2.0   1.0
 e2     0.009820  -0.519605  -0.319574   2.0   1.0 
5
na,ne,core 2 3 a | An all-electron system (H2+)
  H   0.0   0.0   0.0
  H   0.8   0.0   0.0
  e+  0.4   0.0   0.0   1.0  1.2
  e+  0.4   0.5   0.0   1.0  1.2
  e-  0.4  -0.5   0.0  -1.0  1.2
"""
    with open("test_systems.xyz", "w") as f:
        f.write(xyz_content)

    # 1. Initialize
    eff = EFF_OCL()
    
    # 2. Load systems
    eff.load_xyzs("test_systems.xyz")
    
    # 3. Allocate GPU memory
    eff.realloc_buffers()
    
    # 4. Upload initial data
    eff.upload_data()
    
    # 5. Run the relaxation
    pos_out = eff.relax_systems(n_steps=500, dt=0.05, damping=0.2)

    print( "pos_out ", pos_out )
    
    # 6. Get results
    #final_geometries = eff_runner.get_results()
    
    # print("\n--- Final Geometries ---")
    # for i, geom in enumerate(final_geometries):
    #     #print(f"\nSystem {i} (na={eff.systems[i]['na']}, ne={eff.systems[i]['ne']}, mode='{eff.systems[i]['core_mode']}'):")
    #     print("Final positions (x, y, z) and electron size (w):")
    #     print(geom)