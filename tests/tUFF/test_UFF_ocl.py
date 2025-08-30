#!/usr/bin/env python3
"""
UFF Validation Test Script

Tests UFF force field implementation comparing CPU vs GPU/OpenCL versions.
Allows component-wise testing (bonds, angles, dihedrals, inversions) and
detailed force/energy comparison.

Usage:
    python test_uff_validation.py H2O.mol [--gpu] [--components bonds,angles]
"""

import sys
import os
import numpy as np
import argparse
from pathlib import Path

# Add FireCore to path
sys.path.append(str(Path(__file__).parent.parent))

# Import FireCore modules
try:
    from pyBall import atomicUtils as au
    from pyBall import MMFF as mmff_cpu
    from pyBall.OCL import UFF as uff_ocl
    # TODO: Import UFF_multi when implemented
    # from pyBall import UFF_multi as uff_multi
except ImportError as e:
    print(f"Error importing FireCore modules: {e}")
    print("Make sure PYTHONPATH includes FireCore root directory")
    sys.exit(1)

class UFFValidator:
    """Validates UFF implementation between CPU and GPU versions"""
    
    def __init__(self, molecule_file, use_gpu=False, verbosity=1):
        self.molecule_file = molecule_file
        self.use_gpu = use_gpu
        self.verbosity = verbosity
        
        # Component flags
        self.components = {
            'bonds': True,
            'angles': True, 
            'dihedrals': True,
            'inversions': True
        }
        
        # Tolerance for comparisons
        self.energy_tolerance = 1e-6
        self.force_tolerance = 1e-5
        
        # Results storage
        self.cpu_results = {}
        self.gpu_results = {}
        
    def set_components(self, components):
        """Enable/disable specific UFF components
        
        Args:
            components (list): List of component names to enable
                              ['bonds', 'angles', 'dihedrals', 'inversions']
        """
        # Disable all first
        for comp in self.components:
            self.components[comp] = False
            
        # Enable specified components
        for comp in components:
            if comp in self.components:
                self.components[comp] = True
            else:
                print(f"Warning: Unknown component '{comp}'")
                
        if self.verbosity > 0:
            enabled = [k for k, v in self.components.items() if v]
            print(f"Enabled UFF components: {enabled}")
    
    def init_cpu_uff(self):
        """Initialize CPU UFF calculation"""
        if self.verbosity > 0:
            print("Initializing CPU UFF...")
            
        # Set verbosity for MMFF module
        mmff_cpu.setVerbosity(verbosity=self.verbosity, idebug=1)
        
        # Initialize with UFF parameters
        # Note: This uses MMFF interface but with UFF=True flag
        mmff_cpu.init(
            xyz_name=self.molecule_file,
            surf_name=None,
            smile_name=None,
            sElementTypes="data_UFF/ElementTypes.dat",
            sAtomTypes="data_UFF/AtomTypes.dat", 
            sBondTypes="data_UFF/BondTypes.dat",
            sAngleTypes="data_UFF/AngleTypes.dat",
            sDihedralTypes="data_UFF/DihedralTypes.dat",
            bMMFF=False,  # Disable MMFF
            bEpairs=False,  # No electron pairs
            nPBC=(0,0,0),  # No periodic boundary conditions
            gridStep=0.1,
            bUFF=True,  # Enable UFF
            b141=True,
            bSimple=False,
            bConj=False,
            bCumulene=True,
            bNonBonded=False  # Disable non-covalent interactions
        )
        
        # Set component flags
        # TODO: Add interface to control UFF components in CPU version
        # For now, we'll need to modify the C++ interface
        
        if self.verbosity > 0:
            print("CPU UFF initialized successfully")
    
    def init_gpu_uff(self):
        """Initialize GPU/OpenCL UFF calculation"""
        if self.verbosity > 0:
            print("Initializing GPU UFF...")
            
        # Load molecule for UFF_CL
        mol = au.AtomicSystem()
        mol.loadXYZ(self.molecule_file)
        mol.neighs()  # Generate neighbor information
        
        # Initialize UFF_CL
        self.uff_cl = uff_ocl.UFF_CL(nloc=32)
        
        # Convert molecule to UFF representation
        self.uff_data = self.uff_cl.toUFF(mol)
        
        # Set component flags
        self.uff_cl.set_component_flags(
            bBonds=self.components['bonds'],
            bAngles=self.components['angles'], 
            bDihedrals=self.components['dihedrals'],
            bInversions=self.components['inversions']
        )
        
        # Upload topology and parameters
        self.uff_cl.upload_topology_params(self.uff_data)
        self.uff_cl.upload_positions(self.uff_data['apos'])
        
        if self.verbosity > 0:
            print("GPU UFF initialized successfully")
    
    def run_cpu_calculation(self):
        """Run CPU UFF force calculation"""
        if self.verbosity > 0:
            print("Running CPU UFF calculation...")
            
        # TODO: Add component control interface
        # For now, run full UFF calculation
        
        # Get initial energy and forces
        energy = mmff_cpu.eval()
        forces = mmff_cpu.getForces()
        
        # Store results
        self.cpu_results = {
            'total_energy': energy,
            'forces': np.array(forces).reshape(-1, 3),
            'components': {}  # TODO: Add component breakdown
        }
        
        if self.verbosity > 1:
            print(f"CPU Total Energy: {energy:.8f}")
            print(f"CPU Forces shape: {self.cpu_results['forces'].shape}")
    
    def run_gpu_calculation(self):
        """Run GPU/OpenCL UFF force calculation"""
        if self.verbosity > 0:
            print("Running GPU UFF calculation...")
            
        # Run UFF evaluation
        self.uff_cl.run_eval_step(bClearForce=True)
        
        # Get results
        energy = self.uff_cl.get_total_energy()
        forces = self.uff_cl.get_forces()
        
        # Store results
        self.gpu_results = {
            'total_energy': energy,
            'forces': forces,
            'components': {}  # TODO: Add component breakdown
        }
        
        if self.verbosity > 1:
            print(f"GPU Total Energy: {energy:.8f}")
            print(f"GPU Forces shape: {forces.shape}")
    
    def compare_results(self):
        """Compare CPU and GPU results"""
        if not self.cpu_results or not self.gpu_results:
            print("Error: Missing results for comparison")
            return False
            
        print("\n" + "="*60)
        print("UFF VALIDATION RESULTS")
        print("="*60)
        
        # Energy comparison
        cpu_energy = self.cpu_results['total_energy']
        gpu_energy = self.gpu_results['total_energy']
        energy_diff = abs(cpu_energy - gpu_energy)
        energy_rel_diff = energy_diff / abs(cpu_energy) if cpu_energy != 0 else float('inf')
        
        print(f"Energy Comparison:")
        print(f"  CPU Energy:     {cpu_energy:.8f}")
        print(f"  GPU Energy:     {gpu_energy:.8f}")
        print(f"  Absolute Diff:  {energy_diff:.2e}")
        print(f"  Relative Diff:  {energy_rel_diff:.2e}")
        print(f"  Within Tolerance: {'YES' if energy_diff < self.energy_tolerance else 'NO'}")
        
        # Force comparison
        cpu_forces = self.cpu_results['forces']
        gpu_forces = self.gpu_results['forces']
        
        if cpu_forces.shape != gpu_forces.shape:
            print(f"\nError: Force array shape mismatch!")
            print(f"  CPU shape: {cpu_forces.shape}")
            print(f"  GPU shape: {gpu_forces.shape}")
            return False
            
        force_diff = np.abs(cpu_forces - gpu_forces)
        max_force_diff = np.max(force_diff)
        rms_force_diff = np.sqrt(np.mean(force_diff**2))
        
        print(f"\nForce Comparison:")
        print(f"  Max Force Diff:  {max_force_diff:.2e}")
        print(f"  RMS Force Diff:  {rms_force_diff:.2e}")
        print(f"  Within Tolerance: {'YES' if max_force_diff < self.force_tolerance else 'NO'}")
        
        # Per-atom force analysis
        if self.verbosity > 1:
            print(f"\nPer-Atom Force Analysis:")
            natoms = cpu_forces.shape[0]
            for i in range(natoms):
                atom_diff = np.linalg.norm(force_diff[i])
                print(f"  Atom {i:2d}: |ΔF| = {atom_diff:.2e}")
        
        # Component breakdown (if available)
        if self.cpu_results.get('components') and self.gpu_results.get('components'):
            print(f"\nComponent Breakdown:")
            for comp in ['bonds', 'angles', 'dihedrals', 'inversions']:
                if comp in self.cpu_results['components']:
                    cpu_comp = self.cpu_results['components'][comp]
                    gpu_comp = self.gpu_results['components'][comp]
                    comp_diff = abs(cpu_comp - gpu_comp)
                    print(f"  {comp:12s}: CPU={cpu_comp:.6f}, GPU={gpu_comp:.6f}, Diff={comp_diff:.2e}")
        
        # Overall validation result
        validation_passed = (energy_diff < self.energy_tolerance and 
                           max_force_diff < self.force_tolerance)
        
        print(f"\n{'='*60}")
        print(f"VALIDATION: {'PASSED' if validation_passed else 'FAILED'}")
        print(f"{'='*60}")
        
        return validation_passed
    
    def run_validation(self):
        """Run complete validation test"""
        try:
            # Initialize systems
            if self.use_gpu:
                self.init_gpu_uff()
                self.run_gpu_calculation()
                
                # For comparison, also run CPU (if not GPU-only mode)
                self.init_cpu_uff()
                self.run_cpu_calculation()
                
                return self.compare_results()
            else:
                # CPU-only mode
                self.init_cpu_uff()
                self.run_cpu_calculation()
                
                print(f"\nCPU-only mode results:")
                print(f"Total Energy: {self.cpu_results['total_energy']:.8f}")
                print(f"Force array shape: {self.cpu_results['forces'].shape}")
                return True
                
        except Exception as e:
            print(f"Error during validation: {e}")
            if self.verbosity > 1:
                import traceback
                traceback.print_exc()
            return False

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='UFF Validation Test Script')
    parser.add_argument('molecule', help='Molecule file (.mol, .xyz)')
    parser.add_argument('--gpu', action='store_true', help='Use GPU/OpenCL implementation')
    parser.add_argument('--components', default='bonds,angles,dihedrals,inversions', help='Comma-separated list of UFF components to test')
    parser.add_argument('--verbose', '-v', action='count', default=1, help='Increase verbosity level')
    parser.add_argument('--energy-tol', type=float, default=1e-6, help='Energy comparison tolerance')
    parser.add_argument('--force-tol', type=float, default=1e-5, help='Force comparison tolerance')
    
    args = parser.parse_args()
    
    # Check if molecule file exists
    if not os.path.exists(args.molecule):
        print(f"Error: Molecule file '{args.molecule}' not found")
        sys.exit(1)
    
    # Parse components
    components = [comp.strip() for comp in args.components.split(',')]
    
    # Create validator
    validator = UFFValidator(
        molecule_file=args.molecule,
        use_gpu=args.gpu,
        verbosity=args.verbose
    )
    
    # Set tolerances
    validator.energy_tolerance = args.energy_tol
    validator.force_tolerance = args.force_tol
    
    # Set components
    validator.set_components(components)
    
    # Run validation
    success = validator.run_validation()
    
    sys.exit(0 if success else 1)

