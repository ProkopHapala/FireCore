#!/usr/bin/env python3
"""
Surface Potential Fitting Framework
===================================

Efficiently represent non-covalent interactions between molecules and crystal surfaces
by fitting periodic potentials to localized basis functions.

Approach:
1. Calculate periodic Coulomb/Morse potentials on 2D grid (expensive, done once)
2. Fit to compact basis set in unit cell (cheap linear algebra)
3. Fast evaluation using fitted coefficients

Author: Your Name
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize
from scipy.linalg import lstsq
import warnings


class GridManager:
    """Manages 2D grid definition and coordinate transformations."""
    
    def __init__(self, cell_x=10.0, cell_z=10.0, grid_step=0.2, z_offset=2.0):
        """
        Initialize 2D grid within unit cell.
        
        Parameters:
        -----------
        cell_x : float
            Unit cell size in x-direction (Angstrom), periodic
        cell_z : float  
            Unit cell size in z-direction (Angstrom), non-periodic
        grid_step : float
            Grid spacing in both directions (Angstrom)
        z_offset : float
            Offset above surface atoms (Angstrom) - focuses on non-covalent interaction region
        """
        self.cell_x = cell_x
        self.cell_z = cell_z
        self.grid_step = grid_step
        self.z_offset = z_offset
        
        # Create grid points
        self.nx = int(cell_x / grid_step)
        self.nz = int(cell_z / grid_step)
        
        # Grid coordinates (atoms at z=0, grid starts at z_offset above surface)
        self.x_coords = np.linspace(0, cell_x, self.nx, endpoint=False)
        self.z_coords = np.linspace(z_offset, z_offset + cell_z, self.nz)  # Start at z_offset above atoms
        
        # 2D mesh for vectorized operations
        self.X, self.Z = np.meshgrid(self.x_coords, self.z_coords, indexing='ij')
        
        print(f"Grid: {self.nx}x{self.nz} points, total {self.nx*self.nz}")
        print(f"X range: [0, {cell_x:.1f}] Å")
        print(f"Z range: [{z_offset:.1f}, {z_offset + cell_z:.1f}] Å")
    
    def get_grid_shape(self):
        """Returns (nx, nz) grid dimensions."""
        return (self.nx, self.nz)
    
    def get_coordinates(self):
        """Returns X, Z meshgrid coordinates."""
        return self.X, self.Z
    
    def flatten_grid(self, grid_2d):
        """Flatten 2D grid for linear algebra operations."""
        return grid_2d.flatten()
    
    def reshape_grid(self, flat_array):
        """Reshape flattened array back to 2D grid."""
        return flat_array.reshape(self.nx, self.nz)


class PotentialCalculator:
    """Calculate Coulomb and Morse potentials with periodic boundary conditions."""
    
    def __init__(self, grid_manager):
        """
        Initialize potential calculator.
        
        Parameters:
        -----------
        grid_manager : GridManager
            Grid management instance
        """
        self.grid = grid_manager
        self.atoms = []  # List of (x, z, charge, morse_params)
    
    def add_atom(self, x, z, charge, vdw_radius=2.0, vdw_depth=0.01, morse_a=1.6):
        """
        Add atom to the system.
        
        Parameters:
        -----------
        x, z : float
            Atom coordinates (Angstrom)
        charge : float
            Atomic charge (elementary units)
        vdw_radius : float
            Van der Waals radius (Angstrom)
        vdw_depth : float
            Van der Waals energy minimum depth (eV)
        morse_a : float
            Morse potential exponent (1/Angstrom)
        """
        # Morse potential parameters: r0 = vdw_radius, D = vdw_depth, a = morse_a
        self.atoms.append({
            'x': x, 'z': z, 'charge': charge,
            'morse_D': vdw_depth, 'morse_a': morse_a, 'morse_r0': vdw_radius
        })
        print(f"Added atom at ({x:.1f}, {z:.1f}) Å, charge={charge:+.2f}e, vdW_radius={vdw_radius:.1f}Å, vdW_depth={vdw_depth:.3f}eV")
    
    def get_all_atom_images(self, n_images):
        """
        Generates a list of all atoms including their periodic images.

        Parameters:
        -----------
        n_images : int
            Number of periodic images to include (-n to +n) for each original atom.

        Returns:
        --------
        list_of_atom_images : list of dict
            Each dict contains atom properties including 'x_image' for the replicated x-coordinate.
        """
        all_images = []
        for original_idx, atom_spec in enumerate(self.atoms):
            for n in range(-n_images, n_images + 1):
                image_spec = atom_spec.copy()
                image_spec['x_image'] = atom_spec['x'] + n * self.grid.cell_x
                image_spec['original_atom_index'] = original_idx
                all_images.append(image_spec)
        return all_images

    def calculate_coulomb_periodic(self, n_images=5):
        """
        Calculate Coulomb potential 1/r with periodic images in x-direction.
        
        Parameters:
        -----------
        n_images : int
            Number of periodic images to include (-n to +n)
            
        Returns:
        --------
        potential : ndarray
            2D Coulomb potential on grid
        """
        X, Z = self.grid.get_coordinates()
        potential = np.zeros_like(X)
        
        # Coulomb constant in eV*Angstrom/e^2
        k_e = 14.3996448915
        
        all_atom_images = self.get_all_atom_images(n_images)

        for atom_image in all_atom_images:
            # Distance from each grid point to this image
            dx = X - atom_image['x_image']
            dz = Z - atom_image['z'] # Original z of the atom
            r = np.sqrt(dx**2 + dz**2)
            
            # Avoid division by zero (minimum distance cutoff)
            r = np.maximum(r, 0.1)
            
            potential += k_e * atom_image['charge'] / r
        return potential
    
    def calculate_morse_periodic(self, n_images=5):
        """
        Calculate Morse potential with periodic images in x-direction.
        
        Parameters:
        -----------
        n_images : int
            Number of periodic images to include
            
        Returns:
        --------
        potential : ndarray
            2D Morse potential on grid
        """
        X, Z = self.grid.get_coordinates()
        potential = np.zeros_like(X)
        
        all_atom_images = self.get_all_atom_images(n_images)

        for atom_image in all_atom_images:
            D = atom_image['morse_D']
            a = atom_image['morse_a']
            r0 = atom_image['morse_r0']
            
            dx = X - atom_image['x_image']
            dz = Z - atom_image['z'] # Original z of the atom
            r = np.sqrt(dx**2 + dz**2)
            
            exp_term = np.exp(-a * (r - r0))
            potential += D * (exp_term**2 - 2 * exp_term)
        return potential


class BasisFunctions:
    """Generate basis functions for potential fitting."""
    
    def __init__(self, grid_manager):
        """
        Initialize basis function generator.
        
        Parameters:
        -----------
        grid_manager : GridManager
            Grid management instance
        """
        self.grid = grid_manager
        self.basis_grids = []  # List of 2D basis function grids
        self.basis_labels = []  # Descriptive labels for each basis function
        self.n_harmonics_gen = 0  # Stores the max harmonic index used (e.g., n_harmonics input)
        self.n_z_functions_gen = 0 # Stores the number of z-functions used
    
    def generate_plane_wave_exponential_basis(self, n_harmonics=8, n_z_functions=4):
        """
        Generate basis functions: plane waves in x × exponentials in z.
        
        Parameters:
        -----------
        n_harmonics : int
            Number of plane wave harmonics (n=1,2,...,n_harmonics)
        n_z_functions : int
            Number of exponential decay functions in z
        """
        self.n_harmonics_gen = n_harmonics
        self.n_z_functions_gen = n_z_functions

        X, Z = self.grid.get_coordinates()
        self.basis_grids = []
        self.basis_labels = []
        
        # Z-direction exponential decay parameters
        # Generate decay rates like [1.0, 2.0, ..., n_z_functions]
        z_decay_rates = np.linspace(1.0, float(n_z_functions), n_z_functions) # 1/Å

        # Loop for n includes n=0 for the constant term in x
        for n in range(n_harmonics + 1):  # n from 0 to n_harmonics
            if n == 0:
                # Constant term in x-direction
                plane_wave_term = np.ones_like(X)
                x_label_part = "const(x)"
            else:
                # Plane wave frequency for periodic boundary conditions
                k_x = 2 * np.pi * n / self.grid.cell_x
                plane_wave_term = np.cos(k_x * X)
                x_label_part = f"cos({n}πx/L)"
            
            for i, decay_rate in enumerate(z_decay_rates):
                # Basis function: cos(k_x * x) * exp(-decay_rate * z)
                basis_func = plane_wave_term * np.exp(-decay_rate * Z)
                
                self.basis_grids.append(basis_func)
                self.basis_labels.append(f"{x_label_part}·exp(-{decay_rate:.1f}z)")
        
        print(f"Generated {len(self.basis_grids)} basis functions")
        print(f"X harmonics: n=0 to {n_harmonics}")
        print(f"Z decays: {z_decay_rates}")
    
    def get_basis_matrix(self):
        """
        Get basis functions as matrix for linear fitting.
        
        Returns:
        --------
        basis_matrix : ndarray
            Shape (n_grid_points, n_basis_functions)
        """
        n_points = self.grid.nx * self.grid.nz
        n_basis = len(self.basis_grids)
        
        basis_matrix = np.zeros((n_points, n_basis))
        
        for i, basis_grid in enumerate(self.basis_grids):
            basis_matrix[:, i] = self.grid.flatten_grid(basis_grid)
        
        return basis_matrix
    
    def reconstruct_potential(self, coefficients):
        """
        Reconstruct potential from basis function coefficients.
        
        Parameters:
        -----------
        coefficients : ndarray
            Coefficients for each basis function
            
        Returns:
        --------
        potential : ndarray
            Reconstructed 2D potential on grid
        """
        potential = np.zeros_like(self.basis_grids[0])
        
        for coeff, basis_grid in zip(coefficients, self.basis_grids):
            potential += coeff * basis_grid
        
        return potential


class PotentialFitter:
    """Fit grid potentials to basis functions with regularization."""
    
    def __init__(self, basis_functions):
        """
        Initialize potential fitter.
        
        Parameters:
        -----------
        basis_functions : BasisFunctions
            Basis function generator instance
        """
        self.basis = basis_functions
        self.coefficients = None
        self.rmse = None
        self.regularization_strength = 1e-6
    
    def set_regularization(self, strength):
        """Set L2 regularization strength."""
        self.regularization_strength = strength
        print(f"Regularization strength: {strength:.2e}")
    
    def fit(self, target_potential):
        """
        Fit basis functions to target potential using regularized least squares.
        
        Parameters:
        -----------
        target_potential : ndarray
            2D target potential on grid
            
        Returns:
        --------
        coefficients : ndarray
            Fitted coefficients
        rmse : float
            Root mean square error
        """
        # Get basis matrix and flatten target
        A = self.basis.get_basis_matrix()  # (n_points, n_basis)
        b = self.basis.grid.flatten_grid(target_potential)  # (n_points,)
        
        # Regularized least squares: minimize ||Ax - b||² + λ||x||²
        # Solution: x = (A^T A + λI)^(-1) A^T b
        AtA = A.T @ A
        Atb = A.T @ b
        
        # Add regularization
        regularization_matrix = self.regularization_strength * np.eye(AtA.shape[0])
        regularized_AtA = AtA + regularization_matrix
        
        # Solve regularized system
        self.coefficients = np.linalg.solve(regularized_AtA, Atb)
        
        # Calculate RMSE and total error
        fitted_flat = A @ self.coefficients
        residual = fitted_flat - b
        self.rmse = np.sqrt(np.mean(residual**2))
        
        # Calculate total absolute error (sum over all grid points)
        total_abs_error = np.sum(np.abs(residual))
        mean_abs_error = np.mean(np.abs(residual))
        max_abs_error = np.max(np.abs(residual))
        
        # Calculate relative errors
        target_range = np.max(b) - np.min(b)
        relative_rmse = self.rmse / target_range * 100 if target_range > 0 else 0
        
        print(f"Fitting completed:")
        print(f"  RMSE: {self.rmse:.6f} eV")
        print(f"  Relative RMSE: {relative_rmse:.2f}% of potential range")
        print(f"  Total absolute error: {total_abs_error:.3f} eV")
        print(f"  Mean absolute error: {mean_abs_error:.6f} eV")
        print(f"  Max absolute error: {max_abs_error:.6f} eV")
        print(f"  Target potential range: [{np.min(b):.3f}, {np.max(b):.3f}] eV")
        print(f"  Number of basis functions: {len(self.coefficients)}")
        
        # Print coefficients as a matrix
        num_x_funcs = self.basis.n_harmonics_gen + 1
        num_z_funcs = self.basis.n_z_functions_gen
        if len(self.coefficients) == num_x_funcs * num_z_funcs:
            coeffs_matrix = self.coefficients.reshape((num_x_funcs, num_z_funcs))
            print(f"  Coefficients matrix (rows: x-harmonics {0}-{self.basis.n_harmonics_gen}, cols: z-decays {0}-{num_z_funcs-1}):")
            with np.printoptions(precision=3, suppress=True):
                print(coeffs_matrix)
        else:
            print(f"  Coefficients (linear array): {self.coefficients}")
        
        return self.coefficients, self.rmse
    
    def get_fitted_potential(self):
        """Get reconstructed potential using fitted coefficients."""
        if self.coefficients is None:
            raise ValueError("Must fit first!")
        return self.basis.reconstruct_potential(self.coefficients)


class Visualizer:
    """Visualization tools for potentials and fitting quality."""
    
    def __init__(self, grid_manager):
        """Initialize visualizer."""
        self.grid = grid_manager
    
    def plot_potential_comparison(self, original, fitted, title="Potential Comparison"):
        """
        Plot original vs fitted potential side by side.
        
        Parameters:
        -----------
        original, fitted : ndarray
            2D potential grids to compare
        title : str
            Plot title
        """
        fig, axes = plt.subplots(1, 3, figsize=(15, 4))
        
        # Common colorbar limits
        vmin = min(np.min(original), np.min(fitted))
        vmax = max(np.max(original), np.max(fitted))
        
        # Original potential
        im1 = axes[0].imshow(original.T, origin='lower', aspect='auto', 
                            extent=[0, self.grid.cell_x, self.grid.z_offset, self.grid.z_offset + self.grid.cell_z],
                            vmin=vmin, vmax=vmax, cmap='RdBu_r')
        axes[0].set_title('Original Potential')
        axes[0].set_xlabel('x (Å)')
        axes[0].set_ylabel('z (Å)')
        plt.colorbar(im1, ax=axes[0])
        
        # Fitted potential  
        im2 = axes[1].imshow(fitted.T, origin='lower', aspect='auto',
                            extent=[0, self.grid.cell_x, self.grid.z_offset, self.grid.z_offset + self.grid.cell_z], 
                            vmin=vmin, vmax=vmax, cmap='RdBu_r')
        axes[1].set_title('Fitted Potential')
        axes[1].set_xlabel('x (Å)')
        axes[1].set_ylabel('z (Å)')
        plt.colorbar(im2, ax=axes[1])
        
        # Difference
        diff = original - fitted
        im3 = axes[2].imshow(diff.T, origin='lower', aspect='auto',
                            extent=[0, self.grid.cell_x, self.grid.z_offset, self.grid.z_offset + self.grid.cell_z],
                            cmap='RdBu_r')
        axes[2].set_title('Difference (Original - Fitted)')
        axes[2].set_xlabel('x (Å)')
        axes[2].set_ylabel('z (Å)')
        plt.colorbar(im3, ax=axes[2])
        
        plt.suptitle(title)
        plt.tight_layout()
        #plt.savefig(filename)
        #print(f"Potential comparison plot saved to {filename}")
        #plt.show()
    
    def plot_basis_functions(self, basis_object, coefficients=None):
        """
        Plot all generated basis functions in a grid.
        Rows correspond to z-functions, columns to x-functions.

        Parameters:
        -----------
        basis_object : BasisFunctions
            The BasisFunctions instance containing the generated basis.
        coefficients : np.ndarray, optional
            Fitted coefficients. If provided, they are added to the titles.
        """
        num_x_funcs = basis_object.n_harmonics_gen + 1  # n=0 to n_harmonics_gen
        num_z_funcs = basis_object.n_z_functions_gen
        
        if not basis_object.basis_grids:
            print("No basis functions to plot.")
            return

        fig, axes = plt.subplots(num_z_funcs, num_x_funcs,  figsize=(3 * num_x_funcs, 2.8 * num_z_funcs), squeeze=False)
        
        for k, basis_grid in enumerate(basis_object.basis_grids):
            # Order of generation: outer loop x (n), inner loop z (decay_rate)
            # So, to map to (row=z_idx, col=x_idx):
            x_idx = k // num_z_funcs  # Column index (0 to n_harmonics_gen)
            z_idx = k % num_z_funcs   # Row index (0 to n_z_functions_gen - 1)
            
            ax = axes[z_idx, x_idx]
            im = ax.imshow(basis_grid.T, origin='lower', aspect='auto',
                           extent=[0, self.grid.cell_x, self.grid.z_offset, 
                                   self.grid.z_offset + self.grid.cell_z],
                           cmap='RdBu_r')
            title = basis_object.basis_labels[k]
            if coefficients is not None and k < len(coefficients):
                title += f"\nC={coefficients[k]:.2e}"
            ax.set_title(title, fontsize=7)
            # ax.set_xlabel('x (Å)', fontsize=7) # Can be set per column or row if desired
            # ax.set_ylabel('z (Å)', fontsize=7) # Can be set per column or row if desired
            ax.tick_params(axis='both', which='major', labelsize=6)
            plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
            
        plt.suptitle('Basis Functions (Rows: Z-decay, Cols: X-harmonic [n=0 to max])')
        plt.tight_layout(rect=[0, 0.02, 1, 0.96]) # Adjust for suptitle and labels
        #plt.savefig(filename)
        #print(f"Basis functions plot saved to {filename}")
        #plt.show()

    def plot_potential_with_atoms(self, potential_grid, all_atom_images,  title="Potential with Atom Images"):
        """
        Plots the 2D potential grid with all atom images (original and periodic) overlaid.

        Parameters:
        -----------
        potential_grid : np.ndarray
            The 2D potential grid to plot.
        all_atom_images : list of dict
            List of atom image specifications, from PotentialCalculator.get_all_atom_images().
        title : str
            Plot title.
        """
        fig, ax = plt.subplots(figsize=(15, 4)) # Elongated aspect ratio, e.g., 15 width, 4 height

        # Plot potential using imshow
        # The extent should match the unit cell for the potential grid
        im = ax.imshow(potential_grid.T, origin='lower', aspect='auto',
                       extent=[0, self.grid.cell_x, self.grid.z_offset, self.grid.z_offset + self.grid.cell_z],
                       cmap='RdBu_r', interpolation='bicubic',
                       vmin=np.min(potential_grid), vmax=np.max(potential_grid))
        plt.colorbar(im, ax=ax, label="Potential (eV)", shrink=0.8)

        # Plot atoms
        if all_atom_images:
            r0_type_Na = 1.6 # Assuming Na-like has r0=1.6
            for ia,atom_img in enumerate(all_atom_images):
                color = 'blue' if np.isclose(atom_img['morse_r0'], r0_type_Na) else 'red'
                # Scatter plot size 's' is roughly proportional to area
                scatter_size = (atom_img['morse_r0'] * 7)**2 # Adjust scaling factor as needed
                print( f"Atom {ia}: x_image={atom_img['x_image']:.2f}, z={atom_img['z']:.2f}, r0={atom_img['morse_r0']:.2f}, color={color}")
                ax.scatter(atom_img['x_image'], atom_img['z'], marker='o', s=100.0, c=color, edgecolors='k', linewidth=0.0)

            # Set plot limits to show all atom images
            atom_x_coords = [atom['x_image'] for atom in all_atom_images]
            min_x_atom = min(atom_x_coords)
            max_x_atom = max(atom_x_coords)
            # Ensure the central cell [0, cell_x] is clearly visible within the atom image range
            plot_min_x = min(min_x_atom - self.grid.cell_x * 0.1, -0.05 * self.grid.cell_x)
            plot_max_x = max(max_x_atom + self.grid.cell_x * 0.1, self.grid.cell_x * 1.05)
            ax.set_xlim(plot_min_x, plot_max_x)

        ax.set_ylim( -1.0 ,  11.0 )

        # Draw lines for the primary unit cell boundaries
        ax.axvline(0, color='gray', linestyle=':', linewidth=1)
        ax.axvline(self.grid.cell_x, color='gray', linestyle=':', linewidth=1)

        ax.set_xlabel("x (Å)")
        ax.set_ylabel("z (Å)")
        ax.set_title(title, fontsize=10)
        #plt.tight_layout()
        #plt.savefig(filename)
        #print(f"Potential with atoms plot saved to {filename}")
        #plt.show()

def plot_1d_morse_debug(atom_params_list, title="1D Morse Potential Debug"):
    """
    Plots the 1D Morse potential for a list of atom parameters.

    Parameters:
    -----------
    atom_params_list : list of dict
        List of dictionaries, each containing 'morse_D', 'morse_a', 'morse_r0', and 'label'.
    title : str
        Plot title.
    """
    plt.figure(figsize=(10, 6))
    r_values = np.linspace(0.5, 10.0, 500)  # Start r from a small positive value up to 10Å

    max_D_abs = 0.0
    for params in atom_params_list:
        D  = params['morse_D']
        a  = params['morse_a']
        r0 = params['morse_r0']
        label = params.get('label', f'D={D:.2f}, a={a:.1f}, r0={r0:.1f}')
        exp_term = np.exp(-a * (r_values - r0))
        morse_V = D * (exp_term**2 - 2 * exp_term)
        plt.plot(r_values, morse_V, label=label)
        if D > 0: # D is depth, should be positive
            max_D_abs = max(max_D_abs, D)
    plt.ylim(-1.1 * max_D_abs, 1.5 * max_D_abs if max_D_abs > 0 else 0.1) # Ensure well and repulsive part are visible
    plt.axhline(0, color='gray', linestyle='--', linewidth=0.8)
    plt.xlabel("Distance r (Å)")
    plt.ylabel("Morse Potential V(r) (eV)")
    plt.title(title)
    plt.legend()
    plt.grid(True)
    #plt.show()

def main_example():

    
    print("\nExample completed successfully!")
    return grid, calc, basis, fitter, viz


if __name__ == "__main__":
    """Example usage of the surface potential fitting framework."""
    
    print("Surface Potential Fitting Framework")
    print("=" * 40)
    
    # 1. Setup grid with 2.0 Å offset above surface
    print("\n1. Setting up grid...")
    grid = GridManager(cell_x=10.0, cell_z=10.0, grid_step=0.2, z_offset=2.0)
    
    # 2. Setup atoms with realistic parameters (Na+ and Cl- at surface)
    print("\n2. Adding atoms...")
    calc = PotentialCalculator(grid)
    
    # Na-like atom at the edge (x=0) for symmetry: charge 0.0e, vdW radius 1.6Å, vdW depth 0.01eV, Morse exponent 1.6 Å^-1
    calc.add_atom(x=0.0, z=0.0, charge=0.0, vdw_radius=1.6, vdw_depth=0.01, morse_a=1.6)
    
    # Cl-like atom in the center (x=cell_x/2) for symmetry: charge 0.0e, vdW radius 2.3Å, vdW depth 0.01eV, Morse exponent 1.6 Å^-1  
    calc.add_atom(x=grid.cell_x / 2.0, z=0.0, charge=0.0, vdw_radius=2.3, vdw_depth=0.01, morse_a=1.6)

    # --- 1D Morse Potential Debug Plot ---
    # print("\n--- Plotting 1D Morse Potentials for Debugging ---")
    # debug_atom_params = []
    # atom_labels = ["Na-like (r0=1.6Å, D=0.01eV)", "Cl-like (r0=2.3Å, D=0.01eV)"] 
    # for i, atom_p in enumerate(calc.atoms):
    #     debug_atom_params.append({
    #         'label': atom_labels[i] if i < len(atom_labels) else f"Atom {i+1}",
    #         'morse_D': atom_p['morse_D'],
    #         'morse_a': atom_p['morse_a'],
    #         'morse_r0': atom_p['morse_r0']
    #     })
    # plot_1d_morse_debug(debug_atom_params, title="1D Morse Potentials (Individual Atom Parameters)")
    # --- End of 1D Morse Potential Debug Plot ---

    
    # 3. Calculate periodic potentials
    print("\n3. Calculating periodic potentials...")
    coulomb_potential = calc.calculate_coulomb_periodic(n_images=3)
    morse_potential = calc.calculate_morse_periodic(n_images=3)
    
    # Total potential (for realistic modeling)
    # Since charges are 0, coulomb_potential will be 0. We fit/visualize Morse directly.
    total_potential_to_fit = morse_potential 
    
    # Get all atom images for the new visualization
    n_images_for_calc_and_plot = 3 # Should be same as used in calculate_..._periodic
    all_atom_images_for_plot = calc.get_all_atom_images(n_images=n_images_for_calc_and_plot)

    print(f"Coulomb potential range: [{np.min(coulomb_potential):.3f}, {np.max(coulomb_potential):.3f}] eV")
    print(f"Morse potential range: [{np.min(morse_potential):.3f}, {np.max(morse_potential):.3f}] eV")
    print(f"Potential to fit range: [{np.min(total_potential_to_fit):.3f}, {np.max(total_potential_to_fit):.3f}] eV")
    
    # 4. Generate basis functions
    print("\n4. Generating basis functions...")
    basis = BasisFunctions(grid)
    basis.generate_plane_wave_exponential_basis(n_harmonics=3, n_z_functions=4) # n_harmonics for 0,1,2,3
    
    # 5. Fit potential with different regularization strengths
    print("\n5. Fitting potential...")
    fitter = PotentialFitter(basis)
    
    # Try fitting with minimal regularization first
    print("\n--- Fitting with low regularization ---")
    fitter.set_regularization(1e-8)
    coefficients, rmse = fitter.fit(total_potential_to_fit)
    
    # If error is large, try stronger regularization
    if rmse > 0.01:  # If RMSE > 10 meV
        print("\n--- Trying stronger regularization ---")
        fitter.set_regularization(1e-4)
        coefficients, rmse = fitter.fit(total_potential_to_fit)
    
    # 6. Get fitted potential
    fitted_potential = fitter.get_fitted_potential()
    
    # 7. Visualize results
    print("\n6. Visualizing results...")
    viz = Visualizer(grid)
    viz.plot_potential_with_atoms(total_potential_to_fit, all_atom_images_for_plot, title=f"Morse Potential with Atom Images (n_images={n_images_for_calc_and_plot})")
    viz.plot_potential_comparison(total_potential_to_fit, fitted_potential,title="Morse Potential Fitting (Charges = 0)")
    viz.plot_basis_functions(basis, coefficients=coefficients )

    plt.show()