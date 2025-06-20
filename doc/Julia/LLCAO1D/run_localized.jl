# Simple driver for 1-D localized solver in Julia

# --- load local modules ---
include("GaussianBasis1D.jl")
include("LocalizedSolver.jl")

using .GaussianBasis1D
using .LocalizedSolver
using LinearAlgebra

# ---------------- System definition ----------------
num_atoms = 10
spacing    = 1.5          # Bohr
centers = collect(range(-(num_atoms-1)*spacing/2, length=num_atoms, step=spacing))
widths  = fill(1.0, num_atoms)  # Gaussian width w (Bohr)
charges = fill(1.0, num_atoms)  # Hydrogen chain

basis = GaussBasis1D(centers, widths)
H,S = make_HS(basis, charges, centers)  # treat nuclear positions identical to basis centers

# ---------------- Localised solver -----------------
solver = LocalizedSolver1D(H,S; maxiter=400, step=0.05, ortho_damp=0.5,
                           coeff_tol=1e-7, overlap_tol=1e-5,
                           alpha_loc=0.0, hard_cut=false)

println("Running localized solver on $num_atoms-atom chain …")
C = LocalizedSolver.solve_localized(solver)

# energies of each orbital
num_mos = size(C,2)
energies = zeros(num_mos)
for i in 1:num_mos
    num = C[:,i]'*H*C[:,i]
    den = C[:,i]'*S*C[:,i]
    energies[i] = num/den
end
println("Localized orbital energies (Hartree):")
foreach(e->println("  ",e), energies)

# report max off-diag overlap
S_MO = C' * S * C
max_off = maximum(abs.(S_MO .- Diagonal(diag(S_MO))))
println("Max off–diagonal MO overlap: ", max_off)
