# run_localized.jl — driver for 1-D localized solver

using Revise
using Plots
using LinearAlgebra

includet("GaussianBasis1D.jl")
includet("LocalizedSolver.jl")
using .GaussianBasis1D, .LocalizedSolver

function run_localized(; num_atoms=10, spacing=1.5)

    # --- build system ---
    num_atoms = 10
    spacing    = 1.5          # Bohr
    centers = collect(range(-(num_atoms-1)*spacing/2, length=num_atoms, step=spacing))
    widths  = fill(1.0, num_atoms)  # Gaussian width w (Bohr)
    charges = fill(1.0, num_atoms)  # Hydrogen chain

    basis = GaussBasis1D(centers, widths)
    H,S = make_HS(basis, charges, centers)  # treat nuclear positions identical to basis centers

    # ---------------- Localised solver -----------------
    solver = LocalizedSolver1D(H,S, centers; maxiter=400, step=0.05, ortho_damp=0.5,
                            coeff_tol=1e-7, overlap_tol=1e-5,
                            alpha_loc=0.0, hard_cut=false)

    println("Running localized solver on $num_atoms-atom chain …")
    C = LocalizedSolver.solve_localized(solver; loc_centers=centers, cutR=spacing * 8.0)

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



    # ---------------- Plotting -----------------------
    num_plot = min(4, num_mos)
    x_min = minimum(centers) - 5
    x_max = maximum(centers) + 5
    x_grid = range(x_min, x_max, length=1000)

    α(w) = 1.0 / w^2
    norm(a) = (2a/π)^0.25
    gaussian_val(x,X,w) = norm(α(w)) * exp(-α(w)*(x-X)^2)

    function mo_wavefunction(idx)
        psi = zero.(x_grid)
        for μ in eachindex(centers)
            psi .+= C[μ,idx] .* gaussian_val.(x_grid, centers[μ], widths[μ])
        end
        return psi
    end

    plt = plot(x_grid, mo_wavefunction(1) .+ energies[1], label="MO1")
    for i in 2:num_plot
        plot!(plt, x_grid, mo_wavefunction(i) .+ energies[i], label="MO$i")
    end
    xlabel!(plt, "x (Bohr)"); ylabel!(plt, "ψ(x) + energy")
    title!(plt, "First $num_plot Localized Orbitals")
    display(plt)

end

# invoke
run_localized()
