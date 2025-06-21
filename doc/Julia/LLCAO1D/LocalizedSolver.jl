module LocalizedSolver
using LinearAlgebra
export LocalizedSolver1D, solve_localized

struct LocalizedSolver1D
    H::Matrix{Float64}
    S::Matrix{Float64}
    basis_centers::Vector{Float64}
    # parameters
    maxiter::Int
    step::Float64
    ortho_damp::Float64
    ortho_iter::Int
    coeff_tol::Float64
    overlap_tol::Float64
    alpha_loc::Float64
    hard_cut::Bool
end

LocalizedSolver1D(H::AbstractMatrix{T}, S::AbstractMatrix{T}, basis_centers::Vector{Float64};
                  maxiter=200, step=0.05, ortho_damp=0.5, ortho_iter=5,
                  coeff_tol=1e-6, overlap_tol=1e-6, alpha_loc=0.0,
                  hard_cut=true) where {T<:Real} =
    LocalizedSolver1D(Matrix{Float64}(H), Matrix{Float64}(S), basis_centers,
                      maxiter, step, ortho_damp, ortho_iter,
                      coeff_tol, overlap_tol, alpha_loc, hard_cut)

# ---------------- internal helpers ------------------
calc_SMO(S,C) = C' * S * C

function normalize!(S, C)
    for i in axes(C,2)
        n2 = C[:,i]' * S * C[:,i]
        n2 > 1e-12 || continue
        C[:,i] ./= sqrt(n2)
    end
end

function orthogonalize!(solver::LocalizedSolver1D, C)
    for _ in 1:solver.ortho_iter
        SMO = calc_SMO(solver.S, C)
        for i in axes(C,2), j in axes(C,2)
            i==j && continue
            Sij = SMO[i,j]
            C[:,i] .-= solver.ortho_damp * Sij .* C[:,j]
        end
    end
end

function apply_cutoff!(C, solver::LocalizedSolver1D, loc_centers, cutR)
    (cutR === nothing || !solver.hard_cut) && return
    for i in axes(C,2) # for each MO
        for μ in axes(C,1) # for each basis function
            if abs(solver.basis_centers[μ] - loc_centers[i]) > cutR
                C[μ, i] = 0.0
            end
        end
    end
end

function gradient(solver::LocalizedSolver1D, C, loc_centers)
    # electronic part of the gradient G = H*C
    G = solver.H * C

    # localization potential (quadratic in distance)
    if solver.alpha_loc > 0.0 && loc_centers !== nothing
        # dist2 has shape (nbasis, nmo)
        dist2 = (solver.basis_centers .- loc_centers').^2
        G .+= 2.0 * solver.alpha_loc .* C .* dist2
    end

    return -solver.step .* G
end

function solve_localized(solver::LocalizedSolver1D; loc_centers=nothing, cutR=nothing)
    nbas = size(solver.H,1)
    nm   = nbas
    # --- Initialization ---
    C = Matrix{Float64}(I, nbas, nm)
    apply_cutoff!(C, solver, loc_centers, cutR)
    normalize!(solver.S,C)
    for iter in 1:solver.maxiter
        Cold = copy(C)
        # gradient step
        C .+= gradient(solver, C, loc_centers)
        orthogonalize!(solver,C)
        normalize!(solver.S,C)
        apply_cutoff!(C, solver, loc_centers, cutR)
        normalize!(solver.S,C) # re-normalize after cutoff
        # convergence
        Δ = maximum(abs.(C-Cold))
        SMO = calc_SMO(solver.S,C)
        S_off = maximum(abs.(SMO .- Matrix{Float64}(I, size(SMO))))
        (Δ<solver.coeff_tol && S_off<solver.overlap_tol) && break
    end
    return C
end

end # module
