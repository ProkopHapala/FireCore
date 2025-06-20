module LocalizedSolver
using LinearAlgebra
export LocalizedSolver1D, solve_localized

struct LocalizedSolver1D
    H::Matrix{Float64}
    S::Matrix{Float64}
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

LocalizedSolver1D(H::AbstractMatrix{T}, S::AbstractMatrix{T};
                  maxiter=200, step=0.05, ortho_damp=0.5, ortho_iter=5,
                  coeff_tol=1e-6, overlap_tol=1e-6, alpha_loc=0.0,
                  hard_cut=true) where {T<:Real} =
    LocalizedSolver1D(Matrix{Float64}(H), Matrix{Float64}(S),
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

function gradient(solver::LocalizedSolver1D, C)
    return -solver.step .* (solver.H * C)   # simple downhill step for now
end

function solve_localized(solver::LocalizedSolver1D; loc_centers=nothing, cutR=nothing)
    nbas = size(solver.H,1)
    nm   = nbas
    C = Matrix{Float64}(I, nbas, nm)
    normalize!(solver.S,C)
    for iter in 1:solver.maxiter
        Cold = copy(C)
        # gradient step
        C .+= gradient(solver,C)
        orthogonalize!(solver,C)
        normalize!(solver.S,C)
        # convergence
        Δ = maximum(abs.(C-Cold))
        SMO = calc_SMO(solver.S,C)
        S_off = maximum(abs.(SMO .- Matrix{Float64}(I, size(SMO))))
        (Δ<solver.coeff_tol && S_off<solver.overlap_tol) && break
    end
    return C
end

end # module
