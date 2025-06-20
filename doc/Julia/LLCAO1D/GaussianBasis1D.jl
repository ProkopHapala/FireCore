module GaussianBasis1D
"""
GaussianBasis1D – basic 1-D Gaussian basis utilities.
This module is intentionally lightweight; it only knows about
basis-function geometry and analytical integrals.  No localisation
logic is present here.
"""

export GaussBasis1D, overlap, kinetic, coulomb, make_HS

using LinearAlgebra
# minimal erf approximation to avoid external deps
@inline function _erf(x::Float64)
    # Abramowitz & Stegun 7.1.26 approximation (good to 1e-7)
    t = 1.0 / (1.0 + 0.3275911*abs(x))
    τ = (((((1.061405429 * t) - 1.453152027)*t) + 1.421413741)*t - 0.284496736)*t + 0.254829592
    y = 1.0 - τ * exp(-x*x)
    return copysign(y, x)
end

# ---------------------------- Data container -----------------------------

struct GaussBasis1D
    centers::Vector{Float64}   # nuclear / basis centres (Bohr)
    widths ::Vector{Float64}   # Gaussian width w  (phi(x) = N * exp(- (x-X)^2 / w^2 ))
    """Inverse width α = 1/w² is used internally."""
    function GaussBasis1D(centers::Vector{<:Real}, widths::Vector{<:Real})
        length(centers)==length(widths) || error("centers and widths must have same length")
        new(Vector{Float64}(centers), Vector{Float64}(widths))
    end
end

const π = Base.MathConstants.pi

# --------------------------- Helper functions ----------------------------
α(w) = 1.0 / w^2          # convert width -> exponent  (Bohr⁻²)
normfactor(a) = (2a/π)^0.25  # 1-D normalisation of Gaussian with exponent α

# Boys F₀ for 1-D attraction integral (simple implementation)
F0(t) = t < 1e-10 ? 1.0 : 0.5 * sqrt(π / t) * _erf(sqrt(t))

# --------------------------- Integral kernels ----------------------------

"Overlap integral  ⟨ϕ_i|ϕ_j⟩ for normalised 1-D Gaussians"
function overlap(ai, Xi, aj, Xj)
    p = ai + aj
    return normfactor(ai)*normfactor(aj) * exp(-(ai*aj/p)*(Xi-Xj)^2) * sqrt(π/p)
end

"Kinetic energy integral  ⟨ϕ_i| -½∇² |ϕ_j⟩"
function kinetic(ai, Xi, aj, Xj, Sij)
    p = ai + aj
    return Sij*(ai*aj/p - 2*(ai*aj)^2/p^2*(Xi-Xj)^2)
end

"Nuclear attraction  ⟨ϕ_i| -Z/|x-Xk| |ϕ_j⟩  (Gauss-transform route)"
function coulomb(ai, Xi, aj, Xj, Zk, Xk)
    p = ai + aj
    Xp = (ai*Xi + aj*Xj)/p
    Kij = normfactor(ai)*normfactor(aj)*exp(-(ai*aj/p)*(Xi-Xj)^2)
    t = p*(Xp-Xk)^2
    return -Zk * Kij * sqrt(2π/p) * F0(t)
end

# --------------------------- Matrix builder ------------------------------

"""
make_HS(basis, charges, charge_positions) -> H, S
Returns Hamiltonian (kinetic+potential) and overlap matrices using analytical
integrals.
"""
function make_HS(basis::GaussBasis1D, charges::Vector{Float64}, charge_pos::Vector{Float64})
    n = length(basis.centers)
    H = zeros(n,n);  S = zeros(n,n)
    for i in 1:n, j in 1:i
        ai, Xi = α(basis.widths[i]), basis.centers[i]
        aj, Xj = α(basis.widths[j]), basis.centers[j]
        Sij = overlap(ai, Xi, aj, Xj)
        Tij = kinetic(ai, Xi, aj, Xj, Sij)
        Vij = 0.0
        for (Zk,Xk) in zip(charges, charge_pos)
            Vij += coulomb(ai, Xi, aj, Xj, Zk, Xk)
        end
        Hij = Tij + Vij
        H[i,j] = H[j,i] = Hij
        S[i,j] = S[j,i] = Sij
    end
    return H, S
end

end # module
