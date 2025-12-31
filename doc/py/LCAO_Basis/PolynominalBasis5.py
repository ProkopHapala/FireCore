import numpy as np
from numpy.polynomial import Polynomial
from scipy.integrate import fixed_quad

class PolynomialBasis:
    def __init__(self, coeffs, R_cut, fit_degree=12):
        """
        coeffs: list of [c0, c1, c2...] for phi(r)
        R_cut: Cutoff radius
        fit_degree: Degree of polynomial to fit the interactions (10-14 is usually enough)
        """
        self.R = float(R_cut)
        self.poly_phi = Polynomial(coeffs)
        
        # 1. Normalize
        # Norm = 4pi * Int_0^R r^2 phi^2
        integrand = (self.poly_phi**2) * Polynomial([0, 0, 1])
        integral = integrand.integ()(self.R) - integrand.integ()(0)
        norm = 1.0 / np.sqrt(4 * np.pi * integral)
        self.poly_phi *= norm
        
        # 2. Prepare Helper Polynomials
        self.poly_rho = self.poly_phi**2
        
        # Laplacian (for Kinetic)
        # L = 1/r^2 d/dr (r^2 dphi/dr)
        r2_dphi = (self.poly_phi.deriv()) * Polynomial([0, 0, 1])
        d_r2_dphi = r2_dphi.deriv()
        # Divide by r^2 (shift coeffs down by 2)
        c = d_r2_dphi.coef
        lap_coeffs = c[2:] if len(c) > 2 else [0]
        self.poly_laplacian = Polynomial(lap_coeffs)
        
        # Internal Potential V_in (for Coulomb)
        self.poly_V_in = self._solve_poisson_poly(self.poly_rho)
        
        # 3. FIT THE INTERACTION POLYNOMIALS
        # We calculate the exact integrals numerically at Chebyshev nodes
        # to ensure a perfect polynomial fit.
        print(f"--- Fitting Interaction Polynomials (N={fit_degree}) ---")
        self.S_coeffs = self._fit_interaction(self._integrand_overlap, fit_degree)
        self.T_coeffs = self._fit_interaction(self._integrand_kinetic, fit_degree)
        self.J_coeffs = self._fit_interaction(self._integrand_coulomb, fit_degree)
        
    def _solve_poisson_poly(self, rho):
        """ Returns V_in(r) polynomial """
        # Q(r)
        rho_r2 = rho * Polynomial([0, 0, 1])
        Q_prim = rho_r2.integ()
        Q_total = 4 * np.pi * (Q_prim(self.R) - Q_prim(0))
        
        # Indefinite integrals
        rhs = -4 * np.pi * rho_r2
        first_int = rhs.integ()
        
        # Divide by r^2
        c = first_int.coef
        dVdr = Polynomial(c[2:] if len(c) > 2 else [0])
        
        V_shape = dVdr.integ()
        C2 = (Q_total / self.R) - V_shape(self.R)
        return V_shape + C2

    # =========================================================
    #  FITTING ENGINE
    # =========================================================
    def _fit_interaction(self, integ_func, deg):
        """
        Generates 'deg' Chebyshev points, computes exact integral,
        and fits a polynomial to (Integral * d).
        Why multiply by d? Because the form is P(d)/d. 
        Fitting P(d) is smoother and numerically stable.
        """
        # Chebyshev nodes for stability in [0, 2R]
        k = np.arange(deg + 1)
        x = np.cos(np.pi * (2*k + 1) / (2*(deg + 1))) # Nodes in [-1, 1]
        # Map to [0, 2R] (avoid d=0 exactly to prevent div zero in quad)
        d_samples = self.R * (x + 1) 
        d_samples[d_samples < 1e-6] = 1e-6 # Safety
        
        y_samples = []
        for d in d_samples:
            # Compute Exact Integral
            val = self._bipolar_quad(d, integ_func)
            # We fit the function P(d) = Energy(d) * d
            y_samples.append(val * d)
            
        # Fit polynomial
        return Polynomial.fit(d_samples, y_samples, deg, domain=[0, 2*self.R]).convert().coef

    def _bipolar_quad(self, d, func_integrand):
        """ High precision quadrature for derivation """
        if d >= 2*self.R: return 0.0
        
        # Note: Coulomb uses full range [0, R]. Others use [d-R, R] or [0, R] overlap.
        # We handle limits inside the integrand wrappers or here.
        # Simplest: Integrate rA from 0 to R, mask zeros.
        
        val, _ = fixed_quad(func_integrand, 0, self.R, args=(d,), n=60)
        return (2 * np.pi / d) * val

    # =========================================================
    #  INTEGRANDS (Used only during init)
    # =========================================================
    def _integrand_overlap(self, rA, d):
        # Int rA * phi(rA) * Int_rB phi(rB)
        # Antiderivative of r*phi(r)
        r_phi = self.poly_phi * Polynomial([0, 1])
        H = r_phi.integ()
        
        # Vectorized evaluation
        out = np.zeros_like(rA)
        rA = np.atleast_1d(rA)
        
        # Integration limits for inner B
        # rB goes from |d-rA| to min(R, d+rA)
        # Since phi is 0 outside R, H(min(R, ...)) is H(R) effectively if clipped.
        
        rb_upper = np.minimum(d + rA, self.R)
        rb_lower = np.abs(d - rA)
        
        # Mask where intersection is invalid (shouldn't happen in [0, 2R] but good practice)
        valid = rb_upper > rb_lower
        
        val_upper = H(rb_upper[valid])
        val_lower = H(rb_lower[valid])
        
        out[valid] = rA[valid] * self.poly_phi(rA[valid]) * (val_upper - val_lower)
        return out

    def _integrand_kinetic(self, rA, d):
        # Int rA * phi(rA) * Int_rB Laplacian(rB)
        r_lap = self.poly_laplacian * Polynomial([0, 1])
        H = r_lap.integ()
        
        out = np.zeros_like(rA)
        rA = np.atleast_1d(rA)
        
        rb_upper = np.minimum(d + rA, self.R)
        rb_lower = np.abs(d - rA)
        valid = rb_upper > rb_lower
        
        val_upper = H(rb_upper[valid])
        val_lower = H(rb_lower[valid])
        
        # Kinetic factor -0.5
        out[valid] = -0.5 * rA[valid] * self.poly_phi(rA[valid]) * (val_upper - val_lower)
        return out

    def _integrand_coulomb(self, rA, d):
        # Int rA * rho(rA) * Int_rB V_B(rB)
        # V_B is Piecewise: V_in(r) for r<R, 1/r for r>R
        
        # We need primitive of r*V(r).
        # Inside: r*V_in. Outside: r*(1/r) = 1.
        
        r_Vin = self.poly_V_in * Polynomial([0, 1])
        H_in = r_Vin.integ()
        
        # H(r) definition
        def get_H(r_vals):
            res = np.zeros_like(r_vals)
            # Inside R
            mask_in = r_vals <= self.R
            if np.any(mask_in):
                res[mask_in] = H_in(r_vals[mask_in])
            
            # Outside R: Integral_0^r xV = Int_0^R xV_in + Int_R^r 1 dx
            # = H_in(R) + (r - R)
            mask_out = ~mask_in
            if np.any(mask_out):
                res[mask_out] = H_in(self.R) + (r_vals[mask_out] - self.R)
            return res

        out = np.zeros_like(rA)
        rA = np.atleast_1d(rA)
        
        rb_upper = d + rA
        rb_lower = np.abs(d - rA)
        
        val_upper = get_H(rb_upper)
        val_lower = get_H(rb_lower)
        
        out = rA * self.poly_rho(rA) * (val_upper - val_lower)
        return out

    # =========================================================
    #  RUNTIME EVALUATION
    # =========================================================
    def _eval(self, coeffs, d):
        # Evaluate polynomial sum(c_n * d^n)
        return np.polynomial.polynomial.polyval(d, coeffs)

    def overlap(self, d):
        if d >= 2*self.R: return 0.0
        val = self._eval(self.S_coeffs, d)
        return val / (d + 1e-12)

    def kinetic(self, d):
        if d >= 2*self.R: return 0.0
        val = self._eval(self.T_coeffs, d)
        return val / (d + 1e-12)

    def coulomb(self, d):
        if d >= 2*self.R: return 1.0/d
        # We fitted (J * d), so we divide by d
        val = self._eval(self.J_coeffs, d)
        return val / (d + 1e-12)

# =========================================================
#  TEST
# =========================================================
if __name__ == "__main__":
    import matplotlib.pyplot as plt
    
    # Example: Parabolic (1 - r^2)^2 for r < 1.0
    # Coefficients [1, 0, -2, 0, 1] for (1-r^2)^2
    R = 1.0
    coeffs = [1, 0, -2, 0, 1] 
    
    basis = PolynomialBasis(coeffs, R, fit_degree=16)
    
    # Test Points
    d_vals = np.linspace(0.01, 3.0, 100)
    J_vals = [basis.coulomb(d) for d in d_vals]
    S_vals = [basis.overlap(d) for d in d_vals]
    
    # Check limit at d -> 0
    print(f"J(0.01) = {basis.coulomb(0.01):.5f}")
    
    # Check transition at 2R
    print(f"J(1.99) = {basis.coulomb(1.99):.5f}")
    print(f"J(2.01) = {basis.coulomb(2.01):.5f} (Should be 1/2.01 = {1/2.01:.5f})")

    plt.figure(figsize=(10, 5))
    plt.plot(d_vals, J_vals, label="Coulomb J(d)")
    plt.plot(d_vals, 1/d_vals, 'k--', alpha=0.5, label="1/d")
    plt.plot(d_vals, S_vals, label="Overlap S(d)")
    plt.ylim(-0.2, 2.0)
    plt.axvline(2*R, color='r', linestyle=':', label="2R")
    plt.legend()
    plt.grid(True)
    plt.title("Robust Polynomial Fitting Method")
    plt.show()
    
    print("\nFitted Coefficients for Coulomb (J*d):")
    print(basis.J_coeffs)