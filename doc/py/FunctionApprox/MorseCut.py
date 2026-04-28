import numpy             as np
import matplotlib.pyplot as plt

'''

We should use the relation e^x = lim_{n->inf} (1+x/n)^n   to apporximate the exponential function with a polynomial.

slope of 1-x^2 at point x=1 is -2 

'''

def exp_approx( x,  Ri=1.4,  Rcut=6.0, n=8 ):
    x_ = 1-x/Rcut
    k = 1/(1-Ri/Rcut)
    y = (x_*k)**n
    mask = x>Rcut
    y[mask] = 0
    return y

def exp_approx_2( x, Ri=1.4, Rcut=5.0, n=8 ):
    x_ = 1-(x/Rcut)**2
    k  = 1/(1-(Ri/Rcut)**2)
    y = (x_*k)**n
    mask = x>Rcut
    y[mask] = 0
    return y

def numDeriv( x, y ):

    dx = x[2:]-x[:-2]
    dy = y[2:]-y[:-2]
    x_ = x[1:-1]
    return -dy/dx, x_


def Morse( r, Ei=0.1, Ri=1.4, b=1.6 ):
    e = np.exp( -b*( r - Ri ) )
    E    =      Ei*( e*e - 2*e )
    F    =  2*b*Ei*( e*e -   e ) 
    return E,F, e

def LennardJones( r, Ei=0.1, Ri=1.4 ):
    ir6  = (Ri/r)**6
    ir12 = ir6*ir6
    E    = Ei*   ( ir12 - 2*ir6 )
    F    = Ei*12*( ir12 -  ir6  )/r
    return E,F




def MorseCut_1( r, Ei=0.1, Ri=1.4, b=1.6, Rc=6.0, n=8 ):
    x_  = 1-r/Rc
    k   = 1/(1-Ri/Rc)
    e   = (x_*k)**n
    mask    = r>Rc
    e[mask] = 0
    E =    Ei*(e*e - 2*e)
    F =  2*Ei*(e*e -   e)
    return E,F, e

def MorseCut_2( r, Ei=0.1, Ri=1.4, b=1.6, Rc=6.0, n=8 ):
    r2 = r*r
    x_ = 1-(r2/(Rc*Rc))
    k  = 1/(1-(Ri/Rc)**2)
    e = (x_*k)**n
    mask = r2>(Rc*Rc)
    e[mask] = 0
    E =    Ei*(e*e - 2*e)
    F =  2*Ei*(e*e -   e)
    return E,F, e

# def MorseCut_1_( r, Ei=0.1, Ri=1.4, b=1.6, Rc=6.0, n=8 ):
#     """
#     A numerically efficient, physically correct polynomial approximation of the Morse Potential.
    
#     It approximates the exponential term exp(-b*(r-Ri)) using the formula:
    
#         e(r) = (1 - (b/n) * (r-Ri))^n
        
#     This form has three crucial properties:
#     1. It uses the user-provided INTEGER exponent 'n'.
#     2. It correctly uses the physical Morse parameter 'b'.
#     3. It guarantees that at r=Ri, e=1 (correct minimum) and |de/dr|=b (correct curvature).
#     """
#     # Define the base of the polynomial. This correctly incorporates b and n.
#     base = 1 - (b / n) * (r - Ri)
#     Rc    =  
#     e     =  np.maximum(0,base)**n
#     E     =  Ei * (e*e - 2*e)
#     de_dr = -b * np.maximum(0, base)**(n-1)
#     F     = -Ei * (2*e - 2) * de_dr.
#     mask = r >= Rc
#     E[mask] = 0
#     F[mask] = 0
#     e[mask] = 0
    
#     return E, F, e

# def MorseCut_2_( r, Ei=0.1, Ri=1.4, b=1.6, Rc=6.0, n=8 ):
#     """
#     Approximation based on (r^2 - Ri^2). Correctly uses integer n and physical b.
#     e(r) = (1 - (b / (2*n*Ri)) * (r^2 - Ri^2))^n
#     """
#     base = 1 - (b / (2 * n * Ri)) * (r**2 - Ri**2)
#     e = np.power(np.maximum(0, base), n)
#     E = Ei * (e*e - 2*e)
#     de_dr = -(b * r / Ri) * np.power(np.maximum(0, base), n - 1)
#     F = -Ei * (2*e - 2) * de_dr
#     mask = r >= Rc
#     E[mask] = 0
#     F[mask] = 0
#     e[mask] = 0
#     return E, F, e


def MorseCut_1_( r, Ei=0.1, Ri=1.4, b=1.6, n=8 ):
    """
    Approximation based on (r - Ri) with automatic Rc calculation.
    Rc is determined by the condition that the polynomial base is zero at the cutoff.
    """
    Rc    = Ri + n / b
    base  = 1 - (b / n) * (r - Ri)
    base  = np.maximum(0, base)
    e     = base**n
    E     = Ei * (e*e - 2*e)
    de_dr = -b * np.power(base, n - 1)
    F     = -Ei * (2*e - 2) * de_dr
    return E, F, e

def MorseCut_2_( r, Ei=0.1, Ri=1.4, b=1.6, n=8 ):
    """
    Approximation based on (r^2 - Ri^2) with automatic Rc calculation.
    """
    Rc    = np.sqrt(Ri**2 + (2 * n * Ri) / b)
    base  = 1 - (b / (2 * n * Ri)) * (r**2 - Ri**2)
    base  = np.maximum(0, base)
    e     = base**n
    E     = Ei * (e*e - 2*e)
    de_dr = -(b * r / Ri) * np.power(base, n - 1)
    F     = -Ei * (2*e - 2) * de_dr
    return E, F, e


def get_cubic_coeffs(Ei, Rmin, Rcut):
    A = np.array([
        [  Rmin**3 ,   Rmin**2 , Rmin, 1],  # P(Rmin) = -Ei
        [3*Rmin**2 , 2*Rmin    ,    1, 0],  # P'(Rmin) = 0
        [  Rcut**3 ,   Rcut**2 , Rcut, 1],  # P(Rcut) = 0
        [3*Rcut**2 , 2*Rcut    ,    1, 0]   # P'(Rcut) = 0
    ])
    B = np.array([-Ei, 0, 0, 0])
    coeffs = np.linalg.solve(A, B)
    return coeffs

def get_quartic_coeffs(Ei, Rmin, Rcut, b):
    """
    Solves for the 5 coefficients of the BASE quartic polynomial which satisfies
    the 5 hard constraints: P(Rmin), P'(Rmin), P''(Rmin), P(Rcut), P'(Rcut).
    """
    morse_curvature_at_min = 2 * Ei * b**2
    A = np.array([
        [   Rmin**4,   Rmin**3,   Rmin**2,   Rmin, 1],  # P(Rmin) = -Ei
        [   Rcut**4,   Rcut**3,   Rcut**2,   Rcut, 1],  # P(Rcut) = 0
        [ 4*Rmin**3, 3*Rmin**2, 2*Rmin,     1,     0],  # P'(Rmin) = 0
        [ 4*Rcut**3, 3*Rcut**2, 2*Rcut,     1,     0],  # P'(Rcut) = 0
        [12*Rmin**2, 6*Rmin,     2,         0,     0]   # P''(Rmin) = morse_curvature
    ])
    B = np.array([-Ei, 0, 0, 0, morse_curvature_at_min])
    return np.linalg.solve(A, B)

def get_quintic_coeffs(Ei, Rmin, Rcut, b):
    morse_curvature_at_min = 2 * Ei * b**2
    A = np.array([
        [    Rmin**5 ,    Rmin**4,   Rmin**3,   Rmin**2, Rmin, 1 ], # P(Rmin) = -Ei
        [    Rcut**5 ,    Rcut**4,   Rcut**3,   Rcut**2, Rcut, 1 ], # P(Rcut) = 0
        [  5*Rmin**4 ,  4*Rmin**3, 3*Rmin**2, 2*Rmin   , 1,    0 ], # P'(Rmin) = 0
        [  5*Rcut**4 ,  4*Rcut**3, 3*Rcut**2, 2*Rcut   , 1,    0 ], # P'(Rcut) = 0
        [ 20*Rmin**3 , 12*Rmin**2, 6*Rmin,            2, 0,    0 ], # P''(Rmin) = morse_curvature
        [ 20*Rcut**3 , 12*Rcut**2, 6*Rcut,            2, 0,    0 ]  # P''(Rcut) = 0
    ])
    B = np.array([-Ei, 0, 0, 0, morse_curvature_at_min, 0])
    coeffs = np.linalg.solve(A, B)
    return coeffs

def shape_function(r, Rmin, Rcut):
    """
    Calculates the special quintic Hermite basis function and its derivative.
    H(r) = (r - Rmin)^3 * (r - Rc)^2
    This function is zero and has zero derivatives at all the constraint points,
    allowing us to add it without violating the hard constraints.
    """
    term1 = (r - Rmin)
    term2 = (r - Rcut)
    H = (term1**3) * (term2**2)
    # Derivative dH/dr = 3*(r-Rmin)^2*(r-Rc)^2 + 2*(r-Rmin)^3*(r-Rc)
    dH = 3 * (term1**2) * (term2**2) + 2 * (term1**3) * term2
    return H, dH

def Morse_cubic(r, Ei=0.1, Ri=1.4, Rc=6.0, b=1.6):
    coeffs = get_cubic_coeffs(Ei, Ri, Rc)
    deriv_coeffs = coeffs[:-1] * np.array([3, 2, 1])
    E =  np.polyval(coeffs, r)
    F = -np.polyval(deriv_coeffs, r) # Force is -dE/dr
    mask = (r < Ri) | (r > Rc)
    E[mask] = 0
    F[mask] = 0
    return E, F

def Morse_quintic(r, Ei=0.1, Ri=1.4, Rc=6.0, b=1.6):
    """Evaluates the quintic potential defined by the boundary conditions."""
    coeffs = get_quintic_coeffs(Ei, Ri, Rc, b)
    deriv_coeffs = coeffs[:-1] * np.array([5, 4, 3, 2, 1])
    E = np.polyval(coeffs, r)
    F = -np.polyval(deriv_coeffs, r)
    mask = (r < Ri) | (r > Rc)
    E[mask], F[mask] = 0, 0
    return E, F

def Morse_quintic_opt(r, Ei, Ri, Rc, b):
    """
    Constructs an optimized quintic potential using a direct NumPy least-squares fit.
    """
    # 1. Get the base quartic polynomial satisfying the 5 hard constraints
    base_coeffs = get_quartic_coeffs(Ei, Ri, Rc, b)
    base_deriv_coeffs = base_coeffs[:-1] * np.array([4, 3, 2, 1])

    # 2. Define the sample points inside the valid range for the fit
    fit_mask = (r >= Ri) & (r <= Rc)
    r_fit = r[fit_mask]

    # 3. Define the target vector 'b' for the least squares problem.
    #    This is the difference between the true potential and our base polynomial.
    morse_E_target, _, _ = Morse(r_fit, Ei=Ei, Ri=Ri, b=b)
    base_E_fit = np.polyval(base_coeffs, r_fit)
    b_vector = morse_E_target - base_E_fit

    # 4. Define the design matrix 'A' for the least squares problem.
    #    This is our single basis function (the shape function).
    #    It must be a column vector, so we reshape.
    shape_H_fit, _ = shape_function(r_fit, Ri, Rc)
    A_matrix = shape_H_fit.reshape(-1, 1)

    # 5. Solve the linear least squares problem A*x = b for x=[alpha].
    #    rcond=None is recommended for modern NumPy versions.
    result = np.linalg.lstsq(A_matrix, b_vector, rcond=None)
    optimal_alpha = result[0][0] # The solution is the first element of the result tuple
    #optimal_alpha = 0
    print(f"[INFO] NumPy linalg.lstsq found optimal alpha = {optimal_alpha:.4f}")

    # 6. Construct the final potential and force using the optimal alpha.
    #    Evaluate on the original full 'r' array.
    base_E_full = np.polyval(base_coeffs, r)
    base_F_full = -np.polyval(base_deriv_coeffs, r)
    shape_H_full, shape_dH_full = shape_function(r, Ri, Rc)
    
    E = base_E_full + optimal_alpha * shape_H_full
    F = base_F_full - optimal_alpha * shape_dH_full

    # 7. Apply the mask outside the valid region
    final_mask = (r < Ri) | (r > Rc)
    E[final_mask], F[final_mask] = 0, 0
    
    return E, F, None


def Morse_exp_cubic(r, Ei=0.1, Ri=1.4, b=1.6, r1=1.4, r2=6.0):
    """
    Approximates the Morse potential by first creating a cubic polynomial fit
    for the exponential term e(r) = exp(-b*(r - Ri)).
    
    The cubic's coefficients are found by matching the value and derivative of
    the true exponential at two specified ABSOLUTE coordinate points, r1 and r2.
    """
    # 1. Define the target values and derivatives at the two absolute fit points
    e1 = np.exp(-b * (r1 - Ri))
    d1 = -b * e1  # This is de/dr at r1
    e2 = np.exp(-b * (r2 - Ri))
    d2 = -b * e2  # This is de/dr at r2
    
    # 2. Set up the 4x4 linear system Ax = B to solve for the cubic coefficients [c3,c2,c1,c0]
    #    The polynomial is now a function of r: P(r) = c3*r^3 + ...
    A = np.array([
        [r1**3,   r1**2,  r1, 1],
        [3*r1**2, 2*r1,   1,  0],
        [r2**3,   r2**2,  r2, 1],
        [3*r2**2, 2*r2,   1,  0]
    ])
    B = np.array([e1, d1, e2, d2])
    
    # 3. Solve for the polynomial coefficients
    coeffs = np.linalg.solve(A, B)
    deriv_coeffs = coeffs[:-1] * np.array([3, 2, 1])

    # 4. Evaluate the polynomial P(r) to get e_approx
    e_approx = np.polyval(coeffs, r)
    
    # 5. Calculate the final Energy using the standard formula
    E = Ei * (e_approx**2 - 2 * e_approx)
    
    # 6. Calculate the Force using the chain rule: F = -dE/dr = - (dE/de * de/dr)
    dE_de = Ei * (2 * e_approx - 2)
    de_dr = np.polyval(deriv_coeffs, r) # de/dr is the derivative of the polynomial P(r)
    F = -dE_de * de_dr
    
    return E, F, e_approx

# def Morse_pow_cubic(r, Ei, Ri, b, n, r1, r2):
#     """
#     Approximates the Morse potential by first creating a cubic polynomial p(r)
#     that fits the n-th root of the exponential term, g(r) = exp(-b*(r-Ri)/n).
    
#     The final exponential is then reconstructed as e_approx = p(r)^n.
#     This "power-transform" method provides a much more stable and accurate fit.
#     """
#     # 1. Define the target function g(r) and its derivative g'(r)
#     #    g(r) is the n-th root of the true exponential term.
#     def g(r_val):
#         return np.exp(-b * (r_val - Ri) / n)
        
#     def g_prime(r_val):
#         return (-b / n) * g(r_val)

#     # 2. Get target values and derivatives at the two absolute fit points
#     g1 = g(r1)
#     d1 = g_prime(r1)
#     g2 = g(r2)
#     d2 = g_prime(r2)
    
#     # 3. Set up and solve the 4x4 linear system for the cubic coefficients [c3,c2,c1,c0]
#     #    of our polynomial p(r).
#     A = np.array([
#         [r1**3,   r1**2,  r1, 1],
#         [3*r1**2, 2*r1,   1,  0],
#         [r2**3,   r2**2,  r2, 1],
#         [3*r2**2, 2*r2,   1,  0]
#     ])
#     B = np.array([g1, d1, g2, d2])
    
#     coeffs = np.linalg.solve(A, B)
#     deriv_coeffs = coeffs[:-1] * np.array([3, 2, 1])

#     # 4. Evaluate the base polynomial p(r) and its derivative p'(r)
#     p_r = np.polyval(coeffs, r)
#     p_prime_r = np.polyval(deriv_coeffs, r)
    
#     # 5. Reconstruct the full exponential term via integer exponentiation
#     e_approx = p_r**n
    
#     # 6. Calculate the final Energy
#     E = Ei * (e_approx**2 - 2 * e_approx)
    
#     # 7. Calculate the Force using the chain rule: F = - (dE/de * de/dr)
#     #    where de/dr = d/dr(p(r)^n) = n * p(r)^(n-1) * p'(r)
#     dE_de = Ei * (2 * e_approx - 2)
#     de_dr = n * (p_r**(n - 1)) * p_prime_r
#     F = -dE_de * de_dr
    
#     # We return the base polynomial p(r) as the third output for visualization
#     return E, F, p_r


def Morse_pow_cubic(r, Ei, Ri, b, n, r1, r2):
    """
    This is the most advanced version. It creates a cubic polynomial p(r) to
    fit the n-th root of the exponential term.
    
    It then CALCULATES the cutoff Rc by finding the root of p(r)=0. This
    ensures the potential and force go perfectly to zero at the cutoff.
    """
    # 1. Define the target function g(r) = exp(-b*(r-Ri)/n) and its derivative
    def g(r_val): return np.exp(-b * (r_val - Ri) / n)
    def g_prime(r_val): return (-b / n) * g(r_val)

    # 2. Get target values and derivatives at the fit points
    g1, d1 = g(r1), g_prime(r1)
    g2, d2 = g(r2), g_prime(r2)
    
    # 3. Set up and solve the 4x4 linear system for the coefficients of p(r)
    A = np.array([
        [r1**3,   r1**2,  r1, 1], [3*r1**2, 2*r1,   1,  0],
        [r2**3,   r2**2,  r2, 1], [3*r2**2, 2*r2,   1,  0]
    ])
    B = np.array([g1, d1, g2, d2])
    coeffs = np.linalg.solve(A, B)
    
    # 4. === Find the natural cutoff Rc ===
    # Find all roots (real and complex) of the polynomial p(r)
    all_roots = np.roots(coeffs)
    # Filter for only the real roots
    real_roots = all_roots[np.isreal(all_roots)].real
    # Filter for roots that are physically meaningful (must be > Ri)
    physical_roots = real_roots[real_roots > Ri]
    
    if len(physical_roots) == 0:
        raise ValueError("Could not find a physical cutoff root for the given parameters.")
    
    # The cutoff is the largest valid physical root.
    Rc = np.max(physical_roots)

    global Rcs
    Rcs.append(Rc)

    # 5. Evaluate the base polynomial p(r) and its derivative
    deriv_coeffs = coeffs[:-1] * np.array([3, 2, 1])
    p_r = np.polyval(coeffs, r)
    p_prime_r = np.polyval(deriv_coeffs, r)
    
    # Clamp p(r) at zero. This automatically handles the cutoff smoothly.
    p_r_clamped = np.maximum(0, p_r)
    
    # 6. Reconstruct the full exponential term and calculate Energy
    e_approx = p_r_clamped**n
    E = Ei * (e_approx**2 - 2 * e_approx)
    
    # 7. Calculate the Force using the chain rule
    dE_de = Ei * (2 * e_approx - 2)
    de_dr = n * (p_r_clamped**(n - 1)) * p_prime_r
    F = -dE_de * de_dr
    
    # Pass the calculated Rc back out for plotting/information
    return E, F, p_r


def MorseCut_3( r, Ei=0.1, Ri=1.4, b=1.6, Rc=6.0 ):
    r2 = r*r
    x1 = 1-(r2/(Rc*Rc))
    x2 = 1-(r2/(0.7*Ri*Ri))
    x3 = 1-(r2/(2.0*Rc*Rc))
    mask = r2>(Rc*Rc)
    x1[mask] = 0
    #E =  5*Ei* x1*x1*x1*x1*( -1 + x2 )
    E =  35*Ei*x2*(x1**2)*(x3**16)
    #E =  10*Ei*x2*(x1**2)*(x3**4)
    F =  0
    return E,F

def MorseCut_4( r, Ei=0.1, Ri=1.4, b=1.6, Rc=6.0 ):
    r2 = r*r
    x1 = 1-(r2/(Rc*Rc))
    ir2 = (Ri*Ri*1.3)/r2 
    ir4 = ( (Ri*Ri*1.05)/r2 )**2
    #x2 = ir2*ir2 - 2*ir2
    x2 = ir4*ir4 - 2*ir4
    mask = r2>(Rc*Rc)
    x1[mask] = 0
    E =  1.8*Ei*(x1**2)*x2
    #E =  2.0*Ei*(x1**2)*x2*x2
    F =  0
    return E,F

if __name__ == "__main__":

    import func_utils as fu

    global Rcs
    Rcs = []

    xs = np.linspace(2.0,10,1000)

    
    Rc = 5.0
    Ri = 1.4+1.6
    Ei = 1.0
    b  = 1.6

    y_ref = Morse(xs, Ei, Ri, b)

    fig,axs = fu.plot_funcs(
        [
            ('Morse'   ,   Morse     ,'k', {'Ei':Ei, 'Ri':Ri, 'b': b  },   None      ),
            # ('MorseCut_1', MorseCut_1,'b', {'Ei':Ei, 'Ri':Ri, 'b': b, 'Rc': 6.0, 'n': 8} ),
            # ('MorseCut_2', MorseCut_2,'g', {'Ei':Ei, 'Ri':Ri, 'b': b, 'Rc': 5.2, 'n': 8} ),

            #('MorseCut_1', MorseCut_1_,'b', {'Ei':Ei, 'Ri':Ri, 'b': b, 'n': 8} ),
            #('MorseCut_2', MorseCut_2_,'g', {'Ei':Ei, 'Ri':Ri, 'b': b, 'n': 8} ),

            #('Morse_exp_cubic', Morse_exp_cubic,'m', {'Ei':Ei, 'Ri':Ri, 'b': b, 'r1': Ri, 'r2': Ri+2.0 } ),
            ('Morse_pow_cubic', Morse_pow_cubic,'r', {'Ei':Ei, 'Ri':Ri, 'b': b, 'r1': Ri, 'r2': Ri+0.5, 'n': 5 }, y_ref ),
            ('Morse_pow_cubic', Morse_pow_cubic,'c', {'Ei':Ei, 'Ri':Ri, 'b': b, 'r1': Ri, 'r2': Ri+1.0, 'n': 5 }, y_ref ),
            ('Morse_pow_cubic', Morse_pow_cubic,'g', {'Ei':Ei, 'Ri':Ri, 'b': b, 'r1': Ri, 'r2': Ri+2.0, 'n': 5 }, y_ref ),

            # ('Morse_cubic',   Morse_cubic,      'b', {'Ei':Ei, 'Ri':Ri, 'b': b, 'Rc': Rc } ),
            # ('Morse_quintic', Morse_quintic,    'g', {'Ei':Ei, 'Ri':Ri, 'b': b, 'Rc': Rc } ),
            # ('Morse_quintic', Morse_quintic_opt,'m', {'Ei':Ei, 'Ri':Ri, 'b': b, 'Rc': Rc } ),

            #('MorseCut_3', MorseCut_3,'r', {'Ei': Ei, 'Ri': Ri, 'Rc': 6.0}         ),
            #('MorseCut_4', MorseCut_4,'m', {'Ei': Ei, 'Ri': Ri, 'Rc': 6.0}         ),
        ],
        xs,
        bError=True,
        figsize=(10,15),
        bNum=True,
    )

    colors=['r','c','g']
    for i,Rc in enumerate(Rcs):
        axs[0].axvline(Rc, color=colors[i], linestyle='--')
        axs[1].axvline(Rc, color=colors[i], linestyle='--')

    axs[0].set_ylim(-1.0,1.0)
    axs[1].set_ylim(-1.0,1.0)
    #axs[2].set_ylim(-1.0,1.0)
    plt.show()
