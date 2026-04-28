# File: derive_poly_gauss_overlap_separate_cases.py

from pprint import saferepr
import sympy as sp

# infinite line limit, not terminal line limit


# sp.init_printing(use_unicode=True) # Disable pretty printing for standard output

def python_no_junk(expr):
    s = sp.printing.python(expr)
    return s.split("\n")[-1][3:]

def evaluate_integral_for_case(integrand, y_var, lower_limit, upper_limit, bAsPoly=True, bCSE=True ):
    """
    Evaluates the integral for a given case, simplifies it, and prints CSE results.
    """
    print("\n" + "-" * 60)
    print(f"Integral in limits: ({lower_limit}, {upper_limit})")
    
    I = sp.integrate(integrand, (y_var, lower_limit, upper_limit))
    I = sp.expand(I)
    I = sp.simplify(I)

    if bAsPoly:
        I = sp.expand(I)
        # Collect terms with respect to x
        I_collected = sp.collect(I, x)
        # Get the Python code representation as a string
        #s_collected = python_no_junk(I_collected)
        # Print coefficients of x
        print("\n as polynomial of sum_k{ c_k * x^k } :")
        if I_collected.is_polynomial(x):
            poly_in_x = sp.Poly(I_collected, x)
            coeffs = poly_in_x.coeffs()
            all_monoms = poly_in_x.monoms() # gives tuples like (power_of_x,)
            for i, coeff in enumerate(coeffs):
                power = all_monoms[i][0] # Get the power of x
                print(f"  c_{power} = {python_no_junk(coeff)}")
        else:
            print("  Expression is not a simple polynomial in x after collection.")

    if bCSE:
        I = sp.factor(I)
        #I = sp.simplify(I)
        print("\nCommon Subexpression Elimination (CSE):")
        subexprs, cse_result = sp.cse(I) # Perform CSE on the collected form
        for sub_sym, sub_expr_val in subexprs:
            print(f"  {sub_sym.name} = ", end="")
            print(python_no_junk(sub_expr_val))
        print( "I_CSE = ", python_no_junk(cse_result[0]) ) # cse_result is a list
    
    return I

if __name__ == "__main__":
    """
    Symbolically derives the overlap integral S(x) = Int[phi1(y) * phi2(y-x) dy] for
    unnormalized polynomial functions: (w^2 - y^2)^n, by explicitly evaluating
    the integral over distinct overlap regions for x >= 0, assuming w1 > w2.

    phi1 is (w1^2 - y^2)^n, centered at 0. Support: (-w1, w1).
    phi2 is (w2^2 - (y-x)^2)^n, centered at x. Support: (x-w2, x+w2).

    Assumptions: x >= 0 and w1 > w2.

    Args:    
        n_val (int): The exponent 'n' in the polynomial definition.

    Returns:
        dict: A dictionary where keys are the x-ranges (as strings) and values
              are the symbolic expressions for the overlap integral in that range.
    """
    sp.init_printing(use_unicode=True)
    # Define symbols
    x, y   = sp.symbols('x y', real=True)
    w1, w2 = sp.symbols('w1 w2', positive=True)

    n = sp.Integer(2)   # result is 9th degree polynomial 
    #n = sp.Integer(3)  # result is 13th degree polynomial
    #n = sp.Integer(4)  # still fine - result is 17-degree polynomial
    #n = sp.symbols('n', integer=True, positive=True)
    #n = sp.Integer(n_val)

    
    # Core integrand expression
    phi1_core_expr = (w1**2 -  y**2      )**n
    phi2_core_expr = (w2**2 - (y - x)**2 )**n
    
    print( "phi1(x) = ", phi1_core_expr )
    print( "phi2(x) = ", phi2_core_expr )

    bCSE = False

    print("\n" + "=" * 60)
    print( "# Overlap integral of functions: S(x,w1,w2) = int_y{ phi1(y,w1) * phi2(y-x,w2) dy} " )
    #W_diff = w1 - w2  # Boundary for containment vs partial overlap
    # W_sum  = w1 + w2  # Boundary for partial overlap vs no overlap
    phi1_phi2 = phi1_core_expr * phi2_core_expr
    evaluate_integral_for_case( phi1_phi2, y, x - w2, x + w2  , bCSE=bCSE )
    evaluate_integral_for_case( phi1_phi2, y, x - w2, w1      , bCSE=bCSE )

    print("\n" + "=" * 60)
    print( "# Kinetic integral of functions: T(x,w1,w2) = int_y{ phi1(y,w1) * d_x^2 phi2(y-x,w2) dy}" )
    phi1_phi2_d2 = sp.diff(phi1_phi2, x, 2)
    evaluate_integral_for_case( phi1_phi2_d2, y, x - w2, x + w2  , bCSE=bCSE )
    evaluate_integral_for_case( phi1_phi2_d2, y, x - w2, w1      , bCSE=bCSE )
