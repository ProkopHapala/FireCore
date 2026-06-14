# Maxima Tutorial: Symbolic Computation for Scientists

This tutorial introduces Maxima, a powerful open-source computer algebra system, focusing on commands and techniques relevant to scientific computations, particularly those encountered in physics and chemistry. We'll draw examples from typical calculations involving Gaussian functions, special functions, and numerical methods.

## 1. Getting Started with Maxima

### 1.1. Basic Interaction

run script:
`> maxima -b GaussCoulombBoys.mac`

*   **Commands**: You type commands followed by a semicolon (`;`) or a dollar sign (`$`).
    *   Semicolon (`;`): Executes the command and displays the result (e.g., `(%o1)`).
    *   Dollar sign (`$`): Executes the command but suppresses the output. This is useful for assignments or intermediate steps where you don't need to see the immediate result.
*   **Comments**: Use `/* ... */` for multi-line comments.
*   **Clearing Memory**: `kill(all);` removes all user-defined variables and functions from memory. It's good practice to start scripts with this.

```maxima
/* This is a comment */
kill(all)$ /* Clears memory, output suppressed */
a: 5;    /* Assign 5 to 'a' and show output */
b: 10$   /* Assign 10 to 'b', suppress output */
a+b;
```

### 1.2. Assignments

Use the colon (`:`) for assignment.

```maxima
wi: 1.0;
wj: 0.5;
r: 2.0;
```

### 1.3. Assumptions

The `assume` command tells Maxima about properties of variables, which can help simplify expressions or enable certain integrations.

```maxima
assume(wi > 0, wj > 0, r > 0, si > 0, sj > 0);
/* Maxima now knows these variables are positive */
```

## 2. Defining Functions

Use `:=` to define functions.

```maxima
/* Define a Gaussian function fi(x,y,z) centered at (r,0,0) with exponent wi */
fi(x, y, z) := exp(-wi*((x-r)^2 + y^2 + z^2));

/* Define another Gaussian fj(x,y,z) centered at origin with exponent wj */
fj(x, y, z) := exp(-wj*(x^2 + y^2 + z^2));

/* Call the function */
fi(1,0,0);
```

## 3. Symbolic Calculus

### 3.1. Differentiation

The `diff` command performs differentiation.
`diff(expr, var1, n1, var2, n2, ...)` differentiates `expr` `n1` times with respect to `var1`, then `n2` times with respect to `var2`, etc.

```maxima
/* Define the Laplacian of fj */
/* Note: 'diff shows the unevaluated form initially if not forced by ev or ' */
Lfj_expr := diff(fj(x,y,z), x, 2) + diff(fj(x,y,z), y, 2) + diff(fj(x,y,z), z, 2);

/* To evaluate it immediately, you can use ' before diff or ev() */
Lfj(x,y,z) := ev(diff(fj(x,y,z), x, 2) + diff(fj(x,y,z), y, 2) + diff(fj(x,y,z), z, 2));
Lfj(x,y,z);

/* Example from eFF Pauli potential: derivative of T with respect to si */
/* First, define T (assuming si, sj, r are defined and have assumptions) */
si2 : si*si;
sj2 : sj*sj;
si2sj2 : si2 + sj2;
invsi2sj2 :  1/si2sj2;
T : (3/2)*( si2sj2/(si*si*sj*sj)  )   -   ( 6*si2sj2-4*r*r )*invsi2sj2*invsi2sj2;

dT_dsi : diff( T, si );
dT_dsi_simplified : factor(ratsimp(dT_dsi)); /* Simplify the result */
print("dT/dsi = ", dT_dsi_simplified)$
```

### 3.2. Integration

The `integrate` command performs integration.
`integrate(expr, var, low_limit, high_limit)` for definite integrals.
`integrate(expr, var)` for indefinite integrals.

For symbolic infinity, use `inf` and `minf`.

```maxima
/* Product of two Gaussians */
sij_product(x,y,z) := fi(x,y,z) * fj(x,y,z);

/* Overlap integral Sij = Integral(fi*fj dV) */
/* Note: This assumes unnormalized Gaussians as per GaussianIntegrals.md */
Sij_integral_x : integrate(sij_product(x,y,z), x, minf, inf);
Sij_integral_xy : integrate(Sij_integral_x, y, minf, inf);
Sij : integrate(Sij_integral_xy, z, minf, inf);
Sij_simplified : factor(ratsimp(Sij));
print("Sij = ", Sij_simplified)$

/* Kinetic energy integral Tij = Integral(fi * Laplacian(fj) dV) */
tij_integrand(x,y,z) := fi(x,y,z) * Lfj(x,y,z); /* Lfj defined earlier */
Tij_integral_x : integrate(tij_integrand(x,y,z), x, minf, inf);
Tij_integral_xy : integrate(Tij_integral_x, y, minf, inf);
Tij : integrate(Tij_integral_xy, z, minf, inf);
Tij_simplified : factor(ratsimp(Tij));
print("Tij = ", Tij_simplified)$
```

### 3.3. Simplification

Maxima offers several functions to simplify expressions:
*   `ratsimp(expr)`: Simplifies `expr` to a canonical form (ratio of two polynomials).
*   `factor(expr)`: Factors `expr`.
*   `expand(expr)`: Expands `expr`.
*   `trigsimp(expr)`: Simplifies trigonometric functions in `expr`.
*   `radcan(expr)`: Simplifies expressions involving radicals, logarithms, and exponentials.

```maxima
tau_expr : Tij_simplified / Sij_simplified;
tau_simplified : factor(ratsimp(tau_expr));
print("Tij/Sij = ", tau_simplified)$

/* Example with trigsimp (from spherical harmonics) */
/* assoc_legendre_p(l,m,x) gives the associated Legendre polynomial P_l^m(x) */
P_1_1_cos_theta : assoc_legendre_p(1, 1, cos(theta));
print("P_1^1(cos(theta)) before trigsimp: ", P_1_1_cos_theta)$
P_1_1_cos_theta_simplified : trigsimp(P_1_1_cos_theta);
print("P_1^1(cos(theta)) after trigsimp: ", P_1_1_cos_theta_simplified)$ /* Should be -sin(theta) */
```

### 3.4. Substitution

*   `subst(new, old, expr)`: Substitutes `new` for `old` in `expr`.
*   `subst([eq1, eq2, ...], expr)`: Substitutes using a list of equations.
*   `ratsubst(new, old, expr)`: Similar to `subst`, but performs rational substitution. Useful when `old` is not a simple variable but an expression.

```maxima
/* Transformation from width 'w' to size 's' for Gaussians (wi = 1/(2*si^2)) */
eq_wi : wi = 1/(2*si^2);
eq_wj : wj = 1/(2*sj^2);

tau_s : factor(ratsimp(subst([eq_wi, eq_wj], tau_simplified)));
print("tau_s (in terms of s) = ", tau_s)$

Sij_s : factor(ratsimp(subst([eq_wi, eq_wj], Sij_simplified)));
print("Sij_s (in terms of s) = ", Sij_s)$

/* Example of ratsubst from Lorenz-like functions */
Lf : 1/(w*x^2 + 1);
f_expr : (1-x^2)/(w*x^2+1); /* This is Lf*(1-x^2) */
df_expr : diff(f_expr, x);
/* Substitute Lf back into the derivative to express it in terms of Lf */
df_s : factor(ratsubst(Lf_symbolic, Lf, df_expr)); /* Lf_symbolic is just a symbol */
/* print("Derivative in terms of Lf_symbolic: ", df_s)$ */
/* Note: For this to work as intended in the example, Lf_symbolic would be 'L' in the original context. */
```

### 3.5. Limits

The `limit(expr, var, val)` command computes the limit of `expr` as `var` approaches `val`.

```maxima
/* Example: Limit of Euu_sisj as r approaches 0 */
/* Define S22 and T first (simplified for brevity) */
S22_r0_sisj : 1; /* Example value for S22 when r=0, si=sj */
T_r0_sisj : 3/(2*sj^2); /* Example value for T when r=0, si=sj */
rho : 'rho; /* Make rho a symbolic variable */

Euu_sisj_r0_expr : S22_r0_sisj * T_r0_sisj * ( -rho*S22_r0_sisj + (rho - 2)  )/( S22_r0_sisj^2 - 1  );
/* This specific expression has S22^2-1 in denominator, which is 0 if S22=1.
   The original GaussianIntegrals.md uses a more complex Euu_sisj expression
   and then takes the limit. Let's use a placeholder for Euu_sisj for demonstration.
*/
Euu_sisj(r, sj, rho_val) := ( (3/sj^2-(12*sj^2-4*r^2)/(4*sj^4))*exp(-(r^2/sj^2))*(-(rho_val*exp(-(r^2/sj^2)))+rho_val-2) ) / (exp(-(2*r^2/sj^2))-1);

Euu_sisj_r0_limit : limit(Euu_sisj(r, sj, rho), r, 0 );
print("Limit of Euu_sisj as r->0: ", Euu_sisj_r0_limit)$ /* Should be 1/sj^2 */
```

### 3.6. Evaluation with `ev`

`ev(expr, var1=val1, var2=val2, ...)` evaluates `expr` by substituting `val1` for `var1`, etc. It can also take evaluation flags.

```maxima
/* From Lorenz-like functions: evaluate derivative at x=1 */
f_symb : (1-x^2)/(w*x^2+1);
df : diff(f_symb, x);
df1 : ev(df, x=1);
print("Derivative of f_symb at x=1: ", ratsimp(df1))$
```

## 4. Advanced Techniques

### 4.1. `block` for Local Variables and Multi-statement Functions

`block([var1, var2, ...], expr1, expr2, ..., exprN)` creates a block of code. `var1, var2, ...` are local variables within the block. The value of the block is the value of `exprN`.

```maxima
/* Example: Spherical Harmonics function from GaussianIntegrals.md */
/* Note: Maxima's assoc_legendre_p(l,m,x) might have different conventions.
   The original script uses trigsimp(assoc_legendre_p(l, m, cos(theta))).
*/
sph(l, m) := block(
    [P_lm_val, A_norm], /* Local variables */
    assume(theta > 0, sin(theta) > 0), /* Assumption local to this call if not global */
    A_norm: sqrt((2*l + 1) * factorial(l - abs(m)) / (4 * %pi * factorial(l + abs(m)))),
    P_lm_val: ratsimp(trigsimp(assoc_legendre_p(l, abs(m), cos(theta)))),
    /* The original script uses printf for output within the block.
       A function typically returns a value. Let's return the normalized P_lm_val.
    */
    /* printf(true, "Normalization Factor A: ~a~%", A_norm), */
    /* printf(true, "P_lm(cos(theta)): ~a~%", P_lm_val), */
    return(A_norm * P_lm_val * exp(%i*m*phi)) /* Adding phi dependence for Ylm */
);

Y10 : sph(1,0);
print("Y_1,0 proportional to: ", Y10)$
```

### 4.2. Recursive Functions: B-Splines

Maxima can define recursive functions. This is powerfully demonstrated in the B-spline basis function definition.

```maxima
/* Recursive function to define the B-spline basis functions with uniform knots */
/* (Adapted from GaussianIntegrals.md) */
bspline_basis(i, p, t) := block(
    [term1_val, term2_val, numer_val], /* Local variables */

    /* Base case: degree 0 */
    if p = 0 then (
        if i <= t and t < i+1 then 1 else 0
    )
    else (
        /* Recursive case: degree p */
        /* Denominators are p in the Cox-de Boor formula for uniform knots */
        term1_val : (t - i)/p * bspline_basis(i, p-1, t),
        term2_val : (i+p+1 - t)/p * bspline_basis(i+1, p-1, t),

        numer_val : term1_val + term2_val,
        return(numer_val)
    )
);

/* Test the B-spline function */
/* B_0,1(t) for t in [0,1) should be t, for t in [1,2) should be 2-t */
B_0_1_at_0_5 : bspline_basis(0, 1, 0.5); /* (0.5-0)/1 * B(0,0,0.5) + (0+1+1-0.5)/1 * B(1,0,0.5) = 0.5*1 + 1.5*0 = 0.5 */
print("B_0,1(0.5) = ", B_0_1_at_0_5)$

B_0_1_at_1_5 : bspline_basis(0, 1, 1.5); /* (1.5-0)/1 * B(0,0,1.5) + (0+1+1-1.5)/1 * B(1,0,1.5) = 1.5*0 + 0.5*1 = 0.5 */
print("B_0,1(1.5) = ", B_0_1_at_1_5)$

/* Quadratic B-spline B_0,2(t) */
/* This will expand into a piecewise polynomial */
B_0_2_expr : ratsimp(bspline_basis(0, 2, t_var));
/* To see its structure, you might plot it or evaluate at specific points. */
/* plot2d(B_0_2_expr, [t_var, 0, 3]); */
print("B_0,2(t) expression: ", B_0_2_expr)$
```

### 4.3. Processing Lists of Functions with `map` and `lambda`

The `map(function, list)` command applies `function` to each element of `list`.
`lambda([arg1, ...], body)` creates an anonymous function.
This is very useful for applying a complex processing step to multiple items.

```maxima
/* Example: Processing Lorenz-like functions from GaussianIntegrals.md */

/* Define a processing function (simplified from original) */
process_my_function(f_symbol_name, var) := block(
    [f_expr, df_expr, ddf_expr, df_at_1, ddf_at_0], /* Local variables */

    f_expr : ev(f_symbol_name), /* Get the expression associated with the symbol */
    df_expr : diff(f_expr, var),
    ddf_expr : diff(f_expr, var, 2),

    df_at_1 : ratsimp(ev(df_expr, var=1)),
    ddf_at_0 : ratsimp(ev(ddf_expr, var=0)),

    printf(true, "==== Processing: ~a ====", f_symbol_name),
    printf(true, "f(~a): ~a", var, f_expr),
    printf(true, "df/d~a (~a=1): ~a", var, var, df_at_1),
    printf(true, "d2f/d~a2 (~a=0): ~a", var, var, ddf_at_0),
    printf(true, "========================~%"), /* ~% is newline */

    return([f_symbol_name, df_at_1, ddf_at_0])
);

/* Define some Lorenz-like functions (using 'w' as a global variable) */
assume(w>0)$
L1x2 : (1-x^2)/(w*x^2+1);
L0x4 : (1-x^2)^2/(w*x^2+1);

list_of_functions_to_process : ['L1x2, 'L0x4]; /* List of symbols */

/* Use map and lambda to apply process_my_function to each symbol in the list */
processed_results : map(lambda([func_symb], process_my_function(func_symb, x)), list_of_functions_to_process);

print("Processed Results List: ", processed_results)$
```

### 4.4. Plotting

*   `plot2d(expr, [var, low, high], ...options...)`: For 2D plots.
*   `wxplot2d(...)`: Similar, but often uses the WxMaxima interface for interactive plots.
*   You can plot multiple functions: `plot2d([expr1, expr2], [var, low, high])`.

```maxima
/* Plotting B-splines (example from previous section) */
/* Define B1_1 and B1_2 for plotting */
B1_1(t_var) := bspline_basis(0,1,t_var);
B1_2(t_var) := bspline_basis(1,1,t_var);

/* plot2d([B1_1(t_val), B1_2(t_val)], [t_val, 0.0, 3.0], [legend, "B_0,1(t)", "B_1,1(t)"]); */
/* Note: Plotting commands execute and open a plot window or embed in interfaces like WxMaxima.
   In batch mode, they might save to a file if configured or do nothing visible.
   The actual plot command is commented out to prevent issues in non-graphical execution.
   To see plots, run Maxima interactively or use WxMaxima.
*/

/* Example from Lorenz-like functions in GaussianIntegrals.md */
/*
block([w: 10.0],
    wxplot2d([L1x2, L0x4], [x, -0.2, 1.2], [legend, "L1x2", "L0x4"])
);
*/
```

## 5. Controlling Output

*   `display2d: false;` or `display2d: true;`: Toggles between 1D (plain text) and 2D (pretty-printed) output.
*   `linel: 120;`: Sets the line length for 1D output.
*   `echo: false;` or `echo: true;`: Toggles echoing of input commands in batch files (though `maxima -b` often has its own echoing behavior).
*   `printf(true, "format_string", arg1, arg2, ...);`: For formatted C-style printing.
    *   `~a` for general argument, `~f` for float, `~s` for string, `~%` for newline.

```maxima
display2d: false$ /* Switch to plain text output */
linel: 100$     /* Set line length */

my_var : %pi^2 / 6;
printf(true, "The value of my_var is approximately ~6,4f.~%And as an expression: ~a~%", float(my_var), my_var)$

display2d: true$ /* Switch back to pretty-printed output if desired */
```

