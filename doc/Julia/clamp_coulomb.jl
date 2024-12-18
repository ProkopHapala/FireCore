
using Plots
include("plot_utils.jl")


Coulomb_const_eVA = 14.3996448915 

# ========== Functions

# ============ Coulomb potential ============

function getCoulomb( r, Q )
    Q *= Coulomb_const_eVA
    E  = Q/r
    F  = Q/r^2
    return E,F
end

function getClamped( func::Function, clamp_func::Function, r, params, y1, y2)
    E_, F_ = func(r, params)    
    E,  F  = clamp_func(E_, F_, y1, y2)
    return E, F
end

# ========================

# Soft clamp function for scalars
function soft_clamp(y, dy, y1, y2)
    if y > y1
        y12   = y2 - y1
        invdy = 1.0 / y12
        z     = (y - y1) * invdy
        denom = 1 / (1 + z)
        y_    = y1 + y12 * (1 - denom )
        dy_   = dy * denom*denom
        return y_, dy_
    end
    return y, dy
end

function soft_clamp_neg(y, dy, y1, y2)
    if y < y1
        y12   = y2 - y1
        invdy = 1.0 / y12
        z     = (y - y1) * invdy
        denom = 1 / (1 + z)
        y_    = y1 + y12 * (1 - denom )
        dy_   = dy * denom*denom
        return y_, dy_
    end
    return y, dy
end

function smooth_clamp(y, dy, y1, y2)
    if y > y1
        y21   = y2 - y1
        z     = (y - y1) / y21
        denom = 1/(1 + z + z*z)
        y_    = y2 - y21 * denom
        dy_   = dy * (1 + 2*z)  * denom * denom
        return y_, dy_
    end
    return y, dy
end

function smooth_clamp_neg(y, dy, y1, y2)
    if y < y1
        y21   = (y2 - y1)
        z     = (y - y1) / y21
        denom = 1/(1 + z + z*z)
        y_    = y2 - y21 * denom
        dy_   = dy * (1 + 2*z)  * denom * denom
        return y_, dy_
    end
    return y, dy
end

# function smooth_clamp_neg(y, dy, y1, y2)
#     if y < y1
#         z = (y - y1) / (y1 - y2)
#         denominator = 1 + z + z^2
#         y_clamped = y2 - (y2 - y1) / denominator
#         dy_clamped = dy * (1 + 2z) / (denominator^2)
#         return y_clamped, dy_clamped
#     else
#         return y, dy
#     end
# end

# Soft clamp exponential function for scalars
function soft_clamp_exp(y, min_value, max_value, width=10.0)
    if y > max_value && max_value < Inf
        z = (y - max_value) / width
        y_new = max_value + (y - max_value) / (1 + exp(z))
    elseif y < min_value && min_value > -Inf
        z = (y - min_value) / -width
        y_new = min_value + (y - min_value) / (1 + exp(z))
    else
        y_new = y
    end
    return y_new
end



# ========== Body

# eval_forces = (position, velocity) -> eval_force_and_plot(position,velocity, plt, truss.bonds )
Q = -1.0
xs = xrange( 0.001, 0.01, 1000 )

plt = plot( layout = (2, 1), size=(1000, 1000) )
mins = []

xlim = [1.0, 12.0]

#fcut = (x)->smootherstep(x,Rc0,Rc)

#alpha = -1.5

push!( mins, plot_func( plt,  xs, (x)->getCoulomb( x,  Q ),       clr=:black  ,  label="Q/r"       ,  xlim=xlim, dnum=:true ) )
push!( mins, plot_func( plt,  xs, (x)->getCoulomb( x, -Q ),       clr=:black  ,  label="Q/r"       ,  xlim=xlim, dnum=:true ) )
#push!( mins, plot_func( plt,  xs, (x)->getClamped( getCoulomb, soft_clamp, x, -Q, 5.0, 10.0 ),       clr=:green  ,  label="soft_clamp(Q/r)"       ,  xlim=xlim, dnum=:true ) )
push!( mins, plot_func( plt,  xs, (x)->getClamped( getCoulomb, soft_clamp_neg, x, Q, -5.0, -10.0 ),       clr=:red    ,  label="clamp(Q/r)"       ,  xlim=xlim, dnum=:true ) )
# push!( mins, plot_func( plt,  xs, (x)->getClamped( getCoulomb, soft_clamp_neg, x, Q, -6.0, -10.0 ),       clr=:blue   ,  label="clamp(Q/r)"       ,  xlim=xlim, dnum=:true ) )
# push!( mins, plot_func( plt,  xs, (x)->getClamped( getCoulomb, soft_clamp_neg, x, Q, -4.0, -5.0  ),       clr=:magenta,  label="clamp(Q/r)"       ,  xlim=xlim, dnum=:true ) )

#push!( mins, plot_func( plt,  xs, (x)->getClamped( getCoulomb, smooth_clamp, x, -Q, 5.0, 6.0 ),       clr=:green  ,  label="smooth_clamp(Q/r)"       ,  xlim=xlim, dnum=:true ) )
push!( mins, plot_func( plt,  xs, (x)->getClamped( getCoulomb, smooth_clamp_neg, x, Q, -5.0, -6.0 ),       clr=:green  ,  label="smooth_clamp(Q/r)"       ,  xlim=xlim, dnum=:true ) )
#push!( mins, plot_func( plt,  xs, (x)->pexp4( x, alpha ),       clr=:red   ,  label="pexp4"     ,  xlim=xlim, dnum=:true ) )
#push!( mins, plot_func( plt,  xs, (x)->pexp4_( x, alpha ),      clr=:cyan   ,  label="pexp4_"     ,  xlim=xlim, dnum=:true ) )

#Emin = minimum( [min[1] for min in mins] ); ylims!( plt[1], Emin*2.0*0, -Emin*2. )
#Fmin = minimum( [min[2] for min in mins] ); ylims!( plt[2], Fmin*2.0*0, -Fmin*2. )

#ylims!( plt[1], 0.0, 10.0 ); ylims!( plt[2], 0.0, 10.0 )
ylims!( plt[1], -10.0, 0.0 ); ylims!( plt[2], -10.0, 0.0 )


hline!( plt[1], [0.0],  color=:black, label="", linestyle=:dash )
hline!( plt[2], [0.0],  color=:black, label="", linestyle=:dash )

# vline!( plt[1], [Rc],   color=:gray,  label="", linestyle=:dash ) 
# vline!( plt[1], [RvdW], color=:black, label="", linestyle=:dash )
# vline!( plt[1], [RHb],  color=:black, label="", linestyle=:dash )
# vline!( plt[1], [0.0],  color=:black, label="", linestyle=:dash )

# vline!( plt[2], [Rc],   color=:gray,  label="", linestyle=:dash )
# vline!( plt[2], [RvdW], color=:black, label="", linestyle=:dash )
# vline!( plt[2], [RHb],  color=:black, label="", linestyle=:dash )
# vline!( plt[2], [0.0],  color=:black, label="", linestyle=:dash )

display(plt)


