# Load the essential modules
using LsqFit
using LinearAlgebra

using Plots

# =========== Functions

struct Grid2D
    g0::Tuple{Float64, Float64}  # (x0, y0)
    dg::Tuple{Float64, Float64}  # (dx, dy)
    n::Tuple{Int, Int}           # (nx, ny)
    data::Array{Float64, 2}      # Grid data (nx x ny)

    # Constructor
    function Grid2D(g0::Tuple{Float64, Float64}, dg::Tuple{Float64, Float64}, n::Tuple{Int, Int})
        data = zeros(Float64, n...)
        new(g0, dg, n, data)
    end
end

function extent(grid::Grid2D)
    xlims = (grid.g0[1], grid.g0[1] + (grid.n[1] - 1) * grid.dg[1])
    ylims = (grid.g0[2], grid.g0[2] + (grid.n[2] - 1) * grid.dg[2])
    return xlims, ylims
end

function extent2(grid::Grid2D)
    xlims = range(grid.g0[1], length=grid.n[1], step=grid.dg[1])
    ylims = range(grid.g0[2], length=grid.n[2], step=grid.dg[2])
    return xlims, ylims
end



# Function to perform bilinear interpolation projection
function project_charge_to_grid( grid::Grid2D, qs::Vector{Float64})::Array{Float64, 2}
    # Create an empty grid with zeros
    x0, y0 = grid.g0
    dx, dy = grid.dg
    nx, ny = grid.n
    idx::Float64 = 1.0 / dx
    idy::Float64 = 1.0 / dy

    # Iterate over each point and its associated charge
    for i::Int in 1:size(ps, 1)
        x::Float64      = ps[i, 1]
        y::Float64      = ps[i, 2]
        charge::Float64 = qs[i]

        # Calculate the grid cell coordinates
        gx::Float64 = (px - x0) * idx
        gy::Float64 = (py - y0) * idy

        x0::Int = floor(Int, gx);  x1::Int = x0 + 1
        y0::Int = floor(Int, gy);  y1::Int = y0 + 1
        
        # Compute the interpolation weights
        wx1::Float64 = x - x0; wx0::Float64 = 1.0 - wx1
        wy1::Float64 = y - y0; wy0::Float64 = 1.0 - wy1

        # Ensure indices are within grid bounds
        x0 = clamp(x0, 1, nx); x1 = clamp(x1, 1, nx)
        y0 = clamp(y0, 1, ny); y1 = clamp(y1, 1, ny)

        # Distribute the charge using bilinear interpolation
        grid.data[ix0, iy0] += charge * wx0 * wy0
        grid.data[ix1, iy0] += charge * wx1 * wy0
        grid.data[ix0, iy1] += charge * wx0 * wy1
        grid.data[ix1, iy1] += charge * wx1 * wy1
    end
    return grid
end

# Function to calculate total potential at target points given source charges in 2D
function evalV_J( ps::Matrix{Float64}, qs::Vector{Float64}, pvs::Matrix{Float64}, bJ::Bool=:fase )::Tuple{Vector{Float64}, Union{Matrix{Float64},Nothing}}
    nq = size(ps, 1)
    nV = size(pvs,1)
    V  = zeros(Float64, nV)

    #println("evalV_J: nq = ", nq, "  nV = ", nV)
    #println("evalV_J ps:  ", ps )
    #println("evalV_J qs:  ", qs )
    #println("evalV_J pvs: ", pvs )
    if bJ
        J  = zeros(Float64, nV, nq)
    else
        J = nothing
    end
    for i in 1:nV
        pv = pvs[i,:]
        Vi::Float64 = 0.0
        for j in 1:nq
            pi        = ps[j,:]
            r::Float64   = norm(pv-pi)
            Jij::Float64 = 1.0/r
            if bJ
                J[i,j] = Jij
            end
            Vi += qs[j]*Jij
        end
        V[i] = Vi
    end
    return V, J
end

function is_symmetric(H::Matrix{Float64}; atol::Float64=1e-8)
    return norm(H - H', Inf) < atol
end

function lm_fit( model::Function, J::Matrix{Float64}, Vref::Vector{Float64}, params0::Vector{Float64}; regularize::Union{Function,Nothing}=nothing, max_iter::Int=1000, damp::Float64=0.001, tol::Float64=1.e-9 )::Vector{Float64}

    # Initial charges
    params = copy(params0)

    JT = copy(J')
    H = JT*J +  damp*I 

    #println("J: "); display(J)
    #println("H: "); display(H)
    #println( "is_symmetric(H) ", is_symmetric(H) )
    #eigvals_H = eigen(H).values;  println("Eigenvalues of H: ", eigvals_H)

    #Hchol = cholesky(H)
    Hfac = lu(H)
    err=0.0
    errF=0.0
    iter=0
    for iter_ in 1:max_iter
        # Evaluate the potential with the current charges
        Vpred = model(params)

        # Compute the residual (difference between predicted and reference potentials)
        residual = Vpred - Vref

        # Compute the normal equations: (J^T * J + λ * I) Δq = -J^T * residual
        JTresidual = J' * residual  
        force = - ( Hfac \ JTresidual )  #;println( "lm_fit iter=", iter, " JTresidual: ", JTresidual," force ", force )

        # Update the charges
        params += force

        if regularize !== nothing
            params = regularize(params)
        end

        # Check convergence (if residual norm is small enough)
        err  = norm(residual) 
        errF = norm(force)
        #println( "lm_fit iter=", iter, " err=", err, " errF=", errF )
        iter+=1
        if errF < tol
            break
        end
    end
    println( "lm_fit() DONE iter=", iter, " err=", err, " errF=", errF )
    return params
end


function dd_fit( model::Function, J::Matrix{Float64}, Vref::Vector{Float64}, params0::Vector{Float64}; regularize::Union{Function,Nothing}=nothing, max_iter::Int=1000, dt::Float64=0.05, tol::Float64=1.e-9 )::Vector{Float64}

    params = copy(params0)
    vel = zeros(Float64, size(params))

    JT = copy(J')

    err=0.0
    errF=0.0
    iter = 0
    for iter_ in 1:max_iter
        # Evaluate the potential with the current charges
        Vpred = model(params)

        residual = Vpred - Vref

        force   =  - ( JT * residual )
        vel    += force * dt
        params += vel   * dt

        cvf = dot(vel,force)
        if cvf < 0.0
            vel[:] .= 0.0
        end

        # if regularize !== nothing
        #     params = regularize(params)
        # end

        # Check convergence (if residual norm is small enough)
        err   = norm(residual) 
        errF  = norm(force)
        vnorm = norm(vel)
        println( "dd_fit iter=", iter, " err=", err, " errF=", errF, " vnorm ", vnorm, " cos(v,f)", cvf/(vnorm*errF) )
        iter+=1
        if errF < tol
            break
        end
    end
    println( "lm_fit() DONE iter=", iter, " err=", err, " errF=", errF )
    return params
end

function fit_charges( ps::Matrix{Float64}, pVs::Matrix{Float64}, Vref::Vector{Float64}, q0s::Vector{Float64}, Qtot::Float64=1.0 )
    function model(qs)
        Vpred,_ = evalV_J(ps,qs,pVs, :false)
        return Vpred
    end
    function regularize(qs)
        return qs*( Qtot/sum(qs) )
    end
    _,J = evalV_J(ps,q0s,pVs, :true )

    #println("V_ref: ", V_ref );
    #println("J: ", size(J) ); display(J)

    #println("fit_charges ps:   ", size(ps) )
    #println("fit_charges pVs:  ", size(pVs) )
    #println("fit_charges Vref: ", size(Vref) )
    #println("fit_charges q0s:  ", size(q0s) )
    #println("fit_charges J:    ", size(J) )
    #qs = lm_fit( model, J, Vref, q0s, regularize=regularize )
    #qs = lm_fit( model, J, Vref, q0s )

    qs = J \ Vref

    #qs ./= sum(qs)

    #qs = dd_fit( model, J, Vref, qs )
    return qs
end


function bilinear_project( d::Vector{Float64} )::Vector{Float64}
    ws =  Vector{Float64}(undef, 4 )
    wx1::Float64 = d[1]; wx0::Float64 = 1.0 - wx1
    wy1::Float64 = d[2]; wy0::Float64 = 1.0 - wy1
    ws[1] = wx0 * wy0
    ws[2] = wx1 * wy0
    ws[3] = wx0 * wy1
    ws[4] = wx1 * wy1
    return ws
end


# function fit_charges( ps::Matrix{Float64}, pVs::Matrix{Float64}, Vref::Vector{Float64}, q0s::Union{Nothing,Vector{Float64}}=nothing, Qtot::Float64=1.0 )
#     model = (qs, ps) -> evalV(ps, qs, pVs)
#     if q0s === nothing
#         q0s = ones(Float64, size(ps, 1)).*(Qtot/size(ps,1))
#     end
#     println("q0s:  ", size(q0s ) )
#     println("ps:   ", size(ps  ) )
#     println("pVs:  ", size(pVs ) )
#     println("Vref: ", size(Vref) )
#     result = curve_fit( model, q0s, pVs, Vref )
#     qs_out = result.param
#     return qs_out
# end

function generate_circle_points(center::Vector{Float64}, Rs::Vector{Float64}, nPhi::Int)::Matrix{Float64}
    nr = length(Rs)
    ntot = nr * nPhi
    cx = center[1] 
    cy = center[2]
    ps = zeros(Float64, ntot, 2)
    dphi = 2 * π / nPhi
    i::Int = 1
    for ip in 1:nPhi
        phi = dphi * (ip - 1) 
        hx = cos(phi)       
        hy = sin(phi)
        for ir in 1:nr
            r = Rs[ir]       
            ps[i,1] = cx + hx*r
            ps[i,2] = cy + hy*r
            i += 1
        end
    end
    return ps
end

# Function to calculate the electrostatic potential on the grid
function calculate_electrostatic_potential(grid::Grid2D, ps::Matrix{Float64}, qs::Vector{Float64}, mask::Array{Bool,2})::Array{Float64, 2}
    # Unpack grid properties
    x0, y0 = grid.g0
    dx, dy = grid.dg
    nx, ny = grid.n

    # Create a potential grid with zeros
    V::Array{Float64, 2} = zeros(Float64, nx, ny)

    # for k::Int in 1:size(ps, 1)
    #     println( "p[",k,"] ", ps[k,:] )
    # end

    # Iterate over each grid point
    for i::Int in 1:nx
        for j::Int in 1:ny
            # Check if the current grid point is within the masked area
            if mask[i, j]
                # Calculate the physical coordinates of the grid point
                gx::Float64 = x0 + (i - 1) * dx
                gy::Float64 = y0 + (j - 1) * dy

                # Initialize potential at this grid point
                Vi::Float64 = 0.0

                # Sum the contributions from all point charges
                for k::Int in 1:size(ps, 1)
                    
                    px::Float64 = ps[k, 1]
                    py::Float64 = ps[k, 2]
                    q::Float64 = qs[k]

                    # Calculate the distance from the point charge to the grid point
                    r::Float64 = sqrt((gx - px)^2 + (gy - py)^2)

                    # Add the potential contribution if r is non-zero
                    if r > 1.e-32
                        Vi += q/r
                    end
                end
                # Assign the calculated potential to the grid point
                V[i,j] = Vi
            end
        end
    end
    return V
end



# Function to create a mask based on minimum and maximum radius from a certain point
function create_circular_mask(grid::Grid2D, center::Vector{Float64}, r_min::Float64, r_max::Float64)::Array{Bool, 2}
    # Unpack grid properties
    x0, y0 = grid.g0
    dx, dy = grid.dg
    nx, ny = grid.n

    # Initialize the mask with false values
    mask::Array{Bool, 2} = falses(nx, ny)

    # Unpack the center coordinates
    cx, cy = center
    #println("center: ", cx, ", ", cy)
    #println("center: ", x0, ", ", y0)

    # Iterate over each grid point
    for i::Int in 1:nx
        for j::Int in 1:ny
            # Calculate the physical coordinates of the grid point
            gx::Float64 = x0 + (i - 1) * dx
            gy::Float64 = y0 + (j - 1) * dy

            # Calculate the distance from the grid point to the center
            r::Float64 = sqrt((gx-cx)^2 + (gy-cy)^2)

            # Set mask value to true if the distance is within the specified range
            if r_min <= r <= r_max
                mask[i, j] = true
            end
        end
    end
    return mask
end



function get_optimal_projections( pqs, pgs, Rs, nphi )
    Rs   = [1.0, 3.0]
    #pVs  = generate_circle_points( [0.05,0.05], Rs, 32 )
    q0s  = zeros(Float64, size(pVs,1))
    prjs = Matrix{Float64}(undef, size(pqs,1), size(pgs,1) )
    npq = size(pqs,1)
    for ip in 1:npq
        pq = pqs[ip,:]
        #qs_fit = fit_charges( pgs, pVs, V_ref, qgs0, 1.0 )
        pVs   = generate_circle_points( pq, Rs, 32 )
        pq_   = copy(hcat([pq,]...)')
        V_ref,_ = evalV_J( pq_, [1.0], pVs, :false )
        _,J = evalV_J(pgs,q0s,pVs,:true )
        #println("V_ref: ", V_ref );
        #println("J: ", size(J) ); display(J)
        qs_fit = J \ V_ref

        #qs_fit ./= sum(qs_fit)

        prjs[ip,:] = qs_fit
        err = maximum( abs.(V_ref - J*qs_fit) )
        #println( "pq: ", pq," qsum=", sum(qs_fit) ," qs: ", qs_fit )
        println( "pq: ", pq," qsum=", sum(qs_fit) ," err: ", err )
        #return
    end
    return prjs
end

function lerp_points( p1::Vector{Float64}, p2::Vector{Float64}, n::Int, bEnd::Bool=:true )::Matrix{Float64}
    # Initialize the matrix to hold the interpolated points
    m = length(p1)
    ps = Matrix{Float64}(undef, n, m)
    dt = 1.0/n 
    if bEnd
        dt = 1.0/(n-1)
    end
    for i in 1:n
        t = dt*(i-1)
        ps[i,:] = p1 .* (1-t) + p2 .* t
    end
    return ps
end



# =========== Body

grid = Grid2D(  (0.0,0.0), (0.1,0.1), (100,100) )

pg0 = [5.00,5.00]; 
#pq = [5.05,5.05];   
#pq = [5.01,5.05]; 
pq = [5.04,5.06]; 
qq = [1.0,];   
pq_=copy(hcat([pq,]...)')
#ps = [ pq,   [5.0,5.0],[5.0,6.0],[6.0,5.0],[6.0,6.0] ]; ps=copy(hcat(ps...)')
#Qs = [ 1.0,  -0.25, -0.25, -0.25, -0.25 ]
#pgs = [ [0.0,0.0],[0.0,0.1],[0.1,0.0],[0.1,0.1] ]; 

# pgs = [ [0.0,0.0],[0.0,0.1],[0.1,0.0],[0.1,0.1],        ]; 
# qgs = [ 0.25, 0.25, 0.25, 0.25,     ]

pgs = [ [0.0,0.0],[0.0,0.1],[0.1,0.0],[0.1,0.1],    [0.2,0.0],[0.2,0.1], [-0.1,0.0],[-0.1,0.1],    [0.0,0.2],[0.1,0.2], [0.0,-0.1],[0.1,-0.1],       ]; 
qgs = [ 0.25, 0.25, 0.25, 0.25,    0.0, 0.0, 0.0, 0.0,    0.0, 0.0, 0.0, 0.0 ]

pgs = [ 
    [-0.1,-0.1],[-0.1,-0.0],[-0.1,0.1],[-0.1,0.2],
    [ 0.0,-0.1],[ 0.0,-0.0],[ 0.0,0.1],[ 0.0,0.2],
    [ 0.1,-0.1],[ 0.1,-0.0],[ 0.1,0.1],[ 0.1,0.2],
    [ 0.2,-0.1],[ 0.2,-0.0],[ 0.2,0.1],[ 0.2,0.2],
]; 
qgs = [ 0.25, 0.25, 0.25, 0.25,    0.0, 0.0, 0.0, 0.0,    0.0, 0.0, 0.0, 0.0,   0.0, 0.0, 0.0, 0.0  ]

pgs=copy(hcat(pgs...)')




#pgs = [ [5.0,5.0],[5.0,6.0],[6.0,5.0],[6.0,6.0],  [4.0,5.0],[4.0,6.0],    [5.0,4.0],[6.0,4.0],     [7.0,5.0],[4.0,6.0],  [5.0,4.0],[6.0,4.0],     ];  pgs=copy(hcat(pgs...)')
#qgs = [ 0.25, 0.25, 0.25, 0.25 ]

#Rs = [1.0,1.25, 2.0, 3.0]
Rs = [1.5,2.0, 4.0]
p1=[0.0,0.05]; p2=[0.05,0.05];
#pqs = lerp_points( [0.0,0.0], [0.05,0.05], 10 )
#pqs  = lerp_points( [0.0,0.0], [0.00,0.05], 10 )
pqs  = lerp_points( p1, p2, 10 )
prjs = get_optimal_projections( pqs, pgs, Rs, 16 )


#println("prjs "); display(prjs)
plt = plot(title="Charge Projection o Grid Points", xlabel=(string(p1)*"-"*string(p2)), ylabel="Q[i]", legend=:topleft)
npg = size(pgs,1)
for ip in 1:npg 
    #println("ip: ", ip)
    pg = pgs[ip,:]
    label = string(ip) * " " * string(pg)
    plot!( plt, prjs[:,ip], label=label )
end
display(plt)




println("====================")

pgs[:,1].+=pg0[1]
pgs[:,2].+=pg0[2]

# qgs0  = copy(qgs)
# qblin = bilinear_project( pq-pg0 );  
# println("qgs0 ", qgs0)
# println("qblin ", qblin)
# qgs0[1:4]=qblin[:]; 
# println("qgs0 ", qgs0)




ps = vcat(pq_, pgs )
qs = vcat(qq,  qgs )
qs[1] *= -1.0 

#ps = [ pq, ]   ; ps=hcat(ps...)
#Qs = [ 1.0, ]
#println("typeof(ps) ", typeof(ps))

#mask_::Array{Bool, 2} = trues(100,100)
#V = calculate_electrostatic_potential(grid, ps, Qs, mask_)


pVs = generate_circle_points( pq, Rs, 16 )

V_ref,_ = evalV_J( pq_, qq, pVs, :false )
qs_fit = fit_charges( pgs, pVs, V_ref, qgs0, 1.0 )

#qs_fit = qs_fit/sum(qs_fit)

println("qs_fit ", sum(qs_fit), qs_fit)
#println( "size(qs_fit) ", size(qs_fit))
#println( "size(qs)     ", size(qs) )
qsf = copy(qs); qsf[2:end]= qs_fit[:]; #println("qs_fit ", qsf)


mask = create_circular_mask( grid, pq, Rs[1], Rs[end] )
V     = calculate_electrostatic_potential(grid, ps, qs,  mask)
V_fit = calculate_electrostatic_potential(grid, ps, qsf, mask)

#plt_mask = heatmap(mask, color=:blues,   aspect_ratio=1, title="Mask Visualization",      xlabel="X Index", ylabel="Y Index"); display(plt_mask);
#Vmax = max( -minimum(V), maximum(V))
Vmax = max( -minimum(V_fit), maximum(V_fit))

x_sz, y_sz = extent2(grid)
#plt_V      = heatmap( x_sz, y_sz, V,    color=:bwr, aspect_ratio=1, title="Electrostatic Potential", xlabel="X Index", ylabel="Y Index", clims=(-Vmax,Vmax), size=(1000,1000) ); 
plt_V      = heatmap( x_sz, y_sz, V_fit,    color=:bwr, aspect_ratio=1, title="Electrostatic Potential", xlabel="X Index", ylabel="Y Index", clims=(-Vmax,Vmax), size=(1000,1000) ); 

scatter!(plt_V, pVs[:, 1], pVs[:, 2], marker=:circle, color=:black, label="V samples Points", markersize=3)

scatter!(plt_V, pgs[:, 1], pgs[:, 2], marker=:cross, color=:red, label="grid charges", markersize=3)

display(plt_V);


