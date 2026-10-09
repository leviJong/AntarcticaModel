using DynamicalSystemsBase
using CairoMakie
using OrdinaryDiffEq: Tsit5


p_i = Dict(

    :d0 => 200.0, #Starting height by L = 0

    :α_f => 0.5,
    :α_m => 2.0,

    :ϵ => 1.0,
    :δ => 1.127,
    :c => 2.4,
    :β => 0.005,

    :E0 => 100.0,
    :T0 => 265.0,
    :k => 2.0
)

function ice(u,p,t)
    T = p
    L = u[1]
    
    dL = change_ice(L, T)
    return SVector(dL)
end

function ice_cover(L)
    if L > x_c 
        I = 1
    else
        I = L / x_c
    end
    
    return(I)
end

function ground_depth(x, d0; x_s = 7.5e4, s = 1.2e-7, λ = 2.5e2, σ = 1.0e4)
    d = d0 - s * (x)^2 + λ * exp(-( ((x - x_s) / σ) ^2 ) )
    #d = -1.0
    return(d)
end

x_values = 0:1:8e4
x_c = findfirst(ground_depth.(x_values, 240) .< 0)

"""
function ground_depth_lock(x, d0; x_s = 40000.0, s = 0.014, λ = 300.0, σ = 10000.0)
    d = d0 - s*x + λ * ℯ^-((x-x_s)/σ)^2
    return(d)
end"""

function change_ice(L, T)
   if L <= 0.0
        dL = 0.0
   else   
        p = NamedTuple(p_i)

        d = ground_depth(L, p.d0) 

                H_f = max(p.α_f * L^(1/2), -p.ϵ * p.δ * d)
                H_m = p.α_m * L^(1/2)

                h_m = (p.d0 + d + H_m + H_f) / 2

            E = p.E0 + p.k * (T - p.T0)


        F = min(0, p.c * d * H_f) 

        B = p.β * (h_m - E) * L

        dL = L^(-1/2) * 2 * (B + F) / (3 * p.α_m)
   end

   return(dL) 
end
""""

diffeq = (; alg = Tsit5(), abstol=1e-12, reltol=1e-12)

u0 = [1.0]
t0 = 0.0
p = 2003.0
ds = CoupledODEs(ice, u0, p; diffeq)

grid = (
    range(0, 1.0e5; step=5000), # Starting values of L
)
println(range(1, 5; step =1))

bmap = BasinMapRecurrences(
    ds, grid;
    consecutive_recurrences = 100, attractor_locate_steps = 100,
    consecutive_lost_steps = 100, horizon_limit = 1e7,
    sparse = false,
)

basin, attractors = basins_of_attraction(bmap)
println(attractors)
println(attractors[1][1])
println(attractors[2][1])

t_total = 1.0e5
dt = 1.0
X, t = trajectory(ds, t_total, u0; Δt=dt)

X_columns = columns(X)

fig = Figure()
ax = Axis(fig[1, 1]; xlabel = "time", ylabel = "L")
#ax2 = Axis(fig[1, 2]; xlabel = "x", ylabel = "d")
#ax3 = Axis(fig[2,1]; xlabel = "time", ylabel = "ice cover")

lines!(ax, t, X_columns[1])
#lines!(ax2, x_values, ground_depth.(x_values, 200.0))
#lines!(ax3, t, ice_values)
#println(length(ice_values))
#println(length(X_columns[1]))

fig"""