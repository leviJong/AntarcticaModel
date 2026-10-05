using DynamicalSystemsBase
using CairoMakie
using OrdinaryDiffEq: Vern9

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

function ground_depth(x, d0; x_s = 40000.0, s = 0.014, λ = 300.0, σ = 10000.0)
    d = d0 - s*x + λ * ℯ^-((x-x_s)/σ)^2
    #d = d0

    return(d)
end

function change_ice(L, T)
   if L == 0.0
        dL = 0.0
   else   
        p = NamedTuple(p_i)

        d = ground_depth(L, p.d0) 

                H_f = max(p.α_f * L^(1/2), -p.ϵ * p.δ * d)
                H_m = p.α_m * L^(1/2)

                h_m = (p.d0 + d + H_m + H_f) / 2

            E = p.E0 #- p.k * (T - p.T0)


        F = min(0, p.c * d * H_f) 

        B = p.β * (h_m - E) * L

        dL = L^(-1/2) * 2 * (B + F) / (3 * p.α_m)
   end

   return(dL) 
end

diffeq = (; alg = Vern9(), dt=1e-2)

u0 = [1.0e4]
t0 = 0.0
p = 2003.0
ds = CoupledODEs(ice, u0, p; diffeq)

t_total = 10000.0
dt = 1.0
X, t = trajectory(ds, t_total; Δt=dt)

X_columns = columns(X)

fig = Figure()
ax = Axis(fig[1, 1]; xlabel = "time", ylabel = "L")
ax2 = Axis(fig[1, 2]; xlabel = "x", ylabel = "d")

lines!(ax, t, X_columns[1])

x_values = 0:1:5e4
lines!(ax2, x_values, ground_depth.(x_values, 200.0))

fig
