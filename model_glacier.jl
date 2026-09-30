using DynamicalSystemsBase
using CairoMakie

p = Dict(
)

function ice(u,p,t)
    
    return SVector(dI)
end

u0 = [0.63]
t0 = 0.0
ds = CoupledODEs(ice, u0, p)

t_total = 100.0
dt = 1.0
X, t = trajectory(ds, t_total; Δt=dt)


fig = Figure()
ax = Axis(fig[1, 1]; xlabel = "time", ylabel = "ice cover")
for var in columns(X)
    lines!(ax, t, var)
end
fig
