using DynamicalSystemsBase
using CairoMakie

growth_function = 1 #0 is optimal temperature growth, 1 is logistic

p = Dict(
    :T => 270.0,
    :T_opt => 273.0,
    :λ_bopt => 0.2,
    :k => 20.0,     
    :r => 1.0
)

function plant(u,p,t)
    b = u[1]
    p = NamedTuple(p)

    g = 1 - b
    db = growth_plant(b, p.T, p.T_opt, p.k) 
    return SVector(db)
end

function growth_opt(b, T, T_opt, k) # Growth function for plants
    if k > abs(T - T_opt)
        β_b = 1-(k^-2)*(T-T_opt)^2
    else
        β_b = 0
    end
    λ = loss_plant(T, T_opt, k, λ_bopt)

    db = b * (β_b - λ)
    return(db)
end

function loss_b(T, T_opt, k, λ_bopt) # loss function for plants
    if k > abs(T - T_opt) #If temperature is outside deviation, then loss=1
        return(λ_bopt+(1-λ_bopt)*(k^-2)*(T-T_opt)^2) #Ask Dylaan how this works
    else
        return(1)
    end
end

u0 = [0.3]
t0 = 0.0
ds = CoupledODEs(plant, u0, p)

t_total = 100.0
dt = 1.0
X, t = trajectory(ds, t_total; Δt=dt)


fig = Figure()
ax = Axis(fig[1, 1]; xlabel = "time", ylabel = "variable")
for var in columns(X)
    lines!(ax, t, var)
end
fig
