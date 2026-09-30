using DynamicalSystemsBase
using CairoMakie

growth_function = 1 #0 is optimal temperature growth, 1 is logistic

p = Dict(
    :T => 270.0,
    :T_opt => 273.0,
    :λ => 0.2,
    :k => 20.0,     
    :r => 1.0
)

function plant(u,p,t)
    b = u[1]
    p = NamedTuple(p)

    g = 1 - b
    
    if growth_function == 0
        db = growth_opt(b, p.T, p.T_opt, p.k, g, p.λ) 
    elseif growth_function == 1
        db = growth_log(b, g, p.r, p.λ)
    end

    return SVector(db)
end

function growth_opt(b, T, T_opt, k, g) # Growth function for plants
    if k > abs(T - T_opt)
        β_b = 1-(k^-2)*(T-T_opt)^2
    else
        β_b = 0
    end
    
    λ = loss_plant(T, T_opt, k)

    db = b * (g * β_b - λ)
    return(db)
end

function loss_plant(T, T_opt, k)
    if k > abs(T - T_opt)
        λ_b = -1+(k^-2)*(T-T_opt)^2
    else
        λ_b = 1
    end
    return(λ_b)
end

function growth_log(b, g, r, λ)
    db = b*r*(1 - (b/g)) - λ*b
    return(db)
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
