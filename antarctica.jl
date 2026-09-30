using DynamicalSystemsBase
using CairoMakie

include("model_plantGrowth.jl") 
#include("model_glacier.jl")

function antarctica(u,p,t)
    T,b,I = u

    g = 1 - b - I # Vacant space
    p = NamedTuple(p) # Whenever p.n is called, n is defined in the Dict
    
    α_tot = p.α_b * b + p.α_I * I + p.α_g * g #Total albedo is the albedo added up
    ϵ = p.E + b * p.ϕ  # Emmisivity

    β_b = growth_b(T, p.T_opt, p.k) #Defining our growth
        
    dT = (p.S / 4 ) * (1 - α_tot) - ϵ * p.σ * T^4 
   # db = growth_opt(b, T, p.T_opt, p.k, g, p.λ_b)
    db = b * (g * β_b - p.λ_b)
#    dI = differential_ice(I, T, g, p.DFF, p.p_snow)
    dI = 0
    #These are the differential equations

    return SVector(dT,db,dI)
end

function growth_b(T, T_opt, k) # Growth function for plants
    if k > abs(T - T_opt) #If temperature is outside deviation, then growth=0
        return(1-(k^-2)*(T-T_opt)^2) #Ask Dylaan how this works
    else
        return(0)
    end
end

p = Dict(
    :σ => (5.670)*10^(-8), #Stefanboltzman constante
    :S => 1368.0,    #Solar impact (energy from sun)
    
    #Constants we can tweak
    :α_b => 0.4,   # Albedo (reflectivity) value of the plants
    :α_I => 0.9,   # Albedo of the ice
    :α_g => 0.1,   # Albedo of the ground
    :λ_b => 0.05,     # Loss rate of plants
    :E => 0.75,     # Default emmisivity of the atmosphere, without plants
    :ϕ => 0.028,    # Effectivity of the plants on emmisivity 
    :k => 32.5,      # Survivavable deviation in temperature wherein plants can still reproduce 
    :T_opt => 265.65, # Optimal temperature for plant reproducition
    :C => 1.0,        # Time scale for the dT function 
    :DFF => 0.1 * 365 #Ice growth surface area 
    )
 
   
print(p)

u0 = [263.0, 0.3, 0.0] #Starting value of T, b, I
t0 = 0.0 #Starting time
ds = CoupledODEs(antarctica, u0, p) #Runs the function over time

t_total = 100.0
dt = 1.0 
#dt is how much you increment time each calculation, and t_total is when it stops
X, t = trajectory(ds, t_total; Δt=dt)

println("test")
println(X)

println(X[2])

X_columns = columns(X) #X is a matrix of results, columns seperates these
#println(X_columns[1])


fig = Figure()

ax_temp = Axis(fig[1, 1]; xlabel = "time", ylabel = "temperature") 
ax_plant = Axis(fig[1,2]; xlabel = "time", ylabel = "plant cover")

lines!(ax_temp, t, X_columns[1])
lines!(ax_plant, t, X_columns[2], color = :green)

fig


