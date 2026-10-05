using DynamicalSystemsBase
using CairoMakie
using OrdinaryDiffEq: Vern9


include("model_plantGrowth.jl") 
include("model_glacier.jl")

function antarctica(u,p,t)
    T,b,L,E = u

    p = NamedTuple(p) # Whenever p.n is called, n is defined in the Dict
    
    I = 0.2 #Change this to be based on L
    g = 1 - b - I # Vacant space

    α_tot = p.α_b * b * g + p.α_I * I + p.α_g * g #Total albedo is the albedo added up
    dE = p.a - p.ϕ * E * b - p.ϕ_I * E * I # Emmisivity  
    dT = (p.S / 4 ) * (1 - α_tot) - (1 - E) * p.σ * T^4  
    db = change_plant(b, T, g, p.T_opt, p.k, p.λ_bopt)
    dL = change_ice(L, T) 
    #These are the differential equations

    return SVector(dT * p.τ_T, db * p.τ_b, dL * p.τ_L, dE)
end

p = Dict(
    :σ => (5.670)*10^(-8), #Stefanboltzman constante
    :S => 1368.0,    #Solar impact (energy from sun)
    
    #Constants we can tweak
    :α_b => 0.4,   # Albedo (reflectivity) value of the plants
    :α_I => 0.9,   # Albedo of the ice
    :α_g => 0.1,   # Albedo of the ground
    :λ_bopt => 0.05,     # Loss rate of plants 
    :a => 1.0e-4,      #production of emmisivity
    :ϕ_I => 4.0e-5,  # Effectivity of the plants on emmisivity 
    :ϕ => 5.0e-5,  # Effectivity of the plants on emmisivity 
    :k => 32.5,      # Survivavable deviation in temperature wherein plants can still reproduce 
    :T_opt => 265.65, # Optimal temperature for plant reproducition

    #Timescales 
    :τ_T => 1.0e-3, # Temperature
    :τ_L => 1.0, # Glacier
    :τ_b => 1.0 # plants
    )
 
   
print(p)

diffeq = (; alg = Vern9(), dt=1e-4)

u0 = [263.0, 0.3, 1.0e4, 0.05] #Starting value of T, b, L, E
t0 = 0.0 #Starting time
ds = CoupledODEs(antarctica, u0, p; diffeq) #Runs the function over time

t_total = 2000.0
dt = 1.0 
#dt is how much you increment time each calculation, and t_total is when it stops
X, t = trajectory(ds, t_total; Δt=dt)

println("test")
println(X)

println(X[2])

X_columns = columns(X) #X is a matrix of results, columns seperates these
#println(X_columns[1])


fig = Figure(size=(1400,400))

ax_temp = Axis(fig[1, 1]; xlabel = "time", ylabel = "temperature") 
ax_plant = Axis(fig[1,2]; xlabel = "time", ylabel = "plant cover")
ax_ice = Axis(fig[2,1]; xlabel = "time", ylabel = "L")
ax_emmisivity = Axis(fig[2,2]; xlabel = "time", ylabel = "E")

lines!(ax_temp, t, X_columns[1], color = :tomato)
lines!(ax_plant, t, X_columns[2], color = :green)
lines!(ax_ice, t, X_columns[3], color = :blue)
lines!(ax_emmisivity, t, X_columns[4], color = :orange)

fig
