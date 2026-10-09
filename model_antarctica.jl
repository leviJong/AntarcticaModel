using DynamicalSystemsBase
using CairoMakie
using Attractors

using OrdinaryDiffEq: Tsit5


include("model_plantGrowth.jl") 
include("model_glacier.jl")

function antarctica(u,p,t)
    T,b,L = u

    p = NamedTuple(p) # Whenever p.n is called, n is defined in the Dict
    
    I = ice_cover(L) #Change this to be based on L
    g = max(1 - b - I, 0) # Vacant space

    α_tot = p.α_b * b + p.α_I * I + p.α_g * g #Total albedo is the albedo added up
    ϵ = p.E + b * p.ϕ  # Emmisivity
        
    dT = (p.S / 4 ) * (1 - α_tot) - ϵ * p.σ * T^4  
    db = change_plant(b, T, g, p.T_opt, p.k, p.λ_bopt)
    dL = change_ice(L, T) 
    #These are the differential equations

    return SVector(dT * p.τ_T, db * p.τ_b, dL * p.τ_L)
end

p = Dict(
    :σ => 5.670e-8, #Stefanboltzman constante
    :S => 1368.0,    #Solar impact (energy from sun)
    
    #Constants we can tweak
    :α_b => 0.4,   # Albedo (reflectivity) value of the plants
    :α_I => 0.9,   # Albedo of the ice
    :α_g => 0.2,   # Albedo of the ground
    :λ_bopt => 0.05,     # Loss rate of plants
    :E => 0.75,     # Default emmisivity of the atmosphere, without plants
    :ϕ => 2.8e-2,    # Effectivity of the plants on emmisivity 
    :k => 32.5,      # Survivavable deviation in temperature wherein plants can still reproduce 
    :T_opt => 265.65, # Optimal temperature for plant reproducition

    #Timescales 
    :τ_T => 1.0e-2, # Temperature
    :τ_L => 1.0e-3, # Glacier
    :τ_b => 1.0 # plants
    )
    