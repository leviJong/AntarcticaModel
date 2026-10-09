using CriticalTransitions # need to install this if you don't have it!
using CairoMakie

include("antarctica.jl")

u0 = [263.0, 0.3, 1.0e3] # initial condition
diffeq = (; alg = Tsit5(), reltol=1e-10) # solver settings (adjust reltol if needed)

sys = CoupledODEs(antarctica, u0, p) # get p and antarctica from antarctica.jl

# Define time depdendence of emissivity
fp = ForcingProfile(:linear)
# Alternative: sinusoidal ForcingProfile
fp2 = ForcingProfile(sin, (0,10*pi))

# Make a system with time-dependent emissivity (called RateSystem)
# More info on how to do this here: 
# https://juliadynamics.github.io/CriticalTransitions.jl/dev/man/system_construction/#Non-autonomous:-RateSystem

T_tot = 1e6 # total simulation time
times = 0:100.0:T_tot # array of time steps (for plotting later)

rsys = RateSystem(sys, fp, :E;
    forcing_start_time=0.0,
    forcing_duration=T_tot/2,
    forcing_scale=0.1,
    reverse=true)

# Run simulation with changing emissivity
tr = trajectory(rsys, T_tot, u0; Δt=100.0)

# Plot results
begin
    fig = Figure(size=(800,600))
    ax1 = Axis(fig[1,1], xlabel="Time (years)", ylabel="Emissivity")
    ax2 = Axis(fig[2,1], xlabel="Time (years)", ylabel="Temperature (K)")

    lines!(ax1, times, parameter.(rsys, times, :E))
    lines!(ax2, times, tr[1][:,1])

    fig
end