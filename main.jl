using DynamicalSystemsBase, CairoMakie, Attractors
using OrdinaryDiffEq: Tsit5

include("model_antarctica.jl")

function main(runAttractors)

    diffeq = (; alg = Tsit5(), abstol=1e-10, reltol=1e-10)

    u0 = [263.0, 0.3, 1.0e3] #Starting value of T, b, L
    t0 = 0.0 #Starting time
    ds = CoupledODEs(antarctica, u0, p; diffeq) #Runs the function over time

    if runAttractors

        grid = (
            range(164.0, 350.0; length=20), # T
            range(0.0, 1.0; length=2), # b
            range(1.0e3, 1.0e5; length = 20) # L
        )
        println(range(1.0e4, 6.0e4; length = 5))

        println(length(grid[1]))
        println(length(grid[2]))
        println(length(grid[3]))

        println(length(grid[1]) * length(grid[2]) * length(grid[3]))

        bmap = BasinMapRecurrences(
            ds, grid;
            consecutive_recurrences = 100, attractor_locate_steps = 1000,
            consecutive_lost_steps = 100,
            sparse = false,
        )

        sampler = RandomICsSampler(1000, grid)
        #basin, attractors = basins_of_attraction(bmap)


        #println(attractors) 

        ascm = AttractorSeedContinueMatch(bmap)

        prange = range(0.5, 0.8; length = 3)
        pidx = :E

        pcurve = [Dict(pidx => p) for p in prange]
        gco = global_continuation(ascm, pcurve, sampler)

        fractions_cont = gco.attractors
        attractors_cont = gco.attractors

        println(attractors_cont)
        attractors_cont[1][1][1]


        Y = [attractors_cont[n][1][1] for n in 1:length(attractors_cont)]

        Y_T = [Y[n][1] for n in 1:length(Y)]
        Y_b = [Y[n][2] for n in 1:length(Y)]
        Y_I = [Y[n][3] for n in 1:length(Y)]


        fig_attractors = Figure(size=(1400,800))

        ax_T_eq = Axis(fig_attractors[1, 1]; xlabel = "Starting value E", ylabel = "temperature") 
        ax_g_eq = Axis(fig_attractors[1,2]; xlabel = "Starting value E", ylabel = "plant cover")
        ax_L_eq = Axis(fig_attractors[2,1]; xlabel = "Starting value E", ylabel = "L")



        lines!(ax_T_eq, prange, Y_T, color = :tomato)
        lines!(ax_g_eq, prange, Y_b, color = :green)
        lines!(ax_L_eq, prange, Y_I, color = :blue)

        lines!(ax_g_eq, prange, ice_cover.(Y_I), color = :blue)
        ylims!(ax_g_eq, -0.1, 1.1)

        fig_attractors

    else

        t_total = 3.3e4
        dt = 1.0 
        #dt is how much you increment time each calculation, and t_total is when it stops
        X, t = trajectory(ds, t_total; Δt=dt)

        X_columns = columns(X) #X is a matrix of results, columns seperates these
        #println(X_columns[1])


        fig = Figure(size=(1400,800))

        ax_temp = Axis(fig[1, 1]; xlabel = "time", ylabel = "temperature") 
        ax_plant = Axis(fig[1,2]; xlabel = "time", ylabel = "plant cover")
        ax_ice = Axis(fig[2,1]; xlabel = "time", ylabel = "L")

        ax_I = Axis(fig[2,2]; xlabel = "time", ylabel = "I")
        lines!(ax_temp, t, X_columns[1], color = :tomato)
        lines!(ax_plant, t, X_columns[2], color = :green)
        lines!(ax_ice, t, X_columns[3], color = :blue)

        lines!(ax_I, t, ice_cover.(X_columns[3]), color = :blue)

        println(X_columns[3][30000])

        fig

    end
end

runAttractors = false # If true will run the attractors code
main(runAttractors)