function d(x; s=0.00000012, x_s=75000, σ=10000, d_0=200, λ=250) 
    d = d_0-s*(x*40827)^2 + λ * exp(-(((x*40827-x_s)/σ)^2))
    return d
end

using CairoMakie

f = Figure()
ax = Axis(f[1,1]; xlabel = "x", ylabel = "depth")
x_values = 0:0.001:2.5
lines!(ax, x_values, d.(x_values))
ylims!(-1500,500)
f
