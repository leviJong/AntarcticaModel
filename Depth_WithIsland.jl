function d(x; s=0.00000012, x_s=75000, σ=12000, d_0=240, λ=500) 
    d = d_0-s*(x)^2 + λ * exp(-(((x-x_s)/σ)^2))
    return(d)
end
    
using CairoMakie


f = Figure()
ax = Axis(f[1,1])
x_values = 0:1:100000
lines!(ax, x_values, d.(x_values))
ylims!(-2000,500)
f