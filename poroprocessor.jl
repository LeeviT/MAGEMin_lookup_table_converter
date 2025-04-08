# This script reads in a text file containing phase transition look-up table in PerpleX format
# and then scales its density to take account porosity near Earth's surface

using DataFrames
using Pkg
using CSV
using Plots

df = DataFrame(CSV.File("./morb_basalt_pmax27Gpa_ig.csv"))

function porosityscaling(pressure)
    rho0 = 2900.0
    g = 9.81
    phi0 = 0.678
    m = 0.008
    n = 89.53
    # First, convert lithostatic pressure p=\rho*g*h to h=depth in km
    depth = (1.0e2 * pressure[1:129:16641]) / (rho0 * g)
    # 
    phi = phi0 ./ (1 .+ m.*depth).^n

    return phi, depth
end

phi = porosityscaling(df[:,2])[1]
rhoporo = zeros(16641)

for i in 1:129
    rhoporo[(i - 1)*129 + 1:i*129] = df[:,3][(i - 1)*129 + 1:i*129] .* (1.0 - phi[i]) .+ 2700.0 * phi[i]
end

df[:,3] = rhoporo

Plots.CURRENT_PLOT.nullableplot = nothing
# display(plot!(porosityscaling(df[:,2])[1], porosityscaling(df[:,2])[2], yflip=true))
# ylims!(0, 25)
# xlims!(0, 0.1)
# display(heatmap!(df[:,1][1:129], df[:,2][1:129:16641], transpose(reshape(df[:,5], 129, 129)), color=cgrad(:imola, rev=true)))
# display(heatmap!(df[:,1][1:129] .- 273.15, df[:,2][1:129:16641], transpose(reshape(rhoporo, 129, 129)), color=cgrad(:imola, rev=true)))
# gui()

CSV.write("./morb_basalt_pmax27GPa_poro_ig.csv", df, delim='\t')