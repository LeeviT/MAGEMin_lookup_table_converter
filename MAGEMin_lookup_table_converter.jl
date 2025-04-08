# This script reads in a CSV file containing phase transition lookup table in MAGEMin 
# format and then scales its density to take account porosity near Earth's surface

using DataFrames, CSV, Plots, TidierData, DelimitedFiles

# Filename of a lookup table to process, let us assume it is located in a same folder as 
# this processor script. Also define the output filepath.
inputtablefilepath = joinpath(".", "metabasite_mb.csv")

# Read in a table to initialize a dataframe
df = DataFrame(CSV.File(inputtablefilepath))

filtereddf = @chain df begin
    # Select the relevant columns and rename them, as special characters such as square 
    # brackets are messy to handle in column names in Julia
    TidierData.@select(P = var"P[kbar]", T = var"T[°C]", var"phase", 
    density = var"density[kg/m3]", Cp = var"heatCapacity[J/K]", alpha = var"alpha[1/K]", 
    enthalpy = var"Enthalpy[J]", Vp = var"Vp[km/s]", Vs = var"Vs[km/s]")
    # Filter out data of individual minerals and preserve bulk rock (=system) composition data  
    @filter(phase == "system")
    # Convert pressure from kbar to bar and temperature from celsius to kelvin 
    @mutate(P = 1e3 * P, T = 273.15 + T)
    # Drop the 'phase' column and reorder the columns for PerpleX format
    TidierData.@select(T, P, density, alpha, Cp, Vp, Vs, enthalpy)
    @arrange(P, T)
end

filtereddf.Cp = strip.(filtereddf.Cp, '[')
filtereddf.Cp = strip.(filtereddf.Cp, ']')
filtereddf.Cp = parse.(Float64, filtereddf.Cp)


# print(first(filtereddf, 10))

# Write the metadata and headers in PerpleX format
open("./delim_file.txt", "w") do f
    write(f, "|6.6.6\nmetabasite-mb.tab\n\t\t2\nT(K)\n\t273.15\n\t17.1875\n\t\t129\nP(bar)\n\t10\n\t3281.171875\n\t\t129\n\t\t8\nT(K)\tP(bar)\trho,kg/m3\talpha,1/K\tcp,J/K/kg\tvp,km/s\tvs,km/s\th,J/kg")
end


# phi = porosityscaling(df[:,2])[1]
# rhoporo = zeros(16641)
# for i in 1:129
#     rhoporo[(i - 1)*129 + 1:i*129] = df[:,3][(i - 1)*129 + 1:i*129] .* (1.0 - phi[i]) .+ 2700.0 * phi[i]
# end
# df[:,3] = rhoporo

# Plots.CURRENT_PLOT.nullableplot = nothing
# display(plot!(porosityscaling(df[:,2])[1], porosityscaling(df[:,2])[2], yflip=true))
# ylims!(0, 25)
# xlims!(0, 0.1)
# display(heatmap!(df[:,1][1:129], df[:,2][1:129:16641], transpose(reshape(df[:,5], 129, 129)), color=cgrad(:imola, rev=true)))
# display(heatmap!(df[:,1][1:129] .- 273.15, df[:,2][1:129:16641], transpose(reshape(rhoporo, 129, 129)), color=cgrad(:imola, rev=true)))
# gui()
# CSV.write("./morb_basalt_pmax27GPa_poro_ig.csv", df, delim='\t')

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