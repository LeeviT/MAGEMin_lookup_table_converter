# This script reads in a CSV file containing phase transition lookup table in MAGEMin 
# format and then scales its density to take account porosity near Earth's surface

using DataFrames, CSV, CairoMakie, TidierData, DelimitedFiles

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
    @mutate(alphafilt = case_when(alpha < 0.0 => 0.0, alpha > 0.0 => alpha))
    # Convert pressure from kbar to bar and temperature from celsius to kelvin 
    @mutate(P = 1e3 * P, T = 273.15 + T)
    # Drop the 'phase' column and reorder the columns for PerpleX format
    TidierData.@select(T, P, density, alphafilt, Cp, Vp, Vs, enthalpy)
    @arrange(P, T)
end

# print(filtereddf.alphafilt)

# Filter out the square brackets from the heat capacity column, and convert it from string
# to 64-bit float
filtereddf.Cp = strip.(filtereddf.Cp, '[')
filtereddf.Cp = strip.(filtereddf.Cp, ']')
filtereddf.Cp = parse.(Float64, filtereddf.Cp)

filtereddf = @chain filtereddf begin
    @mutate(Cpfilt = case_when(Cp < 750.0 => 750.0, Cp > 6000.0 => 6000.0, Cp > 0.0 => Cp))
end

# Write the metadata and headers in PerpleX format
open("./delim_file.txt", "w") do io
    write(io, "|6.6.6\nmetabasite-mb.tab\n\t\t2\nT(K)\n\t273.15\n\t17.1875\n\t\t129\nP(bar)\n\t10\n\t3281.171875\n\t\t129\n\t\t8\nT(K)\tP(bar)\trho,kg/m3\talpha,1/K\tcp,J/K/kg\tvp,km/s\tvs,km/s\th,J/kg\n")
    writedlm(io, eachrow(filtereddf), '\t')
end

# print(round(minimum(filtereddf.density), sigdigits = 3))

# phi = porosityscaling(df[:,2])[1]
# rhoporo = zeros(16641)
# for i in 1:129
#     rhoporo[(i - 1)*129 + 1:i*129] = df[:,3][(i - 1)*129 + 1:i*129] .* (1.0 - phi[i]) .+ 2700.0 * phi[i]
# end
# df[:,3] = rhoporo

# Plots.CURRENT_PLOT.nullableplot = nothing
# display(plot!(porosityscaling(df[:,2])[1], porosityscaling(df[:,2])[2], yflip=true))
plottingT = filtereddf.T[1:129] .- 273.15
plottingP = filtereddf.P[1:129:16641] ./ 1e4
with_theme(theme_latexfonts()) do
    f = Figure(backgroundcolor = RGBf(0.9, 0.79, 0.78), size = (1000, 400))
    ga = f[1, 1] = GridLayout()
    axleft = Axis(ga[1, 1], ylabel = "pressure (GPa)", tellwidth = false, width = 300, ylabelsize = 12, yticksize = 3, xticksize = 3, xticklabelsize = 10, yticklabelsize = 10, xticklabelpad = 0.5, yticklabelpad = 0.5)
    axcenter = Axis(ga[1, 2], xlabel = L"temperature ($\degree$ C)", tellwidth = false, width = 300, xlabelsize = 12, xticksize = 3, xticklabelsize = 10, xticklabelpad = 0.5)
    axright = Axis(ga[1, 3], tellwidth = false, width = 300, xticksize = 3, xticklabelsize = 10, xticklabelpad = 0.5)
    linkyaxes!(axleft, axcenter, axright)
    linkxaxes!(axleft, axcenter, axright)

    # Plot the heatmaps for each three properties
    hmdensity = CairoMakie.heatmap!(axleft, plottingT, plottingP, reshape(filtereddf.density, 129, 129), colormap = Reverse(:imola))
    hmCp = CairoMakie.heatmap!(axcenter, plottingT, plottingP, reshape(filtereddf.Cpfilt, 129, 129), colormap = Reverse(:glasgow))
    hmalpha = CairoMakie.heatmap!(axright, plottingT, plottingP, reshape(filtereddf.alphafilt, 129, 129), colormap = Reverse(:nuuk)) 
    # Create the colorbars for heatmaps
    cbdensity = Colorbar(ga[0, 1][1, 1], hmdensity, label = L"density (kg/m$^3$)", vertical = false, 
        tellwidth = false, width = 270, spinewidth = 0.1, labelpadding = 1.5, labelsize = 12, 
        ticklabelpad = 0.5, ticklabelsize = 10, ticksize = 3.5, 
        ticks = [round(minimum(filtereddf.density), sigdigits = 4), 2500, 3000, 3500, 4000, round(maximum(filtereddf.density), sigdigits = 4)],
        minorticks = [2250, 2750, 3250, 3750, 4250], minorticksvisible = true, minorticksize = 2.5, minortickwidth = 0.75)
    cbCp = Colorbar(ga[0, 2][1, 1], hmCp, label = L"heat capacity (J/$\degree$C)", vertical = false, 
        tellwidth = false, width = 270, spinewidth = 0.1, labelpadding = 1.5, labelsize = 12, 
        ticklabelpad = 0.5, ticklabelsize = 10, ticksize = 3.5, ticks = [750, 2000, 4000, 6000], 
        minorticks = [1000, 3000, 5000], minorticksvisible = true, minorticksize = 2.5, minortickwidth = 0.75)
    cbalpha = Colorbar(ga[0, 3][1, 1], hmalpha, label = L"thermal expansivity (1/$\degree$C)", vertical = false, 
        tellwidth = false, width = 270, spinewidth = 0.1, labelpadding = 1.5, labelsize = 12, 
        ticklabelpad = 0.5, ticklabelsize = 10, ticksize = 3.5,
        ticks = [round(minimum(filtereddf.alphafilt), sigdigits = 2), 1e-4, 2e-4, round(maximum(filtereddf.alphafilt), sigdigits = 2)],
        minorticks = [0.5e-4, 1.5e-4, 2.5e-4], minorticksvisible = true, minorticksize = 2.5, minortickwidth = 0.75)
    # Bring the grid up to make it visible
    CairoMakie.translate!(hmdensity, 0, 0, -100)
    CairoMakie.translate!(hmCp, 0, 0, -100)
    CairoMakie.translate!(hmalpha, 0, 0, -100)


    # Do some beautification for the plots 
    CairoMakie.ylims!(axright, low = 0)
    CairoMakie.xlims!(axright, low = 0)
    hidespines!(axleft)
    hidespines!(axcenter)
    hidespines!(axright)
    hideydecorations!(axright, grid = false)
    hideydecorations!(axcenter, grid = false)
    colgap!(ga, 10)
    rowgap!(ga, 10)
    display(f)
end
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