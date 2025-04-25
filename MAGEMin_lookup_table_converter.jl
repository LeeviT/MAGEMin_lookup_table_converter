# This script reads in a CSV file containing phase transition lookup table in MAGEMin format and 
# then scales its density to take account porosity near Earth's surface.

# DataFrames for better handling table format, CSV for readind in a CSV file, CairoMakie for 
# plotting, TidierData for filtering/manipulating data, DelimitedFiles for output.
using DataFrames, CSV, CairoMakie, TidierData, DelimitedFiles

function porosityscaling(pressure, resolution)
    rho0 = 2900.0
    g = 9.81
    phi0 = 0.678
    m = 0.008
    n = 89.53
    # First, convert lithostatic pressure p=\rho*g*h to h=depth in km
    depth = (100.0 .* pressure[1:resolution:resolution^2]) / (rho0 * g)
    phi = phi0 ./ (1 .+ m .* depth) .^ n
    return phi
end

# Filename of a lookup table to process, let us assume it is located in a same folder as this 
# processor script. Also define the output filepath.
inputtablefilepath = joinpath(".", "morb_basalt_high_reso_ig.csv")

# Read in a table to initialize a dataframe.
df = DataFrame(CSV.File(inputtablefilepath))
resolution = 128 * 2 + 1

# A macro to filter and modify a dataframe in MAGEMin format
filtereddf = @chain df begin
    # Select the relevant columns and rename them, as special characters such as square brackets are
    # messy to handle in column names in Julia.
    TidierData.@select(P = var"P[kbar]", T = var"T[°C]", var"phase", density = var"density[kg/m3]", 
        Cp = var"heatCapacity[J/K]", alpha = var"alpha[1/K]", enthalpy = var"Enthalpy[J]", 
        Vp = var"Vp[km/s]", Vs = var"Vs[km/s]")
    # Filter out data of individual minerals and preserve bulk rock (=system) composition data.
    @filter(phase == "system")
    # Convert subzero (=unphysical) heat expansivity values to zero.
    @mutate(alphafilt = case_when(alpha < 0.0 => 0.0, alpha > 0.75e-4 => 0.75e-4, alpha > 0.0 => alpha))
    # Convert pressure from kbar to bar and temperature from celsius to kelvin. 
    @mutate(P = 1e3 * P, T = 273.15 + T)
    # Drop the 'phase' column and reorder the columns for PerpleX format.
    TidierData.@select(T, P, density, alphafilt, Cp, Vp, Vs, enthalpy)
    @arrange(P, T)
end

# Remove the square brackets from the heat capacity column, and convert it from string to 64-bit
# float.
filtereddf.Cp = strip.(filtereddf.Cp, '[')
filtereddf.Cp = strip.(filtereddf.Cp, ']')
filtereddf.Cp = parse.(Float64, filtereddf.Cp)

filtereddf = @chain filtereddf begin
    @mutate(Cp = case_when(Cp < 750.0 => 750.0, Cp > 1500.0 => 1500.0, Cp > 0.0 => Cp))
end

# Find the last temperature column by column index, which has temperature < 500 °C.
Tlimitcoli = findlast(x -> x < 600 + 273.15, filtereddf.T[1:resolution])
compdensity = reshape(filtereddf.density[1:resolution] 
    * exp.(1e-11 * 1e5 .* filtereddf.P[1:resolution:resolution^2])', resolution, resolution)
for i in eachindex(filtereddf.T)
    Pscaling = filtereddf.P[1:resolution:resolution^2]
    (i >= 1 || 1 == mod(i, resolution)) && Tlimitcoli > mod(i - 1, resolution) + 1 ? 
    filtereddf.density[i] = compdensity[i] : 0
end

phioceanic = porosityscaling(filtereddf.P, resolution)
filtereddf.density = reshape((reshape(filtereddf.density, resolution, resolution)' .* (1.0 .- phioceanic) .+ 2600.0 .* phioceanic)', resolution^2)

# open("./morb_gabbro-ig.txt", "w") do io
    # Write the metadata and headers in PerpleX format.
#    write(io, "|6.6.6\nmorb_gabbro-ig.tab\n\t\t2\nT(K)\n\t273.15\n\t17.1875\n\t\t129\nP(bar)\n\t10\n"
#        * "\t3281.171875\n\t\t129\n\t\t8\nT(K)\tP(bar)\trho,kg/m3\talpha,1/K\tcp,J/K/kg\tvp,km/s"
#        * "\tvs,km/s\th,J/kg\n")
    # Write the actual phase diagram data.
#    writedlm(io, eachrow(filtereddf), '\t')
# end

# Arrays for plotting in celcius degrees and gigapascals.
plottingT = filtereddf.T[1:resolution] .- 273.15
plottingP = filtereddf.P[1:resolution:resolution^2] ./ 1e4

# Use LaTeX font for plotting.
# with_theme(theme_latexfonts()) do
    # Define light rose-ish background color for the plots.
    f = Figure(backgroundcolor = RGBf(0.9, 0.79, 0.78), size = (1000, 400), fonts = (;regular = "TeX Gyre Pagella Makie"))
    # Layout with square subplots. 
    ga = f[1, 1] = GridLayout()
    # Define axis features for left, center and right subplots.
    axleft = Axis(ga[1, 1], ylabel = "pressure (GPa)", tellwidth = false, width = 300, 
        ylabelsize = 12, yticksize = 3, xticksize = 3, xticklabelsize = 10, yticklabelsize = 10, 
        xticklabelpad = 0.5, yticklabelpad = 0.5)
    axcenter = Axis(ga[1, 2], xlabel = L"temperature ($\degree$C)", tellwidth = false, width = 300,
        xlabelsize = 12, xticksize = 3, xticklabelsize = 10, xticklabelpad = 0.5)
    axright = Axis(ga[1, 3], tellwidth = false, width = 300, xticksize = 3, xticklabelsize = 10, 
        xticklabelpad = 0.5)
    # Align/link axis.
    linkyaxes!(axleft, axcenter, axright)
    linkxaxes!(axleft, axcenter, axright)

    # Plot the heatmaps for each three properties.
    hmdensity = CairoMakie.heatmap!(axleft, plottingT, plottingP, 
        reshape(filtereddf.density, resolution, resolution), colormap = Reverse(:imola))
    hmCp = CairoMakie.heatmap!(axcenter, plottingT, plottingP, 
        reshape(filtereddf.Cp, resolution, resolution), colormap = Reverse(:glasgow))
    hmalpha = CairoMakie.heatmap!(axright, plottingT, plottingP, 
        reshape(filtereddf.alphafilt, resolution, resolution) .* 1e4, colormap = Reverse(:nuuk)) 
    # Create the colorbars for heatmaps.
    cbdensity = Colorbar(ga[0, 1][1, 1], hmdensity, label = L"density (kg/m$^3$)", vertical = false, 
        tellwidth = false, width = 270, spinewidth = 0.1, labelpadding = 1.5, labelsize = 12, 
        ticklabelpad = 0.5, ticklabelsize = 10, ticksize = 3.5, 
        ticks = [round(minimum(filtereddf.density), sigdigits = 4), 2500, 3000, 3500, 4000, 
        round(maximum(filtereddf.density), sigdigits = 4)],
        minorticks = [2250, 2750, 3250, 3750, 4250], minorticksvisible = true, minorticksize = 2.5, 
        minortickwidth = 0.75, labelfont=:regular)
    cbCp = Colorbar(ga[0, 2][1, 1], hmCp, label = L"heat capacity (J/$\degree$C)", vertical = false, 
        tellwidth = false, width = 270, spinewidth = 0.1, labelpadding = 1.5, labelsize = 12, 
        ticklabelpad = 0.5, ticklabelsize = 10, ticksize = 3.5, ticks = [750, 1000, 1250, 1500], 
        minorticks = [1000, 3000, 5000], minorticksvisible = true, minorticksize = 2.5, 
        minortickwidth = 0.75)
    cbalpha = Colorbar(ga[0, 3][1, 1], hmalpha, label = L"thermal expansivity ($10^{-4}/\degree$C)",
        vertical = false, tellwidth = false, width = 270, spinewidth = 0.1, labelpadding = 1.5, 
        labelsize = 12, ticklabelpad = 0.5, ticklabelsize = 10, ticksize = 3.5,
        ticks = [round(minimum(filtereddf.alphafilt), sigdigits = 3) * 1e4, 1, 2, 
        round(maximum(filtereddf.alphafilt), sigdigits = 3) * 1e4],
        minorticks = [0.5, 1.5, 2.5], minorticksvisible = true, minorticksize = 2.5, 
        minortickwidth = 0.75)
    
    # Bring the grid up to make it visible.
    CairoMakie.translate!(hmdensity, 0, 0, -100)
    CairoMakie.translate!(hmCp, 0, 0, -100)
    CairoMakie.translate!(hmalpha, 0, 0, -100)

    # Do some beautification for the plot.
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
# end