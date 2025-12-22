# This script reads in a CSV file containing phase transition lookup table in MAGEMin format and 
# then scales its density to take account porosity near Earth's surface.

# DataFrames for better handling table format, CSV for readind in a CSV file, CairoMakie for 
# plotting, TidierData for filtering/manipulating data, DelimitedFiles for output.
using DataFrames, CSV, CairoMakie, TidierData, DelimitedFiles, MAGEMin_C

# A function which scales the density based on the near-surface porosity. Uses the empirical
# model by Chen et al. (2020), https://doi.org/10.1007/s10040-020-02214-x.
function porosityscaling(ρ, pressure, resolution, crusttype)
    # Select the relevant compositional parameters, depending the crustal type.
    if crusttype == "oceanic"
        ϕ0 = 0.678
        m = 0.008
        n = 89.53
    elseif crusttype == "continental" 
        ϕ0 = 0.474
        m = 0.071
        n = 5.989
    end
    # Other parameters, assumes constant density above the calculated point for simplicity.
    ρ0 = 2700.0
    g = 9.81
    # First, convert lithostatic pressure p=\rho*g*h to h=depth in km.
    depth = (100.0 .* pressure[1:resolution:resolution^2]) / (ρ0 * g)
    ϕ = ϕ0 ./ (1 .+ m .* depth) .^ n
    ρ = reshape((reshape(ρ, resolution, resolution)' .* (1.0 .- ϕ) .+ 1000.0 .* ϕ)', resolution^2)
    return ρ
end

# --- Input parameters start ----
lookuptabletype = "interface"

if lookuptabletype == "file"
    # Filename of a lookup table to process, let us assume it is located in a same folder as this 
    # processor script.
    inputtablefolder = "."
    inputtablefilename = "nmorb-gabbro-reso263k-p15gpa-t2000c-h20wt05_ig.csv"
    inputtablefilepath = joinpath(inputtablefolder, inputtablefilename)
    # Read in a table to initialize a dataframe.
    df = DataFrame(CSV.File(inputtablefilepath))
elseif lookuptabletype == "interface"
    n = 256
    resolution = n
    Prange = repeat(range(0.01, 150, n), outer = n)
    Trange = repeat(range(0, 1700, n), inner = n)
    db = Initialize_MAGEMin("ig", verbose = false)
    # data = use_predefined_bulk_rock(db, 4)
    Xoxides = ["SiO2"; "Al2O3"; "CaO"; "MgO"; "FeO"; "Fe2O3"; "K2O"; "Na2O"; "TiO2"; "Cr2O3"; "H2O"] # Stoner et al. composition
    X = [50.03; 15.01; 11.26; 7.70; 9.73; 0.0; 0.01; 2.81; 1.52; 0.01; 0.5] # Stoner et al. composition
    # Xoxides = ["SiO2"; "Al2O3"; "CaO"; "MgO"; "FeO"; "K2O"; "Na2O"; "TiO2"; "O"; "MnO"; "H2O"] garnet-migmatite
    # X = [69.66; 13.76; 1.77; 1.73; 4.32; 2.61; 2.41; 0.8; 0.03; 0.07; 2.82] garnet-migmatite
    sys_in = "wt"
    df = DataFrame(multi_point_minimization(Prange, Trange, db, X = X, Xoxides = Xoxides, sys_in = sys_in))
    # df = DataFrame(multi_point_minimization(Prange, Trange, db, test = 4))
    Finalize_MAGEMin(db)
end

# Regardless of the lookup table creation method, save it as a .txt file which is the only
# format that ASPECT reads.
outputtablefolder = "."
outputtablefilename = "nmorb-gabbro-h20wt05-delay_ig.txt"
outputtablefilepath = joinpath(outputtablefolder, outputtablefilename)

# ---- Input parameters end ----

if lookuptabletype == "file"
    # A macro to filter and modify a dataframe in MAGEMin format
    filtereddf = @chain df begin
        # Select the relevant columns and rename them, as special characters such as square brackets are
        # messy to handle in column names in Julia.
        TidierData.@select(P = var"P[kbar]", T = var"T[°C]", var"phase", density = var"density[kg/m3]", 
            Cp = var"heatCapacity[J/K]", alpha = var"alpha[1/K]", enthalpy = var"Enthalpy[J]", 
            Vp = var"Vp[km/s]", Vs = var"Vs[km/s]")
        # Filter out data of individual minerals and preserve bulk rock (=system) composition data.
        @filter(phase == "system")
        # Convert pressure from kbar to bar and temperature from celsius to kelvin. 
        @mutate(P = 1e3 * P, T = 273.15 + T)
        # Drop the 'phase' column and reorder the columns for PerpleX format.
        TidierData.@select(T, P, density, alphafilt, Cp, Vp, Vs, enthalpy)
        @arrange(P, T)
    end
    # Resolution of the lookup table
    resolution = Int.(sqrt(length(filtereddf.T)))
    # Remove the square brackets from the heat capacity column, and convert it from string to float.
    filtereddf.Cp = strip.(filtereddf.Cp, '[')
    filtereddf.Cp = strip.(filtereddf.Cp, ']')
    filtereddf.Cp = parse.(Float64, filtereddf.Cp)
elseif lookuptabletype == "interface"
    # A macro to filter and modify a dataframe in MAGEMin format
    filtereddf = @chain df begin
        # Select the relevant columns and rename them, as special characters such as square brackets are
        # messy to handle in column names in Julia.
        TidierData.@select(P = var"P_kbar", T = var"T_C", density = var"rho", Cp = var"s_cp", 
            alpha = var"alpha", enthalpy = var"enthalpy", Vp = var"Vp", Vs = var"Vs", bulk_S_wt = var"bulk_S_wt")
        # Convert pressure from kbar to bar and temperature from celsius to kelvin. 
        @mutate(P = 1e3 * P, T = 273.15 + T)
        # Drop the 'phase' column and reorder the columns for PerpleX format.
        TidierData.@select(T, P, density, alpha, Cp, Vp, Vs, enthalpy, bulk_S_wt)
        @arrange(P, T)
        # Rename the column for later use.
        @rename(H20minebound = bulk_S_wt)
    end
    # Convert from array[array] to array
    for (col, idx) in [(:H20minebound, 11), (:alpha, 1), (:Cp, 1), (:enthalpy, 1)]
        filtereddf[!, col] = getindex.(filtereddf[!, col], idx)
    end
end

#= Find the last temperature column by column index, which has temperature < 450 °C.
Tlimitcoli = findlast(x -> x < 600 + 273.15, filtereddf.T[1:resolution])
compdensity = reshape(filtereddf.density[1:resolution] 
   * exp.(1e-11 * 1e5 .* filtereddf.P[1:resolution:resolution^2])', resolution, resolution)
for i in eachindex(filtereddf.T)
   Pscaling = filtereddf.P[1:resolution:resolution^2]
   (i >= 1 || 1 == mod(i, resolution)) && Tlimitcoli > mod(i - 1, resolution) + 1 ? 
   filtereddf.density[i] = compdensity[i] : 0
=#
scalingfactor = ones(resolution)
for i in range(1, resolution)
    if filtereddf.T[i] < 600 + 273.15
        scalingfactor[i] = 0
    elseif filtereddf.T[i] >= 600 + 273.15 && filtereddf.T[i] < 900 + 273.15
        # println(filtereddf.T[i] - 273.15)
        scalingfactor[i] = ((filtereddf.T[i] - (600 + 273.15)) / 300)^3
    else
        scalingfactor[i] = 1
    end
end
filtereddf.scalingfactor = repeat(scalingfactor, outer = resolution)
filtereddf.density = porosityscaling(filtereddf.density, filtereddf.P, resolution, "oceanic")

filtereddf = @chain filtereddf begin
    # Convert subzero (=unphysical) heat expansivity and specific heat values to zero.
    @mutate(alphafilt = case_when(alpha < 0.0 => 0.0, alpha > 0.75e-4 => 0.75e-4, alpha > 0.0 => alpha))
    @mutate(Cp = case_when(Cp < 750.0 => 750.0, Cp > 1500.0 => 1500.0, Cp > 0.0 => Cp))
    # Apply kinetic delay scaling to density values.
    @mutate(density = 3000 * (1.0 - scalingfactor) + density .* scalingfactor)
    # Filter outlier density values for better stability of a geodynamic model.
    @mutate(density = case_when(density < 2600 => 2600, density > 4500 => 4500, density > 0 => density))
    # @mutate(density = density * (1.0 - 0.04 * 0.225))
end

conversionP = [3.12e4, 3.0e4] # in bars, depl. mantle, lower OC, upper OC
conversionT = 600 + 273.15 # in kelvins
# cutoffresolution = Int.(floor(resolution * conversionP[1] / maximum(filtereddf.P)))
# filtereddf.density[1:cutoffresolution*resolution] = filtereddf.density[1:cutoffresolution*resolution] .- 48.0
# Plimitcoli = [Int.(findlast(x -> x < conversionP[1], filtereddf.P[1:resolution^2])),
#               Int.(findlast(x -> x < conversionP[2], filtereddf.P[1:resolution^2]))]
# Tlimitcoli = findlast(x -> x < conversionT, filtereddf.T[1:resolution])
# for i in eachindex(filtereddf.T)
#    (i >= 1 || 1 == mod(i, resolution)) && Tlimitcoli <= mod(i - 1, resolution) + 1 && Plimitcoli[1] > i ? 
#    filtereddf.density[i] = filtereddf.density[i] - 48.0 : 0
# end
# for i in eachindex(filtereddf.T)
#    (i >= 1 || 1 == mod(i, resolution)) && Tlimitcoli > mod(i - 1, resolution) + 1 && Plimitcoli[2] > i ? 
#    filtereddf.density[i] = filtereddf.density[i] - 48.0 : 0
# end


open(outputtablefilepath, "w") do io
    # Write the metadata and headers in PerpleX format.
    write(io, 
        "|6.6.6\n"                                                  # PerpleX version
        * replace(outputtablefilename, ".txt" => ".tab") * "\n"     # Lookup table name
        * "\t\t2\n"                                                 # No. independent variables (T and P)
        * "T(K)\n"                                                  # Temperature parameters
        * "\t" * string(filtereddf.T[1]) * "\n"                     # The lower limit of temperature 
        * "\t" * string(filtereddf.T[2] - filtereddf.T[1]) * "\n"   # The upper limit of temperature
        * "\t\t" * string(resolution) * "\n"                        # Temperature resolution
        * "P(bar)\n"                                                # Pressure parameters
        * "\t" * string(filtereddf.P[1]) * "\n"                     # The lower limit of pressure
        * "\t" * string(filtereddf.P[resolution + 1] - filtereddf.P[1]) * "\n" # The upper limit of pressure
        * "\t\t" * string(resolution) * "\n"                        # Pressure resolution
        * "\t\t8\n"                                                 # Number of columns
        * "T(K)\tP(bar)\trho,kg/m3\talpha,1/K\tcp,J/K/kg\tvp,km/s\tvs,km/s\th,J/kg\n") # Column names and units
    # Write the actual phase diagram data.
    writedlm(io, eachrow(filtereddf), '\t')
end

# Arrays for plotting in celcius degrees and gigapascals.
plottingT = filtereddf.T[1:resolution] .- 273.15
plottingP = filtereddf.P[1:resolution:resolution^2] ./ 1e4

# Use LaTeX font for plotting.
with_theme(theme_latexfonts()) do
    # Define light rose-ish background color for the plots.
    f = Figure(backgroundcolor = :transparent, size = (2100, 700))
    # Layout with square subplots. 
    ga = f[1, 1] = GridLayout()
    # Define axis features for left, center and right subplots.
    axleft = Axis(ga[1, 1], xlabel = L"temperature ($\degree$C)", ylabel = "pressure (GPa)", tellwidth = false, width = 600, 
        xlabelsize = 20, ylabelsize = 20, yticksize = 5, xticksize = 5, xticklabelsize = 18, yticklabelsize = 18, 
        xticklabelpad = 0.5, yticklabelpad = 0.5)
    axcenter = Axis(ga[1, 2], xlabel = L"temperature ($\degree$C)", tellwidth = false, width = 600,
        xlabelsize = 20, xticksize = 5, xticklabelsize = 18, xticklabelpad = 0.5)
    axright = Axis(ga[1, 3], tellwidth = false, width = 600, xticksize = 5, xticklabelsize = 18, 
        xticklabelpad = 0.5)
    # Align/link axis.
    linkyaxes!(axleft, axcenter, axright)
    linkxaxes!(axleft, axcenter, axright)

    # Plot the heatmaps for each three properties.
    hmdensity = CairoMakie.heatmap!(axleft, plottingT, plottingP, 
        reshape(filtereddf.density, resolution, resolution), colormap = Reverse(:imola10))
    #= hmCp = CairoMakie.heatmap!(axcenter, plottingT, plottingP, 
        reshape(filtereddf.Cp, resolution, resolution), colormap = Reverse(:glasgow)) =#
    hmscaling = CairoMakie.heatmap!(axcenter, plottingT, plottingP, 
        reshape(filtereddf.scalingfactor, resolution, resolution), colormap = Reverse(:glasgow))
    #= hmalpha = CairoMakie.heatmap!(axright, plottingT, plottingP, 
        reshape(filtereddf.alphafilt, resolution, resolution) .* 1e5, colormap = Reverse(:nuuk)) =#
    hmH2Owt = CairoMakie.heatmap!(axright, plottingT, plottingP, 
        reshape(filtereddf.H20minebound .* 100, resolution, resolution), colormap = Reverse(:nuuk))
    # Create the colorbars for heatmaps.
    cbdensity = Colorbar(ga[0, 1][1, 1], hmdensity, label = L"density (kg/m$^3$)", vertical = false, 
        tellwidth = false, width = 350, spinewidth = 0.1, labelpadding = 1.5, labelsize = 16, 
        ticklabelpad = 0.5, ticklabelsize = 14, ticksize = 3.5, 
        ticks = [round(minimum(filtereddf.density), sigdigits = 4), 2500, 3000, 3500, 4000, 4500,
        floor(maximum(filtereddf.density), sigdigits = 4)],
        minorticks = [2250, 2750, 3250, 3750, 4250], minorticksvisible = true, minorticksize = 2.5, 
        minortickwidth = 0.75)
    #= cbCp = Colorbar(ga[0, 2][1, 1], hmCp, label = L"heat capacity (J/$\degree$C)", vertical = false, 
        tellwidth = false, width = 350, spinewidth = 0.1, labelpadding = 1.5, labelsize = 16, 
        ticklabelpad = 0.5, ticklabelsize = 14, ticksize = 3.5, ticks = [750, 1000, 1250, 1500, 
        floor(maximum(filtereddf.Cp), sigdigits = 4)]) =#
    cbCp = Colorbar(ga[0, 2][1, 1], hmscaling, label = "scaling factor", vertical = false, 
        tellwidth = false, width = 350, spinewidth = 0.1, labelpadding = 1.5, labelsize = 16, 
        ticklabelpad = 0.5, ticklabelsize = 14, ticksize = 3.5, ticks = [0.0, 0.25, 0.5, 0.75, 1.0])
    #= cbalpha = Colorbar(ga[0, 3][1, 1], hmalpha, label = L"thermal expansivity ($10^{-5}/\degree$C)",
        vertical = false, tellwidth = false, width = 350, spinewidth = 0.1, labelpadding = 1.5, 
        labelsize = 16, ticklabelpad = 0.5, ticklabelsize = 14, ticksize = 3.5,
        ticks = [round(minimum(filtereddf.alphafilt), sigdigits = 3) * 1e5, 2, 3, 4, 5,
        floor(maximum(filtereddf.alphafilt), sigdigits = 3) * 1e5]) =#
    cbwtH2O = Colorbar(ga[0, 3][1, 1], hmH2Owt, label = "bound water (wt. %)",
        vertical = false, tellwidth = false, width = 350, spinewidth = 0.1, labelpadding = 1.5, 
        labelsize = 16, ticklabelpad = 0.5, ticklabelsize = 14, ticksize = 3.5,
        ticks = [0, 0.5, 1.0, 1.5, 2.0, 2.5])

    # Pline = hlines!(axleft, 2.78, 0, 1400, color = "red", linewidth = 0.5)
    # Tline = vlines!(axleft, 565, 0, 15, color = "red", linewidth = 0.5)

    # Bring the grid up to make it visible.
    CairoMakie.translate!(hmdensity, 0, 0, -100)
    CairoMakie.translate!(hmscaling, 0, 0, -100)
    # CairoMakie.translate!(hmCp, 0, 0, -100)
    # CairoMakie.translate!(hmalpha, 0, 0, -100)
    CairoMakie.translate!(hmH2Owt, 0, 0, -100)
    # CairoMakie.translate!(Pline, 0, 0, 100)
    # CairoMakie.translate!(Tline, 0, 0, 100)

    # Do some beautification of a plot.
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
