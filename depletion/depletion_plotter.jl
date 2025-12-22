using DataFrames, CSV, CairoMakie, TidierData, DelimitedFiles, ColorSchemes

function depletionscaling(depth)
    deplfunction = zeros(length(depth))
    # deplfunction = 1.0
    for i in eachindex(depth)
        deplfunction[i] =         depth[i] ≤ 20.0  ? 0.0                             :
                           20.0 < depth[i] ≤ 90.0  ? (-1 / 350) * depth[i] + 16 / 35 :
                           90.0 < depth[i] ≤ 115.0 ? 0.2 * ((115 - depth[i]) / 25)   :
                           #= otherwise =#           0.0
    end
    return deplfunction
end

inputtablefolder = "."
inputtablefilename = "40Ma_nolookup.csv"
inputtablefilename2 = "40Ma_depllookup.csv"
inputtablefilename3 = "40Ma_discontdepllookup.csv"
inputtablefilepath = joinpath(inputtablefolder, inputtablefilename)
inputtablefilepath2 = joinpath(inputtablefolder, inputtablefilename2)
inputtablefilepath3 = joinpath(inputtablefolder, inputtablefilename3)

# Read in a table to initialize a dataframe.
df_40Ma_nolookup = DataFrame(CSV.File(inputtablefilepath))
df_40Ma_depllookup = DataFrame(CSV.File(inputtablefilepath2))
df_40Ma_discontdepllookup = DataFrame(CSV.File(inputtablefilepath3))

df_40Ma_nolookup = @chain df_40Ma_nolookup begin
    @select(depth = var"arc_length", T = var"T", ρ = var"density")
    @arrange(depth, T, ρ)
    @mutate(depth = depth ./ 1e3)  # convert from m to km
end

df_40Ma_depllookup = @chain df_40Ma_depllookup begin
    @select(depth = var"arc_length", T = var"T", ρ = var"density")
    @arrange(depth, T, ρ)
    @mutate(depth = depth ./ 1e3)  # convert from m to km
end

df_40Ma_discontdepllookup = @chain df_40Ma_discontdepllookup begin
    @select(depth = var"arc_length", T = var"T", ρ = var"density")
    @arrange(depth, T, ρ)
    @mutate(depth = depth ./ 1e3)  # convert from m to km
end

df_40Ma_depllookup.T = reverse(df_40Ma_depllookup.T)
df_40Ma_depllookup.ρ = reverse(df_40Ma_depllookup.ρ)
# df_40Ma_discontdepllookup.T = reverse(df_40Ma_discontdepllookup.T)
# df_40Ma_discontdepllookup.ρ = reverse(df_40Ma_discontdepllookup.ρ)  
df_40Ma_nolookup.T = reverse(df_40Ma_nolookup.T)
df_40Ma_nolookup.ρ = reverse(df_40Ma_nolookup.ρ)    

depl_profile = depletionscaling(df_40Ma_nolookup.depth)
ρ_depl_40Ma = (1 .- 0.04 .* depl_profile) .* df_40Ma_nolookup.ρ

# f = Figure(backgroundcolor = RGBf(0.9, 0.79, 0.78), size = (900, 700))
f = Figure(backgroundcolor = :white, size = (1100, 700))
ax1 = CairoMakie.Axis(f[1, 1], limits = (3000, 3400, 0, 150), xlabel = "density (kg/m³)", ylabel = "depth (km)", yreversed = true, 
        tellwidth = false, width = 750, xlabelsize = 20, ylabelsize = 20, yticksize = 5, xticksize = 5, xticklabelsize = 18, 
        yticklabelsize = 18, xticklabelpad = 0.5, yticklabelpad = 0.5, yticks = [0, 20, 30, 50, 60, 90, 115, 150])
ax2 = CairoMakie.Axis(f[1, 1], limits = (0, 40, 0, 150), xlabel = "melt depletion (%)", xaxisposition = :top, yreversed = true, 
        tellwidth = false, width = 750, xlabelsize = 20, ylabelsize = 20, yticksize = 5, xticksize = 5, xticklabelsize = 18, 
        yticklabelsize = 18, xticklabelpad = 0.5, yticklabelpad = 0.5, xticklabelcolor = colorschemes[:buda10].colors[4],
        xlabelcolor = colorschemes[:buda10].colors[4], xtickcolor = colorschemes[:buda10].colors[4], yticksvisible = false, 
        yticklabelsvisible = false, ygridvisible = false)

crustplot = hspan!(20, 0, alpha = 0.35, color = colorschemes[:hawaii100].colors[40])
depllayerplot = hspan!(115, 20, alpha = 0.35, color = colorschemes[:navia100].colors[65])
asthplot = hspan!(150, 115, alpha = 0.35, color = colorschemes[:navia100].colors[80])

depl_profile_plot = CairoMakie.lines!(ax2, depl_profile .* 100, df_40Ma_nolookup.depth, color = colorschemes[:buda10].colors[4], linewidth = 1.5, alpha = 1.0, label = "variable depletion profile")
# nolookup_nodepl_plot = CairoMakie.lines!(ax1, df_40Ma_nolookup.ρ, df_40Ma_nolookup.depth, color = colorschemes[:managua10].colors[2], linewidth = 3.5, linestyle = (:dash, :dense), label = "no depletion")
lookup_discontdepl_plot = CairoMakie.lines!(ax1, df_40Ma_discontdepllookup.ρ, df_40Ma_discontdepllookup.depth, color = colorschemes[:managua10].colors[2], linewidth = 3.5, linestyle = (:dash, :dense), label = "discontinuous depletion")
lookup_depl_plot = CairoMakie.lines!(ax1, df_40Ma_depllookup.ρ, df_40Ma_depllookup.depth, color = colorschemes[:managua10].colors[5], linewidth = 3.5, linestyle = (:dash, :dense), label = "constant 20% depletion\nusing lookup table")
nolookup_depl_plot = CairoMakie.lines!(ax1, ρ_depl_40Ma, df_40Ma_nolookup.depth, color = colorschemes[:managua10].colors[8], linewidth = 3.5, linestyle = (:dash, :dense), label = "variable depletion from 40% to\n0%, using analytical function\nto calculate density reduction") 

text!(15, 10, text = "oceanic crust", align = (:center, :center), alpha = 0.3)
text!(15, 70, text = "depleted mantle layer", align = (:center, :center), alpha = 0.3)
text!(15, 132.5, text = "chemical asthenosphere/ \n fertile mantle", align = (:center, :center), alpha = 0.3)

f[1, 2] = Legend(f, ax2, framevisible = false, halign = :left, padding = (6.0f0, 6.0f0, 180.0f0, 6.0f0))
f[1, 2] = Legend(f, ax1, "density profiles", framevisible = false, halign = :left) 

CairoMakie.translate!(depl_profile_plot, 0, 0, -50)
CairoMakie.translate!(crustplot, 0, 0, -100)
CairoMakie.translate!(depllayerplot, 0, 0, -100)
CairoMakie.translate!(asthplot, 0, 0, -100)

display(f)