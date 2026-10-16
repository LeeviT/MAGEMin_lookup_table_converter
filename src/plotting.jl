# Figures of the lookup table, saved into the run's output folder `prefix`/:
#   <prefix>_<variable>.png – one figure per mapped variable (see plot_variables)
#   <prefix>_summary.png    – 2×2 figure: density and Cp (top), thermal expansivity and flow law (bottom)
# For the ZOOM_DATABASES, all figures are drawn again with the axes cut to the `zoom` limits and saved with
# a _zoom postfix (e.g. <prefix>_density_zoom.png). The zoomed figures keep the colour ranges and colorbar
# ticks of the full figures, so both can be compared directly.

# Databases that also get zoomed figures.
const ZOOM_DATABASES = ("ig", "mp")

"""Temperature (°C) and pressure (GPa) axes of a (P, T)-sorted table with T varying fastest."""
function plotting_axes(table, resolution)
    plottingT = table.T[1:resolution] .- 273.15
    plottingP = table.P[1:resolution:resolution^2] ./ 1e4
    return plottingT, plottingP
end

"""Reshape a flat table column (T varying fastest) into the (T, P) grid matrix used by heatmap!."""
togrid(v, plottingT, plottingP) = reshape(v, length(plottingT), length(plottingP))

"""Path of the figure of `name` in the output folder `prefix`: prefix/prefix_name<suffix>.png."""
figure_path(prefix, name, suffix = "") = joinpath(prefix, "$(prefix)_$(name)$(suffix).png")

# Tick spacing of the zoomed figures: temperature (°C) and pressure (GPa).
const ZOOM_T_TICK_STEP = 250
const ZOOM_P_TICK_STEP = 1

"""Cut the axis `ax` of a heatmap on `plottingT` × `plottingP` to T ≤ `zoom.T_max_C` and P ≤ `zoom.P_max_GPa`,
with ticks every ZOOM_T_TICK_STEP °C and ZOOM_P_TICK_STEP GPa. The lower limits are the heatmap cell edges,
as in the full figures."""
function zoom_axis!(ax, plottingT, plottingP, zoom)
    T_edge = plottingT[1] - (plottingT[2] - plottingT[1]) / 2
    P_edge = plottingP[1] - (plottingP[2] - plottingP[1]) / 2
    limits!(ax, T_edge, zoom.T_max_C, P_edge, zoom.P_max_GPa)
    ax.xticks = 0:ZOOM_T_TICK_STEP:zoom.T_max_C
    ax.yticks = 0:ZOOM_P_TICK_STEP:zoom.P_max_GPa
    return ax
end

# Tick spacing of the full-range figures: temperature (°C) and pressure (GPa).
const FULL_T_TICK_STEP = 500
const FULL_P_TICK_STEP = 10

"""Ticks every `step` from 0 up to `vmax`, plus `vmax` itself. If `vmax` is off the grid but within a quarter
step of the last grid tick, it replaces that tick so that the labels don't overlap."""
function ticks_with_max(step, vmax)
    ticks = collect(0:step:vmax)
    isapprox(ticks[end], vmax) && return ticks
    vmax - ticks[end] < step / 4 && pop!(ticks)
    return [ticks; vmax]
end

"""Ticks of a full-range axis: temperature every FULL_T_TICK_STEP °C plus the maximum temperature, pressure
every FULL_P_TICK_STEP GPa."""
function full_axis_ticks!(ax, plottingT, plottingP)
    ax.xticks = ticks_with_max(FULL_T_TICK_STEP, plottingT[end])
    ax.yticks = 0:FULL_P_TICK_STEP:plottingP[end]
    return ax
end

"""Map phase names to integers for heatmaps. Returns (values, sorted unique names, number of names)."""
function phase_to_numeric(phase_col)
    unique_ph = sort(unique(phase_col))
    n = length(unique_ph)
    mapping = Dict(p => Float64(i) for (i, p) in enumerate(unique_ph))
    return [mapping[p] for p in phase_col], unique_ph, n
end

# Colour of every flow law label; mp reuses the colours of the ig regimes with the same output index.
const FLOW_LAW_COLOURS = Dict(
    "blueschist" => "#0081a7", "greenschist" => "#588157", "metastable protolith" => "#e9c46a",
    ">15% melt" => "#e07a5f", "eclogite" => "#735d78",
    "quartz" => "#0081a7", "feldspar" => "#588157", "garnet" => "#735d78",
    "olivine" => "#a1cca5", "ring/wad" => "#ad6a6c", "bridgmanite+fp" => "#5d2e46")
# Colorbar text of labels that are too long for it.
const FLOW_LAW_SHORT_LABELS = Dict("metastable protolith" => "metast. protol.")

"""Flow law labels of `database` that are left out of the flow law figure because they cannot occur with the
ig flow law scheme `ig_flow_laws` (IG_FLOW_LAW_LABELS): blueschist and greenschist with "protolith", metastable
protolith with "blueschist"."""
hidden_flow_laws(database, ig_flow_laws) =
    database in IG_DATABASES ? setdiff(keys(FLOW_LAW_PLOT_MAP_IG), IG_FLOW_LAW_LABELS[ig_flow_laws]) : String[]

"""Categorical flow law map (values, heatmap keywords, colorbar keywords) from the flow law labels of `table`.
The categories are ordered by their plot value (flow_law_plot_map), i.e. by the output index of the rheology
table, and the hidden_flow_laws are left out (kept, with a warning, if they occur after all)."""
function flow_law_map(table, database, ig_flow_laws)
    plot_map = flow_law_plot_map(database)
    present = Set(table.flow_law)
    hidden = hidden_flow_laws(database, ig_flow_laws)
    for label in hidden
        label in present && @warn "Flow law \"$label\" occurs although it should not with these settings; shown anyway."
    end
    labels = sort([l for l in keys(plot_map) if !(l in hidden) || l in present]; by = l -> plot_map[l])
    position = Dict(l => Float64(i) for (i, l) in enumerate(labels))
    n = length(labels)
    heatmap = (; colormap = cgrad([FLOW_LAW_COLOURS[l] for l in labels], n, categorical = true), colorrange = (0.5, n + 0.5))
    colorbar = (; label = "flow law", ticks = (collect(1.0:n), [get(FLOW_LAW_SHORT_LABELS, l, l) for l in labels]))
    return [position[l] for l in table.flow_law], heatmap, colorbar
end

# Colour ranges of the density (kg/m3), heat capacity (J/K/kg) and thermal expansivity (1e-5 1/K) figures,
# the same for every database so that the figures can be compared. Cp and α use the clamps of the tables
# (CP_RANGE, ALPHA_RANGE in table_processing.jl). Density uses a narrower range than the tables (DENSITY_MIN
# to DENSITY_MAX, no upper bound for the mantle databases); the plotted values are clamped to it.
const DENSITY_PLOT_RANGE = (2600.0, 4400.0)
const CP_PLOT_RANGE      = CP_RANGE
const ALPHA_PLOT_RANGE   = ALPHA_RANGE .* 1e5

# Fixed density colorbar ticks (kg/m3), and the unlabelled minor ticks between them.
const DENSITY_INNER_TICKS = 2500:500:4500
const DENSITY_MINOR_TICKS = [2250, 2750, 3250, 3750, 4250]
# Minimum distance of an inner density tick from the min/max ticks, as a fraction of the density range
# (about one label width on the colorbar).
const DENSITY_TICK_GAP = 0.1

"""Density colorbar ticks for the colour range `range`: its ends plus DENSITY_INNER_TICKS, leaving out inner
ticks that would overlap the end labels."""
density_ticks(range) = minmax_ticks(range, DENSITY_INNER_TICKS, density_tick_gap(range))
density_tick_gap(range) = DENSITY_TICK_GAP * (maximum(range) - minimum(range))

"""Unlabelled density minor ticks: DENSITY_MINOR_TICKS plus the inner ticks left out by density_ticks."""
density_minor_ticks(range) =
    sort([DENSITY_MINOR_TICKS; dropped_ticks(range, DENSITY_INNER_TICKS, density_tick_gap(range))])

"""Colorbar ticks: the minimum and maximum of `v` (the colour range; to `sigdigits` significant digits) plus
the fixed `inner` ticks, leaving out inner ticks within `gap` of the minimum or maximum so that the labels
don't overlap."""
function minmax_ticks(v, inner, gap; sigdigits = 4)
    lo, hi = minimum(v), maximum(v)
    # The end ticks are rounded inwards (ceil/floor), since a tick just outside the colour range is not drawn.
    # The tiny factors keep unit-conversion results like 7.4999999 from being rounded to 7.49.
    lo_tick = ceil(lo * (1 - 1e-12); sigdigits)
    hi_tick = floor(hi * (1 + 1e-12); sigdigits)
    return [lo_tick; filter(t -> lo + gap < t < hi - gap, inner); hi_tick]
end

"""The `inner` ticks inside the range of `v` that minmax_ticks leaves out for being within `gap` of the
minimum or maximum; they can be drawn as minor ticks without a label."""
function dropped_ticks(v, inner, gap)
    lo, hi = minimum(v), maximum(v)
    return filter(t -> (lo < t <= lo + gap) || (hi - gap <= t < hi), collect(inner))
end

# Heat capacity colorbar ticks (J/K/kg): labelled every 200 from the ends of CP_PLOT_RANGE (600, 800, …, 1600),
# unlabelled halfway between them (700, 900, …, 1500).
const CP_TICKS       = CP_PLOT_RANGE[1]:200:CP_PLOT_RANGE[2]
const CP_MINOR_TICKS = CP_PLOT_RANGE[1] + 100:200:CP_PLOT_RANGE[2]
# Thermal expansivity colorbar ticks (1e-5 1/K): labelled every 2 over ALPHA_PLOT_RANGE (−2, 0, …, 8),
# unlabelled halfway between them (−1, 1, …, 7).
const ALPHA_TICKS       = ALPHA_PLOT_RANGE[1]:2:ALPHA_PLOT_RANGE[2]
const ALPHA_MINOR_TICKS = ALPHA_PLOT_RANGE[1] + 1:2:ALPHA_PLOT_RANGE[2]

# Categorical palette of the phase maps (20 colours, repeated if there are more phases).
const PHASE_PALETTE = cgrad(:tab20, categorical = true)

"""Colour of each phase name in `names` (a sorted list of all phases in the dominant and second dominant
phase maps), so that a phase has the same colour in both maps."""
phase_colors(names) = Dict(p => PHASE_PALETTE[mod1(i, length(PHASE_PALETTE))] for (i, p) in enumerate(names))

"""Categorical heatmap and colorbar settings for a phase-name column. `colors` maps each phase name to its
colour (phase_colors); the colorbar lists only the phases present in this column."""
function phase_map(names_col, label, colors)
    values, names, n = phase_to_numeric(names_col)
    heatmap = (; colormap = cgrad([colors[p] for p in names], n, categorical = true), colorrange = (0.5, n + 0.5))
    colorbar = (; label, ticks = (collect(1.0:n), names))
    return values, heatmap, colorbar
end

"""Every mapped variable as a Dict of name => (values, heatmap keywords, colorbar keywords). The name is
also the figure file suffix. Density, Cp and α are clamped to their fixed colour ranges (DENSITY_PLOT_RANGE,
CP_PLOT_RANGE, ALPHA_PLOT_RANGE). The ig flow law scheme `ig_flow_laws` decides the flow law categories shown
(flow_law_map)."""
function plot_variables(table, database, ig_flow_laws)
    density  = clamp.(table.density, DENSITY_PLOT_RANGE...)
    cp       = clamp.(table.Cp, CP_PLOT_RANGE...)
    alpha_e5 = clamp.(table.alphafilt .* 1e5, ALPHA_PLOT_RANGE...)  # the α written to the tables, in 1e-5 1/K
    colors = phase_colors(sort(unique(vcat(table.dominant_phase, table.second_dominant_phase))))
    dominant = phase_map(table.dominant_phase, "dominant phase", colors)
    second   = phase_map(table.second_dominant_phase, "second dominant phase", colors)
    return Dict(
        "density" => (density, (; colormap = Reverse(:imola10), colorrange = DENSITY_PLOT_RANGE),
            (; label = rich("density (kg/m", superscript("3"), ")"), ticks = density_ticks(DENSITY_PLOT_RANGE),
               minorticks = density_minor_ticks(DENSITY_PLOT_RANGE), minorticksvisible = true, minorticksize = 3, minortickwidth = 0.75)),
        "cp" => (cp, (; colormap = Reverse(:glasgow), colorrange = CP_PLOT_RANGE),
            (; label = rich("heat capacity (J kg", superscript("−1"), " K", superscript("−1"), ")"), ticks = CP_TICKS,
               minorticks = CP_MINOR_TICKS, minorticksvisible = true, minorticksize = 3, minortickwidth = 0.75)),
        "alpha" => (alpha_e5, (; colormap = :bilbao, colorrange = ALPHA_PLOT_RANGE),
            (; label = rich("thermal expansivity (10", superscript("−5"), " K", superscript("−1"), ")"), ticks = ALPHA_TICKS,
               minorticks = ALPHA_MINOR_TICKS, minorticksvisible = true, minorticksize = 3, minortickwidth = 0.75)),
        "dominant_phase" => dominant,
        "second_dominant_phase" => second,
        "melt_fraction" => (table.liq_frac_vol, (; colormap = :thermal), (; label = "liq (melt) vol. fraction")),
        "flow_law" => flow_law_map(table, database, ig_flow_laws),
    )
end

# Figure fonts: Figtree, loaded from the static font files in Figtree/static/ next to src/ (copy that
# folder along with the .jl files).
# Labels are plain or rich-text strings, not LaTeXStrings, since Makie renders LaTeX in its own math font.
# Exponents use rich-text superscripts with U+2212 (−): Figtree has no superscript minus (⁻).
figure_font(style) = joinpath(@__DIR__, "..", "Figtree", "static", "Figtree-$(style).ttf")
const FIGURE_THEME = Theme(fonts = (; regular = figure_font("Regular"), bold = figure_font("Bold"),
    italic = figure_font("Italic"), bold_italic = figure_font("BoldItalic"), medium = figure_font("Medium")))

# Shared axis and colorbar styles of all figures (medium-weight axis and colorbar titles; every colorbar uses the same
# tick label size, COLORBAR_STYLE.ticklabelsize); every panel has a square plot area (the figures grow to fit
# the colorbars and labels).
const AXIS_STYLE = (; xlabel = "temperature (°C)", ylabel = "pressure (GPa)", width = 600, height = 600,
    xlabelsize = 24, ylabelsize = 24, xlabelfont = :medium, ylabelfont = :medium,
    xticksize = 6, yticksize = 6, xticklabelsize = 18, yticklabelsize = 18)
const COLORBAR_STYLE = (; vertical = false, tellwidth = false, width = 450, size = 20, spinewidth = 0.1,
    labelpadding = 1.5, labelsize = 19, labelfont = :medium, ticklabelpad = 0.5, ticklabelsize = 14, ticksize = 4)

# Resolution of the saved PNGs. Makie figure units are CSS pixels (1/96 inch), so px_per_unit = dpi / 96.
const FIGURE_DPI = 300
const PX_PER_UNIT = FIGURE_DPI / 96

# Figure padding (left, right, bottom, top); the extra right padding fits the last temperature tick label,
# which sits on the right edge of the zoomed axes.
const FIGURE_PADDING = (16, 32, 16, 16)

"""Draw one variable `(values, heatmap keywords, colorbar keywords)` into the layout cell `gl`: a P-T
heatmap with temperature and pressure labels and ticks, and a horizontal colorbar above it. With `zoom` =
(P_max_GPa, T_max_C), the axis is cut to those limits (zoom_axis!), otherwise it gets the full-range ticks
(full_axis_ticks!). Returns the axis."""
function draw_panel!(gl, variable, plottingT, plottingP; zoom = nothing)
    values, heatmap_kw, colorbar_kw = variable
    ax = Axis(gl[2, 1]; AXIS_STYLE...)
    hm = CairoMakie.heatmap!(ax, plottingT, plottingP, togrid(values, plottingT, plottingP); heatmap_kw...)
    Colorbar(gl[1, 1], hm; COLORBAR_STYLE..., colorbar_kw...)
    CairoMakie.translate!(hm, 0, 0, -100)  # draw ticks and grid in front of the heatmap
    hidespines!(ax)
    zoom === nothing ? full_axis_ticks!(ax, plottingT, plottingP) : zoom_axis!(ax, plottingT, plottingP, zoom)
    rowgap!(gl, 7.5)
    return ax
end

"""Save one figure per mapped variable, with a square plot area, as `prefix`/`prefix`_<name><suffix>.png."""
function plot_variable_figures(variables, plottingT, plottingP, prefix, suffix, zoom)
    for (name, variable) in variables
        with_theme(FIGURE_THEME) do
            f = Figure(backgroundcolor = :transparent, figure_padding = FIGURE_PADDING)
            gl = f[1, 1] = GridLayout()
            draw_panel!(gl, variable, plottingT, plottingP; zoom)
            resize_to_layout!(f)
            save(figure_path(prefix, name, suffix), f; px_per_unit = PX_PER_UNIT)
        end
    end
end

"""Save the 2×2 summary figure (density and Cp on top, thermal expansivity and flow law at the bottom)
as `prefix`/`prefix`_summary<suffix>.png."""
function plot_summary(variables, plottingT, plottingP, prefix, suffix, zoom)
    with_theme(FIGURE_THEME) do
        f = Figure(backgroundcolor = :transparent, figure_padding = FIGURE_PADDING)
        ga = f[1, 1] = GridLayout()
        layout = [("density", 1, 1), ("cp", 1, 2), ("alpha", 2, 1), ("flow_law", 2, 2)]
        axes = map(layout) do (name, row, col)
            gl = ga[row, col] = GridLayout()
            ax = draw_panel!(gl, variables[name], plottingT, plottingP; zoom)
            # Shared axes: temperature labels and values only on the bottom row, pressure only on the left
            # column (tick marks and grid lines are kept).
            row == 1 && hidexdecorations!(ax, grid = false, ticks = false)
            col == 2 && hideydecorations!(ax, grid = false, ticks = false)
            ax
        end
        linkxaxes!(axes...)
        linkyaxes!(axes...)
        colgap!(ga, 50)  # room for the last temperature tick label of the left column
        rowgap!(ga, 20)
        resize_to_layout!(f)
        display(f)
        save(figure_path(prefix, "summary", suffix), f; px_per_unit = PX_PER_UNIT)
    end
end

"""Draw the single-variable figures and the summary with file-name `suffix`, cut to `zoom` unless it is nothing."""
function plot_figure_set(variables, plottingT, plottingP, prefix, suffix, zoom)
    plot_variable_figures(variables, plottingT, plottingP, prefix, suffix, zoom)
    plot_summary(variables, plottingT, plottingP, prefix, suffix, zoom)
end

"""Draw all figures and save them into the output folder `prefix`: first for the full P-T range, then,
for the ZOOM_DATABASES, with the axes cut to `zoom` = (P_max_GPa, T_max_C) and a _zoom postfix. Both sets
share the colour scales and colorbar ticks: fixed ranges for density, Cp and α, the others fitted to the full table.
`ig_flow_laws` is the ig flow law scheme ("protolith" or "blueschist")."""
function plot_results(table, resolution, database, prefix, zoom, ig_flow_laws)
    plottingT, plottingP = plotting_axes(table, resolution)
    variables = plot_variables(table, database, ig_flow_laws)
    plot_figure_set(variables, plottingT, plottingP, prefix, "", nothing)
    database in ZOOM_DATABASES && plot_figure_set(variables, plottingT, plottingP, prefix, "_zoom", zoom)
end
