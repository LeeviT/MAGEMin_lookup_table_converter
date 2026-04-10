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

# Quadratic kinetic-delay scaling used to blend low-temperature and high-temperature states.
function kinetic_delay_scaling(T)
    if T < 600 + 273.15
        return 0.0
    elseif T < 800 + 273.15
        return ((T - (600 + 273.15)) / 200)^2
    else
        return 1.0
    end
end

# --- Input parameters start ----
lookuptabletype = "interface"
database = "mtl"  # "ig" for igneous, "mtl" for mantle
usescaling = true
outputtablefolder = "."
outputtablefilename = "3200kgm3_mtl.txt"
outputtablefilepath = joinpath(outputtablefolder, outputtablefilename)
# ---- Input parameters end ----

function generate_dataframe(lookuptabletype, database)
    if lookuptabletype == "file"
        inputtablefolder = "."
        inputtablefilename = "nmorb-gabbro-reso263k-p15gpa-t2000c-h20wt05_ig.csv"
        inputtablefilepath = joinpath(inputtablefolder, inputtablefilename)
        df = DataFrame(CSV.File(inputtablefilepath))
        resolution = Int.(sqrt(nrow(df)))
        return df, resolution
    end
    # lookuptabletype == "interface"
    n = 64
    resolution = n
    Prange = repeat(range(0.01, 300, n), outer = n)
    Trange = repeat(range(0, 1800, n), inner = n)
    db = Initialize_MAGEMin(database, verbose = false)
    Xoxides = ["SiO2";    "Al2O3";   "CaO";     "MgO";     "FeO";    "Na2O"] # mantle oxides
    X =       [45.219696; 4.31225;   3.364569;  38.997326; 8.065189; 0.011995] # pyrolite
    # Xoxides = ["SiO2";    "Al2O3";   "CaO";     "MgO";    "FeO";    "K2O";    "Na2O";   "TiO2";   "O";     "Cr2O3";  "H2O"] # OC oxides
    # X =       [50.992811; 15.304009; 10.922183; 7.848858; 9.913537; 0.135232; 2.867050; 1.541733; 0.176099; 0.048487; 0.25] # gabbro
    # X =       [50.376296; 15.118980; 10.790131; 7.753963; 9.793680; 0.133597; 2.832387; 1.523093; 0.17397; 0.047901; 1.456] # basalt
    # Xoxides = ["SiO2"; "Al2O3"; "CaO"; "MgO"; "FeO"; "K2O"; "Na2O"; "TiO2"; "O"; "MnO"; "H2O"] # garnet-migmatite
    # X = [69.66; 13.76; 1.77; 1.73; 4.32; 2.61; 2.41; 0.8; 0.03; 0.07; 2.82] # garnet-migmatite
    # X = [64.84; 14.99; 3.50; 2.41; 4.90; 2.72; 3.18; 0.62; 0.0; 0.1; 2.71] # average upper continental crust
    sys_in = "wt"
    df = DataFrame(multi_point_minimization(Prange, Trange, db, X = X, Xoxides = Xoxides, sys_in = sys_in)) # For user-specified bulk composition
    Finalize_MAGEMin(db)
    return df, resolution
end

"""Helper: safely get a fraction from a vector by index, returning 0.0 if index is nothing."""
get_frac(idx, fracs) = idx === nothing ? 0.0 : fracs[idx]

"""Find the dominant and second-dominant solution phase by volume fractions."""
function find_dominant_phases(ss_names, ss_fracs_adj)
    isempty(ss_names) && return "none", "none"
    dom = ss_names[argmax(ss_fracs_adj)]
    if length(ss_names) < 2
        return dom, "none"
    end
    adj2 = copy(ss_fracs_adj)
    adj2[argmax(ss_fracs_adj)] = -Inf
    return dom, ss_names[argmax(adj2)]
end

# ── Flow law classification for "ig" (igneous) database ──

# Flow law plot value mapping for "ig" labels.
const FLOW_LAW_PLOT_MAP_IG = Dict(
    "blueschist" => 1.0, "eclogite" => 2.0, "greenschist" => 3.0, ">15% melt" => 4.0
)

"""Classify the flow law regime for the igneous database from weight fractions and temperature.
Returns (label, plot_value, phase_indices).
The blueschist→eclogite transition is discrete: eclogite is assigned only when
the kinetic-delay scaling reaches ≥ 0.5; below that it stays blueschist."""
function classify_flow_law_ig(blueschist_wt, greenschist_wt, eclogite_wt, liq_wt, T_C)
    scaling = kinetic_delay_scaling(273.15 + T_C)

    # Step 1: determine raw label from weight fractions
    raw = if liq_wt > 0.15
        ">15% melt"
    elseif eclogite_wt > 0.5 || (eclogite_wt > blueschist_wt && eclogite_wt > greenschist_wt)
        "eclogite"
    elseif (blueschist_wt + greenschist_wt) == 0.0
        "blueschist"
    elseif blueschist_wt >= greenschist_wt
        "blueschist"
    else
        "greenschist"
    end

    # Step 2: apply discrete kinetic delay to eclogite transition (threshold at 0.5)
    label = if raw == "eclogite"
        scaling >= 0.5 ? "eclogite" : "blueschist"
    else
        raw
    end

    # Step 3: plot value (all categorical, no continuous scaling)
    plot_val = FLOW_LAW_PLOT_MAP_IG[label]

    # Step 4: phase indices (blueschist, greenschist, eclogite, melt)
    indices = if label == "eclogite"
        (0.0, 0.0, 1.0, 0.0)
    elseif label == "greenschist"
        (0.0, 1.0, 0.0, 0.0)
    elseif label == ">15% melt"
        (0.0, 0.0, 0.0, 1.0)
    elseif label == "blueschist"
        (1.0, 0.0, 0.0, 0.0)
    else  # blueschist fallback
        (1.0, 0.0, 0.0, 0.0)
    end

    return label, plot_val, indices
end

# ── Flow law classification for "mtl" (mantle) database ──

# Phase abbreviations in the MAGEMin "mtl" database (Holland et al., 2013):
#   ol = olivine, wad = wadsleyite, ring = ringwoodite,
#   mpv = MgSi-perovskite (bridgmanite), fp = ferropericlase
const FLOW_LAW_PLOT_MAP_MTL = Dict(
    "olivine" => 1.0, "ring/wad" => 2.0, "bridgmanite+fp" => 3.0
)

"""Classify the flow law regime for the mantle database from weight fractions.
Returns (label, plot_value, phase_indices)."""
function classify_flow_law_mtl(ol_wt, ring_wad_wt, bm_fp_wt)
    total = ol_wt + ring_wad_wt + bm_fp_wt
    if total == 0.0
        return "bridgmanite+fp", 3.0, (0.0, 0.0, 0.0, 0.0)
    end

    label = if bm_fp_wt >= ol_wt && bm_fp_wt >= ring_wad_wt
        "bridgmanite+fp"
    elseif ring_wad_wt >= ol_wt
        "ring/wad"
    else
        "olivine"
    end

    plot_val = FLOW_LAW_PLOT_MAP_MTL[label]
    # Phase indices: (olivine, ring/wad, bridgmanite+fp, unused)
    indices = if label == "olivine"
        (1.0, 0.0, 0.0, 0.0)
    elseif label == "ring/wad"
        (0.0, 1.0, 0.0, 0.0)
    else
        (0.0, 0.0, 1.0, 0.0)
    end

    return label, plot_val, indices
end

"""Convert a phase-fraction tuple into a single dominant phase index (1-based).
For mtl, defaults to bridgmanite+fp (index 3). For ig, defaults to blueschist (index 1)."""
function dominant_index_from_tuple(indices, database, label)
    # Find the index of the maximum fraction; use only first 3 for mtl (4th is unused)
    n = database == "mtl" ? 3 : 4
    # If all zero (no stable phases), default to bridgmanite+fp for mtl
    if all(indices[1:n] .== 0.0)
        return database == "mtl" ? 3.0 : 1.0
    end
    idx = argmax(indices[1:n])
    return Float64(idx)
end

function compute_phase_fractions!(df, database)
    n = nrow(df)

    # Pre-allocate result columns.
    dominant_phase        = Vector{String}(undef, n)
    second_dominant_phase = Vector{String}(undef, n)
    flow_law_col          = Vector{String}(undef, n)
    flow_law_plot_col     = Vector{Float64}(undef, n)
    dominant_phase_idx    = zeros(n)
    liq_frac_vol          = zeros(n)
    ep_frac_vol           = zeros(n)

    for i in 1:n
        ph, pf, pt, pw = df.ph[i], df.ph_frac_vol[i], df.ph_type[i], df.ph_frac_wt[i]

        # No stable phases at this P-T point — keep zero/default values.
        if isempty(ph)
            dominant_phase[i] = database == "mtl" ? "bridgmanite+fp" : "blueschist"
            second_dominant_phase[i] = database == "mtl" ? "bridgmanite+fp" : "blueschist"
            flow_law_col[i] = database == "mtl" ? "bridgmanite+fp" : "blueschist"
            flow_law_plot_col[i] = database == "mtl" ? 3.0 : 1.0
            dominant_phase_idx[i] = database == "mtl" ? 3.0 : 1.0
            continue
        end

        # ── 1. Separate solution phases from pure phases ──
        ss_mask  = pt .== 1
        ss_names = any(ss_mask) ? ph[ss_mask] : String[]
        ss_fracs = any(ss_mask) ? pf[ss_mask] : Float64[]
        sw       = any(ss_mask) ? pw[ss_mask] : Float64[]

        # ── 2. Dominant and second-dominant phases ──
        dominant_phase[i], second_dominant_phase[i] = find_dominant_phases(ss_names, ss_fracs)

        # ── 3. Compute weight fractions and classify flow law ──
        if database == "ig"
            # Glaucophane weight: amp instances where compVariables[3] > 0.5
            amp_idxs    = findall(==("amp"), ph)
            amp_ss_idxs = findall(==("amp"), ss_names)
            gl_wt = 0.0
            for k in amp_ss_idxs
                cv = df.SS_vec[i][k].compVariables
                if length(cv) >= 3 && cv[3] > 0.5
                    gl_wt += sw[k]
                end
            end
            amp_only_wt = isempty(amp_idxs) ? 0.0 : sum(pw[k] for k in amp_idxs) - gl_wt

            # Na-cpx weight: cpx where compVariables[4] > 0.5
            cpx_ss_idx = findfirst(==("cpx"), ss_names)
            nacpx_wt = if cpx_ss_idx !== nothing
                cv = df.SS_vec[i][cpx_ss_idx].compVariables
                (length(cv) >= 4 && cv[4] > 0.5) ? sw[cpx_ss_idx] : 0.0
            else
                0.0
            end

            blueschist_wt  = gl_wt + get_frac(findfirst(==("ep"), ph), pw) + get_frac(findfirst(==("chl"), ph), pw)
            liq_wt         = get_frac(findfirst(==("liq"), ph), pw)
            greenschist_wt = amp_only_wt + get_frac(findfirst(==("mu"), ph), pw) + get_frac(findfirst(==("fsp"), ph), pw) +
                             get_frac(findfirst(==("chl"), ph), pw) + get_frac(findfirst(==("ep"), ph), pw) +
                             (liq_wt <= 0.15 ? liq_wt : 0.0)
            eclogite_wt    = get_frac(findfirst(==("g"), ph), pw) + nacpx_wt

            label, plot_val, indices = classify_flow_law_ig(blueschist_wt, greenschist_wt, eclogite_wt, liq_wt, df.T_C[i])

        elseif database == "mtl"
            # Mantle phase weight fractions (Holland et al., 2013 phase names)
            ol_wt      = get_frac(findfirst(==("ol"), ph), pw)
            wad_wt     = get_frac(findfirst(==("wad"), ph), pw)
            ring_wt    = get_frac(findfirst(==("ring"), ph), pw)
            mpv_wt     = get_frac(findfirst(==("mpv"), ph), pw)
            fp_wt      = get_frac(findfirst(==("fp"), ph), pw)
            ring_wad_wt = wad_wt + ring_wt
            bm_fp_wt    = mpv_wt + fp_wt

            label, plot_val, indices = classify_flow_law_mtl(ol_wt, ring_wad_wt, bm_fp_wt)
        end

        flow_law_col[i]      = label
        flow_law_plot_col[i] = plot_val
        dominant_phase_idx[i] = dominant_index_from_tuple(indices, database, label)

        # ── 4. Store volume fractions ──
        liq_frac_vol[i] = get_frac(findfirst(==("liq"), ph), pf)
        ep_frac_vol[i]  = get_frac(findfirst(==("ep"),  ph), pf)
    end

    # Attach results as DataFrame columns.
    df.dominant_phase        = dominant_phase
    df.second_dominant_phase = second_dominant_phase
    df.flow_law              = flow_law_col
    df.flow_law_plot_value   = flow_law_plot_col
    df[!, :dominant_phase_index] = dominant_phase_idx
    df.liq_frac_vol = liq_frac_vol
    df.ep_frac_vol  = ep_frac_vol
    return df
end

function filter_and_scale(df, lookuptabletype, resolution, usescaling, database)

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
        # Select and rename columns from the MAGEMin interface output.
        # The first nine entries are standard bulk thermodynamic properties.
        # The last six (dominant_phase … nacpx_frac_vol) are phase-fraction / dominant-phase
        # columns added by Copilot; they are carried through the pipeline for visualisation
        # but excluded from the PerpleX .txt output via the names filter applied later.
        TidierData.@select(P = var"P_kbar", T = var"T_C", density = var"rho", Cp = var"s_cp", 
            alpha = var"alpha", enthalpy = var"enthalpy", Vp = var"Vp", Vs = var"Vs", bulk_S_wt = var"bulk_S_wt",
            dominant_phase = var"dominant_phase", second_dominant_phase = var"second_dominant_phase", flow_law = var"flow_law", liq_frac_vol = var"liq_frac_vol",
            flow_law_plot_value = var"flow_law_plot_value", dominant_phase_index = var"dominant_phase_index", ep_frac_vol = var"ep_frac_vol")
        # Convert pressure from kbar to bar and temperature from celsius to kelvin. 
        @mutate(P = 1e3 * P, T = 273.15 + T)
        # Drop the 'phase' column and reorder the columns for PerpleX format.
        # Reorder columns: standard thermodynamic columns first (required for PerpleX output),
        # then phase-fraction columns last (to be stripped by the names filter before writing .txt).
        TidierData.@select(T, P, density, alpha, Cp, Vp, Vs, enthalpy, bulk_S_wt, dominant_phase, second_dominant_phase, flow_law, flow_law_plot_value, dominant_phase_index, liq_frac_vol, ep_frac_vol)
        @arrange(P, T)
        # Rename the column for later use.
        @rename(H20minebound = bulk_S_wt)
    end
    # Convert from array[array] to array
    for (col, idx) in [(:alpha, 1), (:Cp, 1), (:enthalpy, 1)]
        filtereddf[!, col] = getindex.(filtereddf[!, col], idx)
    end
    # H2O bound water: extract last element only when the oxide vector is long enough (i.e. H2O is present);
    # otherwise fill with zeros.
    filtereddf[!, :H20minebound] = map(filtereddf[!, :H20minebound]) do v
        length(v) >= 11 ? v[11] : 0.0
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
if usescaling == true
    scalingfactor = ones(resolution)
    for i in range(1, resolution)
        scalingfactor[i] = kinetic_delay_scaling(filtereddf.T[i])
    end
    filtereddf.scalingfactor = repeat(scalingfactor, outer = resolution)
    # filtereddf.density = porosityscaling(filtereddf.density, filtereddf.P, resolution, "oceanic")
end

filtereddf = @chain filtereddf begin
    # Convert subzero (=unphysical) heat expansivity and specific heat values to zero.
    @mutate(alphafilt = case_when(alpha < 0.0 => 0.0, alpha > 0.75e-4 => 0.75e-4, alpha > 0.0 => alpha))
    @mutate(Cp = case_when(Cp < 750.0 => 750.0, Cp > 1500.0 => 1500.0, Cp > 0.0 => Cp))
end

# Apply density and Cp scaling only for igneous ("ig") database.
if database == "ig"
    filtereddf = @chain filtereddf begin
        # Apply kinetic delay scaling to density values (only where density > 3000).
        @mutate(density = case_when(density > 3000 => 3000 * (1.0 - scalingfactor) + density .* scalingfactor, density < 3000 => density))
        # Filter outlier density values for better stability of a geodynamic model.
        @mutate(density = case_when(density < 2600 => 2600, density > 4600 => 4600, density > 0 => density))
    end
end
    return filtereddf
end

function write_output(filtereddf, df, lookuptabletype, resolution, outputtablefilepath, outputtablefilename)

open(outputtablefilepath, "w") do io
    # Ensure index column exists also in file-based input mode.
    if !(:dominant_phase_index in propertynames(filtereddf))
        filtereddf[!, :dominant_phase_index] = zeros(nrow(filtereddf))
    end

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
        * "\t\t10\n"                                                # Number of columns
        * "T(K)\tP(bar)\ts,J/K/kg\trho,kg/m3\talpha,1/K\tcp,J/K/kg\tvp,km/s\tvs,km/s\th,J/kg\tdominant_phase_index,1\n") # Column names and units
    # Add dummy entropy column (copy of Cp).
    filtereddf[!, :s] = filtereddf[!, :Cp]
    # Write the actual phase diagram data (exclude non-numeric columns).
    output_df = filtereddf[!, [:T, :P, :s, :density, :alphafilt, :Cp, :Vp, :Vs, :enthalpy, :dominant_phase_index]]
    writedlm(io, eachrow(output_df), '\t')
end

# Write the full phase assemblage (all phases + volume fractions) at the first (minimum) temperature
# to a long-format CSV. In the original df layout (Prange = repeat(P1..Pn, outer=n),
# Trange = repeat(T1..Tn, inner=n)), rows 1..resolution have T=T_min and P varying P_min→P_max —
# i.e. the first pressure column.
if lookuptabletype == "interface"
    first_P_indices = 1:resolution
    first_P_rows = NamedTuple{(:T_C, :P_kbar, :phase, :frac_vol), Tuple{Float64, Float64, String, Float64}}[]
    for idx in first_P_indices
        T_val = df.T_C[idx]
        P_val = df.P_kbar[idx]
        ph_list = df.ph[idx]
        pf_list = df.ph_frac_vol[idx]
        if isempty(ph_list)
            push!(first_P_rows, (T_C = T_val, P_kbar = P_val, phase = "none", frac_vol = 0.0))
        else
            for (ph, pf) in zip(ph_list, pf_list)
                push!(first_P_rows, (T_C = T_val, P_kbar = P_val, phase = ph, frac_vol = pf))
            end
        end
    end
    CSV.write("first_pressure_column_phases.csv", DataFrame(first_P_rows))
end
end

function plot_results(filtereddf, resolution, database)
# Arrays for plotting in celcius degrees and gigapascals.
plottingT = filtereddf.T[1:resolution] .- 273.15
plottingP = filtereddf.P[1:resolution:resolution^2] ./ 1e4

# CairoMakie heatmaps require a numeric matrix; build phase-name → integer mappings.
function phase_to_numeric(phase_col)
    unique_ph = sort(unique(phase_col))
    n = length(unique_ph)
    mapping = Dict(p => Float64(i) for (i, p) in enumerate(unique_ph))
    return [mapping[p] for p in phase_col], unique_ph, n
end

dominant_phase_num, unique_phases, n_phases = phase_to_numeric(filtereddf.dominant_phase)
second_dominant_phase_num, unique_phases2, n_phases2 = phase_to_numeric(filtereddf.second_dominant_phase)

# Continuous plotting scale for flow_law, database-dependent.
if database == "ig"
    flow_law_ticks = ([1.0, 2.0, 3.0, 4.0], ["1: blueschist", "2: eclogite", "3: greenschist", "4: >15% melt"])
    flowlaw_cmap = cgrad([:royalblue3, :firebrick3, :darkolivegreen3, :darkorange2], 4, categorical = true)
    flowlaw_colorrange = (0.5, 4.5)
elseif database == "mtl"
    flow_law_ticks = ([1.0, 2.0, 3.0], ["1: olivine", "2: ring/wad", "3: bridgmanite+fp"])
    flowlaw_cmap = cgrad([:forestgreen, :darkorange2, :mediumpurple3], 3, categorical = true)
    flowlaw_colorrange = (0.5, 3.5)
end

# Use LaTeX font for plotting.
with_theme(theme_latexfonts()) do
    f = Figure(backgroundcolor = :transparent, size = (3000, 1400))
    ga = f[1, 1] = GridLayout()
    # Shared axis keyword arguments.
    ax_base = (; xlabel = L"temperature ($\degree$C)", tellwidth = false, width = 600,
        xlabelsize = 20, xticksize = 5, xticklabelsize = 18, xticklabelpad = 0.5)
    ax_with_y = (; ax_base..., ylabel = "pressure (GPa)", ylabelsize = 20, yticksize = 5,
        yticklabelsize = 18, yticklabelpad = 0.5)

    # Row 1 axes: density, Cp, bound water
    axleft   = Axis(ga[1, 1]; ax_with_y...)
    axcenter = Axis(ga[1, 2]; ax_base...)
    axright  = Axis(ga[1, 3]; ax_base...)
    # Row 2 axes (gridlayout row 3; row 2 reserved for colorbars)
    axphase   = Axis(ga[3, 1]; ax_with_y...)
    axliq     = Axis(ga[3, 2]; ax_base...)
    axflowlaw = Axis(ga[3, 3]; ax_base...)
    axphase2  = Axis(ga[3, 4]; ax_base...)

    all_axes = [axleft, axcenter, axright, axphase, axliq, axflowlaw, axphase2]
    linkyaxes!(all_axes...)
    linkxaxes!(all_axes...)

    # Helper to reshape a flat vector into the P-T grid matrix.
    togrid(v) = reshape(v, resolution, resolution)
    # Helper to create a heatmap on plottingT × plottingP.
    makehm!(ax, data; kw...) = CairoMakie.heatmap!(ax, plottingT, plottingP, togrid(data); kw...)

    # Row 1: density, Cp, bound water
    hmdensity = makehm!(axleft, filtereddf.density, colormap = Reverse(:imola10))
    hmCp      = makehm!(axcenter, filtereddf.Cp, colormap = Reverse(:glasgow))
    hmH2Owt   = makehm!(axright, filtereddf.H20minebound .* 100, colormap = Reverse(:nuuk))
    # Row 2: phase maps, melt fraction, flow law
    phase_cmap = cgrad(:tab20, n_phases, categorical = true)
    hmphase   = makehm!(axphase, dominant_phase_num, colormap = phase_cmap, colorrange = (0.5, n_phases + 0.5))
    hmliq     = makehm!(axliq, filtereddf.liq_frac_vol, colormap = :thermal)
    hmflowlaw = makehm!(axflowlaw, filtereddf.flow_law_plot_value, colormap = flowlaw_cmap, colorrange = flowlaw_colorrange)
    phase2_cmap = cgrad(:tab20, n_phases2, categorical = true)
    hmphase2  = makehm!(axphase2, second_dominant_phase_num, colormap = phase2_cmap, colorrange = (0.5, n_phases2 + 0.5))

    # Shared colorbar keyword arguments.
    cb_base = (; vertical = false, tellwidth = false, width = 350, spinewidth = 0.1,
        labelpadding = 1.5, labelsize = 16, ticklabelpad = 0.5, ticklabelsize = 14, ticksize = 3.5)

    # Row 1 colorbars
    Colorbar(ga[0, 1][1, 1], hmdensity; cb_base..., label = L"density (kg/m$^3$)",
        ticks = [round(minimum(filtereddf.density), sigdigits = 4), 2500, 3000, 3500, 4000, 4500,
        floor(maximum(filtereddf.density), sigdigits = 4)],
        minorticks = [2250, 2750, 3250, 3750, 4250], minorticksvisible = true, minorticksize = 2.5, minortickwidth = 0.75)
    Colorbar(ga[0, 2][1, 1], hmCp; cb_base..., label = "Cp", ticks = [0.0, 0.25, 0.5, 0.75, 1.0])
    Colorbar(ga[0, 3][1, 1], hmH2Owt; cb_base..., label = "bound water (wt. %)", ticks = [0, 0.5, 1.0, 1.5, 2.0, 2.5])
    # Row 2 colorbars
    Colorbar(ga[2, 1][1, 1], hmphase; cb_base..., label = "dominant phase", ticklabelsize = 10,
        ticks = (collect(1.0:n_phases), unique_phases))
    Colorbar(ga[2, 2][1, 1], hmliq; cb_base..., label = "liq (melt) vol. fraction")
    Colorbar(ga[2, 3][1, 1], hmflowlaw; cb_base..., label = "flow law", ticklabelsize = 12, ticks = flow_law_ticks)
    Colorbar(ga[2, 4][1, 1], hmphase2; cb_base..., label = "second dominant phase", ticklabelsize = 10,
        ticks = (collect(1.0:n_phases2), unique_phases2))

    # Push all heatmaps to z = -100 so ticks/contours render in front.
    for hm in [hmdensity, hmCp, hmH2Owt, hmphase, hmliq, hmflowlaw, hmphase2]
        CairoMakie.translate!(hm, 0, 0, -100)
    end

    # Beautification
    CairoMakie.ylims!(axright, low = 0)
    CairoMakie.xlims!(axright, low = 0)
    for ax in all_axes
        hidespines!(ax)
    end
    for ax in [axcenter, axright, axliq, axflowlaw, axphase2]
        hideydecorations!(ax, grid = false)
    end
    colgap!(ga, 10)
    rowgap!(ga, 10)
    display(f)
    save("test.png", f)
end
end

# ---- Main execution ----
df, resolution = generate_dataframe(lookuptabletype, database)
if lookuptabletype == "interface"
    compute_phase_fractions!(df, database)
    CSV.write("magemin_output.csv", df)
end
filtereddf = filter_and_scale(df, lookuptabletype, resolution, usescaling, database)
write_output(filtereddf, df, lookuptabletype, resolution, outputtablefilepath, outputtablefilename)
plot_results(filtereddf, resolution, database)
