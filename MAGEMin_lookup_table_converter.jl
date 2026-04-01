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
usescaling = true

if lookuptabletype == "file"
    # Filename of a lookup table to process, let us assume it is located in a same folder as this processor script.
    inputtablefolder = "."
    inputtablefilename = "nmorb-gabbro-reso263k-p15gpa-t2000c-h20wt05_ig.csv"
    inputtablefilepath = joinpath(inputtablefolder, inputtablefilename)
    # Read in a table to initialize a dataframe.
    df = DataFrame(CSV.File(inputtablefilepath))
elseif lookuptabletype == "interface"
    n = 512
    resolution = n
    Prange = repeat(range(0.01, 200, n), outer = n)
    Trange = repeat(range(0, 1500, n), inner = n)
    db = Initialize_MAGEMin("ig", verbose = false)
    Xoxides = ["SiO2";    "Al2O3";   "CaO";     "MgO";    "FeO";    "K2O";    "Na2O";   "TiO2";   "O";     "Cr2O3";  "H2O"] # Stoner et al. composition
    X =       [50.992811; 15.304009; 10.922183; 7.848858; 9.913537; 0.135232; 2.867050; 1.541733; 0.176099; 0.048487; 0.25] # gabbro
    # X =       [50.376296; 15.118980; 10.790131; 7.753963; 9.793680; 0.133597; 2.832387; 1.523093; 0.17397; 0.047901; 1.456] # basalt
    # Xoxides = ["SiO2"; "Al2O3"; "CaO"; "MgO"; "FeO"; "K2O"; "Na2O"; "TiO2"; "O"; "MnO"; "H2O"] # garnet-migmatite
    # X = [69.66; 13.76; 1.77; 1.73; 4.32; 2.61; 2.41; 0.8; 0.03; 0.07; 2.82] # garnet-migmatite
    # X = [64.84; 14.99; 3.50; 2.41; 4.90; 2.72; 3.18; 0.62; 0.0; 0.1; 2.71] # average upper continental crust
    sys_in = "wt"
    df = DataFrame(multi_point_minimization(Prange, Trange, db, X = X, Xoxides = Xoxides, sys_in = sys_in)) # For user-specified bulk composition
    Finalize_MAGEMin(db)

    # Compute dominant phase of a P-T-point and several phase fractions
    dominant_phase        = String[]     # Acronym of the volumetrically dominant solid solution phase
    second_dominant_phase = String[]     # Acronym of the second-most dominant solution phase
    flow_law              = String[]
    flow_law_plot_value   = Float64[]
    dominant_phase_indices = [Float64[], Float64[], Float64[], Float64[]]  # blueschist, greenschist, eclogite, melt
    dominant_phase_index_1, dominant_phase_index_2, dominant_phase_index_3, dominant_phase_index_4 = dominant_phase_indices
    
    liq_frac_vol     = Float64[]    # (Silicate) melt, i.e. "liq" volume fraction
    amp_frac_vol     = Float64[]    # Amphibole "amp" volume fraction
    amp_is_gl        = Bool[]       # true when amp compVariables[3] > 0.5, i.e. glaucophane
    gl_frac_vol      = Float64[]    # Glaucophane "gl" volume fraction: amp fraction when amp_is_gl, otherwise 0.0
    g_frac_vol       = Float64[]
    g_nacpx_frac_vol = Float64[]    # Garnet "g" + Na-rich cpx "Na-cpx" combined volume fraction
    cpx_frac_vol     = Float64[]    # Clinopyroxene "cpx" volume fraction, regardless of Na content
    nacpx_frac_vol   = Float64[]    # Na-rich cpx fraction, cpx with compVariables[4] > 0.5
    ep_frac_vol      = Float64[]
    blsc_frac_vol    = Float64[]

    for i in 1:nrow(df)
        ph, pf, pt = df.ph[i], df.ph_frac_vol[i], df.ph_type[i]
        pw = df.ph_frac_wt[i]  # weight fractions parallel to ph/pf
        # Check if no stable phases returned at this P-T point due to numerical failure or boundary condition
        if isempty(ph)
            # Push safe defaults so all vectors stay the same length as nrow(df); if a dominant phase, use "none", if phase fraction, use 0.0
            push!(dominant_phase, "none")
            push!(second_dominant_phase, "none")
            push!(flow_law, "none")
            push!(flow_law_plot_value, 1.0)
            for idx_arr in dominant_phase_indices
                push!(idx_arr, 0.0)
            end
            push!(liq_frac_vol, 0.0)
            push!(amp_frac_vol, 0.0)
            push!(amp_is_gl, false)
            push!(gl_frac_vol, 0.0)
            push!(g_frac_vol, 0.0)
            push!(cpx_frac_vol, 0.0)
            push!(nacpx_frac_vol, 0.0)
            push!(ep_frac_vol, 0.0)
            push!(blsc_frac_vol, 0.0)
        else
            # Separate solution phases (ph_type == 1) from pure phases (ph_type == 0)
            # Only solution phases participate in the dominant-phase logic; pure phases such as quartz or rutile are excluded because they don't substitute into one another
            ss_mask  = pt .== 1
            ss_names = any(ss_mask) ? ph[ss_mask]  : String[]   # names of stable solution phases only
            ss_fracs = any(ss_mask) ? pf[ss_mask]  : Float64[]  # their corresponding volume fractions
            sw       = any(ss_mask) ? pw[ss_mask]  : Float64[]  # their corresponding weight fractions

            # Pick up melt and phases fractions, returns 0.0 if the phases are not present in the stable assemblage
            idx_liq  = findfirst(==("liq"), ph)
            amp_idxs = findall(==("amp"), ph)      # all "amp" entries in ph (MAGEMin may produce n instances)
            idx_g    = findfirst(==("g"),   ph)
            idx_cpx  = findfirst(==("cpx"), ph)
            idx_ep   = findfirst(==("ep"),  ph)
            idx_chl  = findfirst(==("chl"), ph)
            idx_mu   = findfirst(==("mu"),  ph)
            idx_fsp  = findfirst(==("fsp"), ph)

            # Glaucophane "gl" identification + adjust amp fractions in single pass
            # amp_ss_idxs indexes into ss_names (solution phases only), which is the same ordering as SS_vec uses.
            # MAGEMin may produce n "amp" solution phases; loop over all of them.
            amp_ss_idxs = findall(==("amp"), ss_names)
            gl_vol = 0.0  # glaucophane volume fraction (for dominant_phase logic and output vectors)
            gl_wt  = 0.0  # glaucophane weight fraction (for flow_law logic)
            ss_fracs_adj = copy(ss_fracs)
            
            for k in amp_ss_idxs
                cv_amp = df.SS_vec[i][k].compVariables
                if length(cv_amp) >= 3 && cv_amp[3] > 0.5
                    gl_vol += ss_fracs[k]
                    gl_wt  += sw[k]
                    ss_fracs_adj[k] = 0.0  # subtract glaucophane from dominant phase calc
                else
                    ss_fracs_adj[k] = ss_fracs[k]  # keep non-glaucophane amp
                end
            end
            # The dominant phase is the solution phase with the largest adjusted volume fraction at this P-T point
            push!(dominant_phase, isempty(ss_names) ? "none" : ss_names[argmax(ss_fracs_adj)])
            # The second dominant phase is the phase with the second-largest adjusted volume fraction
            if isempty(ss_names) || length(ss_names) < 2
                push!(second_dominant_phase, "none")
            else
                ss_fracs_adj2 = copy(ss_fracs_adj)
                ss_fracs_adj2[argmax(ss_fracs_adj)] = -Inf
                push!(second_dominant_phase, ss_names[argmax(ss_fracs_adj2)])
            end

            # Na-cpx identification follows the MAGEMin ig-database solvus naming rule; cpx is called "Na-cpx" when compVariables[4] > 0.5, i.e jadeite content >50%
            # cpx_ss_idx indexes into ss_names (solution phases only), which is the same ordering as SS_vec uses
            cpx_ss_idx = findfirst(==("cpx"), ss_names)
            if cpx_ss_idx !== nothing
                # Retrieve the internal compositional variables of the cpx at this P-T point
                cv = df.SS_vec[i][cpx_ss_idx].compVariables
                nacpx_vol = (length(cv) >= 4 && cv[4] > 0.5) ? ss_fracs[cpx_ss_idx] : 0.0
                nacpx_wt  = (length(cv) >= 4 && cv[4] > 0.5) ? sw[cpx_ss_idx]       : 0.0
            else
                nacpx_vol = 0.0  # "cpx" not present in the stable assemblage at this point
                nacpx_wt  = 0.0
            end
            
            # Extract all phase fractions (volume and weight) at once
            get_frac = (idx, fracs) -> idx === nothing ? 0.0 : fracs[idx]
            
            ep_vol  = get_frac(idx_ep,  pf)
            chl_vol = get_frac(idx_chl, pf)
            g_vol   = get_frac(idx_g,   pf)
            liq_vol = get_frac(idx_liq, pf)
            
            ep_wt   = get_frac(idx_ep,  pw)
            chl_wt  = get_frac(idx_chl, pw)
            mu_wt   = get_frac(idx_mu,  pw)
            fsp_wt  = get_frac(idx_fsp, pw)
            g_wt    = get_frac(idx_g,   pw)
            liq_wt  = get_frac(idx_liq, pw)
            
            amp_only_vol = isempty(amp_idxs) ? 0.0 : sum(pf[k] for k in amp_idxs) - gl_vol
            amp_only_wt  = isempty(amp_idxs) ? 0.0 : sum(pw[k] for k in amp_idxs) - gl_wt
            # flow_law: dominant rock type determined by weight fractions.
            # blueschist  = gl + ep + chl weight fractions
            # greenschist = non-gl amp + mu + fsp + melt (when liq ≤ 0.15) weight fractions
            # eclogite    = Na-cpx + g weight fractions
            flow_law_scaling = kinetic_delay_scaling(273.15 + df.T_C[i])
            blueschist_wt   = gl_wt + ep_wt + chl_wt
            amp_only_wt     = isempty(amp_idxs) ? 0.0 : sum(pw[k] for k in amp_idxs) - gl_wt
            greenschist_wt  = amp_only_wt + mu_wt + fsp_wt + chl_wt + ep_wt + (liq_wt <= 0.15 ? liq_wt : 0.0)
            eclogite_wt     = g_wt + nacpx_wt  # Na-cpx + g combined weight fraction

            raw_flow_law_label = if liq_wt > 0.15
                ">15% melt"
            elseif eclogite_wt > 0.5 || eclogite_wt > blueschist_wt && eclogite_wt > greenschist_wt
                "eclogite"
            elseif (blueschist_wt + greenschist_wt) == 0.0
                "none"
            elseif blueschist_wt >= greenschist_wt
                "blueschist"
            else
                "greenschist"
            end

            flow_law_label = if raw_flow_law_label == "eclogite" && flow_law_scaling == 0.0
                "blueschist"
            elseif raw_flow_law_label == "eclogite" && flow_law_scaling == 1.0
                "eclogite"
            elseif raw_flow_law_label == "eclogite"
                "blueschist-eclogite transition"
            else
                raw_flow_law_label
            end
            push!(flow_law, flow_law_label)

            # Map flow_law_label to plot value; raw eclogite gets scaled blend (2.0-3.0)
            flow_law_plot = if raw_flow_law_label == "eclogite"
                2.0 + flow_law_scaling
            else
                Dict("none" => 1.0, "blueschist" => 2.0, "greenschist" => 4.0, ">15% melt" => 5.0, "blueschist-eclogite transition" => 2.0 + 0.5)[flow_law_label]
            end
            push!(flow_law_plot_value, flow_law_plot)

            # Dominant phase indices from flow_law state: 
            # 1=blueschist, 2=greenschist, 3=eclogite, 4=>15% melt
            # Greenschist & melt are binary 0/1; blueschist-eclogite is continuous (1-scaling, scaling)
            idx_values = (0.0, 0.0, 0.0, 0.0)  # (blueschist, greenschist, eclogite, melt)
            if raw_flow_law_label == "eclogite"
                idx_values = (1.0 - flow_law_scaling, 0.0, flow_law_scaling, 0.0)
            elseif flow_law_label == "greenschist"
                idx_values = (0.0, 1.0, 0.0, 0.0)
            elseif flow_law_label == ">15% melt"
                idx_values = (0.0, 0.0, 0.0, 1.0)
            elseif flow_law_label == "blueschist"
                idx_values = (1.0, 0.0, 0.0, 0.0)
            elseif flow_law_label == "eclogite"
                idx_values = (0.0, 0.0, 1.0, 0.0)
            end
            for (idx_arr, val) in zip(dominant_phase_indices, idx_values)
                push!(idx_arr, val)
            end

            push!(g_frac_vol, g_vol)
            push!(cpx_frac_vol, idx_cpx === nothing ? 0.0 : pf[idx_cpx])
            push!(liq_frac_vol, liq_vol)
            push!(amp_frac_vol, amp_only_vol)  # amp with gl (glaucophane) fraction subtracted
            push!(gl_frac_vol, gl_vol)
            push!(ep_frac_vol, ep_vol)
            push!(nacpx_frac_vol, nacpx_vol) # Na-rich cpx volume fraction
            push!(g_nacpx_frac_vol, g_frac_vol[i] + nacpx_vol) # Eclogite = garnet + Na-cpx

        end
    end

    # Rename certain dominant_phase labels, now all fraction vectors are fully populated and the conditions are well-defined
    for i in eachindex(dominant_phase)
        # Amphibole dominant and glaucophane-bearing (compVariables[3] > 0.5) → "gl".
        # Follows the MAGEMin ig-database solvus rule from julia/name_solvus.jl.
        if dominant_phase[i] == "amp" && gl_frac_vol[i] > 0.0
            dominant_phase[i] = "gl"
        elseif dominant_phase[i] == "cpx" && nacpx_frac_vol[i] > 0.0
            dominant_phase[i] = "Na-cpx"
        end
        # Apply the same renaming rules to second_dominant_phase
        if second_dominant_phase[i] == "amp" && gl_frac_vol[i] > 0.0
            second_dominant_phase[i] = "gl"
        elseif second_dominant_phase[i] == "cpx" && nacpx_frac_vol[i] > 0.0
            second_dominant_phase[i] = "Na-cpx"
        end
    end
    # Attach all per-point vectors as new columns on df so they are accessible
    # in the @chain filtereddf pipeline below and in the visualisation section.
    df.dominant_phase        = dominant_phase         # relabeled dominant solution-phase string
    df.second_dominant_phase = second_dominant_phase  # relabeled second-most dominant solution-phase string
    df.flow_law              = flow_law          # dominant rock type: "gl+ep" (blueschist/greenschist) or "amp" (amphibolite)
    df.flow_law_plot_value   = flow_law_plot_value
    for (i, col_name) in enumerate(["dominant_phase_index_1", "dominant_phase_index_2", "dominant_phase_index_3", "dominant_phase_index_4"])
        df[!, Symbol(col_name)] = dominant_phase_indices[i]
    end
    df.liq_frac_vol      = liq_frac_vol      # melt volume fraction
    df.amp_frac_vol      = amp_frac_vol      # amphibole volume fraction
    df.gl_frac_vol       = gl_frac_vol       # glaucophane volume fraction (amp fraction where compVariables[3] > 0.5)
    df.ep_frac_vol       = ep_frac_vol       # epidote volume fraction
    df.g_nacpx_frac_vol  = g_nacpx_frac_vol  # garnet + Na-cpx combined volume fraction
    df.cpx_frac_vol      = cpx_frac_vol      # raw clinopyroxene volume fraction
    df.nacpx_frac_vol    = nacpx_frac_vol    # Na-rich clinopyroxene volume fraction
end

# Save the complete, unprocessed MAGEMin output DataFrame to CSV file
if lookuptabletype == "interface"
    CSV.write("magemin_output.csv", df)
end

# Regardless of the lookup table creation method, save it as a .txt file which is the only
# format that ASPECT reads.
outputtablefolder = "."
outputtablefilename = "3200kgm3_mtl.txt"
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
        # Select and rename columns from the MAGEMin interface output.
        # The first nine entries are standard bulk thermodynamic properties.
        # The last six (dominant_phase … nacpx_frac_vol) are phase-fraction / dominant-phase
        # columns added by Copilot; they are carried through the pipeline for visualisation
        # but excluded from the PerpleX .txt output via the names filter applied later.
        TidierData.@select(P = var"P_kbar", T = var"T_C", density = var"rho", Cp = var"s_cp", 
            alpha = var"alpha", enthalpy = var"enthalpy", Vp = var"Vp", Vs = var"Vs", bulk_S_wt = var"bulk_S_wt",
            dominant_phase = var"dominant_phase", second_dominant_phase = var"second_dominant_phase", flow_law = var"flow_law", liq_frac_vol = var"liq_frac_vol", amp_frac_vol = var"amp_frac_vol",
            flow_law_plot_value = var"flow_law_plot_value", dominant_phase_index_1 = var"dominant_phase_index_1", dominant_phase_index_2 = var"dominant_phase_index_2", dominant_phase_index_3 = var"dominant_phase_index_3", dominant_phase_index_4 = var"dominant_phase_index_4", gl_frac_vol = var"gl_frac_vol", ep_frac_vol = var"ep_frac_vol", g_nacpx_frac_vol = var"g_nacpx_frac_vol", cpx_frac_vol = var"cpx_frac_vol", nacpx_frac_vol = var"nacpx_frac_vol")
        # Convert pressure from kbar to bar and temperature from celsius to kelvin. 
        @mutate(P = 1e3 * P, T = 273.15 + T)
        # Drop the 'phase' column and reorder the columns for PerpleX format.
        # Reorder columns: standard thermodynamic columns first (required for PerpleX output),
        # then phase-fraction columns last (to be stripped by the names filter before writing .txt).
        TidierData.@select(T, P, density, alpha, Cp, Vp, Vs, enthalpy, bulk_S_wt, dominant_phase, second_dominant_phase, flow_law, flow_law_plot_value, dominant_phase_index_1, dominant_phase_index_2, dominant_phase_index_3, dominant_phase_index_4, liq_frac_vol, amp_frac_vol, gl_frac_vol, ep_frac_vol, g_nacpx_frac_vol, cpx_frac_vol, nacpx_frac_vol)
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
    # Apply kinetic delay scaling to density values (only where density > 3000).
    @mutate(density = case_when(density > 3000 => 3000 * (1.0 - scalingfactor) + density .* scalingfactor, density < 3000 => density))
    # Filter outlier density values for better stability of a geodynamic model.
    @mutate(density = case_when(density < 2600 => 2600, density > 4600 => 4600, density > 0 => density))
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
    # Ensure index columns exist also in file-based input mode.
    for i in 1:4
        col_name = Symbol("dominant_phase_index_$i")
        if !(col_name in propertynames(filtereddf))
            filtereddf[!, col_name] = zeros(nrow(filtereddf))
        end
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
        * "\t\t12\n"                                                # Number of columns
        * "T(K)\tP(bar)\trho,kg/m3\talpha,1/K\tcp,J/K/kg\tvp,km/s\tvs,km/s\th,J/kg\tdominant_phase_index,1\tdominant_phase_index,2\tdominant_phase_index,3\tdominant_phase_index,4\n") # Column names and units
    # Write the actual phase diagram data (exclude non-numeric columns).
    # Drop visualisation-only columns before writing the PerpleX-format .txt file.
    # ASPECT expects a purely numeric table with the eight standard thermodynamic properties;
    # the string column "dominant_phase" and all Float64 phase-fraction columns must be excluded.
    output_df = filtereddf[!, [:T, :P, :density, :alphafilt, :Cp, :Vp, :Vs, :enthalpy, :dominant_phase_index_1, :dominant_phase_index_2, :dominant_phase_index_3, :dominant_phase_index_4]]
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

# Arrays for plotting in celcius degrees and gigapascals.
plottingT = filtereddf.T[1:resolution] .- 273.15
plottingP = filtereddf.P[1:resolution:resolution^2] ./ 1e4

# CairoMakie heatmaps require a numeric matrix; dominant_phase is a String column.
# Build a phase-name → integer mapping so the heatmap can be rendered with a categorical
# colour scale. Sorting alphabetically ensures a reproducible colour assignment across runs,
# even if the stable mineral assemblage changes slightly at different P-T ranges.
unique_phases      = sort(unique(filtereddf.dominant_phase))  # all distinct dominant phase names, sorted
n_phases           = length(unique_phases)                     # total number of distinct phases present
phase_to_int       = Dict(p => Float64(i) for (i, p) in enumerate(unique_phases))  # name → Float64 integer map
dominant_phase_num = [phase_to_int[p] for p in filtereddf.dominant_phase]           # flat numeric vector for heatmap

# Integer mapping for second_dominant_phase; uses the same unique_phases set so colours are consistent.
unique_phases2      = sort(unique(filtereddf.second_dominant_phase))
n_phases2           = length(unique_phases2)
phase2_to_int       = Dict(p => Float64(i) for (i, p) in enumerate(unique_phases2))
second_dominant_phase_num = [phase2_to_int[p] for p in filtereddf.second_dominant_phase]

# Continuous plotting scale for flow_law.
# 1 = none, 2 = blueschist, 3 = eclogite, 4 = greenschist, 5 = >15% melt.
# Values between 2 and 3 visualize the temperature-scaled blueschist → eclogite transition.
flow_law_ticks = ([1.0, 2.0, 3.0, 4.0, 5.0], ["none", "blueschist", "eclogite", "greenschist", ">15% melt"])

# Use LaTeX font for plotting.
with_theme(theme_latexfonts()) do
    f = Figure(backgroundcolor = :transparent, size = (4000, 2100))
    ga = f[1, 1] = GridLayout()
    # Row 1 axes
    axleft = Axis(ga[1, 1], xlabel = L"temperature ($\degree$C)", ylabel = "pressure (GPa)", tellwidth = false, width = 600, 
        xlabelsize = 20, ylabelsize = 20, yticksize = 5, xticksize = 5, xticklabelsize = 18, yticklabelsize = 18, 
        xticklabelpad = 0.5, yticklabelpad = 0.5)
    axcenter = Axis(ga[1, 2], xlabel = L"temperature ($\degree$C)", tellwidth = false, width = 600,
        xlabelsize = 20, xticksize = 5, xticklabelsize = 18, xticklabelpad = 0.5)
    axright = Axis(ga[1, 3], tellwidth = false, width = 600, xticksize = 5, xticklabelsize = 18, 
        xticklabelpad = 0.5)
    # Row 2 axes: dominant-phase map (left) | melt fraction (centre) | amphibole fraction (right).
    # GridLayout row 3 is used; row 2 (ga[2, *]) is reserved for the row-1 colorbars.
    axphase = Axis(ga[3, 1], xlabel = L"temperature ($\degree$C)", ylabel = "pressure (GPa)", tellwidth = false, width = 600,
        xlabelsize = 20, ylabelsize = 20, yticksize = 5, xticksize = 5, xticklabelsize = 18, yticklabelsize = 18,
        xticklabelpad = 0.5, yticklabelpad = 0.5)  # left: categorical dominant-phase heatmap
    axliq = Axis(ga[3, 2], xlabel = L"temperature ($\degree$C)", tellwidth = false, width = 600,
        xlabelsize = 20, xticksize = 5, xticklabelsize = 18, xticklabelpad = 0.5)  # centre: melt volume fraction
    axamp = Axis(ga[3, 3], xlabel = L"temperature ($\degree$C)", tellwidth = false, width = 600,
        xlabelsize = 20, xticksize = 5, xticklabelsize = 18, xticklabelpad = 0.5)  # third: raw amphibole volume fraction
    # Fourth column in row 2: glaucophane volume fraction (amp where compVariables[3] > 0.5).
    axgl = Axis(ga[3, 4], xlabel = L"temperature ($\degree$C)", tellwidth = false, width = 600,
        xlabelsize = 20, xticksize = 5, xticklabelsize = 18, xticklabelpad = 0.5)  # fourth: glaucophane volume fraction
    # Fifth column in row 2: flow_law categorical map ("gl+ep" vs "amp").
    axflowlaw = Axis(ga[3, 5], xlabel = L"temperature ($\degree$C)", tellwidth = false, width = 600,
        xlabelsize = 20, xticksize = 5, xticklabelsize = 18, xticklabelpad = 0.5)  # fifth: flow law rock type
    # Sixth column in row 2: second dominant phase categorical map.
    axphase2 = Axis(ga[3, 6], xlabel = L"temperature ($\degree$C)", tellwidth = false, width = 600,
        xlabelsize = 20, xticksize = 5, xticklabelsize = 18, xticklabelpad = 0.5)  # sixth: second dominant phase
    # Row 3 axes: g + Na-cpx combined (left) | raw cpx (centre) | Na-cpx (right).
    # GridLayout row 5 is used; row 4 (ga[4, *]) is reserved for the row-2 colorbars.
    axgomph = Axis(ga[5, 1], xlabel = L"temperature ($\degree$C)", ylabel = "pressure (GPa)", tellwidth = false, width = 600,
        xlabelsize = 20, ylabelsize = 20, yticksize = 5, xticksize = 5, xticklabelsize = 18, yticklabelsize = 18,
        xticklabelpad = 0.5, yticklabelpad = 0.5)  # left: garnet + Na-cpx combined volume fraction
    axcpx = Axis(ga[5, 2], xlabel = L"temperature ($\degree$C)", tellwidth = false, width = 600,
        xlabelsize = 20, xticksize = 5, xticklabelsize = 18, xticklabelpad = 0.5)  # centre: raw cpx volume fraction
    axnacpx = Axis(ga[5, 3], xlabel = L"temperature ($\degree$C)", tellwidth = false, width = 600,
        xlabelsize = 20, xticksize = 5, xticklabelsize = 18, xticklabelpad = 0.5)  # right: Na-rich cpx volume fraction
    # Link all twelve axes so zooming or panning one panel simultaneously updates all others.
    linkyaxes!(axleft, axcenter, axright, axphase, axliq, axamp, axgl, axflowlaw, axphase2, axgomph, axcpx, axnacpx)
    linkxaxes!(axleft, axcenter, axright, axphase, axliq, axamp, axgl, axflowlaw, axphase2, axgomph, axcpx, axnacpx)

    # Row 1: density, Cp, bound water
    hmdensity = CairoMakie.heatmap!(axleft, plottingT, plottingP, 
        reshape(filtereddf.density, resolution, resolution), colormap = Reverse(:imola10))
    hmCp = CairoMakie.heatmap!(axcenter, plottingT, plottingP, 
        reshape(filtereddf.Cp, resolution, resolution), colormap = Reverse(:glasgow))
    hmH2Owt = CairoMakie.heatmap!(axright, plottingT, plottingP, 
        reshape(filtereddf.H20minebound .* 100, resolution, resolution), colormap = Reverse(:nuuk))
    # --- Row 2 heatmaps ---
    # Categorical colour gradient: :tab20 provides 20 visually distinct colours (one per phase).
    # colorrange = (0.5, n_phases + 0.5) centres each colour bin on its integer tick value.
    phase_cmap = cgrad(:tab20, n_phases, categorical = true)
    hmphase = CairoMakie.heatmap!(axphase, plottingT, plottingP,
        reshape(dominant_phase_num, resolution, resolution),
        colormap = phase_cmap, colorrange = (0.5, n_phases + 0.5))
    # Melt (liq) volume fraction heatmap.
    hmliq = CairoMakie.heatmap!(axliq, plottingT, plottingP,
        reshape(filtereddf.liq_frac_vol, resolution, resolution), colormap = :thermal)
    # Raw amphibole volume fraction heatmap (all amp regardless of composition).
    hmamp = CairoMakie.heatmap!(axamp, plottingT, plottingP,
        reshape(filtereddf.amp_frac_vol, resolution, resolution), colormap = :thermal)
    # Glaucophane volume fraction heatmap: amp fraction at points where compVariables[3] > 0.5.
    hmgl = CairoMakie.heatmap!(axgl, plottingT, plottingP,
        reshape(filtereddf.gl_frac_vol, resolution, resolution), colormap = :thermal)
    # Flow law heatmap: fixed classes with a continuous blue → red blend for blueschist → eclogite.
    flowlaw_cmap = cgrad([:grey85, :royalblue3, :firebrick3, :darkolivegreen3, :darkorange2])
    hmflowlaw = CairoMakie.heatmap!(axflowlaw, plottingT, plottingP,
        reshape(filtereddf.flow_law_plot_value, resolution, resolution),
        colormap = flowlaw_cmap, colorrange = (1.0, 5.0))
    # Second dominant phase categorical heatmap.
    phase2_cmap = cgrad(:tab20, n_phases2, categorical = true)
    hmphase2 = CairoMakie.heatmap!(axphase2, plottingT, plottingP,
        reshape(second_dominant_phase_num, resolution, resolution),
        colormap = phase2_cmap, colorrange = (0.5, n_phases2 + 0.5))
    # --- Row 3 heatmaps ---
    # Combined garnet + Na-cpx volume fraction: a proxy for eclogite-facies assemblages.
    hmgomph = CairoMakie.heatmap!(axgomph, plottingT, plottingP,
        reshape(filtereddf.g_nacpx_frac_vol, resolution, resolution), colormap = :thermal)
    # Dark red dashed contour at the 0.75 vol. fraction level: marks the boundary above which
    # garnet + Na-cpx together exceed 75 vol% of the assemblage.
    # contour! is intentionally not captured in a variable so the translate! call below
    # (which pushes heatmaps to z = -100) does not affect this line — it stays in front.
    CairoMakie.contour!(axgomph, plottingT, plottingP,
        reshape(filtereddf.g_nacpx_frac_vol, resolution, resolution),
        levels = [0.75], color = :darkred, linestyle = :dash, linewidth = 1.5)
    # Raw clinopyroxene volume fraction heatmap (all cpx regardless of Na content).
    hmcpx = CairoMakie.heatmap!(axcpx, plottingT, plottingP,
        reshape(filtereddf.cpx_frac_vol, resolution, resolution), colormap = :thermal)
    # Na-rich clinopyroxene volume fraction heatmap.
    hmnacpx = CairoMakie.heatmap!(axnacpx, plottingT, plottingP,
        reshape(filtereddf.nacpx_frac_vol, resolution, resolution), colormap = :thermal)

    # Colorbars for row 1
    cbdensity = Colorbar(ga[0, 1][1, 1], hmdensity, label = L"density (kg/m$^3$)", vertical = false, 
        tellwidth = false, width = 350, spinewidth = 0.1, labelpadding = 1.5, labelsize = 16, 
        ticklabelpad = 0.5, ticklabelsize = 14, ticksize = 3.5, 
        ticks = [round(minimum(filtereddf.density), sigdigits = 4), 2500, 3000, 3500, 4000, 4500,
        floor(maximum(filtereddf.density), sigdigits = 4)],
        minorticks = [2250, 2750, 3250, 3750, 4250], minorticksvisible = true, minorticksize = 2.5, 
        minortickwidth = 0.75)
    cbCp = Colorbar(ga[0, 2][1, 1], hmCp, label = "Cp", vertical = false, 
        tellwidth = false, width = 350, spinewidth = 0.1, labelpadding = 1.5, labelsize = 16, 
        ticklabelpad = 0.5, ticklabelsize = 14, ticksize = 3.5, ticks = [0.0, 0.25, 0.5, 0.75, 1.0])
    cbwtH2O = Colorbar(ga[0, 3][1, 1], hmH2Owt, label = "bound water (wt. %)",
        vertical = false, tellwidth = false, width = 350, spinewidth = 0.1, labelpadding = 1.5, 
        labelsize = 16, ticklabelpad = 0.5, ticklabelsize = 14, ticksize = 3.5,
        ticks = [0, 0.5, 1.0, 1.5, 2.0, 2.5])
    # Row 2 colorbars (placed in GridLayout row 2, above the row-2 heatmap axes at row 3).
    # Categorical ticks: positions are the integer values 1 … n_phases; labels are the phase name strings.
    cbphase = Colorbar(ga[2, 1][1, 1], hmphase, label = "dominant phase", vertical = false,
        tellwidth = false, width = 350, spinewidth = 0.1, labelpadding = 1.5, labelsize = 16,
        ticklabelpad = 0.5, ticklabelsize = 10, ticksize = 3.5,
        ticks = (collect(1.0:n_phases), unique_phases))  # maps integer bins → phase-name tick labels
    # Melt volume fraction colorbar.
    cbliq = Colorbar(ga[2, 2][1, 1], hmliq, label = "liq (melt) vol. fraction", vertical = false,
        tellwidth = false, width = 350, spinewidth = 0.1, labelpadding = 1.5, labelsize = 16,
        ticklabelpad = 0.5, ticklabelsize = 14, ticksize = 3.5)
    # Raw amphibole volume fraction colorbar.
    cbamp = Colorbar(ga[2, 3][1, 1], hmamp, label = "amp vol. fraction", vertical = false,
        tellwidth = false, width = 350, spinewidth = 0.1, labelpadding = 1.5, labelsize = 16,
        ticklabelpad = 0.5, ticklabelsize = 14, ticksize = 3.5)
    # Glaucophane volume fraction colorbar.
    cbgl = Colorbar(ga[2, 4][1, 1], hmgl, label = "gl (glaucophane) vol. fraction", vertical = false,
        tellwidth = false, width = 350, spinewidth = 0.1, labelpadding = 1.5, labelsize = 16,
        ticklabelpad = 0.5, ticklabelsize = 14, ticksize = 3.5)
    # Flow law colorbar.
    cbflowlaw = Colorbar(ga[2, 5][1, 1], hmflowlaw, label = "flow law", vertical = false,
        tellwidth = false, width = 350, spinewidth = 0.1, labelpadding = 1.5, labelsize = 16,
        ticklabelpad = 0.5, ticklabelsize = 12, ticksize = 3.5,
        ticks = flow_law_ticks)
    # Second dominant phase categorical colorbar.
    cbphase2 = Colorbar(ga[2, 6][1, 1], hmphase2, label = "second dominant phase", vertical = false,
        tellwidth = false, width = 350, spinewidth = 0.1, labelpadding = 1.5, labelsize = 16,
        ticklabelpad = 0.5, ticklabelsize = 10, ticksize = 3.5,
        ticks = (collect(1.0:n_phases2), unique_phases2))
    # Row 3 colorbars (placed in GridLayout row 4, above the row-3 heatmap axes at row 5).
    # Garnet + Na-cpx combined volume fraction colorbar.
    cbgomph = Colorbar(ga[4, 1][1, 1], hmgomph, label = "g + Na-cpx vol. fraction", vertical = false,
        tellwidth = false, width = 350, spinewidth = 0.1, labelpadding = 1.5, labelsize = 16,
        ticklabelpad = 0.5, ticklabelsize = 14, ticksize = 3.5)
    # Raw clinopyroxene volume fraction colorbar.
    cbcpx = Colorbar(ga[4, 2][1, 1], hmcpx, label = "cpx vol. fraction", vertical = false,
        tellwidth = false, width = 350, spinewidth = 0.1, labelpadding = 1.5, labelsize = 16,
        ticklabelpad = 0.5, ticklabelsize = 14, ticksize = 3.5)
    # Na-rich clinopyroxene volume fraction colorbar.
    cbnacpx = Colorbar(ga[4, 3][1, 1], hmnacpx, label = "Na-cpx vol. fraction", vertical = false,
        tellwidth = false, width = 350, spinewidth = 0.1, labelpadding = 1.5, labelsize = 16,
        ticklabelpad = 0.5, ticklabelsize = 14, ticksize = 3.5)

    # Push all heatmap objects to z = -100 so that tick marks, grid lines, and overlaid elements
    # (e.g. the dark red contour on axgomph) render in front of them at the default z = 0.
    # The contour! call above is not captured in a variable and is therefore unaffected.
    CairoMakie.translate!(hmdensity, 0, 0, -100)
    CairoMakie.translate!(hmCp, 0, 0, -100)
    CairoMakie.translate!(hmH2Owt, 0, 0, -100)
    CairoMakie.translate!(hmphase, 0, 0, -100)   # dominant-phase categorical heatmap
    CairoMakie.translate!(hmliq, 0, 0, -100)      # melt volume fraction heatmap
    CairoMakie.translate!(hmamp, 0, 0, -100)      # amphibole volume fraction heatmap
    CairoMakie.translate!(hmgl, 0, 0, -100)       # glaucophane volume fraction heatmap
    CairoMakie.translate!(hmflowlaw, 0, 0, -100)  # flow law categorical heatmap
    CairoMakie.translate!(hmphase2, 0, 0, -100)   # second dominant phase categorical heatmap
    CairoMakie.translate!(hmgomph, 0, 0, -100)    # g + Na-cpx heatmap (contour stays at z = 0)
    CairoMakie.translate!(hmcpx, 0, 0, -100)      # raw cpx volume fraction heatmap
    CairoMakie.translate!(hmnacpx, 0, 0, -100)    # Na-cpx volume fraction heatmap

    # Beautification
    CairoMakie.ylims!(axright, low = 0)
    CairoMakie.xlims!(axright, low = 0)
    hidespines!(axleft)
    hidespines!(axcenter)
    hidespines!(axright)
    hidespines!(axphase)  # remove box frame from dominant-phase axis
    hidespines!(axliq)    # remove box frame from melt-fraction axis
    hidespines!(axamp)    # remove box frame from amphibole-fraction axis
    hidespines!(axgl)       # remove box frame from glaucophane-fraction axis
    hidespines!(axflowlaw)  # remove box frame from flow-law axis
    hidespines!(axphase2)   # remove box frame from second-dominant-phase axis
    hidespines!(axgomph)    # remove box frame from g + Na-cpx axis
    hidespines!(axcpx)    # remove box frame from cpx-fraction axis
    hidespines!(axnacpx)  # remove box frame from Na-cpx-fraction axis
    hideydecorations!(axright, grid = false)
    hideydecorations!(axcenter, grid = false)
    # Hide y-axis tick marks and labels from all non-leftmost axes in rows 2 and 3;
    # y is shared via linkyaxes! with axphase (row 2) and axgomph (row 3) respectively.
    hideydecorations!(axliq, grid = false)    # shares y with axphase
    hideydecorations!(axamp, grid = false)    # shares y with axphase
    hideydecorations!(axgl, grid = false)       # shares y with axphase
    hideydecorations!(axflowlaw, grid = false)  # shares y with axphase
    hideydecorations!(axphase2, grid = false)   # shares y with axphase
    hideydecorations!(axcpx, grid = false)      # shares y with axgomph
    hideydecorations!(axnacpx, grid = false)  # shares y with axgomph
    colgap!(ga, 10)
    rowgap!(ga, 10)
    display(f)
    save("test.png", f)
end
