# Running MAGEMin and turning its output into a clean P-T lookup table.

# MAGEMin output columns kept for every P-T point: those read by select_magemin_columns, including the
# ones added by compute_phase_fractions!. The rest of MAGEMin's output (e.g. the solution phase data in
# SS_vec and mSS_vec, ~50 kB per point for ig/mp) is dropped chunk by chunk in run_magemin.
const MAGEMIN_COLUMNS = [:P_kbar, :T_C, :rho, :s_cp, :alpha, :enthalpy, :Vp, :Vs, :bulk_S_wt, :entropy, :M_sys,
    :dominant_phase, :second_dominant_phase, :flow_law, :flow_law_plot_value, :dominant_phase_index,
    :liq_frac_vol, :ep_frac_vol, :fluid_frac_wt, :fluid_frac_vol]

# Free fluid phases: the H2O fluid solution of ig and the pure H2O of mp. MAGEMin keeps the water released by
# dehydration in the system as a free fluid (closed system), which lowers the density and raises the heat
# capacity and thermal expansivity. With the remove_free_fluid input the fluid is assumed to escape instead.
const FLUID_PHASES = ("fl", "H2O")

"""Add the free fluid weight and volume fractions (`fluid_frac_wt`, `fluid_frac_vol`) of every point of
MAGEMin's output `out` and, if `remove_fluid`, replace the density `rho`, specific heat `s_cp` and thermal
expansivity `alpha` with those of the solid + melt part. MAGEMin computes these system properties from the
phase fractions (normalised to sum to 1) as 1/ρ = Σ w_k/ρ_k, Cp = Σ w_k·Cp_k (weight fractions w_k, Cp_k per kg)
and α = Σ n_k·α_k (mole fractions n_k; src/toolkit.c and dump_function.c of MAGEMin_C 2.3.8), so the terms of
the FLUID_PHASES are subtracted and the rest renormalised. Points without fluid are unchanged. The phase α
has the same factor-10 error as the system α, so alpha_factor still applies afterwards."""
function free_fluid!(out, remove_fluid)
    n = nrow(out)
    fluid_wt, fluid_vol = zeros(n), zeros(n)
    for i in 1:n
        ph, ph_type = out.ph[i], out.ph_type[i]
        w_F = n_F = volume_F = cp_F = alpha_F = 0.0
        # ph lists the solution phases (ph_type 1, in the order of SS_vec) before the pure phases (PP_vec).
        i_ss = i_pp = 0
        for k in eachindex(ph)
            if ph_type[k] == 1
                i_ss += 1
                data = out.SS_vec[i][i_ss]
            else
                i_pp += 1
                data = out.PP_vec[i][i_pp]
            end
            ph[k] in FLUID_PHASES || continue
            w, n_mol = out.ph_frac_wt[i][k], out.ph_frac[i][k]
            w_F      += w
            n_F      += n_mol
            volume_F += w / data.rho
            cp_F     += w * data.cp
            alpha_F  += n_mol * data.alpha
        end
        w_F > 0.0 || continue
        fluid_wt[i]  = w_F
        fluid_vol[i] = volume_F * out.rho[i]  # the system volume per unit mass is 1/ρ
        if remove_fluid
            out.rho[i]   = (1.0 - w_F) / (1.0 / out.rho[i] - volume_F)
            out.s_cp[i]  = [(out.s_cp[i][1] - cp_F) / (1.0 - w_F)]
            out.alpha[i] = [(out.alpha[i][1] - alpha_F) / (1.0 - n_F)]
        end
    end
    out.fluid_frac_wt  = fluid_wt
    out.fluid_frac_vol = fluid_vol
    return out
end

# MAGEMin minimisation status codes of a converged point (MAGEMin.c PrintStatus: 0 success, 1 under-relaxed,
# 2 heavily under-relaxed); 3 (maximum iterations reached) and 4 (slow convergence or divergence) are failures.
const CONVERGED_STATUS = (0, 1, 2)

"""Set the properties of the points of MAGEMin's output `out` whose minimisation did not converge (status
not in CONVERGED_STATUS) to NaN, so that nanfill! fills them from the neighbouring points of the same isobar.
The flow law classification of these points is kept. Returns the number of such points."""
function mask_unconverged!(out)
    failed = findall(s -> !(s in CONVERGED_STATUS), out.status)
    for i in failed
        for col in (:rho, :Vp, :Vs, :fluid_frac_wt, :fluid_frac_vol, :liq_frac_vol, :ep_frac_vol)
            out[!, col][i] = NaN
        end
        for col in (:s_cp, :alpha, :enthalpy, :entropy)  # 1-element vectors
            out[!, col][i] = [NaN]
        end
    end
    return length(failed)
end

"""Check the inputs of the main file before anything is run, with clear error messages. The kinetic-delay window
is checked only with `usescaling` and the deep blend only for the ig databases, where they are used."""
function check_inputs(database, ig_flow_laws, usescaling, kinetic_delay, deep_blend, resolution, P_kbar_limits,
                      T_C_limits, Xoxides, X)
    known = (IG_DATABASES..., MP_DATABASES..., MANTLE_DATABASES...)
    database in known || error("Unknown database \"$database\"; use one of $(join(known, ", ")).")
    database in IG_DATABASES && check_ig_flow_laws(ig_flow_laws)
    if usescaling
        kd = kinetic_delay
        kd.profile in ("quadratic", "linear") ||
            error("Unknown kinetic_delay.profile \"$(kd.profile)\"; use \"quadratic\" or \"linear\".")
        kd.T_high_C > kd.T_low_C ||
            error("kinetic_delay.T_high_C ($(kd.T_high_C)) must be greater than kinetic_delay.T_low_C ($(kd.T_low_C)).")
    end
    database in IG_DATABASES && (deep_blend.P_end_GPa > deep_blend.P_start_GPa ||
        error("deep_blend.P_end_GPa ($(deep_blend.P_end_GPa)) must be greater than deep_blend.P_start_GPa " *
              "($(deep_blend.P_start_GPa))."))
    resolution >= 2 || error("resolution ($(resolution)) must be at least 2.")
    0.0 <= P_kbar_limits[1] < P_kbar_limits[2] ||
        error("P_kbar_limits $(P_kbar_limits) must be (min, max) with 0 ≤ min < max.")
    T_C_limits[1] < T_C_limits[2] || error("T_C_limits $(T_C_limits) must be (min, max) with min < max.")
    length(Xoxides) == length(X) ||
        error("Xoxides and X have different lengths ($(length(Xoxides)) and $(length(X))).")
    all(>=(0.0), X) && sum(X) > 0.0 || error("The bulk composition X must be non-negative and not all zero.")
    return nothing
end

"""Pressures (kbar) and temperatures (°C) of the `resolution` × `resolution` P-T grid between `P_kbar_limits`
and `T_C_limits` (min, max). Returns (P, T)."""
grid_axes(resolution, P_kbar_limits, T_C_limits) =
    collect(Float64, range(P_kbar_limits..., resolution)), collect(Float64, range(T_C_limits..., resolution))

"""Run MAGEMin on the P-T grid of pressures `P` (kbar) × temperatures `T` (°C) for bulk composition `X`
(oxides `Xoxides`, units `sys_in`) and classify the phases and flow laws (compute_phase_fractions!,
kinetic-delay window `kinetic_delay` of the ig eclogite and mp garnet regimes, or nothing for no delay; ig flow
law scheme `ig_flow_laws`). The free fluid fractions are added and, with `remove_free_fluid`, the fluid is
removed from the density, heat capacity and thermal expansivity (free_fluid!). The properties of points where
the minimisation did not converge are set to NaN (mask_unconverged!) and their number is printed, as are the
MAGEMin version, database and dataset.
Returns a DataFrame with the MAGEMIN_COLUMNS and one row per P-T point; P varies fastest.
The full MAGEMin output does not fit in memory for large grids (2048² points need ~200 GB for ig), so
the grid is computed in chunks of Threads.nthreads() isotherms (one isotherm per thread), and each chunk
is classified and reduced to the MAGEMIN_COLUMNS before the next one starts."""
function run_magemin(database, Xoxides, X, sys_in, P, T, kinetic_delay; ig_flow_laws = "blueschist",
                     remove_free_fluid = true)
    isotherms_per_chunk = Threads.nthreads()
    db = Initialize_MAGEMin(database, verbose = false)
    chunks = DataFrame[]
    n_unconverged = 0
    start_time = time()
    for first_T in 1:isotherms_per_chunk:length(T)
        T_chunk = T[first_T:min(first_T + isotherms_per_chunk - 1, length(T))]
        Prange = repeat(P, outer = length(T_chunk))
        Trange = repeat(T_chunk, inner = length(P))
        out = DataFrame(multi_point_minimization(Prange, Trange, db, X = X, Xoxides = Xoxides, sys_in = sys_in,
            progressbar = false))
        first_T == 1 && println("  MAGEMin $(out.MAGEMin_ver[1]) (MAGEMin_C $(pkgversion(MAGEMin_C))), database ",
            "$(out.database[1]), dataset $(out.dataset[1])"); flush(stdout)
        free_fluid!(out, remove_free_fluid)
        compute_phase_fractions!(out, database, kinetic_delay; ig_flow_laws)
        n_unconverged += mask_unconverged!(out)
        push!(chunks, out[:, MAGEMIN_COLUMNS])

        done = first_T + length(T_chunk) - 1
        elapsed = (time() - start_time) / 60
        println("  $(done)/$(length(T)) isotherms done, $(round(elapsed, digits = 1)) min elapsed, ETA ",
            "$(round(elapsed / done * (length(T) - done), digits = 1)) min"); flush(stdout)
    end
    Finalize_MAGEMin(db)
    n_points = length(P) * length(T)
    println("  Minimisation did not converge at $(n_unconverged) of $(n_points) points ",
        "($(round(100 * n_unconverged / n_points, digits = 3)) %); their properties are filled from the neighbouring points")
    flush(stdout)
    return reduce(vcat, chunks)
end

"""Factor to correct MAGEMin's thermal expansivity α. MAGEMin_C 2.0.6 reports α 10× too small: its
src/toolkit.c divides dV/dT by `dGdP*10.0`, where 1.6.5 divides by `dGdP`. This self-check runs on nearly
pure quartz (ig database, 1 kbar, 25 °C) and compares MAGEMin's α with −(1/ρ) dρ/dT from MAGEMin's own
densities. Returns 10.0 if the bug is present, otherwise 1.0, so the correction stops automatically
once MAGEMin is fixed. Verified against MAGEMin_C 1.6.5, which gives identical densities and 10× larger
α (= −(1/ρ) dρ/dT). Valid for the Holland–Powell ("tc") databases used here (ig, mp, mtl). The
Stixrude ("sb") databases compute α by another route without the `*10.0` (still so in MAGEMin_C 2.3.8),
so for `database` sb11/sb21/sb24 the check is skipped and 1.0 is returned."""
function magemin_alpha_correction(database)
    startswith(database, "sb") && return 1.0
    Xoxides = ["SiO2"; "Al2O3"; "CaO"; "MgO"; "FeO"; "K2O"; "Na2O"; "TiO2"; "O"; "Cr2O3"; "H2O"]
    X       = [99.9;   0.01;    0.01;  0.01;  0.01;  0.0;   0.01;   0.01;   0.0; 0.0;     0.0]
    db = Initialize_MAGEMin("ig", verbose = false)
    point(T_C) = single_point_minimization(1.0, T_C, db, X = X, Xoxides = Xoxides, sys_in = "wt")
    o, o_minus, o_plus = point(25.0), point(24.0), point(26.0)
    Finalize_MAGEMin(db)

    alpha_fd = -(o_plus.rho - o_minus.rho) / 2.0 / o.rho
    ratio = alpha_fd / o.alpha[1]
    if isapprox(ratio, 10.0; rtol = 0.05)
        @warn "MAGEMin_C reports thermal expansivity 10× too small (self-check ratio $(round(ratio, digits = 3))); multiplying α by 10."
        return 10.0
    elseif isapprox(ratio, 1.0; rtol = 0.05)
        return 1.0
    end
    @warn "Unexpected MAGEMin α self-check ratio $(round(ratio, digits = 3)) (expected ≈1 or ≈10); α left uncorrected."
    return 1.0
end

"""Select and rename the MAGEMin output columns used downstream: the bulk thermodynamic properties,
followed by the phase-fraction / flow law columns added by compute_phase_fractions!. Converts
P from kbar to bar and T from °C to K, and sorts by (P, T) so that T varies fastest."""
function select_magemin_columns(df)
    return @chain df begin
        TidierData.@select(P = var"P_kbar", T = var"T_C", density = var"rho", Cp = var"s_cp",
            alpha = var"alpha", enthalpy = var"enthalpy", Vp = var"Vp", Vs = var"Vs", bulk_S_wt = var"bulk_S_wt",
            entropy_S = var"entropy", M_sys = var"M_sys",
            dominant_phase = var"dominant_phase", second_dominant_phase = var"second_dominant_phase", flow_law = var"flow_law", liq_frac_vol = var"liq_frac_vol",
            flow_law_plot_value = var"flow_law_plot_value", dominant_phase_index = var"dominant_phase_index", ep_frac_vol = var"ep_frac_vol",
            fluid_frac_wt = var"fluid_frac_wt", fluid_frac_vol = var"fluid_frac_vol")
        @mutate(P = 1e3 * P, T = 273.15 + T)
        # Standard thermodynamic columns first, then phase-fraction columns.
        TidierData.@select(T, P, density, alpha, Cp, Vp, Vs, enthalpy, bulk_S_wt, entropy_S, M_sys, dominant_phase, second_dominant_phase, flow_law, flow_law_plot_value, dominant_phase_index, liq_frac_vol, ep_frac_vol, fluid_frac_wt, fluid_frac_vol)
        @arrange(P, T)
        @rename(H20minebound = bulk_S_wt)
    end
end

"""Unwrap MAGEMin's 1-element vector columns, convert enthalpy and entropy to per-kg units and
extract the bound water content."""
function unwrap_and_convert_units!(table)
    for col in (:alpha, :Cp, :enthalpy, :entropy_S)
        table[!, col] = getindex.(table[!, col], 1)
    end
    # MAGEMin gives enthalpy in kJ and entropy in kJ/K per mole of system, whose mass is M_sys (g).
    # Convert them to J/kg and J/K/kg as required by the ASPECT tables (MAGEMin's s_cp uses the same
    # normalisation, so it is already in J/K/kg).
    table.enthalpy  = table.enthalpy  .* 1e6 ./ table.M_sys
    table.entropy_S = table.entropy_S .* 1e6 ./ table.M_sys
    # Bound water: H2O is the 11th oxide of the solid bulk composition in the ig database; with fewer
    # oxides (e.g. mtl) there is no H2O, so fill with zeros.
    table.H20minebound = map(v -> length(v) >= 11 ? v[11] : 0.0, table.H20minebound)
    return table
end

# ── Deep blend of the ig tables to sb24 ──
# ig has no transition-zone or lower-mantle phases (majorite, Ca-perovskite, bridgmanite, NAL/Ca-ferrite) and
# predicts false melt above 500 °C at high pressure, so its basaltic crust is ~250 kg/m3 lighter than pyrolite
# below 660 km instead of slightly denser. Above the blend pressures (`deep_blend` input in the main file) the ig
# table is therefore replaced by an sb24 (Stixrude & Lithgow-Bertelloni, 2024) table of the same composition.

# Oxides of the sb24 database.
const SB24_OXIDES = ["SiO2", "Al2O3", "CaO", "MgO", "FeO", "Na2O", "O", "Cr2O3"]

"""Bulk composition for sb24 from the oxides `Xoxides` and amounts `X` of another database: the sb24 oxides
only (K2O, TiO2, H2O etc. are dropped), renormalised to 100. Returns (oxides, amounts)."""
function sb24_composition(Xoxides, X)
    keep = findall(in(SB24_OXIDES), Xoxides)
    return Xoxides[keep], X[keep] .* (100 / sum(X[keep]))
end

# Columns blended linearly from the ig to the sb24 table (join_deep_sb24).
const BLENDED_COLUMNS = [:density, :alpha, :Cp, :Vp, :Vs, :enthalpy, :entropy_S, :liq_frac_vol, :ep_frac_vol,
    :H20minebound, :fluid_frac_wt, :fluid_frac_vol]

"""Weight of the deep (sb24) table at pressure `P_GPa`: 0 up to `blend.P_start_GPa`, rising linearly to 1 at
`blend.P_end_GPa`."""
deep_blend_weight(P_GPa, blend) = clamp((P_GPa - blend.P_start_GPa) / (blend.P_end_GPa - blend.P_start_GPa), 0.0, 1.0)

"""Linear blend (1 − w)·a + w·b; returns `a` at w = 0 and `b` at w = 1 unchanged, even if the other is NaN."""
blend_values(a, b, w) = w <= 0.0 ? a : (w >= 1.0 ? b : (1.0 - w) * a + w * b)

"""Join the ig `table` (isobars up to `blend.P_end_GPa`) and the sb24 table `deep` (isobars from
`blend.P_start_GPa` up) into one table of the full grid. In the overlap (the last isobars of `table`) the
BLENDED_COLUMNS are blended linearly from ig at `blend.P_start_GPa` to sb24 at `blend.P_end_GPa`; above it the
rows are sb24's. Both tables are (P, T)-sorted, α-corrected and in per-kg units, and not yet clamped. From
`blend.P_start_GPa` up the flow law is eclogite (the ig eclogite law, not the mantle laws; user's choice), and
the dominant phases are sb24's from the middle of the blend up. Enthalpy and entropy of ig and sb24 have different
reference states, so they are not continuous across the blend (unused while ASPECT's latent heat is off).
Returns the joined table."""
function join_deep_sb24(table, deep, blend)
    in_overlap = deep.P .<= blend.P_end_GPa * 1e4  # P in bar
    overlap = deep[in_overlap, :]
    rows = (nrow(table) - nrow(overlap) + 1):nrow(table)
    @assert table.P[rows] == overlap.P && table.T[rows] == overlap.T "The sb24 table does not match the ig grid."
    w = deep_blend_weight.(overlap.P ./ 1e4, Ref(blend))
    for col in BLENDED_COLUMNS
        table[!, col][rows] .= blend_values.(table[rows, col], overlap[!, col], w)
    end
    upper = rows[w .>= 0.5]
    table.dominant_phase[upper]        .= overlap.dominant_phase[w .>= 0.5]
    table.second_dominant_phase[upper] .= overlap.second_dominant_phase[w .>= 0.5]

    joined = vcat(table, deep[.!in_overlap, :])
    eclogite_rows = first(rows):nrow(joined)
    joined.flow_law[eclogite_rows] .= "eclogite"
    joined.flow_law_plot_value[eclogite_rows] .= FLOW_LAW_PLOT_MAP_IG["eclogite"]
    joined.dominant_phase_index[eclogite_rows] .= FLOW_LAW_OUTPUT_MAP_IG["eclogite"]
    return joined
end

# Density bounds (kg/m3) for a more stable geodynamic model: the lower bound applies to every database,
# the upper bound only to DENSITY_MAX_DATABASES. The mantle databases and ig (blended to sb24 at high pressure,
# see join_deep_sb24) are not clamped from above, since their lower-mantle densities go above it.
const DENSITY_MIN = 2600.0
const DENSITY_MAX = 4900.0
const DENSITY_MAX_DATABASES = ("mp",)
# Clamps of heat capacity (J/K/kg) and thermal expansivity (1/K), for every database. The figures use the
# same colour ranges (plotting.jl).
const CP_RANGE    = (600.0, 1600.0)
const ALPHA_RANGE = (-2.0e-5, 8.0e-5)

"""Clamp unphysical heat expansivity (ALPHA_RANGE, into the new column alphafilt), specific heat
(CP_RANGE) and density (from DENSITY_MIN, and up to DENSITY_MAX for the DENSITY_MAX_DATABASES) values for
`database`. clamp keeps NaN as NaN (not missing), so nanfill! can still fill it."""
function clamp_properties!(table, database)
    table.alphafilt = clamp.(table.alpha, ALPHA_RANGE...)
    table.Cp = clamp.(table.Cp, CP_RANGE...)
    density_max = database in DENSITY_MAX_DATABASES ? DENSITY_MAX : Inf
    table.density = clamp.(table.density, DENSITY_MIN, density_max)
    return table
end

"""Kinetic-delay scaling factor for every row of the table (T in K, P in bar) for the window `kd`."""
kinetic_scaling_column(table, kd) = kinetic_delay_scaling.(table.T, table.P ./ 1e4, Ref(kd))

# Upper caps of density (kg/m3), heat capacity (J/K/kg) and thermal expansivity (1/K) for the
# metastable state below the kinetic-delay window. The density cap is the density of the metastable rock,
# which compresses and expands with constant compressibility 1/METASTABLE_BULK_MODULUS and constant thermal
# expansivity METASTABLE_ALPHA (the α cap, so that density and α agree in the metastable state), from
# METASTABLE_DENSITY at METASTABLE_T_REF and METASTABLE_P_REF (see metastable_density).
const METASTABLE_DENSITY       = 3000.0   # kg/m3 at the reference conditions
const METASTABLE_T_REF         = 298.15   # K
const METASTABLE_P_REF         = 1.0      # bar
const METASTABLE_BULK_MODULUS  = 85.0e9   # Pa, isothermal; gabbro (plagioclase + clinopyroxene)
const METASTABLE_CP            = 1100.0   # J/K/kg
const METASTABLE_ALPHA         = 2.75e-5

"""Density (kg/m3) of the metastable rock at temperature `T` (K) and pressure `P` (bar), with constant
compressibility β = 1/METASTABLE_BULK_MODULUS and thermal expansivity α = METASTABLE_ALPHA:
ρ = ρ_ref·exp(β(P − P_ref) − α(T − T_ref)), the integral of β = (1/ρ)∂ρ/∂P and α = −(1/ρ)∂ρ/∂T."""
metastable_density(T, P) = METASTABLE_DENSITY * exp((P - METASTABLE_P_REF) * 1e5 / METASTABLE_BULK_MODULUS -
                                                    METASTABLE_ALPHA * (T - METASTABLE_T_REF))

"""Metastability cap: values above `cap` are blended from `cap` (kinetic-delay scaling `s` = 0) to
the lookup table value (s = 1); values at or below `cap` are kept as they are. `cap` is a number or, like
the density cap, one value per row."""
metastable_cap(cap, values, s) = ifelse.(values .> cap, cap .* (1.0 .- s) .+ values .* s, values)

"""ig metastability caps of density (metastable_density, P-T dependent), heat capacity and thermal
expansivity (METASTABLE_CP, METASTABLE_ALPHA), applied to the clamped values. The caps keep the values within
their clamps, since they only lower values above the cap and never below it."""
function scale_ig_metastable!(table)
    density_cap     = metastable_density.(table.T, table.P)
    table.density   = metastable_cap(density_cap, table.density, table.scalingfactor)
    table.Cp        = metastable_cap(METASTABLE_CP, table.Cp, table.scalingfactor)
    table.alphafilt = metastable_cap(METASTABLE_ALPHA, table.alphafilt, table.scalingfactor)
    return table
end

"""Fill NaN values in a vector in place. Leading NaNs are back-filled with the first non-NaN value,
later NaNs are forward-filled; an all-NaN vector is left unchanged."""
function nanfill1d!(v)
    first_valid = findfirst(!isnan, v)
    first_valid === nothing && return v
    v[1:first_valid-1] .= v[first_valid]
    prev = v[first_valid]
    for i in first_valid:lastindex(v)
        if isnan(v[i])
            v[i] = prev
        else
            prev = v[i]
        end
    end
    return v
end

"""Fill NaN values in every AbstractFloat column of a DataFrame sorted by (P, T), i.e. T varies
fastest. Each isobar is filled along T separately, so values never leak between isobars.
NaNs left over (isobars that are entirely NaN) are then filled along P at constant T."""
function nanfill!(df::DataFrame, resolution)
    for col in names(df)
        v = df[!, col]
        eltype(v) <: AbstractFloat || continue
        grid = reshape(v, resolution, resolution)  # dim1 = T, dim2 = P; shares memory with v
        for j in axes(grid, 2)
            nanfill1d!(view(grid, :, j))
        end
        for i in axes(grid, 1)
            nanfill1d!(view(grid, i, :))
        end
    end
    return df
end

"""Turn the classified MAGEMin output into the lookup table written by write_output: select and
convert the columns, correct α by `alpha_factor` (see magemin_alpha_correction), join the deep sb24
output `deep_df` with the blend pressures `deep_blend` (if given; join_deep_sb24), clamp properties,
apply the ig metastability caps of density, heat capacity and thermal expansivity (kinetic-delay window
`kinetic_delay`) and fill NaNs."""
function filter_and_scale(df, resolution, usescaling, database, kinetic_delay, alpha_factor;
                          deep_df = nothing, deep_blend = nothing)
    table = select_magemin_columns(df)
    unwrap_and_convert_units!(table)
    table.alpha = table.alpha .* alpha_factor
    if deep_df !== nothing
        deep = unwrap_and_convert_units!(select_magemin_columns(deep_df))  # sb24 α needs no correction
        table = join_deep_sb24(table, deep, deep_blend)
    end
    if usescaling
        table.scalingfactor = kinetic_scaling_column(table, kinetic_delay)
    end
    clamp_properties!(table, database)
    database == "ig" && usescaling && scale_ig_metastable!(table)
    nanfill!(table, resolution)
    return table
end

"""Scale density for near-surface porosity with the empirical model of Chen et al. (2020),
https://doi.org/10.1007/s10040-020-02214-x. `crusttype` is "oceanic" or "continental".
Currently unused."""
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
    else
        error("Unknown crusttype \"$crusttype\"; use \"oceanic\" or \"continental\".")
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
