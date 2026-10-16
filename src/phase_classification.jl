# Dominant phases and flow law regimes at each P-T point of the MAGEMin output.

# ── General helpers ──

# The kinetic-delay (metastability) window is set by the `kinetic_delay` input in the main script:
# a NamedTuple (T_low_C, T_high_C, P_onset_GPa, dTdP_C_per_GPa, profile). The temperature limits are T_low_C
# and T_high_C up to P_onset_GPa; above it both change by dTdP_C_per_GPa per GPa. `profile` ("quadratic" or
# "linear") is the shape of the scaling from 0 to 1 across the window.

"""Lower and upper temperature limits (K) of the kinetic-delay window `kd` at pressure `P_GPa`."""
function kinetic_delay_limits(P_GPa, kd)
    shift = kd.dTdP_C_per_GPa * max(P_GPa - kd.P_onset_GPa, 0.0)
    return kd.T_low_C + shift + 273.15, kd.T_high_C + shift + 273.15
end

"""Kinetic-delay scaling (T in K, P in GPa) for the window `kd`, used to blend low-temperature and
high-temperature states: 0 below the lower limit, rising to 1 at the upper limit (see
kinetic_delay_limits), quadratically or linearly depending on `kd.profile`."""
function kinetic_delay_scaling(T, P_GPa, kd)
    T_low, T_high = kinetic_delay_limits(P_GPa, kd)
    T < T_low && return 0.0
    T >= T_high && return 1.0
    x = (T - T_low) / (kd.T_high_C - kd.T_low_C)
    if kd.profile == "quadratic"
        return x^2
    elseif kd.profile == "linear"
        return x
    end
    error("Unknown kinetic_delay.profile \"$(kd.profile)\"; use \"quadratic\" or \"linear\".")
end

"""Total fraction of all instances of phase `name` (e.g. two coexisting fsp or cpx),
returning 0.0 if the phase is absent."""
phase_frac(name, ph, fracs) = sum((fracs[k] for k in eachindex(ph) if ph[k] == name); init = 0.0)

"""Find the dominant and second-dominant solution phase by volume fractions. Volume fractions of
multiple instances of the same phase (e.g. two coexisting cpx) are summed first, so the two
returned names are always distinct."""
function find_dominant_phases(ss_names, ss_fracs)
    isempty(ss_names) && return "none", "none"
    names = unique(ss_names)
    fracs = [phase_frac(name, ss_names, ss_fracs) for name in names]
    order = sortperm(fracs, rev = true)
    dom = names[order[1]]
    length(names) < 2 && return dom, "none"
    return dom, names[order[2]]
end

"""Total weight fraction of the solution-phase instances named `name` whose composition variable
`cv_index` exceeds `threshold`. `ss_names`, `ss_wt` and `ss_vec` list the solution phases only,
in the same order."""
function na_rich_weight(name, cv_index, threshold, ss_names, ss_wt, ss_vec)
    wt = 0.0
    for k in findall(==(name), ss_names)
        cv = ss_vec[k].compVariables
        if length(cv) >= cv_index && cv[cv_index] > threshold
            wt += ss_wt[k]
        end
    end
    return wt
end

# ── Flow law classification for "ig" (igneous) database ──

# Databases classified with the ig flow law rules (blueschist/greenschist/eclogite/melt).
const IG_DATABASES = ("ig",)

# ig flow law schemes (input ig_flow_laws in the main file) and the labels each can produce:
#   "blueschist" (e.g. the hydrated basalt): blueschist, greenschist, >15% melt, eclogite;
#   "protolith" (e.g. the dry gabbro, which cannot form blueschist or much greenschist): metastable protolith (the
#   igneous rock, wherever it is not eclogite or melt), >15% melt, eclogite.
const IG_FLOW_LAW_LABELS = Dict(
    "blueschist" => ["blueschist", "greenschist", ">15% melt", "eclogite"],
    "protolith"  => ["metastable protolith", ">15% melt", "eclogite"],
)

"""Throw an error unless `scheme` is one of the ig flow law schemes (IG_FLOW_LAW_LABELS)."""
check_ig_flow_laws(scheme) = haskey(IG_FLOW_LAW_LABELS, scheme) ||
    error("Unknown ig_flow_laws \"$scheme\"; use \"protolith\" or \"blueschist\".")

# Output file index mapping: blueschist=1, greenschist=2, >15% melt=2, eclogite=3 ("blueschist" scheme);
# metastable protolith=1, >15% melt=2, eclogite=3 ("protolith" scheme)
const FLOW_LAW_OUTPUT_MAP_IG = Dict(
    "blueschist" => 1.0, "metastable protolith" => 1.0, "greenschist" => 2.0, ">15% melt" => 2.0, "eclogite" => 3.0
)
# Visualization plot values (distinct for coloring), in the order of the output index: blueschist=1,
# metastable protolith=2 (index 1), greenschist=3, >15% melt=4 (index 2), eclogite=5 (index 3). The flow law
# figures order their categories by these values.
const FLOW_LAW_PLOT_MAP_IG = Dict(
    "blueschist" => 1.0, "metastable protolith" => 2.0, "greenschist" => 3.0, ">15% melt" => 4.0, "eclogite" => 5.0
)

# Na-rich phases per database, as (phase name, compVariables index) pairs; checked against the
# MAGEMin_C 2.0.6 sources (tc_gss_function.c CV/SF lists, objective_functions.c site fractions):
#   glaucophane-like amphibole: ig amp z = xNaM4, compVariables[3]
#   Na-clinopyroxene: ig cpx n = xNaM2 [4]
const GLAUCOPHANE_PHASES = Dict("ig" => [("amp", 3)])
const NACPX_PHASES       = Dict("ig" => [("cpx", 4)])
const GLAUCOPHANE_MIN_XNAM4 = 0.5   # sodic amphibole
const NACPX_MIN_XNAM2 = 0.3         # omphacite/jadeite, at least ~30% jadeite component

"""Total weight fraction of the Na-rich instances of the (phase, compVariables index) pairs `phases`."""
na_rich_weight(phases, threshold, ss_names, ss_wt, ss_vec) =
    sum((na_rich_weight(name, i, threshold, ss_names, ss_wt, ss_vec) for (name, i) in phases); init = 0.0)

# A kinetically delayed eclogite becomes greenschist only if the greenschist phases make up at least
# this weight fraction; otherwise it stays blueschist (e.g. trace phengite in eclogite at > 4 GPa).
# Used in the "blueschist" scheme only.
const GREENSCHIST_MIN_WT = 0.1

"""Weight fractions of the ig flow law regimes at one P-T point for `database` (one of IG_DATABASES).
`ph`/`pw` list all phases and their weight fractions; `ss_names`/`ss_wt`/`ss_vec` the solution
phases only. Returns (blueschist_wt, greenschist_wt, eclogite_wt, liq_wt).
Note: ep and chl count towards both blueschist and greenschist (intentional)."""
function ig_regime_weights(database, ph, pw, ss_names, ss_wt, ss_vec)
    gl_wt       = na_rich_weight(GLAUCOPHANE_PHASES[database], GLAUCOPHANE_MIN_XNAM4, ss_names, ss_wt, ss_vec)
    amp_only_wt = phase_frac("amp", ph, pw) - gl_wt
    nacpx_wt    = na_rich_weight(NACPX_PHASES[database], NACPX_MIN_XNAM2, ss_names, ss_wt, ss_vec)

    blueschist_wt  = gl_wt + phase_frac("ep", ph, pw) + phase_frac("chl", ph, pw)
    liq_wt         = phase_frac("liq", ph, pw)
    greenschist_wt = amp_only_wt + phase_frac("mu", ph, pw) + phase_frac("fsp", ph, pw) +
                     phase_frac("chl", ph, pw) + phase_frac("ep", ph, pw) +
                     (liq_wt <= 0.15 ? liq_wt : 0.0)
    eclogite_wt    = phase_frac("g", ph, pw) + nacpx_wt
    return blueschist_wt, greenschist_wt, eclogite_wt, liq_wt
end

"""Classify the flow law regime for the igneous database from weight fractions, temperature (°C)
and pressure (kbar) with the flow law scheme `scheme` (IG_FLOW_LAW_LABELS). Returns (label, plot_value).
Melt above 15 % comes first. Eclogite (garnet + Na-cpx over half of the rock, or more than the blueschist and the
greenschist phases) is assigned only when the kinetic-delay scaling (window `kd`) reaches ≥ 0.5; with `kd` =
nothing (usescaling = false) it is never delayed. Everything else is
  "protolith" scheme: metastable protolith (the igneous rock, also where eclogite is delayed);
  "blueschist" scheme: a delayed eclogite is greenschist if greenschist dominates with at least
  GREENSCHIST_MIN_WT, otherwise blueschist; other points the larger of blueschist and greenschist (blueschist wins
  ties and empty sums)."""
function classify_flow_law_ig(blueschist_wt, greenschist_wt, eclogite_wt, liq_wt, T_C, P_kbar, kd, scheme)
    scaling = kd === nothing ? 1.0 : kinetic_delay_scaling(273.15 + T_C, P_kbar / 10, kd)
    eclogite = eclogite_wt > 0.5 || (eclogite_wt > blueschist_wt && eclogite_wt > greenschist_wt)

    label = if liq_wt > 0.15
        ">15% melt"
    elseif eclogite && scaling >= 0.5
        "eclogite"
    elseif scheme == "protolith"
        "metastable protolith"
    elseif eclogite  # kinetically delayed eclogite
        greenschist_wt > blueschist_wt && greenschist_wt >= GREENSCHIST_MIN_WT ? "greenschist" : "blueschist"
    else
        blueschist_wt >= greenschist_wt ? "blueschist" : "greenschist"
    end

    return label, FLOW_LAW_PLOT_MAP_IG[label]
end

# ── Flow law classification for "mp" (metapelite, White et al. 2014) database ──

# mp has no amphibole or clinopyroxene, so the ig blueschist/greenschist/eclogite rules don't apply. Its regimes
# follow the framework minerals of felsic to intermediate crust instead: quartz (the SiO2 polymorphs, pure
# phases), feldspar (fsp) and garnet (g). Garnet stands for eclogitization: at high pressure the anorthite
# component of plagioclase breaks down to grossular + kyanite + quartz. (mp has no jadeite, so albite stays in
# fsp at high pressure.) Other phases (bi, mu, opx, cd, ky, ilm, ...) are not counted.
const MP_DATABASES = ("mp",)
const SILICA_PHASES = ("q", "crst", "trd", "coe", "stv")

# Output file index mapping: quartz=1, feldspar=2, garnet=3, >15% melt=2 (the same slots as the ig regimes)
const FLOW_LAW_OUTPUT_MAP_MP = Dict(
    "quartz" => 1.0, "feldspar" => 2.0, "garnet" => 3.0, ">15% melt" => 2.0
)
# Visualization plot values (distinct for coloring), in the order of the output index: quartz=1 (index 1),
# feldspar=2, >15% melt=3 (index 2), garnet=4 (index 3)
const FLOW_LAW_PLOT_MAP_MP = Dict(
    "quartz" => 1.0, "feldspar" => 2.0, ">15% melt" => 3.0, "garnet" => 4.0
)

"""Weight fractions of the mp flow law regimes at one P-T point, from all phases `ph` and their weight
fractions `pw`. Returns (quartz_wt, feldspar_wt, garnet_wt, liq_wt)."""
function mp_regime_weights(ph, pw)
    quartz_wt = sum(phase_frac(name, ph, pw) for name in SILICA_PHASES)
    return quartz_wt, phase_frac("fsp", ph, pw), phase_frac("g", ph, pw), phase_frac("liq", ph, pw)
end

"""Classify the flow law regime for the metapelite database from weight fractions, temperature (°C)
and pressure (kbar). Returns (label, plot_value).
Melt above 15 % comes first. Garnet is assigned if it makes up over half of the rock or is at least as abundant
as both feldspar and quartz, but, like the ig eclogite, only once the kinetic-delay scaling (window `kd`)
reaches ≥ 0.5; with `kd` = nothing (usescaling = false) it is never delayed. Otherwise the regime is the larger
of feldspar and quartz (feldspar wins ties and empty sums)."""
function classify_flow_law_mp(quartz_wt, feldspar_wt, garnet_wt, liq_wt, T_C, P_kbar, kd)
    scaling = kd === nothing ? 1.0 : kinetic_delay_scaling(273.15 + T_C, P_kbar / 10, kd)
    garnet = garnet_wt > 0.5 || (garnet_wt > 0.0 && garnet_wt >= feldspar_wt && garnet_wt >= quartz_wt)

    label = if liq_wt > 0.15
        ">15% melt"
    elseif garnet && scaling >= 0.5
        "garnet"
    else
        feldspar_wt >= quartz_wt ? "feldspar" : "quartz"
    end

    return label, FLOW_LAW_PLOT_MAP_MP[label]
end

# ── Flow law classification for the mantle databases "mtl" and "sb24" ──

# Databases classified with the mantle flow law rules (olivine/ring+wad/bridgmanite+fp).
const MANTLE_DATABASES = ("mtl", "sb24")

# Mantle phase names per database (MAGEMin_C 2.3.8 sources):
#   mtl (Holland et al., 2013): ol = olivine, wad = wadsleyite, ring = ringwoodite,
#     mpv = MgSi-perovskite (bridgmanite), fp = ferropericlase
#   sb24 (Stixrude & Lithgow-Bertelloni, 2024): ol, wa = wadsleyite, ri = ringwoodite,
#     pv = bridgmanite, mw = magnesiowüstite (ferropericlase). sb24 post-perovskite (ppv) is not counted,
#     as mtl has none.
const MANTLE_PHASES = Dict(
    "mtl"  => (ol = "ol", wad = "wad", ring = "ring", bm = "mpv", fp = "fp"),
    "sb24" => (ol = "ol", wad = "wa",  ring = "ri",   bm = "pv",  fp = "mw"),
)

# The plot values double as the output file index.
const FLOW_LAW_PLOT_MAP_MTL = Dict(
    "olivine" => 1.0, "ring/wad" => 2.0, "bridgmanite+fp" => 3.0
)

"""Weight fractions of the mantle flow law regimes at one P-T point for `database` (one of
MANTLE_DATABASES). Returns (ol_wt, ring_wad_wt, bm_fp_wt)."""
function mtl_regime_weights(database, ph, pw)
    names = MANTLE_PHASES[database]
    ring_wad_wt = phase_frac(names.wad, ph, pw) + phase_frac(names.ring, ph, pw)
    bm_fp_wt    = phase_frac(names.bm, ph, pw) + phase_frac(names.fp, ph, pw)
    return phase_frac(names.ol, ph, pw), ring_wad_wt, bm_fp_wt
end

"""Classify the flow law regime for the mantle database from weight fractions.
Returns (label, plot_value). Ties go to the deeper phase; no mantle phases gives bridgmanite+fp."""
function classify_flow_law_mtl(ol_wt, ring_wad_wt, bm_fp_wt)
    total = ol_wt + ring_wad_wt + bm_fp_wt
    total == 0.0 && return "bridgmanite+fp", 3.0

    label = if bm_fp_wt >= ol_wt && bm_fp_wt >= ring_wad_wt
        "bridgmanite+fp"
    elseif ring_wad_wt >= ol_wt
        "ring/wad"
    else
        "olivine"
    end

    return label, FLOW_LAW_PLOT_MAP_MTL[label]
end

# ── Classification of the whole table ──

"""Plot value of every flow law label of `database`."""
function flow_law_plot_map(database)
    database in IG_DATABASES && return FLOW_LAW_PLOT_MAP_IG
    database in MP_DATABASES && return FLOW_LAW_PLOT_MAP_MP
    database in MANTLE_DATABASES && return FLOW_LAW_PLOT_MAP_MTL
    error("No flow law classification for database \"$database\"; use \"ig\", \"mp\", \"mtl\" or \"sb24\".")
end

"""Index of the flow law regime written to the rheology table (1-based)."""
function flow_law_output_index(database, label, plot_val)
    database in IG_DATABASES && return FLOW_LAW_OUTPUT_MAP_IG[label]
    database in MP_DATABASES && return FLOW_LAW_OUTPUT_MAP_MP[label]
    return plot_val
end

"""Flow law label used where MAGEMin found no stable phases (ig flow law scheme `ig_flow_laws`)."""
function default_flow_law(database, ig_flow_laws)
    database in MANTLE_DATABASES && return "bridgmanite+fp"
    database in MP_DATABASES && return "feldspar"
    return ig_flow_laws == "protolith" ? "metastable protolith" : "blueschist"
end

"""Classify the flow law regime at MAGEMin output row `i` (ig flow law scheme `ig_flow_laws`). Returns
(label, plot_value)."""
function classify_flow_law(df, i, database, ph, pw, ss_names, ss_wt, kinetic_delay, ig_flow_laws)
    if database in IG_DATABASES
        weights = ig_regime_weights(database, ph, pw, ss_names, ss_wt, df.SS_vec[i])
        return classify_flow_law_ig(weights..., df.T_C[i], df.P_kbar[i], kinetic_delay, ig_flow_laws)
    elseif database in MP_DATABASES
        return classify_flow_law_mp(mp_regime_weights(ph, pw)..., df.T_C[i], df.P_kbar[i], kinetic_delay)
    elseif database in MANTLE_DATABASES
        return classify_flow_law_mtl(mtl_regime_weights(database, ph, pw)...)
    end
    error("No flow law classification for database \"$database\"; use \"ig\", \"mp\", \"mtl\" or \"sb24\".")
end

"""Add dominant phase, flow law and phase fraction columns to the MAGEMin output `df`:
dominant_phase, second_dominant_phase, flow_law, flow_law_plot_value, dominant_phase_index,
liq_frac_vol and ep_frac_vol. `kinetic_delay` is the kinetic-delay window of the ig eclogite and the mp garnet
regime, or nothing for no delay (usescaling = false); `ig_flow_laws` the ig flow law scheme (IG_FLOW_LAW_LABELS)."""
function compute_phase_fractions!(df, database, kinetic_delay; ig_flow_laws = "blueschist")
    check_ig_flow_laws(ig_flow_laws)
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

        # No stable phases at this P-T point: use the database's default regime everywhere and
        # keep zero phase fractions.
        if isempty(ph)
            label = default_flow_law(database, ig_flow_laws)
            plot_val = flow_law_plot_map(database)[label]
            dominant_phase[i] = second_dominant_phase[i] = flow_law_col[i] = label
            flow_law_plot_col[i] = plot_val
            dominant_phase_idx[i] = flow_law_output_index(database, label, plot_val)
            continue
        end

        # Solution phases only (ph_type == 1), in the same order as df.SS_vec[i].
        ss_mask  = pt .== 1
        ss_names = ph[ss_mask]
        ss_fracs = pf[ss_mask]
        ss_wt    = pw[ss_mask]

        dominant_phase[i], second_dominant_phase[i] = find_dominant_phases(ss_names, ss_fracs)

        label, plot_val = classify_flow_law(df, i, database, ph, pw, ss_names, ss_wt, kinetic_delay, ig_flow_laws)
        flow_law_col[i]       = label
        flow_law_plot_col[i]  = plot_val
        dominant_phase_idx[i] = flow_law_output_index(database, label, plot_val)

        liq_frac_vol[i] = phase_frac("liq", ph, pf)
        ep_frac_vol[i]  = phase_frac("ep",  ph, pf)
    end

    # Attach results as DataFrame columns.
    df.dominant_phase        = dominant_phase
    df.second_dominant_phase = second_dominant_phase
    df.flow_law              = flow_law_col
    df.flow_law_plot_value   = flow_law_plot_col
    df.dominant_phase_index  = dominant_phase_idx
    df.liq_frac_vol          = liq_frac_vol
    df.ep_frac_vol           = ep_frac_vol
    return df
end
