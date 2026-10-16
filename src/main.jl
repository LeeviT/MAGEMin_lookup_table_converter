# This script computes a phase transition lookup table with MAGEMin_C, classifies the
# flow law regime at each P-T point, and writes the table in PerpleX format plus a rheology table
# for ASPECT. Edit the input parameters below and run this file; the code is split into
#   table_processing.jl     – runs MAGEMin (in chunks of isotherms) and turns its output into a clean P-T table
#   phase_classification.jl – dominant phases and flow law regimes at each P-T point
#   output.jl               – writes the PerpleX and rheology (_rheol) lookup tables
#   plotting.jl             – figures (single-variable, 2×2 summary and _zoom versions)

# DataFrames for better handling table format, CairoMakie for plotting, TidierData for
# filtering/manipulating data, DelimitedFiles for output, MAGEMin_C for the phase equilibria.
using CairoMakie
using TidierData
using DelimitedFiles
using DataFrames
using MAGEMin_C
println("MAGEMin_C loaded"); flush(stdout)

include(joinpath(@__DIR__, "table_processing.jl"))
include(joinpath(@__DIR__, "phase_classification.jl"))
include(joinpath(@__DIR__, "output.jl"))
include(joinpath(@__DIR__, "plotting.jl"))

# --- Input parameters start ----
database = "ig"  # "ig" igneous (see ig_flow_laws), "mp" metapelite (quartz/feldspar/garnet),
                 # "mtl"/"sb24" mantle (olivine/ring+wad/bridgmanite+fp flow law rules)
# ig flow law scheme, i.e. what the flow law indices 1/2/3 of the rheology table mean (match the prm!):
#   "protolith" (dry gabbro): 1 metastable protolith (the igneous rock wherever it is not eclogite or melt),
#                             2 >15% melt, 3 eclogite
#   "blueschist" (hydrated basalt): 1 blueschist, 2 greenschist or >15% melt, 3 eclogite
ig_flow_laws = "protolith"
# Free fluid (ig "fl", mp "H2O"): MAGEMin keeps the water released by dehydration in the rock as a free fluid
# (closed system). With true the fluid is assumed to escape: the density, heat capacity and thermal expansivity
# are those of the solid + melt part. With false the system values, fluid included, are kept. The log gives the
# maximum free fluid fraction either way.
remove_free_fluid = true
# Kinetic delay (metastability): delays the eclogite (ig) or garnet (mp) flow law and caps density, Cp and α in
# the metastable state (ig only). With false, eclogite/garnet is assigned wherever the phase weights give it.
usescaling = true
# Kinetic-delay window, used only with usescaling: the temperature limits are T_low_C–T_high_C (°C) up to
# P_onset_GPa, and above it both change by dTdP_C_per_GPa per GPa. profile: "quadratic" or "linear" rise of
# the scaling from 0 to 1 across the window. Eclogite (garnet for mp) is assigned once the scaling reaches 0.5:
# T_low + 0.71·(T_high − T_low) (quadratic, e.g. + 141 °C for a 200 °C window) or the window midpoint (linear).
# outputprefix gets the postfix _<database>, plus _linear with usescaling and "linear" (e.g. gabbro_ig_linear).
kinetic_delay = (T_low_C = 600.0, T_high_C = 800.0, P_onset_GPa = 3.0, dTdP_C_per_GPa = -100.0, profile = "quadratic")
# ig only: ig has no transition-zone or lower-mantle phases, so from P_start_GPa up MAGEMin is also run with sb24
# for the same composition (sb24 oxides only, renormalised), and the table properties are blended linearly to
# the sb24 values between P_start_GPa and P_end_GPa. The flow law stays eclogite from P_start_GPa up.
deep_blend = (P_start_GPa = 12.0, P_end_GPa = 13.0)
# Output prefix: all outputs go into the folder <outputprefix>/ in the run directory, as
# <outputprefix>.txt (PerpleX), <outputprefix>_rheol.txt, one <outputprefix>_<variable>.png per
# plotted variable and the 2×2 <outputprefix>_summary.png.
# For ig and mp the figures are also drawn for P ≤ P_max_GPa and T ≤ T_max_C, with a _zoom postfix.
zoom = (P_max_GPa = 7.0, T_max_C = 1500.0)
outputprefix = "4_oc_LC_gabbro"

# P-T grid: resolution × resolution points between these limits (P in kbar, T in °C).
resolution = 1024
P_kbar_limits = (0.01, 500)
T_C_limits = (0, 1900)

# Bulk composition; uncomment one Xoxides/X pair matching the database.
sys_in = "wt"

# Mantle compositions are after Stixrude & Lithgow-Bertelloni (2024), "sb24"
# Xoxides = ["SiO2";    "Al2O3";   "CaO";     "MgO";     "FeO";    "Na2O";   "O";      "Cr2O3"]     # mantle oxides
# X =       [45.14165;  4.005032;  3.190329;  39.00577;  7.910777; 0.131924; 0.055728; 0.55879]     # pyrolite (asthenosphere+lower mantle)
# X =       [43.458859; 0.999736;  0.900553;  45.958103; 8.343993; 0.010546; 0.027565; 0.300644]    # harzburgite (lith. mantle)

# Igneous compositions are after Green et al., (2025), corrected after Holland et al., (2018), "ig"
Xoxides = ["SiO2";    "Al2O3";   "CaO";     "MgO";    "FeO";    "K2O";    "Na2O";   "TiO2";   "O";     "Cr2O3";  "H2O"]     # ig oxides
X =       [50.992811; 15.304009; 10.922183; 7.848858; 9.913537; 0.135232; 2.867050; 1.541733; 0.176099; 0.048487; 0.25]     # oc LC gabbro
# X =       [50.376296; 15.118980; 10.790131; 7.753963; 9.793680; 0.133597; 2.832387; 1.523093; 0.17397; 0.047901; 1.456]   # oc UC basalt

# Continental crust compositions are the average UCC and LCC presets of MAGEMinApp after Rudnick & Gao (2003),
# converted from mol to wt%, "mp"
# Xoxides = ["SiO2";    "Al2O3";   "CaO";    "MgO";    "FeO";    "K2O";    "Na2O";   "TiO2";   "O"; "MnO";    "H2O"]       # mp oxides
# X =       [64.844565; 14.994035; 3.495369; 2.414637; 4.907109; 2.726215; 3.183776; 0.62319;  0.0; 0.097348; 2.713756]    # cc UC average UCC
# X =       [53.464222; 16.920327; 9.601513; 7.248689; 8.580332; 0.61066;  2.65316;  0.820996; 0.0; 0.100101; 0.0]         # cc LC average LCC
# ---- Input parameters end ----

# ---- Main execution ----
check_inputs(database, ig_flow_laws, usescaling, kinetic_delay, deep_blend, resolution, P_kbar_limits, T_C_limits,
    Xoxides, X)
# Postfix of the output prefix (folder and files): the database, then _linear for the linear kinetic-delay
# profile (only with usescaling, since the profile is unused otherwise), e.g. gabbro_ig or gabbro_ig_linear.
outputprefix *= "_$(database)"
usescaling && kinetic_delay.profile == "linear" && (outputprefix *= "_linear")
# Create the output folder (in the run directory, output-production/ in the batch job) right away. If the
# batch script sets MAGEMIN_LOG_FILE, the job log is copied into it when Julia exits, also after an error.
mkpath(outputprefix)
log_file = get(ENV, "MAGEMIN_LOG_FILE", "")
if !isempty(log_file)
    atexit() do
        flush(stdout); flush(stderr)
        isfile(log_file) && cp(log_file, joinpath(outputprefix, basename(log_file)); force = true)
    end
end
# Record the Julia environment of the run (exact package versions, e.g. MAGEMin_C) in the output folder.
project_file = Base.active_project()
if project_file !== nothing
    project_dir = dirname(project_file)
    for file in readdir(project_dir)
        (file == "Project.toml" || (startswith(file, "Manifest") && endswith(file, ".toml"))) &&
            cp(joinpath(project_dir, file), joinpath(outputprefix, file); force = true)
    end
end
println("Julia $(VERSION), environment $(project_file === nothing ? "none" : dirname(project_file)) ",
    "(its Project.toml and Manifest.toml are copied into $(outputprefix)/)")
println("Output prefix: $(outputprefix)")
println("Database: $(database)")
database in IG_DATABASES && println("ig flow laws: $(ig_flow_laws) (", join(["$(l) → $(Int(FLOW_LAW_OUTPUT_MAP_IG[l]))"
    for l in sort(IG_FLOW_LAW_LABELS[ig_flow_laws]; by = l -> FLOW_LAW_PLOT_MAP_IG[l])], ", "), ")")
println("Kinetic delay (usescaling): $(usescaling ? kinetic_delay : "off")")
println("Remove free fluid: $(remove_free_fluid)")
println("Bulk composition ($(sys_in)):")
for (oxide, x) in zip(Xoxides, X)
    println("  $(rpad(oxide, 6)) $(x)")
end
flush(stdout)
alpha_factor = magemin_alpha_correction(database)
println("Thermal expansivity α correction factor: $(alpha_factor)"); flush(stdout)
P_kbar, T_C = grid_axes(resolution, P_kbar_limits, T_C_limits)
# ig with the deep blend: ig only up to the end of the blend, sb24 only from its start (the overlap is blended).
P_deep = P_kbar[P_kbar .>= 10 * deep_blend.P_start_GPa]
use_deep_blend = database in IG_DATABASES && !isempty(P_deep)
P_main = use_deep_blend ? P_kbar[P_kbar .<= 10 * deep_blend.P_end_GPa] : P_kbar
println("Running MAGEMin and classifying phases and flow laws ($database, $(length(P_main)) isobars × ",
    "$(length(T_C)) isotherms of the $(resolution)×$(resolution) grid, $(Threads.nthreads()) threads)"); flush(stdout)
# Without usescaling the flow law gets no kinetic-delay window (nothing), i.e. no eclogite (ig) or garnet (mp) delay.
df = run_magemin(database, Xoxides, X, sys_in, P_main, T_C, usescaling ? kinetic_delay : nothing;
    ig_flow_laws, remove_free_fluid)
deep_df = nothing
if use_deep_blend
    Xoxides_sb24, X_sb24 = sb24_composition(Xoxides, X)
    println("Running MAGEMin with sb24 for the $(length(P_deep)) isobars from $(deep_blend.P_start_GPa) GPa up ",
        "(blend $(deep_blend.P_start_GPa)–$(deep_blend.P_end_GPa) GPa); sb24 bulk composition ($(sys_in)):")
    for (oxide, x) in zip(Xoxides_sb24, X_sb24)
        println("  $(rpad(oxide, 6)) $(round(x, digits = 6))")
    end
    flush(stdout)
    deep_df = run_magemin("sb24", Xoxides_sb24, X_sb24, sys_in, P_deep, T_C, nothing; remove_free_fluid)
end
println("Building the lookup table"); flush(stdout)
table = filter_and_scale(df, resolution, usescaling, database, kinetic_delay, alpha_factor; deep_df, deep_blend)
println("Free fluid at $(count(>(0.0), table.fluid_frac_wt)) of $(nrow(table)) points, max ",
    "$(round(100 * maximum(table.fluid_frac_wt), digits = 2)) wt% / $(round(100 * maximum(table.fluid_frac_vol), digits = 2)) vol% (",
    remove_free_fluid ? "removed from" : "kept in", " density, Cp and α)"); flush(stdout)
println("Writing $(perplex_table_path(outputprefix)) and $(rheology_table_path(outputprefix))"); flush(stdout)
write_output(table, resolution, outputprefix)
println("Plotting into $(outputprefix)/"); flush(stdout)
plot_results(table, resolution, database, outputprefix, zoom, ig_flow_laws)
