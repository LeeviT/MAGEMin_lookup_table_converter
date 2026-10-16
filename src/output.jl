# Writing the lookup tables read by ASPECT. Keep both formats and header strings exactly as they are:
# ASPECT reads the column names. The PerpleX header must contain `h,J/kg`; the `_rheol` file follows ASPECT's
# data/material-model/entropy-table/pyrtable/material_table_temperature_pressure_small.txt (no `#` before the
# column-name line, T varies fastest, 1-based dominant_phase_index).

# All outputs of a run go into a folder named `prefix` in the run directory; file names are built
# from the same prefix.
perplex_table_path(prefix)  = joinpath(prefix, prefix * ".txt")
rheology_table_path(prefix) = joinpath(prefix, prefix * "_rheol.txt")

"""Write the PerpleX-format table (T, P, rho, alpha, cp, vp, vs, h) to `filepath`.
`tablename` is the lookup table name written in the header."""
function write_perplex_table(table, resolution, filepath, tablename)
    open(filepath, "w") do io
        write(io,
            "|6.6.6\n"                                                  # PerpleX version
            * tablename * "\n"                                          # Lookup table name
            * "\t\t2\n"                                                 # No. independent variables (T and P)
            * "T(K)\n"                                                  # Temperature parameters
            * "\t" * string(table.T[1]) * "\n"                          # The lower limit of temperature
            * "\t" * string(table.T[2] - table.T[1]) * "\n"             # The temperature step
            * "\t\t" * string(resolution) * "\n"                        # Temperature resolution
            * "P(bar)\n"                                                # Pressure parameters
            * "\t" * string(table.P[1]) * "\n"                          # The lower limit of pressure
            * "\t" * string(table.P[resolution + 1] - table.P[1]) * "\n" # The pressure step
            * "\t\t" * string(resolution) * "\n"                        # Pressure resolution
            * "\t\t8\n"                                                 # Number of columns
            * "T(K)\tP(bar)\trho,kg/m3\talpha,1/K\tcp,J/K/kg\tvp,km/s\tvs,km/s\th,J/kg\n") # Column names and units
        writedlm(io, eachrow(table[!, [:T, :P, :density, :alphafilt, :Cp, :Vp, :Vs, :enthalpy]]), '\t')
    end
end

"""Write the rheology table (T, P, s, rho, alpha, cp, vp, vs, h, dominant_phase_index) in the format of
ASPECT's data/material-model/entropy-table/pyrtable/material_table_temperature_pressure_small.txt."""
function write_rheology_table(table, resolution, filepath)
    open(filepath, "w") do io
        write(io, "# POINTS: $(resolution) $(resolution)\n")
        write(io, "T(K)\tP(bar)\ts,J/K/kg\trho,kg/m3\talpha,1/K\tcp,J/K/kg\tvp,km/s\tvs,km/s\th,J/kg\tdominant_phase_index,1\n")
        columns = [:T, :P, :entropy_S, :density, :alphafilt, :Cp, :Vp, :Vs, :enthalpy, :dominant_phase_index]
        writedlm(io, eachrow(table[!, columns]), '\t')
    end
end

"""Write both lookup tables into the folder `prefix` (created if needed):
`prefix/prefix.txt` (PerpleX, table name `prefix.tab`) and `prefix/prefix_rheol.txt`."""
function write_output(table, resolution, prefix)
    mkpath(prefix)
    write_perplex_table(table, resolution, perplex_table_path(prefix), prefix * ".tab")
    write_rheology_table(table, resolution, rheology_table_path(prefix))
end
