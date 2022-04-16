
# Installation
# julia --project "using Pkg; Pkg.instantiate()"
# usage example:
# julia --project process_results.jl "solved_network.nc"

using PowerSystems
using PowerModelsInterface
using PyPSA2PowerSystems
using PowerModels
include("matpower_m_to_mat.jl")

nc_filename = ARGS[1]
nc_dir = dirname(nc_filename)
base_name = first(splitext(basename(nc_filename)))
mp_filename = joinpath(nc_dir, base_name * ".m")

@info "processing PyPSA result" nc_filename
out_path = PyPSA2PowerSystems.format_pypsa(nc_filename, cleanup = true)
sys = PyPSA2PowerSystems.System(out_path)

pm_data = PowerModelsInterface.get_pm_data(sys)

@info "writing MATPOWER file" mp_filename
export_matpower(mp_filename, pm_data)
convert_matpower(mp_filename)

tsd = joinpath(out_path, "timeseries")
if !isempty(readdir(tsd))
    new_tsd = mkpath(joinpath(nc_dir, base_name, "timeseries"))
    @info "writing timeseries files" new_tsd
    cp(tsd, new_tsd, force = true)
    cp(joinpath(out_path, "timeseries_pointers.json"), joinpath(nc_dir, base_name, "timeseries_pointers.json"), force = true)
end
