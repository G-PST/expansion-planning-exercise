
# Installation
# julia --project "using Pkg; Pkg.instantiate()"
# usage example:
# julia --project process_results.jl "solved_network.nc"

using PowerSystems
using PowerModelsInterface
using PyPSA2PowerSystems
using PowerModels

nc_filename = ARGS[1]
mp_filename = first(splitext(nc_filename)) * ".m"

@info "processing PyPSA result" nc_filename
out_path = PyPSA2PowerSystems.format_pypsa(nc_filename, cleanup = true)
sys = PyPSA2PowerSystems.System(out_path)

pm_data = PowerModelsInterface.get_pm_data(sys)

@info "writing MATPOWER file" mp_filename
export_matpower("solved_network.m", pm_data)

tsd = joinpath(out_path, "timeseries")
if !isempty(readir(tsd))
    @info "writing timeseries files" tsd
    cp(tsd, dirname(mp_filename))
    cp(joinpath(out_path, "timeseries_pointers.json"), joinpath(dirname(mp_filename), "timeseries_pointers.json"))
