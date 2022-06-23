
using PowerSystems
using PyPSA2PowerSystems
using PowerSimulations
using PowerGraphics
plotlyjs()
using PowerSystemCaseBuilder
using Dates
using Logging

logger = configure_logging(console_level = Logging.Info, filename = "log.txt")

#using HiGHS
#solver = optimizer_with_attributes(HiGHS.Optimizer, "mip_rel_gap" => 0.01)
using Xpress
solver = optimizer_with_attributes(Xpress.Optimizer, "MIPRELSTOP" => 0.5, "OUTPUTLOG" => 1)

reference_sys = build_system(PSITestSystems, "RTS_GMLC_sys")

# read data from PyPSA results
pypsa_path = "../example-pypsa/"
nc_filename = "solved_network_RTS_GMLS_base+gen_expansion.nc"
nc_path = joinpath(pypsa_path, nc_filename)

function prep_pypsa_sys(nc_path, ref_sys)
    sys = PyPSA2PowerSystems.System(nc_path)

    set_units_base_system!(sys, "natural_units")
    set_units_base_system!(ref_sys, "natural_units")

    for ref_g in get_components(ThermalStandard, ref_sys, get_available)
        gen = get_component(ThermalStandard, sys, "generators_" * get_name(ref_g))
        if !isnothing(gen)
            set_active_power_limits!(gen, get_active_power_limits(ref_g))
            set_time_limits!(gen, get_time_limits(ref_g))
            set_operation_cost!(gen, get_operation_cost(ref_g))
            set_ramp_limits!(gen, get_ramp_limits(ref_g))
        end
    end
    transform_single_time_series!(sys, 48, Hour(24))

    sim_dir = mkpath(replace(basename(nc_path), ".nc"=>""))
    to_json(sys, joinpath(sim_dir, replace(basename(nc_path), ".nc"=>".json")), force = true)
end

function run_pypsa_sim(sys, sim_dir)

    ###### Defining the problem template for DA problem
    template = ProblemTemplate()
    #set_network_model!(template, NetworkModel(DCPPowerModel, duals = [NodalBalanceActiveConstraint], use_slacks = true))
    set_network_model!(template, NetworkModel(CopperPlatePowerModel, duals = [CopperPlateBalanceConstraint], use_slacks = true))
    set_device_model!(template, Line, StaticBranchUnbounded)
    set_device_model!(template, Transformer2W, StaticBranchUnbounded)
    set_device_model!(template, MonitoredLine, StaticBranch)
    set_device_model!(template, RenewableDispatch, RenewableFullDispatch)
    set_device_model!(template, HydroDispatch, HydroDispatchRunOfRiver)
    set_device_model!(template, ThermalStandard, ThermalStandardUnitCommitment)
    #set_device_model!(template, GenericBattery, BookKeeping)
    set_device_model!(template, PowerLoad, StaticPowerLoad)

    ###### start with a single step model

    UC_model = DecisionModel(
        template,
        sys,
        name = "UC",
        optimizer = solver,
        optimizer_solve_log_print = true,
        initialize_model = false,
    )

    # build!(UC_model, output_dir = tmp_dir)

    # solve!(UC_model)

    # res = ProblemResults(UC_model);
    # plot_fuel(res);

    ###### Defining the simultions models
    # DecisionModel() takes the problem template, system, a name for the model, and an solver.
    models = SimulationModels(
        decision_models =[
            UC_model,
        ]
    )

    # The feedforward defines how information passes between the UC and ED models
    # Here we specify that the commiment varible (OnVariable in SIIP) is pass forward
    # to the ED problem, for Thermal Devices restricting the power output (ActivePowerVariable in SIIP) varible in the ED problem
    # feedforward = Dict(
    #     "ED" => [
    #         SemiContinuousFeedforward(
    #             component_type = ThermalStandard,
    #             source = OnVariable,
    #             affected_values = [ActivePowerVariable],
    #         ),
    #     ],
    # )

    sequence = SimulationSequence(
        models = models,
        ini_cond_chronology = InterProblemChronology(),
        # feedforwards = feedforward,
    )

    # Specify the simulation setup
    sim = Simulation(
        name = "simulation",
        # initial_time = Dates.DateTime("2012-07-01T00:00:00"),
        steps = 364,
        models = models,
        sequence = sequence,
        simulation_folder = sim_dir,
    )

    # This builds the simulation and creates the necesary file structure for saving simulaiton files
    build!(sim)
    # This starts executing the simulation we have setup.
    execute!(sim, enable_progress_bar = true)

    # Create the PowerSimulation results object
    results = SimulationResults(sim);
    results_uc = get_decision_problem_results(results, "UC");
    return results_uc
end

nc_files = [
    "solved_network_RTS_GMLS_1p5xload+0emission+gen_and_line_expansion.nc",
    "solved_network_RTS_GMLC_base+line_expansion.nc",
    "solved_network_RTS_GMLS_base+gen_and_line_expansion.nc",
    "solved_network_RTS_GMLS_base+gen_expansion.nc",
]

for nc_filename in nc_files
    prep_pypsa_sys(joinpath(pypsa_path, nc_filename), reference_sys)
end

results_dict = Dict()
for nc_filename in nc_files
    nc_basename = replace(nc_filename, ".nc"=>"")
    sim_sys = PowerSystems.System(joinpath(nc_basename, nc_basename * ".json"))
    res = run_pypsa_sim(sim_sys, nc_basename)
    results_dict[nc_basename] = res
end

#plot_fuel(results_uc);
plots_dict = Dict()
for (k,v) in results_dict
    export_realized_results(v)
    plots_dict[k] = plot_fuel(v, bar = true, title = k);
end

ts_plots_dict = Dict()
for (k,v) in results_dict
    ts_plots_dict[k] = plot_fuel(v, title = k, len = 96, save = nc_basename);
end

prices = Dict()
for (k,v) in results_dict
    prices[k] = sort(read_realized_dual(v, "CopperPlateBalanceConstraint__System").CopperPlateBalanceConstraint__System, rev = true)
end
prices = DataFrame(prices)
plot(Matrix(prices), lab = names)

prices_ts = Dict()
for (k,v) in results_dict
    prices_ts[k] = read_realized_dual(v, "CopperPlateBalanceConstraint__System").CopperPlateBalanceConstraint__System
end
prices_ts = DataFrame(prices_ts)
insertcols!(prices_ts, 1, :DateTime =>  read_realized_dual(first(values(results_dict)), "CopperPlateBalanceConstraint__System").DateTime)

#=
plots_dict = Dict()
for nc_filename in nc_files
    nc_basename = replace(nc_filename, ".nc"=>"")
    sim_res = SimulationResults(joinpath(nc_basename, "simulation"))
    uc_res = get_decision_problem_results(sim_res, "UC")
    res_sys = PowerSystems.System(joinpath(nc_basename, nc_basename * ".json"))
    set_system!(uc_res, res_sys)

    plots_dict[nc_basename] = plot_fuel(uc_res, bar = true, title = nc_basename);
end
=#