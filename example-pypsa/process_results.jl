
# Installation
# julia --project "using Pkg; Pkg.instantiate()"
# usage example:
# julia --project process_results.jl "solved_network.nc"

using PowerSystems
using PowerModelsInterface
using PyPSA2PowerSystems
using PowerModels
using MAT

function main(nc_filename)

    nc_dir = dirname(nc_filename)
    base_name = first(splitext(basename(nc_filename)))
    mp_filename = joinpath(nc_dir, base_name * ".m")
    mat_filename = joinpath(nc_dir, base_name * ".mat")

    @info "processing PyPSA result" nc_filename
    out_path = PyPSA2PowerSystems.format_pypsa(nc_filename, cleanup = true)
    sys = PyPSA2PowerSystems.System(out_path)

    pm_data = PowerModelsInterface.get_pm_data(sys)

    @info "writing MATPOWER .mat file" mat_filename
    export_matpower_mat(mat_filename, pm_data)

    @info "writing MATPOWER .m file" mp_filename
    export_matpower_m(mp_filename, pm_data)
end

function get_bus(data)
    buses = Float64[]
    bus_indices = sort(parse.(Int, keys(data["bus"])))
    @assert first(bus_indices) == 1 && last(bus_indices) == length(bus_indices)
    for bi in bus_indices
        bus_data = data["bus"][string(bi)]
        pd = 0 # TODO: get pd from load
        qd = 0 # TODO: get qd from load
        gs = 0 # TODO: get better default for gs
        bs = 0 # TODO: get better default for bs
        push!(buses, bus_data["bus_i"])
        push!(buses, bus_data["index"])
        push!(buses, pd)
        push!(buses, qd)
        push!(buses, gs)
        push!(buses, bs)
        push!(buses, bus_data["area"])
        push!(buses, bus_data["vm"])
        push!(buses, bus_data["va"])
        push!(buses, bus_data["base_kv"])
        push!(buses, parse(Float64, bus_data["zone"]))
        push!(buses, bus_data["vmax"])
        push!(buses, bus_data["vmin"])
    end
    permutedims(reshape(buses, 13, length(buses) ÷ 13))
end

function get_bus_name(data)
    buses = String[]
    bus_indices = sort(parse.(Int, keys(data["bus"])))
    @assert first(bus_indices) == 1 && last(bus_indices) == length(bus_indices)
    for bi in bus_indices
        bus_data = data["bus"][string(bi)]
        push!(buses, bus_data["name"])
    end
    buses
end

function get_gen(data)
    gens = Float64[]
    indices = sort(parse.(Int, keys(data["gen"])))
    @assert first(indices) == 1 && last(indices) == length(indices)
    for i in indices
        gen_data = data["gen"][string(i)]
        # TODO: get better defaults
        pc1      = gen_data["pmin"]
        pc2      = gen_data["pmax"]
        qc1min   = gen_data["qmin"]
        qc1max   = gen_data["qmax"]
        qc2min   = gen_data["qmin"]
        qc2max   = gen_data["qmax"]
        ramp_10  = get(gen_data, "ramp_10", gen_data["pmax"])
        ramp_agc = ramp_10
        ramp_30  = ramp_10
        ramp_q   = gen_data["qmax"]
        apf   = 1
        push!(gens, gen_data["gen_bus"])
        push!(gens, gen_data["pg"])
        push!(gens, gen_data["qg"])
        push!(gens, gen_data["qmax"])
        push!(gens, gen_data["qmin"])
        push!(gens, gen_data["vg"])
        push!(gens, gen_data["mbase"])
        push!(gens, gen_data["gen_status"])
        push!(gens, gen_data["pmax"])
        push!(gens, gen_data["pmin"])
        push!(gens, pc1)
        push!(gens, pc2)
        push!(gens, qc1min)
        push!(gens, qc1max)
        push!(gens, qc2min)
        push!(gens, qc1max)
        push!(gens, ramp_agc)
        push!(gens, ramp_10)
        push!(gens, ramp_30)
        push!(gens, ramp_q)
        push!(gens, apf)
    end
    permutedims(reshape(gens, 21, length(gens) ÷ 21))
end

function get_branch(data)
    branches = Float64[]
    indices = sort(parse.(Int, keys(data["branch"])))
    @assert first(indices) == 1 && last(indices) == length(indices)
    for i in indices
        branch_data = data["branch"][string(i)]
        br_b = 0.0
        push!(branches, branch_data["f_bus"])
        push!(branches, branch_data["t_bus"])
        push!(branches, branch_data["br_r"])
        push!(branches, branch_data["br_x"])
        push!(branches, br_b)
        push!(branches, branch_data["rate_a"])
        push!(branches, branch_data["rate_b"])
        push!(branches, branch_data["rate_c"])
        push!(branches, branch_data["tap"])
        push!(branches, branch_data["shift"])
        push!(branches, branch_data["br_status"])
        push!(branches, branch_data["angmin"])
        push!(branches, branch_data["angmax"])
    end
    permutedims(reshape(branches, 13, length(branches) ÷ 13))
end

function get_gencost(data)
    gencost = Float64[]
    indices = sort(parse.(Int, keys(data["gen"])))
    @assert first(indices) == 1 && last(indices) == length(indices)
    costlength = 0
    for i in indices
        gen_data = data["gen"][string(i)]
        model = gen_data["model"]
        startup = gen_data["startup"] # TODO: double check that this is cost
        shutdown = gen_data["shutdown"] # TODO: double check that this is cost
        ncost = gen_data["ncost"]
        cost = gen_data["cost"]
        costlength = length(cost) # TODO: assert all lengths are the same
        push!(gencost, model)
        push!(gencost, startup)
        push!(gencost, shutdown)
        push!(gencost, ncost)
        for item in cost
            push!(gencost, item)
        end
    end
    permutedims(reshape(gencost, 4+costlength, length(gencost) ÷ (4+costlength)))
end

function get_dcline(data)
    dcline = Float64[]
    indices = sort(parse.(Int, keys(data["dcline"])))
    # TODO: implement getting the data from dcline
    for i in indices
        # TODO: this is not implemented.
        # This code will error if this for loop is entered
        push!(dcline, dcline_data["f_bus"])
        push!(dcline, dcline_data["t_bus"])
        push!(dcline, dcline_data["br_status"])
        push!(dcline, dcline_data["pf"])
        push!(dcline, dcline_data["pt"])
        push!(dcline, dcline_data["qf"])
        push!(dcline, dcline_data["qt"])
        push!(dcline, dcline_data["vf"])
        push!(dcline, dcline_data["vt"])
        push!(dcline, dcline_data["pmin"])
        push!(dcline, dcline_data["pmax"])
        push!(dcline, dcline_data["qminf"])
        push!(dcline, dcline_data["qmaxf"])
        push!(dcline, dcline_data["qmint"])
        push!(dcline, dcline_data["qmaxt"])
        push!(dcline, dcline_data["loss0"])
        push!(dcline, dcline_data["loss1"])
        push!(dcline, dcline_data["loss1"])
    end
    if length(dcline) > 0
        permutedims(reshape(dcline, 17, length(gencost) ÷ 17))
    else
        permutedims(reshape(Float64[], 0, 0))
    end
end

function export_matpower_mat(filename, data)
    out = Dict()
    out["bus"] = get_bus(data)
    out["gen"] = get_gen(data)
    out["branch"] = get_branch(data)
    out["gencost"] = get_gencost(data)
    out["dcline"] = get_dcline(data)
    out["bus_name"] = get_bus_name(data)
    out["baseMVA"] = data["baseMVA"]
    MAT.matwrite(filename, out)
end

writeln(f, s) = begin
    write(f, s)
    write(f, '\n')
end

function export_matpower_m(filename, data)
    out = Dict()
    out["bus"] = get_bus(data)
    out["gen"] = get_gen(data)
    out["branch"] = get_branch(data)
    out["gencost"] = get_gencost(data)
    out["dcline"] = get_dcline(data)
    out["bus_name"] = get_bus_name(data)
    out["baseMVA"] = data["baseMVA"]
    name = replace(filename, ".m" => "", "./" => "")
    open(filename, "w") do f
        writeln(f, "function mpc = $name")
        writeln(f, "")
        writeln(f, "%% MATPOWER Case Format : Version 2")
        writeln(f, "mpc.version = '2';")
        writeln(f, "")
        writeln(f, "%% system MVA base")
        writeln(f, "mpc.baseMVA = $(out["baseMVA"]);")
        writeln(f, "")
        writeln(f, "%% bus data")
        writeln(f, "%	bus_i	type	Pd	Qd	Gs	Bs	area	Vm	Va	baseKV	zone	Vmax	Vmin")
        writeln(f, "mpc.bus = [")
        rows, cols = size(out["bus"])
        for r in 1:rows
            row = out["bus"][r, :]
            for e in row
                write(f, "	$(round(e, sigdigits=4))")
            end
            writeln(f, ";")
        end
        writeln(f, "];")
        writeln(f, "")
        writeln(f, "%% generator data")
        writeln(f, "%	bus	Pg	Qg	Qmax	Qmin	Vg	mBase	status	Pmax	Pmin	Pc1	Pc2	Qc1min	Qc1max	Qc2min	Qc2max	ramp_agc	ramp_10	ramp_30	ramp_q	apf")
        writeln(f, "mpc.gen = [")
        rows, cols = size(out["gen"])
        for r in 1:rows
            row = out["gen"][r, :]
            for e in row
                write(f, "	$(round(e, sigdigits=4))")
            end
            writeln(f, ";")
        end
        writeln(f, "];")
        writeln(f, "")
        writeln(f, "%% branch data")
        writeln(f, "%	fbus	tbus	r	x	b	rateA	rateB	rateC	ratio	angle	status	angmin	angmax")
        writeln(f, "mpc.branch = [")
        rows, cols = size(out["branch"])
        for r in 1:rows
            row = out["branch"][r, :]
            for e in row
                write(f, "	$(round(e, sigdigits=4))")
            end
            writeln(f, ";")
        end
        writeln(f, "];")
        writeln(f, "")
        writeln(f, "%%-----  OPF Data  -----%%")
        writeln(f, "%% generator cost data")
        writeln(f, "%	1	startup	shutdown	n	x1	y1	...	xn	yn")
        writeln(f, "%	2	startup	shutdown	n	c(n-1)	...	c0")
        writeln(f, "mpc.gencost = [")
        rows, cols = size(out["gencost"])
        for r in 1:rows
            row = out["gencost"][r, :]
            for e in row
                write(f, "	$(round(e, sigdigits=4))")
            end
            writeln(f, ";")
        end
        writeln(f, "];")
        writeln(f, "")
        writeln(f, "%% bus names")
        writeln(f, "mpc.bus_name = {")
        for r in out["bus_name"]
            writeln(f, "	\"$(r)\";")
        end
        writeln(f, "};")
    end
    nothing
end

# main("solved_network.nc")
# main("solved_network_RTS_GMLC_base.nc")
# main("solved_network_RTS_GMLC_base+line_expansion.nc")
# main("solved_network_RTS_GMLS_1p5xload+0emission+gen_and_line_expansion.nc")
# main("solved_network_RTS_GMLS_base+gen_and_line_expansion.nc")
# main("solved_network_RTS_GMLS_base+gen_expansion.nc")
