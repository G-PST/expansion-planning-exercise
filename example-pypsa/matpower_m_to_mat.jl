import MAT

function parse_file(filename)
    data = readlines(filename)
    data = map(line -> strip(line), data)
    data2 = Dict()
    PARSE = false
    PARSECELL = false
    element = nothing
    for line in data
        if startswith(line, "]")
            PARSE = false
        end
        if startswith(line, "}")
            PARSECELL = false
        end

        if PARSE
            push!(data2[element], replace(x -> isnothing(x) ? 0.0 : x, tryparse.(Float64, split(line, "\t"))))
        elseif PARSECELL
            push!(data2[element], replace.(split(line, "\t"), r"[;\']"=>""))
        end

        if startswith(line, "mpc") && endswith(line, "[")
            PARSE = true
            element = replace(first(split(line, " = ")), "mpc."=>"")
            data2[element] = []
        end
        if startswith(line, "mpc") && endswith(line, "{")
            PARSECELL = true
            element = replace(first(split(line, " = ")), "mpc."=>"")
            data2[element] = []
        end
        if startswith(line, "mpc") && endswith(line, ";") && !occursin(r"[\[{]", line)
            expr = split(line, " = ")
            data2[replace(first(expr), "mpc."=>"")] = replace(last(expr), r"[;\']"=>"")
        end
    end
    return data2
end

function convert_matpower(filename)
    @info "parsing" filename
    data2 = parse_file(filename)
    mat_filename = first(splitext(filename)) * ".mat"
    @info "writing mpc to" mat_filename
    MAT.matwrite(mat_filename, data2)
end