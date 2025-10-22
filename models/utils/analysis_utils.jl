
"""
find_start(totals, Nstart

Finds the turn on of the biochemical circuit
"""
function find_start(totals, Nstart)
    for i in eachindex(totals)
        if totals[i] >= Nstart
            return i
        end
    end
end

function var_distribution(results::SimResults, varname::String, cellfate::Int64, Nstart::Int64)
    points = Vector{Float64}()
    start  = find_start(results.totals, Nstart)
    var = results.vars[findfirst(results.varsnames .== varname)][:, start:end]
    cfates = results.CFATES[:, start:end]
    for cell=1:size(var)[1]
        pps = Vector{Float64}()
        for t=1:size(var)[2]
            if var[cell,t] != 0
                if cfates[cell,t] == cellfate
                    push!(pps, var[cell,t])
                end
            end
        end
        if length(pps)!=0
            push!(points, mean(pps))
        end
    end
    return points
end
