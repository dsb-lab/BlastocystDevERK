
mutable struct Positions
    X::Matrix{Float64}
    Y::Matrix{Float64}
    Z::Matrix{Float64}
end

mutable struct SimResults
    times::Vector{Float64}
    fDP::Vector{Float64}
    fEPI::Vector{Float64}
    fPRE::Vector{Float64}
    totals::Vector{Int64}
    vars::Vector{Matrix{Float64}}
    comvar::Matrix{Float64}
    varsnames::Vector{String}
    comvarname::String
    positions::Positions
    R::Matrix{Float64}
    CFATES::Matrix{Int64}
    tdivs::Vector{Vector{Float64}}
    tdivs_idx::Vector{Vector{Int64}}
end

mutable struct test
    a::Vector{Float64}
end