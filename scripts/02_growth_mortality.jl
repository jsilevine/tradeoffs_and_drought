##---------------------------------------------------------------
## 02_growth_mortality.jl: Calculate basal area growth and mortality rates
##
## author: jacob levine; jacob.levine@utah.edu
##---------------------------------------------------------------

include(joinpath(@__DIR__, "..", "src", "julia", "utils.jl"))

using CSV, DataFrames, DataFramesMeta
using Statistics
using LightGraphs, ProgressBars
using BenchmarkTools, Infiltrator

include(srcdir("julia", "FIAgm", "FIAgm.jl"))

using .FIAgm

states = CSV.read(datadir("state_list.csv"), DataFrame)

for s in states.state[1:nrow(states)]
    println("Processing: $s")
    process_state(convert(String, s), true)
end

## assemble into single file
gm = CSV.read(datadir("gm", "gm_" * states.state[1] * ".csv"), DataFrame)
gm = gm[gm.ba_init .> 0.0, :]

for s in states.state[2:nrow(states)]
    println("Processing: $s")
    local new_gm = CSV.read(datadir("gm", "gm_" * s * ".csv"), DataFrame)
    local new_gm = new_gm[new_gm.ba_init .> 0.0, :]
    global gm = vcat(gm, new_gm)
end

CSV.write(datadir("gm", "gm_full.csv"), gm)

