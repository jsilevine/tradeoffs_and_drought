## module to process fia data into mortality and growth estimates
module FIAgm

include(joinpath(@__DIR__, "..", "..", "src", "julia", "utils.jl"))

using DataFrames, CSV, Infiltrator, ProgressBars

include(srcdir("julia", "FIAutils.jl"))
using .FIAutils

include(srcdir("julia", "FIAgm", "utils.jl"))
include(srcdir("julia", "FIAgm", "trees.jl"))
include(srcdir("julia", "FIAgm", "stand_properties.jl"))
include(srcdir("julia", "FIAgm", "mortality.jl"))
include(srcdir("julia", "FIAgm", "process_data.jl"))

## exports
export pull_plot, match_condids, TreeList, pull_trees
export TPA, calc_tpa, BA, calc_ba
export Mortality, calc_mortality
export process_growth_and_mortality, process_state

end

