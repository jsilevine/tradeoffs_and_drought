##---------------------------------------------------------------
## 01_prep_data.jl: Prepare data for modeling and trend analyses
##
## author: jacob levine; jacob.levine@utah.edu
##---------------------------------------------------------------

##---------------------------------------------------------------
## 00. load data and libraries
##---------------------------------------------------------------

include(joinpath(@__DIR__, "..", "..", "src", "julia", "utils.jl"))

using .Utils

using CSV, DataFrames
using Statistics

##-------------------------------------------------------------------
## 01. Prep data for modeling routines
##-------------------------------------------------------------------

data = CSV.read(datadir("trait_sensitivity", "full_data.csv"), DataFrame, missingstring = "NA");

filtered_data = data[abs.(data.beta_pca) .< 5 .&& abs.(data.growth_pct) .< 0.15, :];

filtered_data.range_drought_pc2 = filtered_data.max_drought_pc2 .- filtered_data.min_drought_pc2;
filtered_data.range_growth_pc1 = filtered_data.max_growth_pc1 .- filtered_data.min_growth_pc1;

variables_to_scale = [
  "ba", "beta_pca", "tpa", "td_pca", "se_pca", "r2_pca", "drought_strength", "prop_drought",
  "nondrought_pdsi", "growth_pct", "map", "mat",
  "cwm_drought_pc2", "cwm_growth_pc1", "elev",
  "richness", "max_drought_pc2", "min_drought_pc2",
  "max_growth_pc1", "min_growth_pc1", "range_drought_pc2",
  "range_growth_pc1", "stand_age"
]

function scale_with_missing(x)
    μ = mean(skipmissing(x))
    σ = std(skipmissing(x))
    return [(ismissing(xi) ? missing : (xi - μ) / σ) for xi in x]
end

for var in variables_to_scale
    println("Scaling variable: $var")
    filtered_data[!, var * "_scaled"] = scale_with_missing(filtered_data[:, Symbol(var)])
end

rename!(filtered_data, :measyear => :year);
filtered_data[!, :year_centered] = filtered_data.year .- mean(filtered_data.year);

CSV.write(datadir("trait_sensitivity", "data_for_modeling.csv"), filtered_data, missingstring = "NA");
