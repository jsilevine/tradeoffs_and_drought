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

data = CSV.read(datadir("full_data.csv"), DataFrame, missingstring = "NA");

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

CSV.write(datadir("data_for_modeling.csv"), filtered_data, missingstring = "NA");

##-------------------------------------------------------------------
## 02. Prep data for trend analyses
##-------------------------------------------------------------------

trend_data = select(filtered_data, [:grid_cell_id, :lon, :lat,
                                    :year, :prev_measyear,
                                    :map, :mat, :stand_age, :ba, :richness,
                                    :map_scaled, :mat_scaled, :stand_age_scaled, :ba_scaled, :richness_scaled,
                                    :beta_pca, :se_pca, :cwm_drought_pc2, :cwm_growth_pc1,
                                    :beta_pca_current, :se_pca_current,
                                    :cwm_drought_pc2_current, :cwm_growth_pc1_current]);

results = DataFrame[];

for b in unique(trend_data.grid_cell_id)

    subdata = trend_data[trend_data.grid_cell_id .== b, :]
    sort!(subdata, :year)

    final_row = DataFrame(subdata[end, :])
    select!(final_row, Not([:beta_pca, :se_pca, :cwm_drought_pc2, :cwm_growth_pc1]))
    rename!(final_row, :beta_pca_current => :beta_pca, :se_pca_current => :se_pca,
            :cwm_drought_pc2_current => :cwm_drought_pc2, :cwm_growth_pc1_current => :cwm_growth_pc1)

    subdata = subdata[1:end-1, :]
    subdata[!, :year] = subdata[:, :prev_measyear]
    select!(subdata, Not([:beta_pca_current, :se_pca_current, :cwm_drought_pc2_current, :cwm_growth_pc1_current]))

    out = vcat(subdata, final_row)
    push!(results, out)
    
end

trend_data_long = vcat(results...);
trend_data_long[!, :year_centered] = trend_data_long.year .- mean(trend_data_long.year);

CSV.write(datadir("trend_data.csv"), trend_data_long, missingstring = "NA");
