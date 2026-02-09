##-------------------------------------------------------------------
## Test sensitivity of growth model to presence of trait range variables
##
## author: jacob levine; jacob.levine@utah.edu
##-------------------------------------------------------------------

##-------------------------------------------------------------------
## 00. load libraries and includes
##-------------------------------------------------------------------

## basic file utilites
include(joinpath(@__DIR__, "..", "..", "src", "julia", "utils.jl"))

using .Utils

## libraries
using CSV, DataFrames, PrettyTables
using RCall

R"""
library(here)
"""

## includes
include(srcdir("julia", "ModelFitting.jl"))
include(srcdir("julia", "PlottingUtils.jl"))

using .ModelFitting
using .PlottingUtils

##-------------------------------------------------------------------
## 01. Fit model
##-------------------------------------------------------------------

data = CSV.read(datadir("data_for_modeling.csv"), DataFrame, missingstring = "NA");

select!(data, [:year, :lon, :lat,
               :richness_scaled,
               :ba, :ba_scaled,
               :drought_strength, :drought_strength_scaled,
               :prop_drought, :prop_drought_scaled,
               :beta_pca, :beta_pca_scaled,
               :td_pca, :td_pca_scaled,
               :growth_pct, :growth_pct_scaled,
               :map, :map_scaled,
               :mat, :mat_scaled,
               :cwm_drought_pc2, :cwm_drought_pc2_scaled,
               :cwm_growth_pc1, :cwm_growth_pc1_scaled,
               :range_drought_pc2, :range_drought_pc2_scaled,
               :range_growth_pc1, :range_growth_pc1_scaled,
               :elev, :elev_scaled,
               :stand_age, :stand_age_scaled]); 

dropmissing!(data);

covariates =  ["s(year)",
               "s(lon, k = 20)",
               "s(lat, k = 20)",
               "ba_scaled",
               "drought_strength_scaled", "prop_drought_scaled",
               "beta_pca_scaled",
               "td_pca_scaled",
               "drought_strength_scaled:ba_scaled", "ba_scaled:beta_pca_scaled", "ba_scaled:td_pca_scaled",
               "map_scaled", "mat_scaled",
               "cwm_drought_pc2_scaled", "cwm_growth_pc1_scaled",
               "elev_scaled", "stand_age_scaled"]

growth_model_norange = ModelFitting.bam("growth_pct_scaled", covariates, data; method = "fREML", discrete = true)
growth_model_summary_norange = ModelFitting.summarize_gam(growth_model_norange)

ptable = rcopy(growth_model_summary_norange[Symbol("p.table")])
cnames = rcopy(R"colnames"(growth_model_summary_norange[Symbol("p.table")]))
rnames = rcopy(R"rownames"(growth_model_summary_norange[Symbol("p.table")]))

summary_table = DataFrame(hcat(rnames, ptable), vcat(:Coefficient, Symbol.(cnames)))
CSV.write(tabledir("growth_model_norange.csv"), summary_table)

## no tradeoffs

covariates =  ["s(year)",
               "s(lon, k = 20)",
               "s(lat, k = 20)",
               "ba_scaled",
               "drought_strength_scaled", "prop_drought_scaled",
               "drought_strength_scaled:ba_scaled",
               "map_scaled", "mat_scaled",
               "cwm_drought_pc2_scaled", "cwm_growth_pc1_scaled",
               "range_drought_pc2_scaled", "range_growth_pc1_scaled",
               "elev_scaled", "stand_age_scaled"]

growth_model_notradeoff = ModelFitting.bam("growth_pct_scaled", covariates, data; method = "fREML", discrete = true)
growth_model_summary_notradeoff = ModelFitting.summarize_gam(growth_model_notradeoff)

ptable = rcopy(growth_model_summary_notradeoff[Symbol("p.table")])
cnames = rcopy(R"colnames"(growth_model_summary_notradeoff[Symbol("p.table")]))
rnames = rcopy(R"rownames"(growth_model_summary_notradeoff[Symbol("p.table")]))

summary_table = DataFrame(hcat(rnames, ptable), vcat(:Coefficient, Symbol.(cnames)))
CSV.write(tabledir("growth_model_notradeoff.csv"), summary_table)
