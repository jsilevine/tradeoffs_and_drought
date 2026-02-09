##-------------------------------------------------------------------
## 03_growth_models.jl: Fit growth models and visualize output
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
               #"s(lon, lat, k = 750)",
               "s(lon, k = 20)",
               "s(lat, k = 20)",
               "ba_scaled",
               "drought_strength_scaled", "prop_drought_scaled",
               "beta_pca_scaled",
               "td_pca_scaled",
               "drought_strength_scaled:ba_scaled", "ba_scaled:beta_pca_scaled", "ba_scaled:td_pca_scaled",
               "map_scaled", "mat_scaled",
               "cwm_drought_pc2_scaled", "cwm_growth_pc1_scaled",
               "range_drought_pc2_scaled", "range_growth_pc1_scaled",
               "elev_scaled", "stand_age_scaled"]

growth_model = ModelFitting.bam("growth_pct_scaled", covariates, data; method = "fREML", discrete = true)

R"saveRDS($growth_model, here('data', 'model_objects', 'growth_model.rds'))"

growth_model_summary = ModelFitting.summarize_gam(growth_model)

##-------------------------------------------------------------------
## 02. Generate and save output table
##-------------------------------------------------------------------

## grab model output
ptable = rcopy(growth_model_summary[Symbol("p.table")])
cnames = rcopy(R"colnames"(growth_model_summary[Symbol("p.table")]))
rnames = rcopy(R"rownames"(growth_model_summary[Symbol("p.table")]))

summary_table = DataFrame(hcat(rnames, ptable), vcat(:Coefficient, Symbol.(cnames)))
CSV.write(tabledir("growth_model.csv"), summary_table)

##-------------------------------------------------------------------
## 03. Visualize model predictions for non-smooth covariates
##-------------------------------------------------------------------

## tradeoff characteristics
PlottingUtils.plot_growth_model_predictions(growth_model, data,
    interaction_vars = ["td_pca"],
    vars_to_plot = ["beta_pca", "td_pca", "ba"],
    ncol = 3, align_axes = true,
    save_output = true, output_path = figdir("model_predictions", "growth", "growth_tradeoff_new.pdf"),
    axis_limits = [1, 3.4])

## for drought characteristics
PlottingUtils.plot_growth_model_predictions(growth_model, data,
    vars_to_plot = ["drought_strength", "prop_drought"],
    ncol = 2, align_axes = true,
    save_output = true, output_path = figdir("model_predictions", "growth", "growth_drought.pdf"),
    axis_limits = [1, 3.4])

## for stand characteristics
PlottingUtils.plot_growth_model_predictions(growth_model, data,
    vars_to_plot = ["stand_age", "map", "mat", "elev"],
    ncol = 4, align_axes = true,
    save_output = true, output_path = figdir("model_predictions", "growth", "growth_standchar.pdf"),
    axis_limits = [1, 3.4])

## for traits
PlottingUtils.plot_growth_model_predictions(growth_model, data,
    vars_to_plot = ["cwm_drought_pc2", "cwm_growth_pc1", "range_drought_pc2", "range_growth_pc1"],
    ncol = 4, align_axes = true,
    save_output = true, output_path = figdir("model_predictions", "growth", "growth_traits.pdf"),
    axis_limits = [1, 3.4])

##-------------------------------------------------------------------
## 04. Plot predictions for spatial smooth
##-------------------------------------------------------------------

PlottingUtils.plot_spatial_smooth(growth_model, data,
    save_output = true, output_path = figdir("smooths", "growth_spatial_smooth.pdf"))

##-------------------------------------------------------------------
## 05. Evaluate spatial and temporal autocorrelation,
##     generate diagnostic figures. 
##-------------------------------------------------------------------

## examine temporal autocorrelation
PlottingUtils.plot_acf(growth_model, true, figdir("autocorrelation", "growth_model_acf.pdf"))

## examine residual spatial autocorrelation
PlottingUtils.plot_residual_variogram(growth_model, true, figdir("autocorrelation", "growth_model_variogram.pdf"))



