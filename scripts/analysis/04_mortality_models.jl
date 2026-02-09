##-------------------------------------------------------------------
## 04_mortality_models.jl: Fit mortality models to data
##
## author: jacob levine; email: jacob.levine@utah.edu
##-------------------------------------------------------------------

##-------------------------------------------------------------------
## 00. Load dependencies and data
##-------------------------------------------------------------------

include(joinpath(@__DIR__, "..", "..", "src", "julia", "utils.jl"))

using .Utils

using DataFrames, CSV
using RCall
using Revise

include(srcdir("julia", "ModelFitting.jl"))
include(srcdir("julia", "PlottingUtils.jl"))

using .ModelFitting
using .PlottingUtils

data = CSV.read(datadir("data_for_modeling.csv"), DataFrame, missingstring = "NA");

select!(data, [:year, :year_centered,
               :lon, :lat,
               :richness_scaled,
               :ba, :ba_scaled,
               :tpa, :tpa_scaled,
               :drought_strength, :drought_strength_scaled,
               :prop_drought, :prop_drought_scaled,
               :beta_pca, :beta_pca_scaled,
               :td_pca, :td_pca_scaled,
               :growth_pct, :growth_pct_scaled,
               :mort_tpa_pct,
               :map, :map_scaled,
               :mat, :mat_scaled,
               :cwm_drought_pc2, :cwm_drought_pc2_scaled,
               :cwm_growth_pc1, :cwm_growth_pc1_scaled,
               :range_drought_pc2, :range_drought_pc2_scaled,
               :range_growth_pc1, :range_growth_pc1_scaled,
               :elev, :elev_scaled,
               :stand_age, :stand_age_scaled]); 

dropmissing!(data);

##-------------------------------------------------------------------
## 01. Fit mortality GAM
##-------------------------------------------------------------------

mu_fixed = [
    "pb(year_centered)",
    "pb(lon, df = 18)",
    "pb(lat, df = 18)",
    ## "pb(year_centered):NA_L2NAME", ## can include these if interested in accounting for by-group autocorrelation
    ## however, temporal autocorrelation is small already (lag-1 = 0.06)
    "ba_scaled",
    "drought_strength_scaled",
    "drought_strength_scaled:ba_scaled",
    "prop_drought_scaled",
    "beta_pca_scaled",
    "td_pca_scaled",
    "ba_scaled:beta_pca_scaled",
    "ba_scaled:td_pca_scaled",
    "map_scaled", "mat_scaled",
    "cwm_drought_pc2_scaled", "cwm_growth_pc1_scaled",
    "range_drought_pc2_scaled", "range_growth_pc1_scaled",
    "elev_scaled", "stand_age_scaled"
];

nu_fixed = [
    "pb(year_centered, df = 10)",
    "pb(lon, df = 18)",
    "pb(lat, df = 18)",
    ## "pb(year_centered):NA_L2NAME", ## see above
    "ba_scaled",
    "drought_strength_scaled:ba_scaled",
    "drought_strength_scaled",
    "prop_drought_scaled",
    "map_scaled", "mat_scaled",
    "elev_scaled", "stand_age_scaled"
];

## this takes forever, so only run if you need to refit the model
mortality_model = ModelFitting.bezi_gam("mort_tpa_pct", mu_fixed, nu_fixed, data; n_cyc = 30)

R"saveRDS($mortality_model, file = here::here('data', 'model_objects', 'mortality_model.rds'))"

##-------------------------------------------------------------------
## 02. generate and save mortality model summary
##-------------------------------------------------------------------

if isdefined(Base,:mortality_model) == false
    mortality_model = R"readRDS(here::here('data', 'model_objects', 'mortality_model.rds'))"
end

R"mm_s <- summary($mortality_model)"
mortality_summary = DataFrame(rcopy(R"mm_s"), :auto)
rename!(mortality_summary, rcopy(R"colnames(mm_s)"))
rename!(mortality_summary, Symbol("Pr(>|t|)") => :p_value,
    Symbol("Std. Error") => :se,
    Symbol("t value") => :t_value,
    Symbol("Estimate") => :estimate)

mortality_summary.variable = rcopy(R"rownames(mm_s)")

mortality_summary.component = Vector{String}(undef, nrow(mortality_summary))
flag = 0
for row in eachrow(mortality_summary)
    if row.variable == "(Intercept)"
        flag += 1
    end

    if flag == 1
        row.component = "mu"
    elseif flag == 2
        row.component = "link"
    else
        row.component = "nu"
    end
end

mortality_summary = mortality_summary[:, [:variable, :component, :estimate, :se, :t_value, :p_value]]

CSV.write(tabledir("mortality_model_summary.csv"), mortality_summary)


##-------------------------------------------------------------------
## 03. Plot mortality model results
##-------------------------------------------------------------------

## tradeoff characteristics
PlottingUtils.plot_mortality_model_predictions(mortality_model, data,
    interaction_vars = [],
    vars_to_plot = ["beta_pca", "td_pca", "ba"],
    ncol = 3, align_axes = true,
    save_output = true, output_path = figdir("model_predictions", "mortality", "mortality_tradeoff.pdf"),
    axis_limits = [0.0058, 0.00725])

## for drought characteristics
PlottingUtils.plot_mortality_model_predictions(mortality_model, data,
    vars_to_plot = ["drought_strength", "prop_drought"],
    ncol = 2, align_axes = true,
    save_output = true, output_path = figdir("model_predictions", "mortality", "mortality_drought.pdf"),
    axis_limits = [0.0058, 0.00725])

## for stand characteristics
PlottingUtils.plot_mortality_model_predictions(mortality_model, data,
    vars_to_plot = ["stand_age", "map", "mat", "elev"],
    ncol = 4, align_axes = true,
    save_output = true, output_path = figdir("model_predictions", "mortality", "mortality_standchar.pdf"),
    axis_limits = [0.0058, 0.00725])

## for traits
PlottingUtils.plot_mortality_model_predictions(mortality_model, data,
    vars_to_plot = ["cwm_drought_pc2", "cwm_growth_pc1", "range_drought_pc2", "range_growth_pc1"],
    ncol = 4, align_axes = true,
    save_output = true, output_path = figdir("model_predictions", "mortality", "mortality_traits.pdf"),
    axis_limits = [0.0058, 0.00725])

## will include the following plots in text figure: beta_pca, se_pca, ba, drought_strength, cwm_drought_pc2, cwm_growth_pc1

##-------------------------------------------------------------------
## 04. Plot predictions for spatial smooth
##-------------------------------------------------------------------

PlottingUtils.plot_spatial_smooth_bezi(mortality_model, data,
    save_output = true, output_path = figdir("smooths", "mortality_spatial_smooth.pdf"))

##-------------------------------------------------------------------
## 05. Evaluate spatial and temporal autocorrelation,
##     generate diagnostic figures. 
##-------------------------------------------------------------------

## examine temporal autocorrelation
PlottingUtils.plot_acf(mortality_model, true, figdir("autocorrelation", "mortality_model_acf.pdf"))

## only if fit interaction between year_centered and ecoregion
## PlottingUtils.plot_group_acf(mortality_model, group_var = "NA_L2NAME", true, figdir("autocorrelation", "mortality_model_group_acf.pdf"))

## examine residual spatial autocorrelation
PlottingUtils.plot_residual_variogram(mortality_model, true, figdir("autocorrelation", "mortality_model_variogram.pdf"))
