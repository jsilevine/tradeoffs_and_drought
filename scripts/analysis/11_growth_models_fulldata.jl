##-------------------------------------------------------------------
## 03_growth_models_fulldata.jl: Sensitivity analysis (reviewer 2, round 2).
## Refits the growth model on ALL grid cell-intervals, including those
## with no drought year, by replacing drought_strength (undefined without
## drought) with min_pdsi (minimum annual PDSI over the interval, defined
## for every interval; more negative = more severe drought, same sign
## convention as drought_strength).
##
## author: jacob levine; jacob.levine@utah.edu
##-------------------------------------------------------------------

include(joinpath(@__DIR__, "..", "..", "src", "julia", "utils.jl"))
using .Utils

using CSV, DataFrames, PrettyTables, Statistics
using RCall

R"""
library(here)
"""
R"install.packages('gamlss')"

include(srcdir("julia", "ModelFitting.jl"))
using .ModelFitting

##-------------------------------------------------------------------
## 01. Prepare data
##-------------------------------------------------------------------

data = CSV.read(datadir("data_for_modeling.csv"), DataFrame, missingstring = "NA");

## drought_strength is intentionally omitted so dropmissing! no longer
## removes non-drought intervals
select!(data, [:year, :lon, :lat,
               :richness_scaled,
               :ba, :ba_scaled,
               :min_pdsi,
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

## guard against NaN values written by 09_assemble_data.jl
filter!(r -> isfinite(r.min_pdsi), data);

data.min_pdsi_scaled = (data.min_pdsi .- mean(data.min_pdsi)) ./ std(data.min_pdsi);

println("N intervals: ", nrow(data),
        "; without drought: ", sum(data.prop_drought .== 0))

##-------------------------------------------------------------------
## 02. Fit model (same structure as 03_growth_models.jl,
##     drought_strength -> min_pdsi)
##-------------------------------------------------------------------

covariates =  ["s(year)",
               "s(lon, k = 20)",
               "s(lat, k = 20)",
               "ba_scaled",
               "min_pdsi_scaled", "prop_drought_scaled",
               "beta_pca_scaled",
               "td_pca_scaled",
               "min_pdsi_scaled:ba_scaled", "ba_scaled:beta_pca_scaled", "ba_scaled:td_pca_scaled",
               "min_pdsi_scaled:td_pca_scaled", "min_pdsi_scaled:beta_pca_scaled",
               "map_scaled", "mat_scaled",
               "cwm_drought_pc2_scaled", "cwm_growth_pc1_scaled",
               "range_drought_pc2_scaled", "range_growth_pc1_scaled",
               "elev_scaled", "stand_age_scaled"]

growth_model = ModelFitting.bam("growth_pct_scaled", covariates, data; method = "fREML", discrete = true)

R"saveRDS($growth_model, here('data', 'model_objects', 'growth_model_fulldata.rds'))"

growth_model_summary = ModelFitting.summarize_gam(growth_model)

##-------------------------------------------------------------------
## 03. Save output table
##-------------------------------------------------------------------

ptable = rcopy(growth_model_summary[Symbol("p.table")])
cnames = rcopy(R"colnames"(growth_model_summary[Symbol("p.table")]))
rnames = rcopy(R"rownames"(growth_model_summary[Symbol("p.table")]))

summary_table = DataFrame(hcat(rnames, ptable), vcat(:Coefficient, Symbol.(cnames)))
summary_table.ci_lower = summary_table.Estimate .- 1.96 .* summary_table[:,"Std. Error"]
summary_table.ci_upper = summary_table.Estimate .+ 1.96 .* summary_table[:,"Std. Error"]
CSV.write(tabledir("growth_model_fulldata.csv"), summary_table)

pretty_table(summary_table)
