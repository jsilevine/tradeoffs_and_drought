##-------------------------------------------------------------------
## 04_mortality_models_fulldata.jl: Sensitivity analysis (reviewer 2, round 2).
## Refits the mortality BEZI model on ALL grid cell-intervals, including
## those with no drought year, by replacing drought_strength (undefined
## without drought) with min_pdsi (minimum annual PDSI over the interval;
## more negative = more severe drought, same sign as drought_strength).
##
## author: jacob levine; email: jacob.levine@utah.edu
##-------------------------------------------------------------------

include(joinpath(@__DIR__, "..", "..", "src", "julia", "utils.jl"))
using .Utils

using DataFrames, CSV, Statistics
using RCall

include(srcdir("julia", "ModelFitting.jl"))
using .ModelFitting

##-------------------------------------------------------------------
## 00. Prepare data
##-------------------------------------------------------------------

data = CSV.read(datadir("data_for_modeling.csv"), DataFrame, missingstring = "NA");

## drought_strength is intentionally omitted so dropmissing! no longer
## removes non-drought intervals
select!(data, [:year, :year_centered,
               :lon, :lat,
               :richness_scaled,
               :ba, :ba_scaled,
               :tpa, :tpa_scaled,
               :min_pdsi,
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

## guard against NaN values written by 09_assemble_data.jl
filter!(r -> isfinite(r.min_pdsi), data);

data.min_pdsi_scaled = (data.min_pdsi .- mean(data.min_pdsi)) ./ std(data.min_pdsi);

println("N intervals: ", nrow(data),
        "; without drought: ", sum(data.prop_drought .== 0))

##-------------------------------------------------------------------
## 01. Fit mortality GAM (same structure as 04_mortality_models.jl,
##     drought_strength -> min_pdsi)
##-------------------------------------------------------------------

mu_fixed = [
    "pb(year_centered)",
    "pb(lon, df = 18)",
    "pb(lat, df = 18)",
    "ba_scaled",
    "min_pdsi_scaled",
    "min_pdsi_scaled:ba_scaled",
    "prop_drought_scaled",
    "beta_pca_scaled",
    "td_pca_scaled",
    "ba_scaled:beta_pca_scaled",
    "ba_scaled:td_pca_scaled",
    "min_pdsi_scaled:beta_pca_scaled",
    "min_pdsi_scaled:td_pca_scaled",
    "map_scaled", "mat_scaled",
    "cwm_drought_pc2_scaled", "cwm_growth_pc1_scaled",
    "range_drought_pc2_scaled", "range_growth_pc1_scaled",
    "elev_scaled", "stand_age_scaled"
];

nu_fixed = [
    "pb(year_centered, df = 10)",
    "pb(lon, df = 18)",
    "pb(lat, df = 18)",
    "ba_scaled",
    "min_pdsi_scaled:ba_scaled",
    "min_pdsi_scaled",
    "prop_drought_scaled",
    "map_scaled", "mat_scaled",
    "elev_scaled", "stand_age_scaled"
];

## slow: ~1.8x the rows of the original fit
mortality_model = ModelFitting.bezi_gam("mort_tpa_pct", mu_fixed, nu_fixed, data; n_cyc = 30)

R"saveRDS($mortality_model, file = here::here('data', 'model_objects', 'mortality_model_fulldata.rds'))"

##-------------------------------------------------------------------
## 02. Save summary table
##-------------------------------------------------------------------

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
    global flag
    if row.variable == "(Intercept)"
        flag += 1
    end
    row.component = flag == 1 ? "mu" : (flag == 2 ? "link" : "nu")
end

mortality_summary = mortality_summary[:, [:variable, :component, :estimate, :se, :t_value, :p_value]]
mortality_summary.ci_lower = mortality_summary.estimate .- 1.96 .* mortality_summary.se
mortality_summary.ci_upper = mortality_summary.estimate .+ 1.96 .* mortality_summary.se

CSV.write(tabledir("mortality_model_summary_fulldata.csv"), mortality_summary)

show(mortality_summary[occursin.("min_pdsi", mortality_summary.variable), :], allrows = true)
