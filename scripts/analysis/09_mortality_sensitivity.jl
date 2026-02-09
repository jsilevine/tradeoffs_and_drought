##-------------------------------------------------------------------
## Test sensitivity of mortality model to presence of trait range variables
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
## 01. NoRange
##-------------------------------------------------------------------

mu_fixed = [
    "pb(year_centered)",
    "pb(lon, df = 18)",
    "pb(lat, df = 18)",
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
mortality_model_norange = ModelFitting.bezi_gam("mort_tpa_pct", mu_fixed, nu_fixed, data; n_cyc = 30)

R"mm_s <- summary($mortality_model_norange)"
mortality_summary_norange = DataFrame(rcopy(R"mm_s"), :auto)
rename!(mortality_summary_norange, rcopy(R"colnames(mm_s)"))
rename!(mortality_summary_norange, Symbol("Pr(>|t|)") => :p_value,
    Symbol("Std. Error") => :se,
    Symbol("t value") => :t_value,
    Symbol("Estimate") => :estimate)

mortality_summary_norange.variable = rcopy(R"rownames(mm_s)")

mortality_summary_norange.component = Vector{String}(undef, nrow(mortality_summary_norange))
flag = 0
for row in eachrow(mortality_summary_norange)
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

mortality_summary_norange = mortality_summary_norange[:, [:variable, :component, :estimate, :se, :t_value, :p_value]]

CSV.write(tabledir("mortality_model_summary_norange.csv"), mortality_summary_norange)


##-------------------------------------------------------------------
## 01. NoTradeoff
##-------------------------------------------------------------------

mu_fixed = [
    "pb(year_centered)",
    "pb(lon, df = 18)",
    "pb(lat, df = 18)",
    "ba_scaled",
    "drought_strength_scaled",
    "drought_strength_scaled:ba_scaled",
    "prop_drought_scaled",
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
mortality_model_notradeoff = ModelFitting.bezi_gam("mort_tpa_pct", mu_fixed, nu_fixed, data; n_cyc = 30)

R"mm_s <- summary($mortality_model_notradeoff)"
mortality_summary_notradeoff = DataFrame(rcopy(R"mm_s"), :auto)
rename!(mortality_summary_notradeoff, rcopy(R"colnames(mm_s)"))
rename!(mortality_summary_notradeoff, Symbol("Pr(>|t|)") => :p_value,
    Symbol("Std. Error") => :se,
    Symbol("t value") => :t_value,
    Symbol("Estimate") => :estimate)

mortality_summary_notradeoff.variable = rcopy(R"rownames(mm_s)")

mortality_summary_notradeoff.component = Vector{String}(undef, nrow(mortality_summary_notradeoff))
flag = 0
for row in eachrow(mortality_summary_notradeoff)
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

mortality_summary_notradeoff = mortality_summary_notradeoff[:, [:variable, :component, :estimate, :se, :t_value, :p_value]]

CSV.write(tabledir("mortality_model_summary_notradeoff.csv"), mortality_summary_notradeoff)
