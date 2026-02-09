##-------------------------------------------------------------------
## 09_assemble_data.jl: Assemble complete dataset
##
## author: jacob levine; jacob.levine@utah.edu
##-------------------------------------------------------------------

include(joinpath(@__DIR__, "..", "src", "julia", "utils.jl"))

using CSV, DataFrames
using NCDatasets
using ProgressBars, BenchmarkTools
using Plots
using Statistics
using Downloads, Printf, Serialization

using .Utils

##---------------------------------------------------------------
## 00. Load modules and data
##---------------------------------------------------------------

include(srcdir("julia", "FIAutils.jl"))

using .FIAutils

states = CSV.read(datadir("state_list.csv"), DataFrame);

tradeoffs = CSV.read(datadir("tradeoffs", "tradeoffs_full.csv"), DataFrame);

## merge in tradeoff metrics for current measyear, will use to understand trends in these metrics over time
tradeoffs_current = CSV.read(datadir("tradeoffs", "tradeoffs_current_full.csv"), DataFrame);
select!(tradeoffs_current, [:PLOT_ID, :MEASYEAR, :beta_pca, :se_pca, :cwm_drought_pc2, :cwm_growth_pc1]);

exclude = ["PLOT_ID", "MEASYEAR"];
rn_dict = Dict(name => Symbol(string(name, "_current")) for name in names(tradeoffs_current) if !(name in exclude));
rename!(tradeoffs_current, rn_dict);

leftjoin!(tradeoffs, unique(tradeoffs_current), on = [:PLOT_ID, :MEASYEAR]);

gm_grid = CSV.read(datadir("gm", "gm_full_grid.csv"), DataFrame);

filter!(row -> isfinite(row.growth_pct) && isfinite(row.mort_ba_pct), gm_grid);

grids = CSV.read(datadir("templates", "grid_cells.csv"), DataFrame);
plots_to_cells = CSV.read(datadir("templates", "plot_to_cell.csv"), DataFrame);
plots_to_cells = unique(plots_to_cells[:,2:3], :PLOT_ID); ## remove duplicates

## join pull grid information into tradeoffs
leftjoin!(tradeoffs, plots_to_cells, on = :PLOT_ID);

grid_years = unique(select(gm_grid, [:grid_cell_id, :MEASYEAR]), [:grid_cell_id, :MEASYEAR]);

leftjoin!(grid_years, grids[:,[:LAT, :LON, :grid_cell_id]], on = :grid_cell_id);

select!(grid_years, :LON, :LAT, :grid_cell_id, :MEASYEAR);

##---------------------------------------------------------------
## 01. Assign previous measurement year to each plot/grid
##---------------------------------------------------------------

## first for plots
tradeoffs.PREV_MEASYEAR = Vector{Int64}(undef, nrow(tradeoffs));
for s in eachrow(states)

    println("processing: " * s.state)

    plot_data = FIAutils.read_fia(String(s.state), "PLOT")

    tradeoff_sub = tradeoffs[tradeoffs.STATECD .== s.STATECD,:]

    for i in ProgressBar(1:nrow(tradeoff_sub))
        tradeoff_sub[i, :PREV_MEASYEAR] = plot_data[plot_data.CN .== tradeoff_sub[i,:PREV_PLT_CN], :MEASYEAR][1]
    end

    tradeoffs[tradeoffs.STATECD .== s.STATECD, :PREV_MEASYEAR] .= tradeoff_sub.PREV_MEASYEAR

end

## remove rows that have no tie to grid cell
tradeoffs = tradeoffs[.!ismissing.(tradeoffs.grid_cell_id), :];

## then for grids, using plot-level information
grid_years.PREV_MEASYEAR = Vector{Int64}(undef, nrow(grid_years));
mismatch = 0;
for r in ProgressBar(eachrow(grid_years))
    plot_subset = tradeoffs[tradeoffs.grid_cell_id .== r.grid_cell_id .&& tradeoffs.MEASYEAR .== r.MEASYEAR, :]

    if length(unique(plot_subset.PREV_MEASYEAR)) > 1
        mismatch += 1;
    end
    r.PREV_MEASYEAR = round(mean(plot_subset.PREV_MEASYEAR))
    
end

##---------------------------------------------------------------
## 02. Average across plots within each grid cell
##---------------------------------------------------------------

function avg_slope_se(beta, se, ba)

    # keep valid entries
    keep = .!ismissing.(beta) .& .!ismissing.(se) .& .!ismissing.(ba) .&
           isfinite.(beta) .& isfinite.(se) .& isfinite.(ba) .&
           (se .> 0) .& (ba .> 0)

    beta = beta[keep]
    se   = se[keep]
    ba   = ba[keep]

    n = length(beta)

    ## no valid plots
    if n == 0
        return missing, missing
    elseif n == 1
        return disallowmissing(beta)[1], disallowmissing(se)[1]
    end

    ## sampling variances and inverse-variance weights
    vi = se .^ 2
    w_iv = 1 ./ vi

    ## IV-weighted mean used only for heterogeneity estimation
    beta_iv_mean = sum(w_iv .* beta) / sum(w_iv)

    ## Cochran's Q and DerSimonian-Laird estimator for tau2
    Q = sum(w_iv .* (beta .- beta_iv_mean).^2)
    c = sum(w_iv) - sum(w_iv .^ 2) / sum(w_iv)
    tau2 = (c > 0) ? max(0, (Q - (length(beta)-1)) / c) : 0
    
    ## BA adjusted weights, including heterogeneity
    w_final = ba ./ (vi .+ tau2)
    w_final_notau = ba ./ vi
    beta_hat = sum(w_final_notau .* beta) / sum(w_final_notau)
    se_hat = sqrt(1 / sum(w_final))

    ## variance of weighted mean using normalized weights
    w_norm = w_final_notau ./ sum(w_final_notau)
    vi_tau = vi .+ tau2
    var_beta = sum((w_norm .^ 2) .* vi_tau)
    se_hat = sqrt(sum((w_norm .^ 2) .* vi_tau))
    
    return beta_hat, se_hat
    
end

function ba_weighted_avg(var, ba)
    # Create a mask for complete cases (non-missing entries in both inputs)
    keep = .!(ismissing.(var) .| ismissing.(ba))

    var = var[keep]
    ba = ba[keep]

    ## no valid data
    if isempty(var)
        return missing
    end

    # Compute weighted mean
    return sum(ba .* var) / sum(ba)
end

## first average covariates for grid_years
to_avg = [:tpa_init_subset, :ba_init, :growth_ba, :growth_pct, :mort_ba_pct, :mort_tpa_pct, :ELEV, :STDAGE]

## initialize columns in grid_years for averaged values
for col in to_avg
    grid_years[!, col] = Vector{Union{Missing, Float64}}(undef, nrow(grid_years))
end

## average across plots within each grid cell
for r in ProgressBar(eachrow(grid_years))
    plots_subset = tradeoffs[tradeoffs.grid_cell_id .== r.grid_cell_id .&&
        tradeoffs.MEASYEAR .== r.MEASYEAR, :]
    for col in to_avg
       r[col] = ba_weighted_avg(plots_subset[:, col], plots_subset[:, :ba_init])
    end
end

## next average key variables for tradeoffs to grid level

## create template for averaged slopes and standard errors
tradeoffs_grid = select(grid_years, [:LON, :LAT, :grid_cell_id, :MEASYEAR, :PREV_MEASYEAR, :tpa_init_subset, :ba_init,
                                           :growth_ba, :growth_pct, :mort_ba_pct, :mort_tpa_pct, :ELEV, :STDAGE]);

to_avg = [:richness, :cwm_Amax, :cwm_p50, :cwm_growth_pc1, :cwm_drought_pc2,
          :min_Amax, :min_p50, :min_growth_pc1, :min_drought_pc2,
          :max_Amax, :max_p50, :max_growth_pc1, :max_drought_pc2,
          :td_base, :td_pca, :td_pca_p50,
          :r2_base, :r2_pca, :r2_pca_p50,
          :os_r2_base, :os_r2_pca, :os_r2_pca_p50,
          :cwm_drought_pc2_current, :cwm_growth_pc1_current];

## initialize averaged fields
for col in to_avg
    tradeoffs_grid[!, col] = Vector{Union{Missing, Float64}}(undef, nrow(grid_years))
end

## initialize slope and SE by model type
for type in [:base, :pca, :pca_p50]
    tradeoffs_grid[!, Symbol(string("beta_", type))] = Vector{Union{Missing, Float64}}(undef, nrow(grid_years))
    tradeoffs_grid[!, Symbol(string("se_", type))] = Vector{Union{Missing, Float64}}(undef, nrow(grid_years))
    tradeoffs_grid[!, Symbol(string("td_se_", type))] = Vector{Union{Missing, Float64}}(undef, nrow(grid_years))
end

tradeoffs_grid[!, :beta_pca_current] = Vector{Union{Missing, Float64}}(undef, nrow(grid_years))
tradeoffs_grid[!, :se_pca_current] = Vector{Union{Missing, Float64}}(undef, nrow(grid_years))

for r in ProgressBar(eachrow(tradeoffs_grid))
    plots_subset = tradeoffs[tradeoffs.grid_cell_id .== r.grid_cell_id .&& tradeoffs.MEASYEAR .== r.MEASYEAR, :]

    ## average beta and se for each model
    for type in [:base, :pca, :pca_p50]
        r[Symbol(string("beta_", type))], r[Symbol(string("se_", type))] =
            test = avg_slope_se(
                plots_subset[:, Symbol(string("beta_", type))],
                plots_subset[:, Symbol(string("se_", type))],
                plots_subset[:, :ba_init]
        )
    end

    ## current-year PCA slope and se
    r[:beta_pca_current], r[:se_pca_current] =
        avg_slope_se(
            plots_subset[:, :beta_pca_current],
            plots_subset[:, :se_pca_current],
            plots_subset[:, :ba_init]
        )

    ## compute grid-level adherence (TD) standard errors
    for type in (:base, :pca, :pca_p50)
        td_col = Symbol("td_" * string(type))
        se_col = Symbol("td_se_" * string(type))

        td = plots_subset[:, td_col]
        se = plots_subset[:, se_col]
        ba = plots_subset[:, :ba_init]

        ## keep valid entries only
        mask = .!ismissing.(td) .& .!ismissing.(se) .& .!ismissing.(ba) .&
               (se .> 0) .& (ba .> 0)

        ## no valid data
        if count(mask) == 0
            r[se_col] = missing
            continue
        elseif count(mask) == 1
            r[se_col] = se[findfirst(mask)]
            continue
        end
        
        td = td[mask]
        se = se[mask]
        ba = ba[mask]

        ## total basal area
        A = sum(ba)

        ## BA-weighted mean
        μ = sum(ba .* td) / A

        ## variance components
        VarA = sum((ba ./ A).^2 .* se.^2)

        ## between-plot variance component
        τ2 = max(0, sum(ba .* (td .- μ).^2) / A - sum(ba .* se.^2) / A)

        ## total variance
        VarB = sum((ba ./ A).^2 .* τ2)

        ## combined SE
        r[se_col] = sqrt(VarA + VarB)
        
    end
    
    for col in to_avg
        r[col] = ba_weighted_avg(plots_subset[:, col], plots_subset[:, :ba_init])

    end
end

nonmissing = findall(!ismissing, tradeoffs_grid.beta_pca);
tradeoffs_grid = tradeoffs_grid[nonmissing, :];
grid_years = grid_years[nonmissing, :];

##---------------------------------------------------------------
## 03. Join climate data
##---------------------------------------------------------------

## read climate data
clim_data = CSV.read(datadir("terraclim", "terraclim_processed_grid.csv"), DataFrame);
climatologies = CSV.read(datadir("terraclim", "climatologies_grid.csv"), DataFrame);

grid_years[!, :drought_strength] = Vector{Union{Missing, Float64}}(undef, nrow(grid_years));
grid_years[!, :nondrought_pdsi] = Vector{Union{Missing, Float64}}(undef, nrow(grid_years));
grid_years[!, :prop_drought] = Vector{Union{Missing, Float64}}(undef, nrow(grid_years));
grid_years[!, :min_pdsi] = Vector{Union{Missing, Float64}}(undef, nrow(grid_years));

grid_ids = unique(grid_years.grid_cell_id)

clim_dict = Dict(g => clim_data[clim_data.grid_cell_id .== g, :] for g in grid_ids)
years_dict = Dict(g => unique(grid_years[grid_years.grid_cell_id .== g, :MEASYEAR]) for g in grid_ids)

prev_measyear_dict =
    Dict((g, y) => grid_years[(grid_years.grid_cell_id .== g) .&&
                              (grid_years.MEASYEAR .== y), :PREV_MEASYEAR][1] for g in grid_ids for y in years_dict[g])

for g in ProgressBar(grid_ids)

    clim_sub = clim_dict[g]
    for y in years_dict[g]
        prev_measyear = prev_measyear_dict[(g, y)]
        idx = (grid_years.grid_cell_id .== g) .&& (grid_years.MEASYEAR .== y)

        pdsi = clim_sub[clim_sub.year .<= y .&& clim_sub.year .> prev_measyear, :pdsi]
        droughts = pdsi[pdsi .<= -2.0]
        drought_strength = isempty(droughts) ? NaN : mean(droughts)
        nondrought_pdsi = isempty(pdsi[pdsi .> -2.0]) ? NaN : mean(pdsi[pdsi .> -2.0])
        prop_drought = length(droughts) / max(y - prev_measyear, 1)
        min_pdsi = isempty(pdsi) ? NaN : minimum(pdsi)

        grid_years[idx, :drought_strength] .= drought_strength
        grid_years[idx, :nondrought_pdsi] .= nondrought_pdsi
        grid_years[idx, :prop_drought] .= prop_drought
        grid_years[idx, :min_pdsi] .= min_pdsi
    end
end

## join climatology data
leftjoin!(grid_years, climatologies[:,[:grid_cell_id, :MAP, :MAT]], on = :grid_cell_id)

## join climate data to tradeoffs_grid
leftjoin!(tradeoffs_grid, grid_years[:,[:grid_cell_id, :MEASYEAR, :drought_strength,
                                        :nondrought_pdsi, :prop_drought, :min_pdsi, :MAP, :MAT]], on = [:grid_cell_id, :MEASYEAR])

## clean up column names
rename!(grid_years, lowercase.(names(grid_years)))
rename!(tradeoffs_grid, lowercase.(names(tradeoffs_grid)))

rename!(grid_years, :tpa_init_subset => :tpa, :ba_init => :ba, :stdage => :stand_age)
rename!(tradeoffs_grid, :tpa_init_subset => :tpa, :ba_init => :ba, :stdage => :stand_age)

CSV.write(datadir("grid_years_clim.csv"), grid_years)
CSV.write(datadir("full_data.csv"), tradeoffs_grid)
