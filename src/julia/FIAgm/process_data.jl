
function process_growth_and_mortality(gmd::DataFrame, td::DataFrame, cond_data::DataFrame,
                                      maxiter::Int64 = 500)

    ## create output array, faster if we later rejoin this to dataframe
    ret = Array{Float64}(undef, nrow(gmd), 10)
    spp = Vector{String}(undef, nrow(gmd))

    ## this takes forever
    for r in ProgressBar(1:minimum([nrow(ret), maxiter]))
   
        # println("Iteration: $r")

        ## pull tree data corresponding to plot
        tl = pull_trees(gmd[r,:], td)

        ## check if there are at least five live trees in measurement
        if nrow(tl.remeas) == 0 || sum(tl.init.live) < 5

            ## write data to output array
            ret[r,:] = repeat([-9999.9], 10)

        else

            cond = gmd[r,:CONDID]

            ## record initial TPA, both full treelist and subset of those contained in remeas.
            ## if cond_ids don't align, can't calculate tpa_full, bc no-way to separate full tree list
            init_tpa = calc_tpa(tl, cond)

            ## calculate basal area
            basal_area = calc_ba(tl, cond_data, cond)

            ## calculate Δ
            Δt = tl.remeas.MEASYEAR[1] - tl.init.MEASYEAR[1]
            growth_ba = (basal_area.remeas - basal_area.init) / Δt
            growth_pct = (basal_area.remeas / basal_area.init)^(1/Δt) - 1

            ## calculate mortality
            mort = calc_mortality(tl, cond_data, cond, basal_area, init_tpa)

            #@infiltrate r == 3490
            ## write data to output array
            if mort.pct_ba == 0
                mort_pct_ba_ann = 1 - (1 - 1e-20)^(1 / Δt)
            else
                mort_pct_ba_ann = 1 - (1 - mort.pct_ba)^(1 / Δt)
            end

            if mort.pct_tpa == 0
                mort_pct_tpa_ann = 1 - (1 - 1e-20)^(1 / Δt)
            else
                mort_pct_tpa_ann = 1 - (1 - mort.pct_tpa)^(1 / Δt)
            end

            ret[r,:] = [init_tpa.full, init_tpa.subset, basal_area.init, basal_area.remeas,
                        growth_ba, growth_pct, mort.ba / Δt, mort.tpa / Δt,
                        mort_pct_ba_ann, mort_pct_tpa_ann]

        end

        cond_align = true

    end

    return hcat(gmd, DataFrame(ret, [:tpa_init_full, :tpa_init_subset, :ba_init,
                                     :ba_remeas, :growth_ba, :growth_pct,
                                     :mort_ba, :mort_tpa, :mort_ba_pct, :mort_tpa_pct]))

end


function process_state(state::String, write_file::Bool = false, dir::String = "./data/gm/")

    cond_data = FIAutils.read_fia(state, "COND")
    plot_data = FIAutils.read_fia(state, "PLOT")
    tree_data = FIAutils.read_fia(state, "TREE")
    grm_data = FIAutils.read_fia(state, "TREE_GRM_COMPONENT")

    ## takes a minute
    fill_tpa_unadj!(tree_data, grm_data)

    ## make tree id variable more explicit then join GRM data to tree data
    rename!(tree_data, :CN => :TRE_CN)
    leftjoin!(tree_data, select(grm_data, [:TRE_CN, :SUBP_TPAMORT_UNADJ_AL_FOREST]), on = :TRE_CN)

    ## remove unnecessary columns and rename
    select!(tree_data, [:TRE_CN, :PREV_TRE_CN, :PLT_CN, :DIA, :HT,
                        :ACTUALHT, :SPCD, :STATUSCD, :CONDID, :PREVCOND,
                        :TPA_UNADJ, :SUBP_TPAMORT_UNADJ_AL_FOREST,
                        :AGENTCD, :DIACHECK, :CARBON_AG, :CARBON_BG]);
    rename!(tree_data, :SUBP_TPAMORT_UNADJ_AL_FOREST => :TPAMORT_UNADJ)

    ## calculate tree-level basal area per acre
    tree_data.UNADJ_BA = π .* (tree_data.DIA ./ (2 * 12)) .^ 2 .* tree_data.TPA_UNADJ

    assign_plot_id!(plot_data)
    assign_sample_number!(plot_data)

    rename!(plot_data, :CN => :PLT_CN)
    leftjoin!(tree_data, plot_data[:,[:PLT_CN, :PREV_PLT_CN, :PLOT_ID, :INVYR, :SAMPLE_NUMBER,
                   :STATECD, :UNITCD, :COUNTYCD, :LAT, :LON, :ELEV, :MEASYEAR]], on = [:PLT_CN])

    growth_mortality_data = combine(groupby(
        tree_data[tree_data.SAMPLE_NUMBER .> 1,
                  [:PLT_CN, :CONDID, :PREVCOND, :PREV_PLT_CN, :PLOT_ID, :INVYR, :MEASYEAR, :SAMPLE_NUMBER,
                   :STATECD, :UNITCD, :COUNTYCD, :LAT, :LON, :ELEV]], [:PLOT_ID, :SAMPLE_NUMBER, :CONDID])) do sdf
                       first(sdf)
                    end

    leftjoin!(growth_mortality_data, cond_data[:,[:PLT_CN, :CONDID, :COND_STATUS_CD, :OWNGRPCD, :FORTYPCD,
                                                 :STDAGE,:CONDPROP_UNADJ, :DSTRBCD1, :DSTRBCD2, :DSTRBCD3]],
             on = [:PLT_CN, :CONDID])

    sort!(growth_mortality_data, [:PLOT_ID, :SAMPLE_NUMBER, :CONDID])

    ## remove non-forest and unsampled conditions from data as well as records without a PREV_PLT_CN
    growth_mortality_data = growth_mortality_data[growth_mortality_data.COND_STATUS_CD .== 1 .&&
        .!ismissing.(growth_mortality_data.PREV_PLT_CN),:]

    ## remove all conditions that represent less than 30% of plot area
    growth_mortality_data = growth_mortality_data[growth_mortality_data.CONDPROP_UNADJ .> 0.3, :]

    println("Starting...")

    out = process_growth_and_mortality(growth_mortality_data, tree_data, cond_data, Int(1e5))

    ## remove all plots with straight 0s:
    out = out[out.tpa_init_full .!= -9999.9, :]

    if write_file
        CSV.write(dir * "gm_" * state * ".csv", out)
    end

    return out

end
