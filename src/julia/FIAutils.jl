##-------------------------------------------------------------------
## Basic utility for reading and working with raw FIA data
##-------------------------------------------------------------------

module FIAutils

include(joinpath(@__DIR__, "utils.jl"))

using .Utils

using DataFrames, CSV, Distributed

export read_fia, fill_tpa_unadj!
export assign_plot_id!, assign_sample_number!

function read_fia(state::String, tabtyp::String = "TREE")
    return CSV.read(datadir("fia", state, state * "_" * tabtyp * ".csv"), DataFrame)
end

function fill_tpa_unadj!(td::DataFrame, grmd::DataFrame)

    ## taking this from Grayson's code from Ecology Letters paper on climate risks, 'approved by JS'
    ## Fills in gaps in the unadjusted trees per acre using the growth TPA factor from the GRM data
    ## This code is quite slow so we parallelize it using @distributed
    Threads.@threads for i in findall(ismissing.(td.TPA_UNADJ))
        if (td[i,:CN] ∈ grmd[:,:TRE_CN])
            td[i,:TPA_UNADJ] = grmd[grmd.TRE_CN .== td[i,:CN], :SUBP_TPAGROW_UNADJ_AL_FOREST][1]
        end
    end

end

function assign_plot_id!(pd::DataFrame)

    ## create variable with which to uniquely identify a plot
    pd.PLOT_ID = Vector{String}(undef, nrow(pd))
    for i in 1:nrow(pd)
        pd[i,:PLOT_ID] = string(pd[i,:PLOT]) * "_" *
            string(pd[i,:STATECD]) * "_" * string(pd[i,:UNITCD]) * "_" *
            string(pd[i,:COUNTYCD])
    end

end

function assign_sample_number!(pd::DataFrame)

    ## determine which sample/resample number each entry in plot_data is
    ## group by unique plot identifier
    group = groupby(pd, :PLOT_ID)

    ## assign prev_plt ids for plots which are missing them
    ## also denote which sample number each plot record is
    for g in group
        sort!(g, :INVYR)
        g.SAMPLE_NUMBER = [1:1:nrow(g);]
        ## some plots are missing PREV_PLT_CN even when they are a resample
        if nrow(g) > 1 && sum(ismissing.(g.PREV_PLT_CN)) > 1
            for i in 2:nrow(g)
                g[i,:PREV_PLT_CN] = g[i-1,:CN]
            end
        end
    end

    ## paste data back together
    pd = DataFrame(group)

end

end
