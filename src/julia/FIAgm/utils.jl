##-------------------------------------------------------------------
## Utility functions for data processing and analysis in FIAgm
##-------------------------------------------------------------------

function pull_plot(df::DataFrame, plot::String)
    df[df.PLOT_ID .== plot .&& .!ismissing.(df.PREV_PLT_CN) .&& df.SAMPLE_NUMBER .> 1,:]
end

function match_condids!(init::DataFrame, remeas::DataFrame)
    # Extract relevant subsets, sorted by tree identifiers
    init_condids   = sort(init[init.live_all, [:TRE_CN, :CONDID]], :TRE_CN)
    remeas_condids = sort(remeas[remeas.PREV_TRE_CN .∈ init_condids.TRE_CN, [:PREV_TRE_CN, :CONDID]], :PREV_TRE_CN)

    if init_condids.CONDID == remeas_condids.CONDID
        return true   # aligned already
    else
        # Overwrite init.CONDID with remeas.CONDID for matching trees
        mapping = Dict(remeas_condids.PREV_TRE_CN .=> remeas_condids.CONDID)
        init[init.live_all, :CONDID] .= get.(Ref(mapping), init[init.live_all, :TRE_CN])
        return false  # had to realign
    end
end
