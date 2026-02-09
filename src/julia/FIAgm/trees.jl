
struct TreeList

    init::DataFrame
    remeas::DataFrame
    cond_align::Bool

    function TreeList(init, remeas, cond_align)
        new(init, remeas, cond_align)
    end

end

function pull_trees(ps::DataFrameRow, td::DataFrame)

    ## pull trees from remeasurement sample
    remeas = td[td.PLT_CN .== ps.PLT_CN .&&
        .!ismissing.(td.PREV_TRE_CN) .&&
        .!ismissing.(td.TPA_UNADJ) .&&
        td.STATUSCD .!= 0,:]

    ## pull trees from initial (prev.) sample
    init = td[td.PLT_CN .== ps.PREV_PLT_CN,:]

    ## delineate dead trees from initial sample
    init.dead .= init.STATUSCD .== 2

    ## delineate trees which were alive in initial sample, but are dead (killed by non-fire natural causes)
    ## in remeasurement (mortality)

    if nrow(remeas) == 0
        init.mort = repeat([0], nrow(init))
    else
        remeas.mort .= (remeas.STATUSCD .== 2 .&&
            .!ismissing.(remeas.UNADJ_BA) .&&
            .!(remeas.PREV_TRE_CN .∈ [init[init.dead, :TRE_CN]]) .&&
            (ismissing.(remeas.AGENTCD) .||
            div.(remeas.AGENTCD, 10) .∈ [[1, 2, 7]]))
            # div.(remeas.AGENTCD, 10) .∈ [[1, 2, 4, 5, 6, 7]]))
        if length(remeas.mort) == 0 | sum(remeas.mort) == 0
            init.mort = repeat([0], nrow(init))
        else
            init.mort .= init.TRE_CN .∈ [remeas[remeas.mort, :PREV_TRE_CN]]
        end

    end

    ## delineate trees which are alive in remeasurement
    remeas.live .= remeas.STATUSCD .== 1 .&&
        .!ismissing.(remeas.UNADJ_BA)

    ## delineate trees which are alive in initial measurement, and are still alive and found in remeasurement
    init.live .= init.TRE_CN .∈ [remeas[remeas.live, :PREV_TRE_CN]]

    ## delineate trees which are alive in inital measurement, whether or not they are alive/found in remeasurement
    init.live_all .= init.STATUSCD .== 1 .&& .!ismissing.(init.UNADJ_BA) .&&
        init.TRE_CN .∈ [remeas.PREV_TRE_CN]

    ## align condition IDs, sometimes these change from measurement to measurement.
    ## If they do, force trees from initial measurement to reflect condition delineation in remeasurement.
    ## This allows us to compare from one to the other.
    cond_align = match_condids!(init, remeas)

    return TreeList(init, remeas, cond_align)

end
