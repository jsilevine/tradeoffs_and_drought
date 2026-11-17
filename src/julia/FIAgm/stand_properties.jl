
struct TPA

    full::Float64
    subset::Float64
    mort::Float64

    function TPA(full, subset, mort)
        new(full, subset, mort)
    end
end

function calc_tpa(tl::TreeList, cond::Int64)

    full = Float64(-99.0)

    if tl.cond_align
        full = sum(skipmissing(tl.init[tl.init.STATUSCD .== 1 .&& tl.init.CONDID .== cond, :TPA_UNADJ]))
    end

    subset::Float64 = sum(skipmissing(tl.init[tl.init.live .&& tl.init.CONDID .== cond, :TPA_UNADJ]))
    mort::Float64 = sum(skipmissing(tl.init[tl.init.live_all .&& tl.init.CONDID .== cond, :TPA_UNADJ]))

    return TPA(full, subset, mort)

end

struct BA
    init::Float64
    init_mort::Float64
    remeas::Float64

    function BA(init, init_mort, remeas)
        new(init, init_mort, remeas)
    end

end

"""
function calc_ba(tl::TreeList, cond_data::DataFrame, cond::Int64)

Calculate basal area for initial measurement and remeasurement. Returns object of
immutable struct 'ba', which has three components:
          1. 'ba_init'
          2. 'ba_init_mort'
          3. 'ba_remeas'
"""
function calc_ba(tl::TreeList, cond_data::DataFrame, cond::Int64)

    ## TODO Figure out best way to deal with mismatches in conditions between plot meas.
    ## Currently, use condprop_unadj from remeasurement after forcing match for trees shared
    ## between measurements. I don't know if this is accurate, but it at least makes comparing
    ## BA numbers reasonable.

    ba_init = Float64
    ba_init_mort = Float64
    condprop_unadj_remeas::Float64 = cond_data[cond_data.PLT_CN .== unique(tl.remeas.PLT_CN) .&&
            cond_data.CONDID .== cond, :CONDPROP_UNADJ][1]

    ba_init = sum(tl.init[tl.init.live .&& tl.init.CONDID .== cond .&&
        .!ismissing.(tl.init.UNADJ_BA), :UNADJ_BA] ./
        condprop_unadj_remeas)
    ba_init_mort = sum(tl.init[tl.init.live_all .&& tl.init.CONDID .== cond .&&
        .!ismissing.(tl.init.UNADJ_BA), :UNADJ_BA] ./
        condprop_unadj_remeas)
    ba_remeas = sum(tl.remeas[tl.remeas.live .&& tl.remeas.CONDID .== cond .&&
        tl.remeas.PREV_TRE_CN .∈ [tl.init[tl.init.live_all .&& tl.init.CONDID .== cond .&&
        .!ismissing.(tl.init.UNADJ_BA), :TRE_CN]], :UNADJ_BA] ./
        condprop_unadj_remeas)

    return BA(ba_init, ba_init_mort, ba_remeas)

end
