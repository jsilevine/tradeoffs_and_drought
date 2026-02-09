
struct Mortality

    ba::Float64
    tpa::Float64
    pct_ba::Float64
    pct_tpa::Float64

    function Mortality(ba, tpa, pct_ba, pct_tpa)
        new(ba, tpa, pct_ba, pct_tpa)
    end

end

function calc_mortality(tl::TreeList, cond_data::DataFrame, cond::Int64, basal_area::ba, tpa_init::tpa)

    ## TODO Same as basal area calculations, figure out how to handle discrepancies in
    ## condition ID between measurements.

    condprop_unadj::Float64 = cond_data[cond_data.PLT_CN .== unique(tl.remeas.PLT_CN) .&&
            cond_data.CONDID .== cond, :CONDPROP_UNADJ][1]

    if sum(tl.init.mort) == 0 || all(ismissing.(tl.init[tl.init.mort .&& tl.init.CONDID .== cond, :UNADJ_BA]))
        return Mortality(0.0, 0.0, 0.0, 0.0)
    else

        ba::Float64 =
            sum(
            tl.init[tl.init.mort .&& tl.init.CONDID .== cond, :UNADJ_BA] ./
                condprop_unadj
        )
        tpa::Float64 =
            sum(
            tl.init[tl.init.mort .&& tl.init.CONDID .== cond, :TPA_UNADJ] ./
                condprop_unadj)

        return Mortality(ba, tpa,
                         minimum([1.0, ba / basal_area.init_mort]),
                         minimum([1.0, tpa / (tpa_init.mort / condprop_unadj)]))

    end
end
