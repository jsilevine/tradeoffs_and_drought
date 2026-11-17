##-------------------------------------------------------------------
## Function to fit GAMs in julia by calling R::mgcv::bam
##-------------------------------------------------------------------


__revise_mode__ = :eval

module ModelFitting

## prerequisites
using DataFrames, RCall
 
R"library(mgcv)"
R"library(gamlss)"

## for fitting mgcv::bam (GAM) models
"""
bam(response::String, covariates::Vector{String}, data::DataFrame; kwargs...)

Fit a GAM using mgcv::bam() in R, with the forumla constructed from Julia input.

Example:
```julia
m = bam(
    "growth_pct_scaled",
    ["s(year)", "s(lon, lat, k = 100)", "ba_scaled"], 
    df;
    method = "fREML", discrete = true
)

"""
function bam(response::String, covariates::Vector{String}, data::DataFrame; kwargs...)

    formula_str = response * " ~ " * join(covariates, " + ")
    r_kwargs = RObject(Dict(string(k) => v for (k,v) in kwargs))

    # Call mgcv::bam via do.call so kwargs expand properly
    R"""
    m <- do.call(
    bam,
    c(list(formula = as.formula($formula_str), data = $data), $r_kwargs)
    )
    """
    return RCall.reval("m")
end

## for fitting gamlss::gamlss (GAM) models, supports beta-binomial hurdle models
"""
bezi_gam(response::String, mu_fixed::Vector{String}, nu_fixed::Vector{String}, data::DataFrame; n.cyc::Int = 30)

Fit a Beta-Binomial GAM using gamlss::gamlss() in R, with the formula constructed from Julia input. The user should supply two
separate Vector{String} arguments, one including all the desired covariates and interactions for the Beta component of the
model (mu_fixed), and the second including all the desired covariates and interactions for the Binomial component.
Additional keyword arguments will be passed directly to the gamlss:gamlss() call.

Example:
```julia
mu_fixed = ["pb(year)", "drought_strength", "drought_strength:beta_pca"]

nu_fixed = ["pb(year)"]

m = bezi_gam("mortality", mu_fixed, nu_fixed, data, family = "BEZI")
```
"""
function bezi_gam(response::String, mu_fixed::Vector{String}, nu_fixed::Vector{String}, data::DataFrame; n_cyc::Int = 30)

    mu_formula_str = response * " ~ " * join(mu_fixed, " + ")
    nu_formula_str = "~ " * join(nu_fixed, " + ")

    @rput data
    R"""
    m <- gamlss::gamlss(formula = as.formula($mu_formula_str),
                nu.formula = as.formula($nu_formula_str),
                data = data,
                family = BEZI,
                control = gamlss::gamlss.control(n.cyc = $n_cyc, trace = TRUE))
    """

    return(RCall.reval("m"))

end


"""
    summarize_gam(model::RObject)

Summarize a GAM model fitted with `bam()` or `beta_binom()`, printing the summary to the console.

"""
function summarize_gam(model)
    R"""
    summary($model)
    """
end


end
