##-------------------------------------------------------------------
## 02_overall_tradeoff.jl: Analyze overall tradeoff between drought
##                         tolerance and growth acquisitiveness
##
## author: jacob levine; jacob.levine@utah.edu
##-------------------------------------------------------------------

##-------------------------------------------------------------------
## 00. load data and libraries
##-------------------------------------------------------------------

include(joinpath(@__DIR__, "..", "..", "src", "julia", "utils.jl"))

using CSV, DataFrames
using MultivariateStats
using Statistics
using RCall
using GLM
using StatsPlots, StatsBase
using Random

using .Utils

## read imputed trait data
traits = CSV.read(datadir("traits_imputed_merged.csv"), DataFrame, missingstring = "NA");

## make p50 positive so larger values indicate greater drought tolerance
traits.p50 .= -1 .* traits.p50;

function project_pca(;data::DataFrame, pca_center::Vector, pca_scale::Vector,
    pca_rotation::Matrix, pca_vars::Vector{Symbol}, pc_names::Vector{Symbol})

    ## extract variables used in PCA
    data_subset = Matrix(data[:, pca_vars]);

    ## apply PCA centering and scaling
    data_scaled = (data_subset .- pca_center') ./ pca_scale';

    ## project into PCA space
    pc_scores = data_scaled * pca_rotation;

    ## store first two PC scores in data frame
    data[:, pc_names[1]] = pc_scores[:, 1];
    data[:, pc_names[2]] = pc_scores[:, 2];

    return data
end

## load PCA objects from R
R"""
library(here)
drought_pca <- readRDS(here("data", "model_data", "drought_pca.rds"))  
growth_pca  <- readRDS(here("data", "model_data", "growth_pca.rds"))
"""

## extract PCA parameters from R objects
drought_center   = rcopy(R"drought_pca$call$centre")
drought_scale    = rcopy(R"drought_pca$call$ecart.type")
drought_rotation = rcopy(R"drought_pca$var$coord")

growth_center   = rcopy(R"growth_pca$call$centre")
growth_scale    = rcopy(R"growth_pca$call$ecart.type")
growth_rotation = rcopy(R"growth_pca$var$coord")

## variables used in each PCA
drought_pca_vars = [:p50, :rdmax];
growth_pca_vars  = [:Amax, :gsmax, :LeafN, :SLA];

## project drought and growth PCA scores into traits table
traits = project_pca(
    data = traits,
    pca_center = drought_center,
    pca_scale = drought_scale,
    pca_rotation = drought_rotation,
    pca_vars = drought_pca_vars,
    pc_names = [:drought_pc1, :drought_pc2]
);

traits = project_pca(
    data = traits,
    pca_center = growth_center,
    pca_scale = growth_scale,
    pca_rotation = growth_rotation,
    pca_vars = growth_pca_vars,
    pc_names = [:growth_pc1, :growth_pc2]
);

## visualize overall tradeoff in PCA space
R"""
library(ggplot2)
ggplot($traits, aes(x = drought_pc2, y = growth_pc1)) +
    geom_point(size = 5, alpha = 0.5) +
    geom_smooth(method = "lm", fill = "lightgray", color = "#e6550d", size = 2) +
    scale_x_continuous(expand = c(0,0)) +
    xlab("Drought tolerance") +
    ylab("Resource Acquisitiveness") +
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          axis.text = element_text(size = 14),
          axis.title = element_text(size = 16))
"""

## save figure
R"ggsave(here('figures', 'overall_tradeoff.pdf'), width = 5, height = 4, dpi = 300)"

## plot raw trait relationship
R"""
ggplot($traits, aes(x = p50, y = Amax)) +
    geom_point(size = 3) +
    geom_smooth(method = "lm", fill = "lightgray") +
    scale_x_continuous(expand = c(0,0)) +
    theme_bw()
"""

## fit linear model in PCA space
model = lm(@formula(growth_pc1 ~ drought_pc2), traits);
println(coeftable(model))
coef(model)

##-------------------------------------------------------------------
## 02. How does expected tradeoff strength vary by diversity?
##-------------------------------------------------------------------

## simulate subsampling by richness
output = Array{Float64}(undef, 500, 40)
for i in eachindex(output[:,1])
    for j in eachindex(output[1,:])
        nd = traits[sample(eachindex(traits[:,1]), j+3, replace = false), :]
        model = lm(@formula(growth_pc1 ~ drought_pc2), nd)
        output[i,j] = coef(model)[2]
    end
end

## summarize coefficient distribution by sample size
cm    = Vector{Float64}(undef, size(output,2))
lower = Vector{Float64}(undef, size(output,2))
upper = Vector{Float64}(undef, size(output,2))

for i in eachindex(output[1,:])
    cm[i]    = mean(skipmissing(output[:,i]))
    lower[i] = quantile(skipmissing(output[:,i]), 0.025)
    upper[i] = quantile(skipmissing(output[:,i]), 0.975)
end

## write null distribution to disk
CSV.write(datadir("null_tradeoffs.csv"), DataFrame(output, :auto))

## visualize null distribution envelope
plot(collect(3:size(output, 2)+2), cm)
plot!(collect(3:size(output, 2)+2), lower, linestyle = :dash)
plot!(collect(3:size(output, 2)+2), upper, linestyle = :dash)

## compare with empirical grid-level tradeoffs
data = CSV.read(datadir("trend_data.csv"), DataFrame, missingstring = "NA");
plot(data.richness, data.beta_pca, seriestype = :scatter)

## mean slope at fixed richness
mean(skipmissing(data.beta_pca[data.richness .== 15]))

## test richness dependence of tradeoff strength
model2 = lm(@formula(beta_pca ~ richness), data);
println(coeftable(model2))
