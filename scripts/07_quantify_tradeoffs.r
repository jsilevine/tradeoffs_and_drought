##--------------------------------------------------------------
## 07_quantify_tradeoffs.r: Quantify plot-level tradeoffs
##
## author: jacob levine; jacob.levine@utah.edu
##--------------------------------------------------------------

library(data.table)
library(FactoMineR)
library(here)

source(here("src", "r", "trait_processing.r"))

gm_data <- fread(here("data", "gm", "gm_full.csv"))
traits <- fread(here("data", "traits_imputed_merged.csv"))
colnames(traits)[2] <- "SPCD"
states <- read.csv(here("data", "state_list.csv"))

## remove plots with Inf values, due to no separation between measurements
gm_data <- gm_data[is.finite(gm_data$growth_pct) & is.finite(gm_data$mort_ba_pct),]

##---------------------------------------------------------------------------
## 01. Fit PCA to traits, use to classify species by growth and drought strategy
##---------------------------------------------------------------------------

drought_pca_vars <- c("p50", "rdmax")

traits$p50 <- traits$p50 * -1
traits$P50 <- traits$P50 * -1

pca_data <- traits[complete.cases(traits[, ..drought_pca_vars]),
                   ..drought_pca_vars]

drought_pca <- PCA(pca_data, scale.unit = TRUE)

growth_pca_vars <- c("Amax", "gsmax", "LeafN", "SLA")

pca_data <- traits[complete.cases(traits[, ..growth_pca_vars]),
                   ..growth_pca_vars]

growth_pca <- PCA(pca_data, scale.unit = TRUE)

saveRDS(drought_pca, here("data", "model_data", "drought_pca.rds"))
saveRDS(growth_pca, here("data", "model_data", "growth_pca.rds"))

## project pca in advance, these functions modify 'traits' in place,
## appending the projected PC axes. 
project_pca(traits, drought_pca, drought_pca_vars, c("drought_pc1", "drought_pc2"))
project_pca(traits, growth_pca, growth_pca_vars, c("growth_pc1", "growth_pc2"))

##---------------------------------------------------------------
## 02. Fit linear models to each plot to estimate tradeoff slope
##     and uncertainty
##---------------------------------------------------------------

trait_vars <- c("species", "p50", "Amax", "gsmax", "rdmax", "SLA", "LeafN", "growth_pc1", "drought_pc2")

output_cols <- c("richness",
                 "intercept_base", "beta_base", "se_base", "td_base", "td_se_base", "r2_base", "os_r2_base",
                 "intercept_pca", "beta_pca", "se_pca", "td_pca", "td_se_pca", "r2_pca", "os_r2_pca",
                 "intercept_pca_p50", "beta_pca_p50", "se_pca_p50", "td_pca_p50", "td_se_pca_p50", "r2_pca_p50", "os_r2_pca_p50",
                 "cwm_Amax", "cwm_p50", "cwm_growth_pc1", "cwm_drought_pc2",
                 "min_Amax", "min_p50", "min_growth_pc1", "min_drought_pc2",
                 "max_Amax", "max_p50", "max_growth_pc1", "max_drought_pc2")

## runs full estimation pipeline, see src/r/trait_processing.r for details. 
estimate_tradeoffs(states, gm_data, traits, trait_vars, output_cols, TRUE, "tradeoffs")

##-----------------------------------------------------------------
## 03. Fit additional tradeoffs for current tree list, to estimate
##     trends over time
##-----------------------------------------------------------------

estimate_tradeoffs(states, gm_data, traits, trait_vars, output_cols, FALSE, "tradeoffs_current")

