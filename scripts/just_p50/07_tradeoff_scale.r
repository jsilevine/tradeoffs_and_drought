##--------------------------------------------------------------
## 07_tradeoff_scale.r -- compare local tradeoffs to regional tradeoffs
##
## author: jacob levine; jacob.levine@utah.edu
##--------------------------------------------------------------

## core libraries for data, spatial ops, PCA, parallel apply, plotting, and models
library(data.table)
library(FactoMineR)
library(here)
library(sf)
library(future.apply)
library(ggplot2)
library(sf)
library(metafor)
library(brms)

## project helper for trait processing (adds project-specific functions)
source(here("src", "r", "trait_processing.r"))
source(here("src", "r", "tradeoff_scale.r"))

## read data
gm_data <- fread(here("data", "gm", "gm_full.csv"))
traits  <- fread(here("data", "traits_imputed_merged.csv"))
colnames(traits)[2] <- "SPCD"
states  <- read.csv(here("data", "state_list.csv"))

## drop plots with infinite values from percent-change calculations
gm_data <- gm_data[is.finite(gm_data$growth_pct) & is.finite(gm_data$mort_ba_pct),]

##---------------------------------------------------------------------------
## 01. Load PCA models and attach PC scores to species trait table
##---------------------------------------------------------------------------

## variables used for each PCA
drought_pca_vars <- c("p50", "rdmax")
growth_pca_vars  <- c("Amax", "gsmax", "LeafN", "SLA")

## make p50 positive so larger = more drought tolerant (consistent sign convention)
traits$p50 <- traits$p50 * -1
traits$P50 <- traits$P50 * -1

## load precomputed PCA objects (estimated once and saved to disk)
drought_pca <- readRDS(here("data", "model_data", "drought_pca.rds"))
growth_pca  <- readRDS(here("data", "model_data", "growth_pca.rds"))

## project PCA scores into the traits table (functions modify `traits` in place)
project_pca(traits, drought_pca, drought_pca_vars, c("drought_pc1", "drought_pc2"))
project_pca(traits, growth_pca,  growth_pca_vars,  c("growth_pc1", "growth_pc2"))

##---------------------------------------------------------------
## 02. Fit per-plot linear models to estimate tradeoff slope + uncertainty
##     (overall vs community comparisons at multiple spatial scales)
##---------------------------------------------------------------

## tradeoff estimates (may be used by check_scale_diff)
tradeoffs <- fread(here("data", "tradeoffs", "tradeoffs_full.csv"))

## trait vars used when joining/aggregating within check_scale_diff
trait_vars <- c("species", "p50", "Amax", "gsmax", "rdmax", "SLA", "LeafN",
                "growth_pc1", "drought_pc2")

## restrict to continental US plots (drop AK/HI by STATECD)
gm_data <- gm_data[!(STATECD %in% c(2, 15))]

## build simple sf object of plot centroids for spatial queries
plots_spatial <- st_as_sf(gm_data[, .(PLOT_ID, LAT, LON, STATECD, PLT_CN, PREV_PLT_CN, tpa_init_full)],
                          coords = c("LON", "LAT"), crs = 4326)

## radii to test (meters), from large to small (so early iterations sample broad neighborhoods)
radii <- seq(10000, 500000, length.out = 15)
radii <- radii[order(radii, decreasing = TRUE)]

## prepare output frame to store results for each radius
output <- data.frame(radius = radii,
                     alpha_diff = numeric(length(radii)),
                     se_diff = numeric(length(radii)),
                     ci_lower = numeric(length(radii)),
                     ci_upper = numeric(length(radii)),
                     prop_pos = numeric(length(radii)),
                     alpha_diff_lm = numeric(length(radii)),
                     se_diff_lm = numeric(length(radii)),
                     stringsAsFactors = FALSE)

## loop over radii and compute scale-dependent difference using helper check_scale_diff
for (i in seq_along(radii)) {
  res <- check_scale_diff(
    plots_spatial,
    radii[i],
    traits,
    tradeoffs,
    trait_vars,
    states,
    min_n = 1000,    ## minimum number of plots to consider reliable
    max_n = 1000,    ## cap sample size so comparisons are consistent
    seed = 1211      ## reproducible sampling inside helper
  )

  ## store results
  output[i, "alpha_diff"]      <- res$alpha
  output[i, "se_diff"]         <- res$se
  output[i, "ci_lower"]        <- res$ci[1]
  output[i, "ci_upper"]        <- res$ci[2]
  output[i, "prop_pos"]        <- res$prop_pos
  output[i, "alpha_diff_lm"]   <- res$alpha_diff_lm
  output[i, "se_diff_lm"]      <- res$se_diff_lm

  ## persist intermediate results to disk so long runs can be resumed / inspected
  write.csv(output, here("data", "tradeoffs", "scale_analysis", "tradeoff_scale_analysis_new.csv"),
            row.names = FALSE)

  message("Completed radius ", radii[i], " m")
  print(output[i, ])
}

## quick diagnostic plot of alpha_diff vs radius
plot(output$radius, output$alpha_diff, type = "b",
     xlab = "Radius (m)", ylab = "Difference in Tradeoff Slope (Overall - Community)")

##---------------------------------------------------------------------------
## 03. Meta-analytic and Bayesian trend of scale effect across radii
##---------------------------------------------------------------------------

## reload saved output (safe if script restarted)
output <- read.csv(here("data", "tradeoffs", "scale_analysis", "tradeoff_scale_analysis_new.csv"))

## scale radius for modeling (center = FALSE keeps origin at zero)
output$radius_scaled <- scale(output$radius, center = FALSE)

## mixed-effects meta-regression (metafor) using reported se_diff
trend_model <- rma(yi = output$alpha_diff, sei = output$se_diff, mods = output$radius_scaled, method = "REML")
summary(trend_model)

## build prediction grid (radius in km) and transform to model's scaled covariate
radius_km_seq   <- seq(0, 510, length.out = 250)
radius_m_seq    <- radius_km_seq * 1000
radius_scaled_seq <- (radius_m_seq) / attr(output$radius_scaled, "scaled:scale")

## design matrix for meta-regression prediction
X_new <- model.matrix(~ radius_scaled_seq)

## predict from metafor model (returns pred, ci.lb, ci.ub)
pred <- predict(trend_model, newmods = X_new[, -1, drop = FALSE])

pred_df <- data.frame(
  radius_km = radius_km_seq,
  fit       = pred$pred,
  lower     = pred$ci.lb,
  upper     = pred$ci.ub
)

## plot fitted trend with uncertainty band and observed points
ggplot() +
  geom_ribbon(data = pred_df, aes(x = radius_km, ymin = lower, ymax = upper),
              fill = "#b3cde3", alpha = 0.6) +
  geom_line(data = pred_df, aes(x = radius_km, y = fit), color = "#333333", size = 1.2) +
  geom_hline(yintercept = 0, color = "darkgray", size = 1) +
  geom_point(data = output, aes(x = radius / 1000, y = alpha_diff), size = 3.5, color = "#2c7fb8") +
  geom_errorbar(data = output, aes(x = radius / 1000, ymin = ci_lower, ymax = ci_upper),
                width = 10, size = 0.8, color = "#2c7fb8") +
  xlab("Radius (km)") +
  ylab("Difference in Tradeoff Slope (Overall - Community)") +
  scale_x_continuous(limits = c(0, 510), expand = c(0, 0)) +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        text = element_text(size = 18))

ggsave(here("figures", "tradeoff_scale", "tradeoff_scale.pdf"),
       width = 9, height = 6)

##---------------------------------------------------------------------------
## 04. Map example radii around three cities (visualize scale choices)
##---------------------------------------------------------------------------

## centers (lon, lat) and radii (meters) to display on map
centers <- data.frame(
  city = c("Los Angeles", "Kansas City", "Washington, DC"),
  lon  = c(-118.2437, -94.5786, -77.0369),
  lat  = c(34.0522,   39.0997,  38.9072),
  r_m  = c(500000, 150000, 50000)   ## radii used in analysis (m)
)

## convert to sf and project to Albers (5070) for buffering in meters
centers_sf <- st_as_sf(centers, coords = c("lon","lat"), crs = 4326) |>
  st_transform(5070)

## build circular polygons (buffers) representing each radius
circles_sf <- st_buffer(centers_sf, centers_sf$r_m)

## colors for map legend
cols <- c(
  "Los Angeles"    = "#1f78b4",
  "Kansas City"    = "#e66101",
  "Washington, DC" = "#238b45"
)

## plot filled circles on base map (semi-transparent)
base_map_filled(crs = 5070) +
  geom_sf(data = circles_sf, aes(fill = city), alpha = 0.5, color = NA) +
  scale_fill_manual(values = cols) +
  scale_color_manual(values = cols) +
  theme(legend.position = "right", panel.grid.major = element_blank())

ggsave(here("figures", "tradeoff_scale", "tradeoff_scale_map.pdf"),
       width = 8, height = 6)
