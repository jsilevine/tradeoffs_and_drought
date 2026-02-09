##---------------------------------------------------------------
## 05_tradeoff_drivers.r: Fit tradeoff models on grid data
##
## author: jacob levine; jacob.levine@utah.edu
##---------------------------------------------------------------

##--------------------------------------------------------------
## 00. Load libraries and data                                                      
##--------------------------------------------------------------

## Load required libraries
invisible(lapply(c("mgcv", "ggplot2", "dplyr", "sf", "gstat", "sp",
                   "terra", "patchwork", "here", "scales", "RColorBrewer",
                   "parallel"), library, character.only = TRUE))

## dependencies
source(here("src", "r", "maps.r"))
source(here("src", "r", "spatial_block_bootstrap.r"))
source(here("src", "r", "plot_tradeoff_drivers.r"))

full_data <- read.csv(here("data", "data_for_modeling.csv"))

full_data <- full_data[complete.cases(full_data[,c("beta_pca", "mat_scaled", "elev",
                                              "map_scaled", "richness_scaled",
                                              "stand_age_scaled", "lon", "lat",
                                              "se_pca")]),]

## drop water rows and NA ecoregion labelsg
full_data <- full_data[full_data$NA_L1NAME != "WATER" & !is.na(full_data$NA_L1NAME), ]

##--------------------------------------------------------------
## 01. Relationship between species richness and tradeoff slope                                                       
##--------------------------------------------------------------

## sample one row per grid cell to avoid within-cell pseudoreplication
set.seed(42)
ids <- unique(full_data$grid_cell_id)
out <- vector("list", length(ids))
for (i in seq_along(ids)) {
  rows_i <- which(full_data$grid_cell_id == ids[i])
  out[[i]] <- full_data[sample(rows_i, 1), , drop = FALSE]
}
sampled_data <- do.call(rbind, out)
rownames(sampled_data) <- NULL

## plot null distribution vs. empirically sampled data
ggplot() +
  geom_smooth(data = null_data_long[sample(1:nrow(null_data_long), 7500), ],
              aes(x = richness, y = beta_pca), method = "lm",
              color = "gray80", fill = "gray90", linewidth = 1.5, alpha = 0.5) +
  geom_smooth(data = sampled_data, aes(x = richness, y = beta_pca), method = "lm",
              color = "red", fill = "pink", linewidth = 1.5, alpha = 0.8) +
  theme_bw() +
  scale_x_continuous(expand = c(0,0), limits = c(3, 19)) +
  scale_y_continuous(expand = c(0,0))

ggsave(here("figures", "tradeoff_richness.png"))

##--------------------------------------------------------------
## 02. Set up data for analysis                                                       
##--------------------------------------------------------------

## use L1 or L2 name as ecoregion label; promote L2 for eastern temperate forests
full_data$ecoregion <- full_data$NA_L1NAME
full_data[full_data$ecoregion == "EASTERN TEMPERATE FORESTS", "ecoregion"] <- full_data[full_data$ecoregion == "EASTERN TEMPERATE FORESTS", "NA_L2NAME"]
full_data$ecoregion <- as.factor(full_data$ecoregion)

## intiialize empty reduced_data frame
reduced_data <- full_data[0, c("grid_cell_id", "lon", "lat", "beta_pca", "td_pca", "map_scaled",
                               "mat_scaled", "stand_age_scaled", "elev_scaled", "ecoregion",
                               "richness_scaled",
                               "se_pca")]

unique_ids <- unique(full_data$grid_cell_id)

## aggregate grid values into a single row per grid cell (aggregation across years using precision-weighted pooling)
for (i in seq_along(unique_ids)) {

  id <- unique_ids[i]
  subset <- full_data[full_data$grid_cell_id == id, ]

  reduced_data[i, "grid_cell_id"] <- id

  beta <- subset$beta_pca
  se   <- subset$se_pca
  td   <- subset$td_pca
  td_se <- subset$td_se_pca

  # keep only entries with valid beta and se (we're doing precision weighting)
  mask <- !is.na(beta) & !is.na(se) & se > 0
  beta <- beta[mask]
  se   <- se[mask]
  td   <- td[mask]
  td_se <- td_se[mask]

  k <- length(beta)

  if (k == 0) {
    ## no valid observations in cell
    beta_mean <- NA_real_
    se_mean   <- NA_real_
    td_mean   <- NA_real_
    td_se_mean <- NA_real_

  } else if (k == 1) {
    # single observation: use its values directly
    beta_mean <- beta[1]
    se_mean   <- se[1]
    td_mean   <- td[1]
    td_se_mean <- td_se[1]

  } else {
    # inverse-variance pooling (fixed effect IV estimator)
    vi <- se^2
    w_iv <- 1 / vi
    beta_mean <- sum(w_iv * beta) / sum(w_iv)
    se_mean <- sqrt(1 / sum(w_iv))

    vi_td <- td_se^2
    w_iv_td <- 1 / vi_td
    td_mean <- sum(w_iv_td * td) / sum(w_iv_td)
    td_se_mean <- sqrt(1 / sum(w_iv_td))
  }

  ## store aggregated results
  reduced_data[i, "beta_pca"] <- beta_mean
  reduced_data[i, "se_pca"]   <- se_mean
  reduced_data[i, "td_pca"]   <- td_mean
  reduced_data[i, "td_se_pca"] <- td_se_mean

  ## other covariates: simple means or first observation, as appropriate
  reduced_data[i, "map_scaled"]       <- mean(subset$map_scaled, na.rm = TRUE)
  reduced_data[i, "mat_scaled"]       <- mean(subset$mat_scaled, na.rm = TRUE)
  reduced_data[i, "stand_age_scaled"] <- mean(subset$stand_age_scaled, na.rm = TRUE)
  reduced_data[i, "richness_scaled"]  <- mean(subset$richness_scaled, na.rm = TRUE)
  reduced_data[i, "lon"]              <- first(subset$lon)
  reduced_data[i, "lat"]              <- first(subset$lat)
  reduced_data[i, "elev_scaled"]      <- mean(subset$elev_scaled, na.rm = TRUE)
  reduced_data[i, "ecoregion"]        <- first(subset$ecoregion)
}

## compute per-ecoregion summary (n, mean, 95% CI)
summary_df <- do.call(rbind, lapply(split(reduced_data, reduced_data$ecoregion), function(df) {
  vec <- df$beta_pca
  n <- sum(!is.na(vec))
  m <- mean(vec, na.rm = TRUE)
  se <- if(n > 1) sd(vec, na.rm = TRUE) / sqrt(n) else NA
  tval <- if(n > 1) qt(0.975, df = n - 1) else NA
  data.frame(ecoregion = df$ecoregion[1],
             n = n, mean = m,
             lower = if(!is.na(se)) m - tval * se else NA,
             upper = if(!is.na(se)) m + tval * se else NA,
             stringsAsFactors = FALSE)
}))

##--------------------------------------------------------------
## 03. Investigate tradeoff strength by ecoregion                                                        
##--------------------------------------------------------------

## fields used for block bootstrap by ecoregion
fields <- c("x", "y", "beta_pca", "ecoregion", "se_pca")

## run block bootstrap over ecoregion model
boot_ecoregion <- bootstrap_model_blocks(
  full_data = reduced_data,
  formula = beta_pca ~ 0 + ecoregion,
  fitfun = lm,
  filepath = here("data", "bootstraps", "bootstrap_coefs_ecoregion.csv"),
  weight_var = "se_pca",
  fields = fields,
  N = 400,
  nbox_per_draw = 10,
  mc.cores = 7
)

## fit baseline (uncorrected) ecoregion model and check results
ecoregion_model <- lm(beta_pca ~ 0 + ecoregion, data = reduced_data, weights = 1 / (se_pca^2))
summary(ecoregion_model)

## compute bootstrap medians and quantiles per ecoregion
for (i in 1:nrow(summary_df)) {
  summary_df$median[i] <- median(boot_ecoregion[,paste0("ecoregion", summary_df$ecoregion[i])], na.rm = TRUE)
  summary_df$lower[i] <- quantile(boot_ecoregion[,paste0("ecoregion", summary_df$ecoregion[i])], 0.025, na.rm = TRUE)
  summary_df$upper[i] <- quantile(boot_ecoregion[,paste0("ecoregion", summary_df$ecoregion[i])], 0.975, na.rm = TRUE)
}

## replace summary means with model coefficients
summary_df$mean <- coef(ecoregion_model)

## drop ecoregions with fewer than 30 observations
summary_df <- summary_df[summary_df$ecoregion %in% c("TROPICAL WET FORESTS", "SOUTHERN SEMIARID HIGHLANDS") == FALSE, ]

## order data by mean tradeoff strenght for plotting
ordered_levels <- row.names(summary_df[order(summary_df$mean), ])
summary_df$ecoregion <- factor(summary_df$ecoregion, levels = ordered_levels)

## save summary table
write.csv(summary_df, here("tables", "ecoregion_model.csv"))


## Plot results and ecoregion distribution

## build color palette
pal <- c(
  "#E64B35", "#4DBBD5", "#00A087", "#3C5488",
  "#F39B7F", "#8491B4", "#91D1C2", "#DC0000",
  "#7E6148", "#B09C85", "#0099CC", "#EEA236"   
)
pal <- rep(pal, length.out = length(ordered_levels))
names(pal) <- ordered_levels

## plot ecoregions across contiguous US
ecoplot_data <- reduced_data[reduced_data$ecoregion %in% c("TROPICAL WET FORESTS", "SOUTHERN SEMIARID HIGHLANDS") == FALSE, ]
ecoplot_data$ecoregion <- factor(ecoplot_data$ecoregion, levels = ordered_levels)
ecoplot_spatial <- st_as_sf(ecoplot_data, coords = c("lon", "lat"), crs = 4326)
ecoplot_spatial <- st_transform(ecoplot_spatial, crs = 5070)

base_map() +
  geom_sf(data = ecoplot_spatial, aes(color = ecoregion), size = 1.25, stroke = 0, shape = 15, alpha = 0.9) +
  scale_color_manual(values = pal, name = "EPA Ecoregion") +
  theme_bw() +
  guides(color = guide_legend(ncol = 1, byrow = FALSE, nrow = 12,
                              override.aes = list(size = 4))) +
  theme(legend.title = element_text(size = 12),
        legend.text = element_text(size = 10),
        legend.key.width = unit(2, "lines"))

ggsave(here("figures", "tradeoff_drivers", "tradeoff_ecoregion_map.pdf"),
       width = 10, height = 6)


## plot average tradeoff strength by ecoregion

## compute label y position below the lowest CI
y_min <- min(summary_df$lower, na.rm = TRUE)
summary_df$ypos <- y_min - 0.2   # adjust offset if needed

## create horizontal reference bands and errorbars for each ecoregion
p <- ggplot() +
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = -0.37 - 1.97*0.122, ymax = -0.37 + 1.97*0.122),
            fill = "red", alpha = 0.08) +
  geom_abline(slope = 0, intercept = -0.37, color = "red", linetype = "dashed", size = 1) +
  geom_abline(slope = 0, intercept = 0, color = "gray50", size = 1.5) +
  geom_errorbar(data = summary_df, aes(x = ecoregion, ymin = lower, ymax = upper, color = ecoregion),
                width = 0.25, size = 0.8) +
  geom_point(data = summary_df, aes(x = ecoregion, y = mean, fill = ecoregion),
             shape = 21, size = 3.5, color = "black") +
  scale_fill_manual(values = pal) +
  scale_color_manual(values = pal) +
  scale_y_continuous(limits = c(NA, NA)) +
  coord_flip() +
  labs(title = "Tradeoff Slope by Ecoregion",
       y = "Growth–Drought Tradeoff (beta_pca)",
       x = "") +
  theme_bw(base_size = 13) +
  theme(legend.position = "none",
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(size = 16, face = "bold", hjust = 0),
        axis.text.y = element_text(size = 11),
        axis.text.x = element_text(size = 10),
        panel.background = element_rect(fill = "white", colour = NA),
        plot.margin = margin(10, 12, 10, 12))

print(p)

ggsave(here("figures", "tradeoff_drivers", "tradeoff_by_ecoregion.pdf"),
       width = 12, height = 6)

##--------------------------------------------------------------
## 04. Quantify the drivers of tradeoff strength and adherence                                                        
##--------------------------------------------------------------

## log transform TD (adherence) and approximate SE on log scale (delta method)
reduced_data$log_td <- log(reduced_data$td_pca)
reduced_data$log_td_se <- reduced_data$td_se_pca / reduced_data$td_pca

## weighted linear model for log(TD) ~ predictors (weights = 1 / var)
td_model <- lm(log_td ~
                  map_scaled +
                  mat_scaled +
                  stand_age_scaled,
                data = reduced_data,
                weights = 1 / (log_td_se^2))
summary(td_model)

## block bootstrap for TD drivers
boot_td <- bootstrap_model_blocks(
  full_data = reduced_data,
  formula = log_td ~ map_scaled + mat_scaled + stand_age_scaled,
  fitfun = lm,
  filepath = here("data", "bootstraps", "bootstrap_coefs_tddrivers.csv"),
  weight_var = "log_td_se",
  fields = c("x", "y", "log_td", "map_scaled", "mat_scaled", "stand_age_scaled", "log_td_se"),
  N = 400,
  nbox_per_draw = 10,
  mc.cores = 7
)

## read bootstrap output
boot_td <- read.csv(here("data", "bootstraps", "bootstrap_coefs_tddrivers.csv"))
colnames(boot_td)[1] <- "(Intercept)"

## summarize bootstrap distribution
td_summary <- data.frame(coef = c("(Intercept)", "map_scaled", "mat_scaled", "stand_age_scaled"),
                         mean = NA_real_,
                         lower = NA_real_,
                         upper = NA_real_)

for (i in 1:nrow(td_summary)) {
  td_summary$mean[i] <- median(boot_td[,td_summary$coef[i]])
  td_summary$lower[i] <- quantile(boot_td[,td_summary$coef[i]], 0.025)
  td_summary$upper[i] <- quantile(boot_td[,td_summary$coef[i]], 0.975)
}

write.csv(td_summary, here("tables", "td_drivers_model.csv"))


## fit model for tradeoff strength (beta)
beta_model <- lm(beta_pca ~
                  map_scaled +
                  mat_scaled +
                  stand_age_scaled,
                 data = reduced_data,
                 weights = 1 / (se_pca^2))
summary(beta_model)

## block bootstrap for beta drivers
boot_beta <- bootstrap_model_blocks(
  full_data = reduced_data,
  formula = beta_pca ~ map_scaled + mat_scaled + stand_age_scaled,
  fitfun = lm,
  filepath = here("data", "bootstraps", "bootstrap_coefs_betadrivers.csv"),
  weight_var = "se_pca",
  fields = c("x", "y", "beta_pca", "map_scaled", "mat_scaled", "stand_age_scaled", "se_pca"),
  N = 400,
  nbox_per_draw = 10,
  mc.cores = 7
)

## read bootstrap output
boot_beta <- read.csv(here("data", "bootstraps", "bootstrap_coefs_betadrivers.csv"))
colnames(boot_beta)[1] <- "(Intercept)"

## summarize bootstrap distribution
beta_summary <- data.frame(coef = c("(Intercept)", "map_scaled", "mat_scaled", "stand_age_scaled"),
                           mean = NA_real_,
                           lower = NA_real_,
                           upper = NA_real_)

for (i in 1:nrow(beta_summary)) {
  beta_summary$mean[i] <- median(boot_beta[,beta_summary$coef[i]])
  beta_summary$lower[i] <- quantile(boot_beta[,beta_summary$coef[i]], 0.025)
  beta_summary$upper[i] <- quantile(boot_beta[,beta_summary$coef[i]], 0.975)
}

write.csv(beta_summary, here("tables", "beta_drivers_model_summary.csv"))


## Quantify residual spatial autocorrelation in baseline models

## project lat/lon to equal-area coords (5070) for variogram distances
df <- reduced_data
df_coords <- data.frame(x = df$lon, y = df$lat)
df_sf <- st_as_sf(df_coords, coords = c("x", "y"), crs = 4326)
df_sf <- st_transform(df_sf, crs = 5070)
coords_proj <- st_coordinates(df_sf)
df$x <- coords_proj[,1]
df$y <- coords_proj[,2]

## empirical variograms
resid_df <- data.frame(z = residuals(beta_model), x = df$x, y = df$y)
v1 <- variogram(z ~ 1, data = resid_df, locations = ~x + y, cutoff = 500000, width = 500000/20)

resid_df <- data.frame(z = residuals(td_model), x = df$x, y = df$y)
v2 <- variogram(z ~ 1, data = resid_df, locations = ~x + y, cutoff = 500000, width = 500000/20)

resid_df <- data.frame(z = residuals(ecoregion_model), x = df$x, y = df$y)
v3 <- variogram(z ~ 1, data = resid_df, locations = ~x + y, cutoff = 500000, width = 500000/20)

v1_df <- transform(v1, model = "Beta model")
v2_df <- transform(v2, model = "TD model")
v3_df <- transform(v3, model = "Ecoregion model")

## plot variograms side-by-side
p1 <- ggplot(v1_df, aes(x = dist / 1000, y = gamma)) +
  geom_point(aes(size = np), alpha = 0.7) +
  geom_line(alpha = 0.7) +
  scale_y_continuous(limits = c(0, 1.3), expand = c(0,0)) +   # ← left-panel limit
  scale_size_continuous(name = "Pair count") +
  labs(
    x = "Distance (km)",
    y = "Semivariance",
    title = "Strength model"
  ) +
  theme_bw()

p2 <- ggplot(v2_df, aes(x = dist / 1000, y = gamma)) +
  geom_point(aes(size = np), alpha = 0.7) +
  geom_line(alpha = 0.7) +
  scale_y_continuous(limits = c(0, 1.1), expand = c(0,0)) +   # ← right-panel limit
  scale_size_continuous(name = "Pair count") +
  labs(
    x = "Distance (km)",
    y = "Semivariance",
    title = "Adherence model"
  ) +
  theme_bw()

p3 <- ggplot(v3_df, aes(x = dist / 1000, y = gamma)) +
  geom_point(aes(size = np), alpha = 0.7) +
  geom_line(alpha = 0.7) +
  scale_y_continuous(limits = c(0, 1.8), expand = c(0,0)) +   # ← right-panel limit
  scale_size_continuous(name = "Pair count") +
  labs(
    x = "Distance (km)",
    y = "Semivariance",
    title = "Ecoregion model"
  ) +
  theme_bw()

p1 + p2 + p3 + plot_layout(ncol = 2)

ggsave(here("figures", "tradeoff_drivers", "variograms_tradeoff_drivers.png"),
       width = 9, height = 7)


## plot estimated effects, incorporating bootstrap uncertainties
p1 <- plot_effect_simple("map_scaled", beta_model, full_data, boot_beta,
                  orig_var = "map",
                  xlab = "",
                  ylab = "Tradeoff Slope",
                  ylim = c(-2.1, 1.1),
                  xlim = c(250, 2000),
                  n = 300,
                  span = 0.5,
                  fill = "#94d2bd",
                  line_color = "#015f73")
p1

p2 <- plot_effect_simple("mat_scaled", beta_model, full_data, boot_beta,
                  orig_var = "mat",
            xlab = "",
            ylab = "",
            n = 300,
            ylim = c(-2.1, 1.1),
            span = 1,
            fill = "#94d2bd",
            line_color = "#015f73")
p2

p3 <- plot_effect_simple("stand_age_scaled", beta_model, full_data, boot_beta,
                  orig_var = "stand_age",
                  xlab = "",
                  ylab = "",
                  ylim = c(-2.2, 1.1),
                  xlim = c(0, 300),
                  span = 1,
                  n = 300,
                  fill = "#94d2bd",
                  line_color = "#015f73")
p3

p4 <- plot_effect_simple("map_scaled", td_model, full_data, boot_td,
                         orig_var = "map",
                         xlab = "Mean Annual Precipitation (mm)",
                         ylab = "Tradeoff Adherence",
                         ylim = c(0, 0.7),
                         xlim = c(250, 2000),
                         n = 100,
                         span = 1,
                         delog = TRUE,
                         fill = "#c7b5d8",
                         line_color = "#6a4d94")
p4

p5 <- plot_effect_simple("mat_scaled", td_model, full_data, boot_td,
                  orig_var = "mat",
                  xlab = "Mean Annual Temperature (C)",
                  ylab = "",
                  ylim = c(0, 0.7),
                  delog = TRUE,
                  span = 1,
                  fill = "#c7b5d8",
                  line_color = "#6a4d94")
p5

p6 <- plot_effect_simple("stand_age_scaled", td_model, full_data, boot_td,
                  orig_var = "stand_age",
                  xlab = "Stand Age (years)",
            ylab = "",
            ylim = c(0, 0.7),
            n = 300,
            span = 0.1,
            xlim = c(0, 300),
            delog = TRUE,
            fill = "#c7b5d8",
            line_color = "#6a4d94")
p6

p1 + p2 + p3 + p4 + p5 + p6 + plot_layout(ncol = 3)

ggsave(here("figures", "tradeoff_drivers", "tradeoff_drivers_effects.pdf"),
       width = 11, height = 7, dpi = 300)

