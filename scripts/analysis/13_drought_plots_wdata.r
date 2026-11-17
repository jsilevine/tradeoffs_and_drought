library(here)
library(patchwork)

source(here("src", "r", "plot_growth_predictions.r"))

growth_model <- readRDS(here("data", "model_objects", "growth_model.rds"))
mortality_model <- readRDS(here("data", "model_objects", "mortality_model.rds"))

data <- read.csv(here("data", "data_for_modeling.csv"))
cols <- c("year", "lon", "lat", "year_centered",
          "mort_tpa_pct", 
          "richness_scaled",
          "ba", "ba_scaled",
          "drought_strength", "drought_strength_scaled",
          "prop_drought", "prop_drought_scaled",
          "beta_pca", "beta_pca_scaled",
          "td_pca", "td_pca_scaled",
          "growth_pct", "growth_pct_scaled",
          "map", "map_scaled",
          "mat", "mat_scaled",
          "cwm_drought_pc2", "cwm_drought_pc2_scaled",
          "cwm_growth_pc1", "cwm_growth_pc1_scaled",
          "range_drought_pc2", "range_drought_pc2_scaled",
          "range_growth_pc1", "range_growth_pc1_scaled",
          "elev", "elev_scaled",
          "stand_age", "stand_age_scaled")

data <- data[cols]
data <- na.omit(data)

data$growth_pct_pct <- data$growth_pct * 100
p1 <- plot_one_gam_prediction(growth_model, data, "beta_pca", lower_bound_qnt = 0.05)
p1$layers <- c(
  list(geom_point(data = data,
                  aes(x = beta_pca, y = growth_pct_pct),
                  alpha = 0.03, size = 1.5)),
  p1$layers
)
p1$layers$geom_line$aes_params$colour <- "#2ecc71"
p1$layers$geom_ribbon$aes_params$fill <- "#2ecc71"
p1$layers$geom_ribbon$aes_params$alpha <- 0.4
p1 <- p1 + scale_y_continuous(limits = c(-1, 8), expand = c(0,0)) + scale_x_continuous(limits = c(-3.7, 0.7), expand = c(0,0))
p1 <- p1 + theme(legend.position = "none")
p1

data[data$ba > 200, "ba"] <- 200
data[data$ba < 50, "ba"] <- 50
p2 <- plot_one_gam_prediction(growth_model, data, "td_pca", interaction = "ba", lower_bound_qnt = 0.05)
p2$layers <- c(
  list(geom_point(data = data,
                  aes(x = td_pca, y = growth_pct_pct, color = ba),
                  alpha = 0.03, size = 1.5, show.legend = FALSE)),
  p2$layers
)
p2 <- p2 + 
  scale_color_gradientn(
    colors = met.brewer("Archambault", n = 5, type = "continuous"),
    name   = "Basal area (m²/ha)",
    limits = range(data$ba, na.rm = TRUE)
  ) +
  scale_fill_gradientn(
    colors = met.brewer("Archambault", n = 5, type = "continuous"),
    name   = "Basal area (m²/ha)",
    limits = range(data$ba, na.rm = TRUE)
  )
p2$layers$geom_ribbon$aes_params$alpha <- 0.4
p2$layers$geom_ribbon$aes_params$color <- NA
p2 <- p2 + scale_y_continuous(limits = c(-1, 8), expand = c(0,0)) + scale_x_continuous(limits = c(0.05, 0.55), expand = c(0,0))
p2 <- p2 + theme(legend.position = "none") + xlab("Tradeoff Deviation")
p2

p3 <- plot_one_gam_prediction(growth_model, data, "drought_strength", interaction = "ba", lower_bound_qnt = 0.05)
p3$layers <- c(
  list(geom_point(data = data,
                  aes(x = drought_strength, y = growth_pct_pct, color = ba),
                  alpha = 0.03, size = 1.5, show.legend = FALSE)),
  p3$layers
)
p3 <- p3 + 
  scale_color_gradientn(
    colors = met.brewer("Archambault", n = 5, type = "continuous"),
    name   = "Basal area (m²/ha)",
    limits = range(data$ba, na.rm = TRUE)
  ) +
  scale_fill_gradientn(
    colors = met.brewer("Archambault", n = 5, type = "continuous"),
    name   = "Basal area (m²/ha)",
    limits = range(data$ba, na.rm = TRUE)
  )
p3$layers$geom_ribbon$aes_params$alpha <- 0.4
p3$layers$geom_ribbon$aes_params$color <- NA
p3 <- p3 + scale_y_continuous(limits = c(-1, 8), expand = c(0,0)) + scale_x_continuous(limits = c(-5, -2.25), expand = c(0,0))
p3 <- p3 + theme(legend.position = "none")
p3

p4 <- plot_one_gam_prediction(growth_model, data, "cwm_drought_pc2", lower_bound_qnt = 0.05)
p4$layers <- c(
  list(geom_point(data = data,
                  aes(x = cwm_drought_pc2, y = growth_pct_pct),
                  alpha = 0.03, size = 1.5)),
  p4$layers
)
p4$layers$geom_line$aes_params$colour <- "#2ecc71"
p4$layers$geom_ribbon$aes_params$fill <- "#2ecc71"
p4$layers$geom_ribbon$aes_params$alpha <- 0.4
p4 <- p4 + scale_y_continuous(limits = c(-1, 8), expand = c(0,0)) + scale_x_continuous(limits = c(-0.5, 0.5), expand = c(0,0))
p4 <- p4 + theme(legend.position = "none")
p4


p5 <- plot_one_gam_prediction(growth_model, data, "cwm_growth_pc1", lower_bound_qnt = 0.05)
p5$layers <- c(
  list(geom_point(data = data,
                  aes(x = cwm_growth_pc1, y = growth_pct_pct),
                  alpha = 0.03, size = 1.5)),
  p5$layers
)
p5$layers$geom_line$aes_params$colour <- "#2ecc71"
p5$layers$geom_ribbon$aes_params$fill <- "#2ecc71"
p5$layers$geom_ribbon$aes_params$alpha <- 0.4
p5 <- p5 + scale_y_continuous(limits = c(-1, 8), expand = c(0,0)) + scale_x_continuous(limits = c(-2.2, 1.6), expand = c(0,0))
p5 <- p5 + theme(legend.position = "none")
p5

p6 <- plot_one_gam_prediction(growth_model, data, "prop_drought", lower_bound_qnt = 0.01, upper_bound_qnt = 0.99)
p6$layers <- c(
  list(geom_point(data = data,
                  aes(x = prop_drought, y = growth_pct_pct),
                  alpha = 0.03, size = 1.5)),
  p6$layers
)
p6$layers$geom_line$aes_params$colour <- "#2ecc71"
p6$layers$geom_ribbon$aes_params$fill <- "#2ecc71"
p6$layers$geom_ribbon$aes_params$alpha <- 0.4
p6 <- p6 + scale_y_continuous(limits = c(-1, 8), expand = c(0,0)) + scale_x_continuous(limits = c(0.25, 0.8), expand = c(0,0))
p6 <- p6 + theme(legend.position = "none") + xlab("Drought Burden")
p6

p1 + p2 + p3 + p4 + p5 + p6 + plot_layout(ncol = 3) +
  plot_annotation(tag_levels = "A") +
  theme(legend.position = "none")

ggsave(here("figures","model_predictions", "growth", "growth_predictions_wdata.png"), width = 16, height = 10)




p1 <- plot_one_bezi_prediction(mortality_model, data, "beta_pca")

p1$layers <- c(
  list(geom_point(data = data,
                  aes(x = beta_pca, y = mort_tpa_pct),
                  alpha = 0.03, size = 1.5)),
  p1$layers
)
p1$layers$geom_line$aes_params$colour <- "#2ecc71"
p1$layers$geom_ribbon$aes_params$fill <- "#2ecc71"
p1$layers$geom_ribbon$aes_params$alpha <- 0.4
p1 <- p1 + scale_y_continuous(limits = c(0, 0.02), expand = c(0,0)) + scale_x_continuous(limits = c(-4, 0.85), expand = c(0,0))
p1 <- p1 + theme(legend.position = "none")
p1

p2 <- plot_one_bezi_prediction(mortality_model, data, "td_pca")

p2$layers <- c(
  list(geom_point(data = data,
                  aes(x = td_pca, y = mort_tpa_pct),
                  alpha = 0.03, size = 1.5, show.legend = FALSE)),
  p2$layers
)
p2$layers$geom_ribbon$aes_params$alpha <- 0.4
p2$layers$geom_ribbon$aes_params$color <- NA
p2$layers$geom_line$aes_params$colour <- "#2ecc71"
p2$layers$geom_ribbon$aes_params$fill <- "#2ecc71"
p2 <- p2 + scale_y_continuous(limits = c(0, 0.02), expand = c(0,0)) + scale_x_continuous(limits = c(0.05, 0.55), expand = c(0,0))
p2 <- p2 + theme(legend.position = "none") + xlab("Tradeoff Deviation")
p2

p3 <- plot_one_bezi_prediction(mortality_model, data, "drought_strength", interaction = "ba")

p3$layers <- c(
  list(geom_point(data = data,
                  aes(x = drought_strength, y = mort_tpa_pct, color = ba),
                  alpha = 0.03, size = 1.5, show.legend = FALSE)),
  p3$layers
)
p3 <- p3 + 
  scale_color_gradientn(
    colors = met.brewer("Archambault", n = 5, type = "continuous"),
    name   = "Basal area (m²/ha)",
    limits = range(data$ba, na.rm = TRUE)
  ) +
  scale_fill_gradientn(
    colors = met.brewer("Archambault", n = 5, type = "continuous"),
    name   = "Basal area (m²/ha)",
    limits = range(data$ba, na.rm = TRUE)
  )
p3$layers$geom_ribbon$aes_params$alpha <- 0.4
p3$layers$geom_ribbon$aes_params$color <- NA
p3 <- p3 + scale_y_continuous(limits = c(0, 0.02), expand = c(0,0)) + scale_x_continuous(limits = c(-5, -2.25), expand = c(0,0))
p3 <- p3 + theme(legend.position = "none")
p3

p4 <- plot_one_bezi_prediction(mortality_model, data, "cwm_drought_pc2")

p4$layers <- c(
  list(geom_point(data = data,
                  aes(x = cwm_drought_pc2, y = mort_tpa_pct),
                  alpha = 0.03, size = 1.5)),
  p4$layers
)
p4$layers$geom_line$aes_params$colour <- "#2ecc71"
p4$layers$geom_ribbon$aes_params$fill <- "#2ecc71"
p4$layers$geom_ribbon$aes_params$alpha <- 0.4
p4 <- p4 + scale_y_continuous(limits = c(0, 0.02), expand = c(0,0)) + scale_x_continuous(limits = c(-0.5, 0.5), expand = c(0,0))
p4 <- p4 + theme(legend.position = "none")
p4


p5 <- plot_one_bezi_prediction(mortality_model, data, "cwm_growth_pc1")

p5$layers <- c(
  list(geom_point(data = data,
                  aes(x = cwm_growth_pc1, y = mort_tpa_pct),
                  alpha = 0.03, size = 1.5)),
  p5$layers
)
p5$layers$geom_line$aes_params$colour <- "#2ecc71"
p5$layers$geom_ribbon$aes_params$fill <- "#2ecc71"
p5$layers$geom_ribbon$aes_params$alpha <- 0.4
p5 <- p5 + scale_y_continuous(limits = c(0, 0.02), expand = c(0,0)) + scale_x_continuous(limits = c(-2.2, 1.6), expand = c(0,0))
p5 <- p5 + theme(legend.position = "none")
p5

p6 <- plot_one_bezi_prediction(mortality_model, data, "prop_drought")

p6$layers <- c(
  list(geom_point(data = data,
                  aes(x = prop_drought, y = mort_tpa_pct),
                  alpha = 0.03, size = 1.5)),
  p6$layers
)
p6$layers$geom_line$aes_params$colour <- "#2ecc71"
p6$layers$geom_ribbon$aes_params$fill <- "#2ecc71"
p6$layers$geom_ribbon$aes_params$alpha <- 0.4
p6 <- p6 + scale_y_continuous(limits = c(0, 0.02), expand = c(0,0)) + scale_x_continuous(limits = c(0.25, 0.6), expand = c(0,0))
p6 <- p6 + theme(legend.position = "none") + xlab("Drought Burden")
p6

p1 + p2 + p3 + p4 + p5 + p6 + plot_layout(ncol = 3) +
  plot_annotation(tag_levels = "A") +
  theme(legend.position = "none")

ggsave(here("figures","model_predictions", "mortality", "mortality_predictions_wdata.png"), width = 16, height = 10)
