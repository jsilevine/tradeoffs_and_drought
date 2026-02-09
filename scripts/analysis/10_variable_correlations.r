##--------------------------------------------------------------
## Quantify correlations between variables used in drought models
##
## author: jacob levine; jacob.levine@utah.edu
##--------------------------------------------------------------

library(ggplot2)
library(here)

data <- fread(here("data", "data_for_modeling.csv"))

data <- data[,c("ba", "elev", "map", "mat", "stand_age", "beta_pca", "td_pca", "cwm_growth_pc1", "cwm_drought_pc2", "range_growth_pc1", "range_drought_pc2",
                "drought_strength", "prop_drought")]

vars <- colnames(data)
corr_mat <- cor(data, use = "pairwise.complete.obs")

corr_df <- data.frame(
  var1 = rep(vars, each = length(vars)),
  var2 = rep(vars, times = length(vars)),
  r    = as.vector(corr_mat),
  stringsAsFactors = FALSE
)

idx <- which(lower.tri(corr_mat, diag = TRUE), arr.ind = TRUE)
keep_idx <- (match(corr_df$var1, vars) - 1) * length(vars) + match(corr_df$var2, vars)
corr_df <- corr_df[keep_idx %in% ( (idx[,1]-1)*length(vars) + idx[,2] ), ]
corr_df[corr_df$var1 == corr_df$var2, "r"] <- NA  ## remove diagonal
corr_df <- na.omit(corr_df)

corr_df$var1 <- factor(corr_df$var1, levels = vars)
corr_df$var2 <- factor(corr_df$var2, levels = rev(vars))

ggplot(corr_df, aes(x = var1, y = var2, fill = r)) +
  geom_tile(color = "white", size = 0.4) +
  geom_text(aes(label = sprintf("%.2f", r)), size = 3) +
  scale_fill_gradient2(low = "#d73027", mid = "white", high = "#4575b4", limits = c(-1, 1), name = "r") +
  coord_fixed() +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank(),
    axis.title = element_blank()
  )
