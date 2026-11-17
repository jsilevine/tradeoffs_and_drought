##--------------------------------------------------------------
## 06_tradeoff_overview.r: Overview of growth-drought tradeoff (Figure 1)
##
## author: Jacob Levine; email: jacob.levine@utah.edu
##--------------------------------------------------------------

##--------------------------------------------------------------
## 00. Load libraries and helper functions
##--------------------------------------------------------------

## core dependencies for plotting, PCA, and spatial data
library(here)
library(data.table)
library(ggplot2)
library(ggfortify)
library(patchwork)
library(ggrepel)
library(FactoMineR)
library(sf)

## project-specific helpers
source(here("src", "r", "trait_processing.r"))
source(here("src", "r", "maps.r"))

##--------------------------------------------------------------
## 01. Load and preprocess data
##--------------------------------------------------------------

## tradeoff estimates at plot level
data <- fread(here("data", "tradeoffs", "tradeoffs_full.csv"))
data <- data[!is.na(beta_pca)]

## FIA metadata and trait tables
states  <- fread(here("data", "state_list.csv"))
gm_data <- fread(here("data", "gm", "gm_full.csv"))
traits  <- fread(here("data", "traits_imputed_merged.csv"))

## standardize species code and flip p50 sign
setnames(traits, 2, "SPCD")
traits[, p50 := p50 * -1]
traits[, P50 := P50 * -1]

##--------------------------------------------------------------
## 02. Select representative example plots
##--------------------------------------------------------------

## select plots spanning weak to strong tradeoffs and different regions
id1 <- which(data$beta_pca > -0.2 & data$beta_pca < 0 &
             data$LON < -98 & data$richness > 7)[5]

id2 <- which(data$beta_pca > -2 & data$beta_pca < -1 &
             data$LON > -98 & data$richness > 10)[200]

id3 <- which(data$beta_pca > -1 & data$beta_pca < -0.2 &
             data$LON < -98 & data$richness > 9)[2]

plot_ids <- data$PLT_CN[c(id1, id2, id3)]
gm_data_sub <- gm_data[PLT_CN %in% plot_ids]

##--------------------------------------------------------------
## 03. Project species traits into PCA space
##--------------------------------------------------------------

## PCA objects estimated elsewhere
drought_pca <- readRDS(here("data", "model_data", "drought_pca.rds"))
growth_pca  <- readRDS(here("data", "model_data", "growth_pca.rds"))

## attach drought and growth PC scores to trait table
project_pca(traits, drought_pca,
            c("p50", "rdmax"),
            c("drought_pc1", "drought_pc2"))

project_pca(traits, growth_pca,
            c("Amax", "gsmax", "LeafN", "SLA"),
            c("growth_pc1", "growth_pc2"))

##--------------------------------------------------------------
## 04. Aggregate traits to plot level
##--------------------------------------------------------------

## variables retained for plot-level summaries
trait_vars <- c("species", "p50", "Amax", "gsmax", "rdmax",
                "SLA", "LeafN", "growth_pc1", "drought_pc2")

data_list <- list()

for (i in seq_len(nrow(gm_data_sub))) {

  ## read state-specific FIA tree table
  state_name <- states[STATECD == gm_data_sub$STATECD[i], state]
  tree_path  <- here("data", "fia", state_name,
                     paste0(state_name, "_TREE.csv"))

  tree_data <- fread(tree_path)
  colnames(tree_data)[1] <- "TRE_CN"

  ## restrict to fields needed for trait aggregation
  tree_data <- tree_data[, .(TRE_CN, PLT_CN, PREV_TRE_CN,
                             STATUSCD, TPA_UNADJ, SPCD)]
  
  ## extract plot-level tree data and merge with species traits
  tl <- pull_tree_data(gm_data_sub[i, ], tree_data, init = TRUE)
  traits_trees <- merge(tl, traits, by = "SPCD", all.x = TRUE)

  ## aggregate to plot using TPA-weighted means
  traits_sum <- aggregate_traits(traits_trees, trait_vars)
  traits_sum$plot_number <- gm_data_sub$PLT_CN[i]
  
  data_list[[i]] <- traits_sum
}

plot_data <- do.call(rbind, data_list)

##--------------------------------------------------------------
## 05. Plot tradeoff structure for example plots
##--------------------------------------------------------------

## color palette for highlighted plots
highlight_cols <- c("#0072B2", "#D55E00", "#009E73")

## plot growth vs drought PCA, weighted by basal area
make_tradeoff_plot <- function(plot_id, color) {
  ggplot(plot_data[plot_data$plot_number == plot_id, ],
         aes(x = drought_pc2, y = growth_pc1, size = tpa)) +
    geom_smooth(method = "lm", aes(weight = tpa),
                color = "black", se = TRUE, linewidth = 2) +
    geom_point(alpha = 1.0, color = color) +
    scale_size_area(max_size = 20) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(limits = c(-4.6, 3)) +
    theme_bw() +
    labs(x = "Drought Tolerance (PC2)",
         y = "Growth Capacity (PC1)") +
    theme(legend.position = "none")
}

p1 <- make_tradeoff_plot(plot_ids[2], highlight_cols[1])
p2 <- make_tradeoff_plot(plot_ids[3], highlight_cols[3])
p3 <- make_tradeoff_plot(plot_ids[1], highlight_cols[2])

combined <- p1 + p2 + p3 + plot_layout(ncol = 3)

ggsave(here("figures", "overview", "overview_tradeoff_examples.pdf"),
       combined, width = 12, height = 4)

##--------------------------------------------------------------
## 06. Map spatial locations of example plots
##--------------------------------------------------------------

## build sf object of all plots (excluding AK/HI)
plot_sf <- st_as_sf(
  unique(gm_data[!(gm_data$STATECD %in% c(2, 15)),
                 c("PLT_CN", "LON", "LAT")]),
  coords = c("LON", "LAT"), crs = 4326
)

## subset highlighted plots
plot_sf_sub <- plot_sf[plot_sf$PLT_CN %in% plot_ids, ]
plot_sf <- plot_sf[!duplicated(plot_sf[, c("geometry")]), ]

plot_sf_sub$PLT_CN <- factor(plot_sf_sub$PLT_CN,
                             levels = plot_ids)

## project to Albers Equal Area
plot_sf <- st_transform(plot_sf, 5070)
plot_sf_sub <- st_transform(plot_sf_sub, 5070)

base_map() +
  geom_sf(data = plot_sf,
          alpha = 0.2, color = "grey40", size = 0.6) +
  geom_sf(data = plot_sf_sub,
          aes(fill = PLT_CN),
          shape = 21, size = 9,
          color = "black", stroke = 1) +
  scale_fill_manual(values = c(highlight_cols[2],
                               highlight_cols[1],
                               highlight_cols[3])) +
  theme()

ggsave(here("figures", "overview", "overview_tradeoff_map.pdf"),
       width = 8, height = 5)

##--------------------------------------------------------------
## 07. PCA biplots for drought and growth axes
##--------------------------------------------------------------

## reload traits to ensure full species coverage
species_traits <- fread(here("data", "traits_imputed_merged.csv"))
setnames(species_traits, 2, "SPCD")

## flip p50 sign for interpretability
species_traits$p50 <- species_traits$p50 * -1
species_traits$P50 <- species_traits$P50 * -1

## drought PCA biplot
p1 <- plot_pca(drought_pca, species_traits,
               c("p50", "rdmax"),
               c("P50", "Max. Rooting Depth"))
ggsave(here("figures", "overview", "overview_drought_pca.pdf"),
       width = 7, height = 5)

## growth PCA biplot
p2 <- plot_pca(growth_pca, species_traits,
               c("Amax", "gsmax", "LeafN", "SLA"),
               c("Amax", "gs max", "Leaf N", "SLA"))
ggsave(here("figures", "overview", "overview_growth_pca.pdf"),
       width = 7, height = 5)
