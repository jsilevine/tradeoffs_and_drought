##---------------------------------------------------------------
## 06_impute_traits.r: Impute mising trait data using phylogenetic relationships
##
## author: jacob levine; jacob.levine@utah.edu
##---------------------------------------------------------------

##--------------------------------------------------------------
## 00. Load libraries                                                       
##--------------------------------------------------------------

library(devtools)

## install BiocManager, if not already installed
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.19")
BiocManager::install(c("ggtree"), force = TRUE)

## install TreVol package from github
devtools::install_github("pablosanchezmart/TrEvol")

## load libraries
invisible(lapply(c("here", "TrEvol", "ape", "rotl", "ggtree", "stringr", "dplyr",
                   "tidytree", "ggplot2", "cowplot"), library, character.only = TRUE))

source(here("src", "r", "plot_trait_data.r"))

##--------------------------------------------------------------
## 01. Load and prep data                                                       
##-------------------------------------------------------------

traits <- read.csv(here("data", "full_traits.csv"))
species_list <- read.csv(here("data", "fulldata_specieslist.csv"))
clade <- read.csv(here("data", "clade.csv"))

sum(duplicated(traits$species))
## rename columns for clarity
colnames(clade) <- c("species", "clade")

clade <- clade[!duplicated(clade$species),]
species_list <- species_list[!duplicated(species_list$species),]

## merge clade and species list information with trait values
traits <- merge(traits, clade, by = "species", all.x = TRUE, all.y = FALSE)
traits <- merge(traits, species_list[, c("species", "fia_id")], all.x = TRUE, all.y = FALSE)

## match species names with open tree of life
ott_info <- rotl::tnrs_match_names(names = species_list$species)

## remove taxa not found in OTL
ott_info <- ott_info[!(ott_info$ott_id %in% c(3915205, 943924)),]

## construct a tree
open_tree <- rotl::tol_induced_subtree(ott_ids = ott_info[,"ott_id"])

## clean tip labels for better matching with trait data
for (i in 1:length(open_tree$tip.label)) {
  genus <- strsplit(open_tree$tip.label[i], "_")[[1]][1]
  species <- strsplit(open_tree$tip.label[i], "_")[[1]][2]
  open_tree$tip.label[i] <- paste0(genus, "_", species)
}

## reformat trait data to match tree labels
for (i in 1:length(traits$species)) {
  genus <- strsplit(traits$species[i], " ")[[1]][1]
  species <- strsplit(traits$species[i], " ")[[1]][2]
  traits[i, "species"] <- paste0(genus, "_", species)
}

## generate plots for each trait
p1 <- plot_trait_coverage(open_tree, traits, "p50")
p2 <- plot_trait_coverage(open_tree, traits, "Amax")
p3 <- plot_trait_coverage(open_tree, traits, "ssd")

plot_grid(p1, p2, p3, nrow = 1)

ggsave(here("figures", "supplemental", "trait_coverage.pdf"))
ggsave(here("figures", "supplemental", "trait_coverage.png"))

## trait correlations
p1a <- ggplot(traits, aes(x = ssd_mean, y = p50_mean, color = clade)) +
  geom_point(size = 2.5) +
  scale_color_brewer(type = "qual", palette = 6) +
  geom_smooth(method = "lm", linewidth = 2, fullrange = TRUE) +
  scale_x_continuous(expand = c(0,0), limits = c(0.3,NA)) +
  theme_bw() +
  theme(legend.position = "none")
p1a

p2a <- ggplot(traits, aes(x = Amax_mean, y = p50_mean, color = clade)) +
  geom_point(size = 2.5) +
  scale_color_brewer(type = "qual", palette = 6) +
  geom_smooth(method = "lm", linewidth = 2, fullrange = TRUE) +
  scale_x_continuous(expand = c(0,0), limits = c(3,NA)) +
  theme_bw() +
  theme(legend.position = "none")
p2a

p3a <- ggplot(traits, aes(x = ssd_mean, y = Amax_mean, color = clade)) +
  geom_point(size = 2.5) +
  scale_color_brewer(type = "qual", palette = 6) +
  geom_smooth(method = "lm", linewidth = 2, fullrange = TRUE) +
  scale_x_continuous(expand = c(0,0), limits = c(NA,NA)) +
  theme_bw() +
  theme(legend.position = "none")
p3a

plot_grid(p1a, p3a, p2a, nrow = 1)

ggsave(here("figures", "supplemental", "trait_correlations.png"))
ggsave(here("figures", "supplemental", "trait_correlations.pdf"))


##---------------------------------------------------------------
## 02. Impute missing trait data
##---------------------------------------------------------------

## prune tree to include only species with trait data
pruned_tree <- drop.tip(open_tree, open_tree$tip.label[!(open_tree$tip.label %in% traits[!is.na(traits$Amax_mean) | !is.na(traits$p50_mean) | !is.na(traits$ssd_mean), "species"])])

## prep pruned trait data
pruned_traits <- traits[traits$species %in% pruned_tree$tip.label, ]
pruned_traits <- pruned_traits[order(pruned_traits$species),]

## compute branch lengths and make tree ultrametric
pruned_tree <- compute.brlen(pruned_tree)
pruned_tree <- chronos(pruned_tree)

## impute missing p50 and ssd data
imp_p50ssd <- imputeTraits(variables_to_impute = c("p50_mean", "ssd_mean"),
                             dataset = pruned_traits,
                             phylogeny = pruned_tree,
                             terminal_taxa = "species",
                             number_clusters = 4)

## merge imputed data
fd <- merge(pruned_traits, imp_p50ssd$round3$ximp,
            suffixes = c("", "_imputed"),
            by = "species",
            all.x = TRUE, all.y = FALSE)

## impute missing Amax and ssd data
imp_Amaxssd <- imputeTraits(variables_to_impute = c("Amax_mean", "ssd_mean"),
                             dataset = pruned_traits,
                             phylogeny = pruned_tree,
                             terminal_taxa = "species",
                             number_clusters = 4)

## merge additional imputed data
fd <- merge(fd, imp_Amaxssd$round3$ximp[,c("species", "Amax_mean")],
            suffixes = c("", "_imputed"),
            by = "species",
            all.x = T)

## clean up column naming and structure
colnames(fd) <- c("species", "p50_mean", "p50_database", "Amax_mean",
                  "Amax_database", "ssd_mean", "ssd_database", "clade", "fia_id", "p50", "ssd", "Amax")

fd <- fd[,c("species", "fia_id", "clade", "p50", "p50_database", "Amax",
                  "Amax_database", "ssd", "ssd_database")]

## write imputed data to CSV
write.csv(fd, here("data", "traits_imputed.csv"), row.names = FALSE)

##--------------------------------------------------------------
## 03. Visualize correlations with imputed traits                                                       
##--------------------------------------------------------------

p1a <- ggplot(fd, aes(x = ssd, y = p50, color = clade)) +
  geom_point(size = 2.5) +
  scale_color_brewer(type = "qual", palette = 6) +
  geom_smooth(method = "lm", linewidth = 2, fullrange = TRUE) +
  scale_x_continuous(expand = c(0,0), limits = c(0.3,NA)) +
  theme_bw() +
  theme(legend.position = "none")
p1a

p2a <- ggplot(fd, aes(x = Amax, y = p50, color = clade)) +
  geom_point(size = 2.5) +
  scale_color_brewer(type = "qual", palette = 6) +
  geom_smooth(method = "lm", linewidth = 2, fullrange = TRUE) +
  scale_x_continuous(expand = c(0,0), limits = c(3,NA)) +
  theme_bw() +
  theme(legend.position = "none")
p2a

p3a <- ggplot(fd, aes(x = ssd, y = Amax, color = clade)) +
  geom_point(size = 2.5) +
  scale_color_brewer(type = "qual", palette = 6) +
  geom_smooth(method = "lm", linewidth = 2, fullrange = TRUE) +
  scale_x_continuous(expand = c(0,0), limits = c(NA,NA)) +
  theme_bw() +
  theme(legend.position = "none")
p3a

plot_grid(p1a, p3a, p2a, nrow = 1)

ggsave(here("figures", "supplemental", "trait_correlations_imputed.png"))
ggsave(here("figures", "supplemental", "trait_correlations_imputed.pdf"))

##---------------------------------------------------------------------------
## 04. check alignment with knighton et al dataset and merge into imputed data
##---------------------------------------------------------------------------

knighton_traits <- read.csv(here("data", "knighton_traits_database", "knighton_traits.csv"))[,1:13]
knighton_traits$spec.name <- gsub(" ", "_", knighton_traits$spec.name)

colnames(knighton_traits)[1] <- "species"

merged_traits <- merge(fd, knighton_traits, by = "species", all.x = TRUE)

ggplot(merged_traits, aes(x = p50, y = P50, color = clade)) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", size = 1, color = "black") +
  geom_smooth(method = "lm", linewidth = 2, fullrange = TRUE, color = "gray") +
  geom_point(size = 2.5) +
  xlab("Imputed P50 (MPa)") +
  ylab("Knighton et al. P50 (MPa)") +
  scale_color_brewer(type = "qual", palette = 6) +
  scale_x_continuous(expand = c(0,0), limits = c(-15,0)) +
  scale_y_continuous(expand = c(0,0), limits = c(-15,0)) +
  #geom_smooth(method = "lm", linewidth = 2, fullrange = TRUE) +
  theme_bw()

ggsave(here("figures", "supplemental", "p50_imputed_vs_knighton.png"))
ggsave(here("figures", "supplemental", "p50_imputed_vs_knighton.pdf"))

summary(lm(merged_traits$P50 ~ merged_traits$p50))

merged_traits <- merged_traits[,c("species", "fia_id", "clade", "p50", "p50_database", "Amax",
                  "Amax_database", "ssd", "ssd_database", "P50", "P12", "P88", "gsmax", "rdmax", "WUE", "height", "SLA", "LeafN")]

write.csv(merged_traits, here("data", "traits_imputed_merged.csv"), row.names = FALSE)
