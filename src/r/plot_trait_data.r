##-----------------------------------------------------------------
## Function to visualize trait data coverage on a phylogenetic tree                                                       
##-----------------------------------------------------------------

library(ggplot2)
library(ggtree)


## function to visualize trait coverage
plot_trait_coverage <- function(data, traitdata, trait = "p50") {

  spp_list <- traitdata[!is.na(traitdata[,paste0(trait, "_mean")]), "species"]

  pt_tib <- as_tibble(data)
  pt_tib$group <- NA
  for (i in 1:nrow(pt_tib)) {
    if (pt_tib$label[i] %in% spp_list) {
      pt_tib$group[i] <- "yes"
    } else {
      pt_tib$group[i] <- "no"
    }
  }

  pruned_tree <- as.treedata(pt_tib)

  return(ggtree(pruned_tree, branch.length='none', layout='circular') +
         ggtitle(paste0("Coverage for ", trait)) +
         geom_tippoint(aes(color = group), alpha = 0.8) +
         scale_color_manual(values = c("gray", "#c51b8a")) +
         theme(legend.position = "none"))
}
