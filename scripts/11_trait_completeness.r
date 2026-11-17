library(here)
library(data.table)
library(ggplot2)
library(patchwork)

source(here("src", "r", "trait_processing.r"))

traits <- read.csv("data/full_traits.csv")

traits_imputed <- read.csv(here("data", "traits_imputed_merged.csv"))
colnames(traits_imputed)[2] <- "SPCD"

traits$species <- gsub(" ", "_", traits$species)

traits <- merge(traits, traits_imputed[, c("species", "SPCD")], by = "species", all.x = TRUE)

aggregate_traits <- function(tree_data, trait_vars) {

  trait_vars <- unique(c("SPCD", trait_vars))
  trait_vars <- intersect(trait_vars, names(tree_data))

  tree_data[
    ,
    .(tpa = sum(TPA_UNADJ, na.rm = TRUE)),
    by = trait_vars
  ]
}


sum(is.na(traits$p50_mean)) / nrow(traits)
sum(is.na(traits$Amax_mean)) / nrow(traits)

trait_vars <- c("p50_mean", "Amax_mean")



gm_data <- fread("data/gm/gm_full.csv")

plots <- unique(gm_data$PLOT_ID)

states <- read.csv(here("data", "state_list.csv"))



total_tpa <- 0
total_missing_p50 <- 0
total_missing_Amax <- 0

for (s in unique(states$state)) {
  
  state_code <- states[states$state == s, "STATECD"]
  gm_data_sub <- gm_data[STATECD == state_code,]
  
  tree_data <- fread(here("data", "fia", s, paste0(s, "_TREE.csv")))
  colnames(tree_data)[1] <- "TRE_CN"
  tree_data <- tree_data[, .(TRE_CN, PLT_CN, PREV_TRE_CN, STATUSCD, TPA_UNADJ, SPCD)]

  missing <- data.frame(
    plot = gm_data_sub$PLOT_ID,
    sample_number = gm_data_sub$SAMPLE_NUMBER,
    p50 = NA_real_,
    Amax = NA_real_,
    tpa = NA_real_
  )
  
  for (i in 1:nrow(gm_data_sub)) {
    
    td <- pull_tree_data(gm_data_sub[i, ], tree_data, TRUE)
    td <- merge(td, traits, by = "SPCD", all.x = TRUE)
    agg <- aggregate_traits(td, trait_vars)

    missing_tpa_p50 <- sum(agg[is.na(agg$p50_mean), "tpa"])
    missing_tpa_Amax <- sum(agg[is.na(agg$Amax_mean), "tpa"])

    total_missing_p50 <- total_missing_p50 + missing_tpa_p50
    total_missing_Amax <- total_missing_Amax + missing_tpa_Amax
    total_tpa <- total_tpa + sum(agg$tpa)

    missing[i, "p50"] <- missing_tpa_p50 / sum(agg$tpa)
    missing[i, "Amax"] <- missing_tpa_Amax / sum(agg$tpa)
    missing[i, "tpa"] <- sum(agg$tpa)
    
  }

  write.csv(missing, here("data", "trait_completeness", paste0(s, "_trait_completeness.csv")), row.names = FALSE)

  print(paste0("State: ", s, " - Missing p50: ", round(mean(missing$p50) * 100, 2), "% - Missing Amax: ", round(mean(missing$Amax) * 100, 2), "%"))
  
}

total_missing_p50 / total_tpa ## 31.313%
total_missing_Amax / total_tpa ## 20.063%

modeling_data <- read.csv(here("data", "data_for_modeling.csv"))

p50_missing_all <- numeric()
Amax_missing_all <- numeric()

state_missing_summary <- data.frame(
  state = character(),
  prop_missing_p50 = numeric(),
  prop_missing_Amax = numeric(),
  prop_over_30_p50 = numeric(),
  prop_over_30_Amax = numeric(),
  prop_over_50_p50 = numeric(),
  prop_over_50_Amax = numeric()
)

for (s in unique(states$state)) {
  missing <- read.csv(here("data", "trait_completeness", paste0(s, "_trait_completeness.csv")))
  p50_missing_all <- c(p50_missing_all, missing$p50)
  Amax_missing_all <- c(Amax_missing_all, missing$Amax)

  state_missing_summary <- rbind(state_missing_summary, data.frame(
    state = s,
    prop_missing_p50 = mean(missing$p50),
    prop_missing_Amax = mean(missing$Amax),
    prop_over_30_p50 = mean(missing$p50 > 0.3),
    prop_over_30_Amax = mean(missing$Amax > 0.3),
    prop_over_50_p50 = mean(missing$p50 > 0.5),
    prop_over_50_Amax = mean(missing$Amax > 0.5)
  ))
  
}

breaks <- seq(0, 1, by = 0.05)

max_count <- max(
  hist(p50_missing_all, breaks = breaks, plot = FALSE)$counts,
  hist(Amax_missing_all, breaks = breaks, plot = FALSE)$counts
)

h1 <- ggplot(aes(x = p50_missing_all), data = data.frame(p50_missing_all = p50_missing_all)) +
  geom_histogram(breaks = breaks, fill = "#0057D9", color = NA) +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 1), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, max_count), expand = c(0, 0)) +
  labs(title = expression("Proportion of Missing Trait Values"), x = expression("Proportion Missing " * P[50]), y = "Count") +
  theme_bw()

h2 <- ggplot(aes(x = Amax_missing_all), data = data.frame(Amax_missing_all = Amax_missing_all)) +
  geom_histogram(breaks = breaks, fill = "#009B3A", color = NA) +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 1), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, max_count), expand = c(0, 0)) +
  labs(x = expression("Proportion Missing " * A[max]), y = "Count") +
  theme_bw() +
  theme(axis.title.y = element_blank(), axis.text.y = element_blank())

h1 + h2 + plot_layout(ncol = 2)

ggsave(here("figures", "trait_completeness_histograms.png"), width = 10, height = 5)

write.csv(state_missing_summary, here("data", "trait_completeness", "state_missing_summary.csv"), row.names = FALSE)


plot_flags <- data.frame(
  state = character(0),
  plot = character(0),
  sample_number = character(0),
  flag_p50 = logical(0),
  flag_Amax = logical(0)
)

for (s in unique(states$state)) {
  missing <- read.csv(here("data", "trait_completeness", paste0(s, "_trait_completeness.csv")))
  missing$state <- s
  missing$flag_p50 <- missing$p50 > 0.30
  missing$flag_Amax <- missing$Amax > 0.30
  plot_flags <- rbind(plot_flags, missing[, c("state", "plot", "sample_number", "flag_p50", "flag_Amax")])
}

plot_flags$flag_both <- plot_flags$flag_p50 & plot_flags$flag_Amax

sum(plot_flags$flag_both) / nrow(plot_flags) ## 17.8% of plots have >30% missing for either trait

write.csv(plot_flags, here("data", "trait_completeness", "plot_flags.csv"), row.names = FALSE)


plot(traits_imputed$p50, traits_imputed$rdmax)

cor.test(traits_imputed$p50, traits_imputed$rdmax, use = "complete.obs", method = 'pearson')

plot(log(abs(traits_imputed$p50)), log(traits_imputed$rdmax))

p1 <- ggplot(traits_imputed, aes(x = p50, y = rdmax)) +
  geom_point(size = 3) +
  xlab(expression(P[50])) +
  ylab(expression("Max. Rooting Depth")) +
  theme_bw()

p2 <- ggplot(traits_imputed, aes(x = log(abs(p50)), y = log(rdmax))) +
  geom_point(size = 3) +
  xlab(expression("Log |" * P[50] * "|")) +
  ylab(expression("Log Max. Rooting Depth")) +
  theme_bw()

p1 + p2 + plot_layout(ncol = 2)

ggsave(here("figures", "p50_vs_rdmax.png"), width = 10, height = 5)
