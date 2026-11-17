
library(data.table)
library(sf)
library(metafor)

make_circles <- function(points, radius = 100000, min_n = 1, seed = 1) {
  
  bbox <- st_as_sfc(st_bbox(points))

  bbox_proj <- st_transform(bbox, crs = st_crs(5070))
  points <- st_transform(points, crs = st_crs(5070))

  bbox_area <- as.numeric(st_area(bbox_proj))

  N <- floor(bbox_area / (pi * radius^2))

  if (N < min_n) {
    N <- min_n
  }

  set.seed(seed)
  centers <- st_sample(bbox_proj, size = N, type = "random")
  centers_sf <- st_sf(id = seq_len(length(centers)), geometry = centers)

  circles <- st_buffer(centers_sf, dist = radius)
  circles <- st_transform(circles, crs = st_crs(4326))

  return(circles)
  
}

read_tree <- function(state) {
  f <- here::here("data", "fia", state, paste0(state, "_TREE.csv"))
  if (!file.exists(f)) return(NULL)
  dt <- fread(f)
  colnames(dt)[1] <- "TRE_CN"
  dt <- dt[, .(TRE_CN, PLT_CN, PREV_TRE_CN, STATUSCD, TPA_UNADJ, SPCD)]
  dt
}

process_circle <- function(i, intersects_list, tree_data, plots_dt, tradeoffs_dt, trait_vars) {

  pidx <- intersects_list[[i]]
  if (length(pidx) < 3) return(list(i = i, beta_overall = NA_real_, se_overall = NA_real_,
                                    beta_community = NA_real_, se_community = NA_real_))
  p <- plots_dt[PLT_IDX %in% pidx]
  plt_cns <- unique(p$PLT_CN)

  remeas <- tree_data[PLT_CN %in% plt_cns & !is.na(TPA_UNADJ) & STATUSCD != 0]
  trees_prev <- tree_data[PLT_CN %in% unique(p$PREV_PLT_CN)]

  remeas[, mort := STATUSCD == 2 & !is.na(TPA_UNADJ) &
               !(PREV_TRE_CN %in% trees_prev[STATUSCD == 2, TRE_CN])]
  remeas[, live := STATUSCD == 1 & !is.na(TPA_UNADJ)]
  trees_prev[, live := TRE_CN %in% remeas[live == TRUE, PREV_TRE_CN]]
  trees_prev[, live_all := STATUSCD == 1 & !is.na(TPA_UNADJ) & TRE_CN %in% remeas$PREV_TRE_CN]

  traits_sum <- aggregate_traits(trees_prev, trait_vars)
  if (sum(!is.na(traits_sum$species) & traits_sum$tpa > 0) < 3) {
    return(list(i = i, beta_overall = NA_real_, se_overall = NA_real_,
                beta_community = NA_real_, se_community = NA_real_))
  }

  tr <- fit_tradeoff_models(traits_sum)$out
  comm <- tradeoffs_dt[PLOT_ID %in% unique(p$PLOT_ID)]
  cb <- comm[["beta_pca"]]; cs <- comm[["se_pca"]]
  if (sum(!is.na(cb)) < 3) {
    return(list(i = i, beta_overall = tr$beta_pca, se_overall = tr$se_pca,
                beta_community = NA_real_, se_community = NA_real_))
  }

  comm_sub <- comm[!is.na(comm$se_pca),]
  v <- comm_sub$se_pca^2
  a <- comm_sub$tpa_init_full / sum(comm_sub$tpa_init_full, na.rm = TRUE)
  b <- (1 / v) / sum(1 / v, na.rm = TRUE)

  cb <- comm_sub$beta_pca; vi <- comm_sub$se_pca^2
  w <- 1/vi
  beta_hat <- sum(w * cb) / sum(w)
  Q <- sum(w * (cb - beta_hat)^2)
  dfm <- length(cb) - 1
  C <- sum(w) - sum(w^2)/sum(w)
  tau2 <- max(0, (Q - dfm) / C)
  v <- vi + tau2
  b <- (1/v); b <- b / sum(b)

  var_over <- sum(a^2 * v, na.rm = TRUE)
  var_comm <- sum(b^2 * v, na.rm = TRUE)
  cov_ac   <- sum(a * b * v, na.rm = TRUE)

  diff <- tr$beta_pca - beta_hat
  se_diff <- sqrt(pmax(0, var_over + var_comm - 2 * cov_ac))

  list(i = i, beta_overall = tr$beta_pca, se_overall = tr$se_pca,
       beta_community = beta_hat, se_community = sqrt(var_comm),
       diff = diff, se_diff = se_diff )
}

## pool_boot <- function(diff, se, B = 2000) {
##   w <- 1 / (se^2)
##   boot <- numeric(B)
##   for (b in 1:B) {
##     sim <- rnorm(length(diff), mean = diff, sd = se)
##     boot[b] <- sum(w * sim, na.rm = TRUE) / sum(w[!is.na(sim)])
##   }
##   boot
## }

check_scale_diff <- function(plots, radius, traits, tradeoffs, trait_vars, states, min_n = 100, max_n = 1000, seed = 1) {
  
  circles <- make_circles(plots, radius, min_n = min_n, seed = seed)

  plots_dt <- as.data.table(st_drop_geometry(plots))
  plots_dt[, PLT_IDX := .I]
  intersects_list <- st_intersects(circles, plots, sparse = TRUE)
  plot_counts <- lengths(intersects_list)
  circles <- circles[plot_counts >= 3,]
  intersects_list <- intersects_list[plot_counts >= 3]
  plot_counts <- plot_counts[plot_counts >= 3]

  if (nrow(circles) > max_n) {
    tokeep <- sample(1:nrow(circles), max_n)
    circles <- circles[tokeep, ]
    intersects_list <- intersects_list[tokeep]
    plot_counts <- plot_counts[tokeep]
  }

  needed_plot_rows <- unique(unlist(intersects_list[plot_counts >= 3]))
  states_needed_code <- unique(plots_dt[PLT_IDX %in% needed_plot_rows, STATECD])
  states_needed <- states[states$STATECD %in% states_needed_code, "state"]

  tree_list <- lapply(states_needed, read_tree)
  tree_data_all <- rbindlist(tree_list[!sapply(tree_list, is.null)], use.names = TRUE)
  setkey(tree_data_all, PLT_CN)

  traits_dt <- as.data.table(traits); setkey(traits_dt, SPCD)
  tree_data_all <- traits_dt[tree_data_all, on = "SPCD"]
  tradeoffs_dt <- as.data.table(tradeoffs); setkey(tradeoffs_dt, PLOT_ID)

  res <- list()
  for (i in seq_len(nrow(circles))) {
    res[[i]] <- process_circle(i, intersects_list, tree_data_all, plots_dt, tradeoffs_dt, trait_vars)
    print(paste0("Completed circle ", i, " of ", nrow(circles)))
  }

  res_dt <- rbindlist(lapply(res, as.data.table), fill = TRUE)

  res_dt <- res_dt[!is.na(beta_community) & beta_community != 0]

  ## res_dt$diff <- res_dt$beta_overall - res_dt$beta_community
  ## res_dt$se_diff <- sqrt(res_dt$se_overall^2 + res_dt$se_community^2)

  ## boot <- pool_boot(res_dt$diff, res_dt$se_diff, B = 2000)

  diff_mod <- rma(yi = diff, sei = se_diff, data = res_dt, method = "REML")
  diff_lm <- lm(diff ~ 1, weights = 1 / se_diff^2, data = res_dt)

  return(list(alpha = diff_mod$b,
              se = diff_mod$se,
              ci = c(diff_mod$ci.lb,
                     diff_mod$ci.ub),
              prop_pos = mean(res_dt$diff > 0, na.rm = TRUE),
              alpha_diff_lm = summary(diff_lm)$coefficients[1,"Estimate"],
              se_diff_lm = summary(diff_lm)$coefficients[1,"Std. Error"]))

}
