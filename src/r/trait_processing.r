##--------------------------------------------------------------
## Functions for processing tree trait data                                                       
##--------------------------------------------------------------

library(data.table)
library(FactoMineR)
library(Hmisc)
library(progress)
library(MASS)
library(here)

## pull initial census data for a plot
pull_tree_data <- function(ps, tree_data, init = TRUE) {
  
  remeas <- tree_data[PLT_CN == ps$PLT_CN & !is.na(PREV_TRE_CN) &
                      !is.na(TPA_UNADJ) & STATUSCD != 0]

  if (init) {
    trees <- tree_data[PLT_CN == ps$PREV_PLT_CN]
  } else {
    trees <- tree_data[PLT_CN == ps$PLT_CN]
  }
  
  trees$dead <- trees$STATUSCD == 2
  
  remeas$mort <- (remeas$STATUSCD == 2 & !is.na(remeas$TPA_UNADJ) &
                  !(remeas$PREV_TRE_CN %in% trees[trees$dead, TRE_CN]))
  remeas$live <- remeas$STATUSCD == 1 & !is.na(remeas$TPA_UNADJ)

  trees$live <- trees$TRE_CN %in% remeas[remeas$live, PREV_TRE_CN]
  trees$live_all <- trees$STATUSCD == 1 & !is.na(trees$TPA_UNADJ) & trees$TRE_CN %in% remeas$PREV_TRE_CN

  return(trees)
  
}

## compute distance from a point to a line defined by slope and intercept
dist_point_line <- function(point, slope, intercept) {
  abs(slope * point[1] - point[2] + intercept) / sqrt(slope^2 + 1)
}

## compute average weighted distance from points to a line
avg_dist <- function(traits_data, xvar, yvar, slope, intercept) {
  dists <- mapply(function(x, y) dist_point_line(c(x, y), slope, intercept),
                  traits_data[[xvar]], traits_data[[yvar]])
  weighted.mean(dists, traits_data$tpa)
}

se_avg_dist <- function(traits_data, beta, intercept, v, xvar, yvar, n = 50) {
  td_list <- numeric(n)
  for (i in 1:n) {
    sample <- MASS::mvrnorm(1, mu = c(intercept = intercept, slope = beta), Sigma = v[1:2, 1:2])
    td_list[i] <- avg_dist(traits_data, xvar, yvar, sample["slope"], sample["intercept"])
  }
  return(sd(td_list))
}

signed_dist <- function(x, y, slope, intercept) {
   (y - (intercept + slope * x)) / sqrt(1 + slope^2)
}

one_sided_r2 <- function(traits_data, xvar, yvar, slope, intercept) {

  x <- traits_data[[xvar]]
  y <- traits_data[[yvar]]
  w <- traits_data$tpa
  p <- signed_dist(x, y, slope, intercept)

  ## select only the 'below' points (negative p)
  below_idx <- which(!is.na(p) & p < 0 & !is.na(w) & w > 0)
  if (length(below_idx) < 2) return(NA_real_)

  p_below <- p[below_idx]
  x_below <- x[below_idx]
  y_below <- y[below_idx]
  w_below <- w[below_idx]

  r_var_below <- Hmisc::wtd.var(p_below, weights = w_below, normwt = FALSE, na.rm = TRUE)

  ## need to calculate variance of data projected onto line perpendicular to tradeoff
  theta <- atan(slope)
  perp_coord <- -x_below * sin(theta) + y_below * cos(theta)
  y_proj_var <- Hmisc::wtd.var(perp_coord, weights = w_below, normwt = FALSE, na.rm = TRUE)
  
  if (is.na(r_var_below) || is.na(y_proj_var) || y_proj_var == 0) return(NA_real_)
  r2_below <- 1 - (r_var_below / y_proj_var)
  return(r2_below)
}

## aggregate traits by specified grouping variables
aggregate_traits <- function(tree_data, trait_vars) {
  traits_sum <- aggregate(tree_data[, TPA_UNADJ],
                          by = lapply(trait_vars,
                                      function(var) as.data.frame(tree_data)[, var]),
                          FUN = sum,
                          na.rm = TRUE)
  colnames(traits_sum) <- c(trait_vars, "tpa")
  return(traits_sum)
}

## standardize specified traits to standard units (z-scores)
standardize_traits <- function(traits_data, traits_to_standardize) {
  for (trait in traits_to_standardize) {
    traits_data[[paste0(trait, "_scaled")]] <- scale(traits_data[[trait]], center = TRUE, scale = TRUE)
  }
  return(traits_data)
}

## Define a generic function for standardizing and projecting data onto PCA space
project_pca <- function(data, pca_object, pca_vars, pc_names = c("pc1", "pc2")) {
  
  if ("data.table" %in% class(data)) {
    data_subset <- data[, pca_vars, with = FALSE]
  } else {
    data_subset <- data[, pca_vars, drop = FALSE]
  }

  data_scaled <- scale(data_subset,
                       center = pca_object$call$centre,
                       scale = pca_object$call$ecart.type)

    pc_scores <- as.data.frame(as.matrix(data_scaled) %*% pca_object$var$coord)

  if ("data.table" %in% class(data)) {
    data[, (pc_names[1]) := pc_scores[, 1]]
    data[, (pc_names[2]) := pc_scores[, 2]]
    return(data)
  } else {
    data[[pc_names[1]]] <- pc_scores[, 1]
    data[[pc_names[2]]] <- pc_scores[, 2]
    return(data)
  }
}

## compute community-weighted mean
cwm <- function(trait, tpa) {
  sum(trait * tpa, na.rm = TRUE) / sum(tpa, na.rm = TRUE)
}

init_template <- function(output_cols) {
  setNames(data.frame(matrix(NA, nrow = 0, ncol = length(output_cols))), output_cols)
}

fit_tradeoff_models <- function(traits_sum) {
  out <- list()
  tryCatch({
    lm_base <- lm(Amax ~ p50, weights = tpa, data = traits_sum)
    s_base  <- summary(lm_base)
    out$intercept_base <- coef(s_base)[1, "Estimate"]
    out$intercept_se_base <- coef(s_base)[1, "Std. Error"]
    out$beta_base      <- coef(s_base)[2, "Estimate"]
    out$se_base        <- coef(s_base)[2, "Std. Error"]
    out$r2_base        <- s_base$r.squared

    lm_pca <- lm(growth_pc1 ~ drought_pc2, weights = tpa, data = traits_sum)
    s_pca  <- summary(lm_pca)
    out$intercept_pca <- coef(s_pca)[1, "Estimate"]
    out$intercept_se_pca <- coef(s_pca)[1, "Std. Error"]
    out$beta_pca      <- coef(s_pca)[2, "Estimate"]
    out$se_pca        <- coef(s_pca)[2, "Std. Error"]
    out$r2_pca        <- s_pca$r.squared

    lm_pca_p50 <- lm(growth_pc1 ~ p50, weights = tpa, data = traits_sum)
    s_pca_p50  <- summary(lm_pca_p50)
    out$intercept_pca_p50 <- coef(s_pca_p50)[1, "Estimate"]
    out$intercept_se_pca_p50 <- coef(s_pca_p50)[1, "Std. Error"]
    out$beta_pca_p50      <- coef(s_pca_p50)[2, "Estimate"]
    out$se_pca_p50        <- coef(s_pca_p50)[2, "Std. Error"]
    out$r2_pca_p50        <- s_pca_p50$r.squared

    ret <- list(out = out, v_base = vcov(lm_base), v_pca = vcov(lm_pca), v_pca_p50 = vcov(lm_pca_p50))
    
  }, error = function(e) ret <<- NULL)
  return(ret)
}

## computes distances, CWMs, and quantiles
compute_functional_metrics <- function(traits_sum, coefs, v_base, v_pca, v_pca_p50) {
  td <- c(
    avg_dist(traits_sum, "p50", "Amax", coefs$beta_base, coefs$intercept_base),
    avg_dist(traits_sum, "drought_pc2", "growth_pc1", coefs$beta_pca, coefs$intercept_pca),
    avg_dist(traits_sum, "p50", "growth_pc1", coefs$beta_pca_p50, coefs$intercept_pca_p50)
  )

  td_se <- c(
    se_avg_dist(traits_sum, coefs$beta_base, coefs$intercept_base, v_base,
                "p50", "Amax", n = 100),
    se_avg_dist(traits_sum, coefs$beta_pca, coefs$intercept_pca, v_pca,
                "drought_pc2", "growth_pc1", n = 100),
    se_avg_dist(traits_sum, coefs$beta_pca_p50, coefs$intercept_pca_p50, v_pca_p50,
                "p50", "growth_pc1", n = 100)
  )

  os_r2 <- c(
    one_sided_r2(traits_sum, "p50", "Amax", coefs$beta_base, coefs$intercept_base),
    one_sided_r2(traits_sum, "drought_pc2", "growth_pc1", coefs$beta_pca, coefs$intercept_pca),
    one_sided_r2(traits_sum, "p50", "growth_pc1", coefs$beta_pca_p50, coefs$intercept_pca_p50)
  )

  cwm_vals <- c(
    cwm(traits_sum$Amax, traits_sum$tpa),
    cwm(traits_sum$p50, traits_sum$tpa),
    cwm(traits_sum$growth_pc1, traits_sum$tpa),
    cwm(traits_sum$drought_pc2, traits_sum$tpa)
  )

  min_vals <- unname(c(
    wtd.quantile(traits_sum$Amax, weights = traits_sum$tpa, probs = 0.05),
    wtd.quantile(traits_sum$p50, weights = traits_sum$tpa, probs = 0.05),
    wtd.quantile(traits_sum$growth_pc1, weights = traits_sum$tpa, probs = 0.05),
    wtd.quantile(traits_sum$drought_pc2, weights = traits_sum$tpa, probs = 0.05)
  ))

  max_vals <- unname(c(
    wtd.quantile(traits_sum$Amax, weights = traits_sum$tpa, probs = 0.95),
    wtd.quantile(traits_sum$p50, weights = traits_sum$tpa, probs = 0.95),
    wtd.quantile(traits_sum$growth_pc1, weights = traits_sum$tpa, probs = 0.95),
    wtd.quantile(traits_sum$drought_pc2, weights = traits_sum$tpa, probs = 0.95)
  ))
  return({
    list(td = td, td_se = td_se, os_r2 = os_r2, cwm = cwm_vals, min = min_vals, max = max_vals)
  })
}


## processes one plot row
process_plot <- function(row, tree_data, traits, trait_vars, init = TRUE) {
  tl <- pull_tree_data(row, tree_data, init)
  traits_trees <- merge(tl, traits, by = "SPCD", all.x = TRUE)
  traits_sum <- aggregate_traits(traits_trees, trait_vars)
  traits_sum <- standardize_traits(traits_sum, c("p50", "Amax", "growth_pc1", "drought_pc2"))
  richness <- sum(!is.na(traits_sum$species) & traits_sum$tpa > 0)

  if (richness <= 2) return(NULL)

  m_res <- fit_tradeoff_models(traits_sum)
  coefs <- m_res$out
  
  if (is.null(coefs)) return(NULL)

  metrics <- compute_functional_metrics(traits_sum, coefs, m_res$v_base, m_res$v_pca, m_res$v_pca_p50)

  c(
    richness          = richness,
    intercept_base    = coefs$intercept_base,
    beta_base         = coefs$beta_base,
    se_base           = coefs$se_base,
    intercept_pca     = coefs$intercept_pca,
    beta_pca          = coefs$beta_pca,
    se_pca            = coefs$se_pca,
    intercept_pca_p50 = coefs$intercept_pca_p50,
    beta_pca_p50      = coefs$beta_pca_p50,
    se_pca_p50        = coefs$se_pca_p50,
    td_base           = metrics$td[1],
    td_pca            = metrics$td[2],
    td_pca_p50        = metrics$td[3],
    td_se_base        = metrics$td_se[1],
    td_se_pca         = metrics$td_se[2],
    td_se_pca_p50     = metrics$td_se[3],
    r2_base           = coefs$r2_base,
    r2_pca            = coefs$r2_pca,
    r2_pca_p50        = coefs$r2_pca_p50,
    os_r2_base        = metrics$os_r2[1],
    os_r2_pca         = metrics$os_r2[2],
    os_r2_pca_p50     = metrics$os_r2[3],
    cwm_Amax          = metrics$cwm[1],
    cwm_p50           = metrics$cwm[2],
    cwm_growth_pc1    = metrics$cwm[3],
    cwm_drought_pc2   = metrics$cwm[4],
    min_Amax          = metrics$min[1],
    min_p50           = metrics$min[2],
    min_growth_pc1    = metrics$min[3],
    min_drought_pc2   = metrics$min[4],
    max_Amax          = metrics$max[1],
    max_p50           = metrics$max[2],
    max_growth_pc1    = metrics$max[3],
    max_drought_pc2   = metrics$max[4]
  )
}

## processes one state
process_state <- function(state_row, gm_data, traits, trait_vars, output_cols, init = TRUE, filename = "tradeoffs") {
  gm_data_sub <- gm_data[STATECD == state_row$STATECD]
  tree_data <- fread(here("data", "fia", state_row$state, paste0(state_row$state, "_TREE.csv")))
  colnames(tree_data)[1] <- "TRE_CN"
  tree_data <- tree_data[, .(TRE_CN, PLT_CN, PREV_TRE_CN, STATUSCD, TPA_UNADJ, SPCD)]

  N <- nrow(gm_data_sub)
  pb <- progress_bar$new(total = N)
  res <- setNames(data.frame(matrix(NA, nrow = N, ncol = length(output_cols))), output_cols)

  for (i in 1:N) {
    vals <- process_plot(gm_data_sub[i,], tree_data, traits, trait_vars, init)
    if (!is.null(vals)) res[i, names(vals)] <- vals
    pb$tick()
  }

  res <- res[, output_cols]
  out <- cbind(gm_data_sub, res)
  
  fwrite(out, here("data", "tradeoffs", paste0(filename, "_", state_row$state, ".csv")))

  return(out)
}

## master loop
estimate_tradeoffs <- function(states, gm_data, traits, trait_vars, output_cols,
                               init = TRUE, filename = "tradeoffs") {
  
  template <- init_template(output_cols)
  out <- cbind(gm_data[0,], template)

  for (s in seq_len(nrow(states))) {
    message("Processing: ", states[s, "state"], " [", s, "/", nrow(states), "]")
    state_out <- process_state(states[s,], gm_data, traits, trait_vars, output_cols, init, filename)
    out <- rbind(out, state_out)
    fwrite(out, here("data", "tradeoffs", paste0(filename, "_full.csv")))
  }

  return(out)

}
