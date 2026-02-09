##--------------------------------------------------------------
## Functions to plot predictions from GAM models
##
## author: jacob levine; jacob.levine@utah.edu
##--------------------------------------------------------------

library(ggplot2)
library(MetBrewer)
library(patchwork)

plot_one_gam_prediction <- function(model, data, variable, interaction = NULL,
                                    n_points = 100, n_interaction_levels = 3,
                                    ci_level = 0.95) {

  variable_name_dict <- list(
    "beta_pca" = "Tradeoff slope",
    "se_pca" = "Tradeoff uncertainty",
    "drought_strength" = "Mean drought strength",
    "prop_drought" ="Proportion of drought years",
    "ba" = "Basal area (m²/ha)",
    "stand_age" = "Stand age (years)",
    "map" = "Mean annual precipitation (mm)",
    "mat" = "Mean annual temperature (°C)",
    "range_drought_pc2" = "Range in drought traits (PC2)",
    "range_growth_pc1" = "Range in growth traits (PC1)",
    "cwm_drought_pc2" = "Community-weighted mean drought traits (PC2)",
    "cwm_growth_pc1" = "Community-weighted mean growth traits (PC1)",
    "elev" = "Elevation (m)")
  
  model_vars <- all.vars(formula(model))
  model_vars <- model_vars[2:length(model_vars)]

  pred_data <- data.frame(matrix(ncol = length(model_vars), nrow = n_points * max(1, n_interaction_levels)))
  colnames(pred_data) <- model_vars

  vars_to_fix <- setdiff(model_vars, c(variable, interaction))
  for (var in vars_to_fix) {
    if (!is.null(data[[var]])) {
      pred_data[[var]] <- median(data[[var]], na.rm = TRUE)
    }
  }

  var_s <- paste0(variable, "_scaled")
  int_s <- if (!is.null(interaction)) paste0(interaction, "_scaled") else NULL

  # Sequence for main variable
  pred_data[[var_s]] <- rep(seq(quantile(data[[var_s]], 0.25, na.rm = TRUE),
                                quantile(data[[var_s]], 0.95, na.rm = TRUE),
                                length.out = n_points),
                            times = ifelse(is.null(interaction), 1, n_interaction_levels))

  if (!is.null(interaction)) {
    interaction_levels <- seq(quantile(data[[int_s]], 0.25, na.rm = TRUE),
                              quantile(data[[int_s]], 0.95, na.rm = TRUE),
                              length.out = n_interaction_levels)
    pred_data[[int_s]] <- rep(interaction_levels, each = n_points)
  }

  exclude_terms <- if (!is.null(interaction)) {
                     setdiff(names(model$smooth), c(variable, interaction))
                   } else {
                     setdiff(names(model$smooth), variable)
                   }

  pred <- predict(model, newdata = pred_data, type = "lpmatrix")
  b <- coef(model)
  v <- vcov(model)

  pred_no_int <- pred
  pred_no_int[, c("(Intercept)")] <- 0

  var_param_only <- rowSums((pred_no_int %*% v) * pred_no_int)

  pred_data$se.fit <- sqrt(pmax(0, var_param_only))
  pred_data$fit <- as.numeric(pred %*% b)

  crit_val <- qnorm(1 - (1 - ci_level) / 2)
  pred_data$lower <- pred_data$fit - crit_val * pred_data$se.fit
  pred_data$upper <- pred_data$fit + crit_val * pred_data$se.fit
  
  ## inverse transform to prediction & CI
  pred_data$fit <- (pred_data$fit * sd(data$growth_pct, na.rm = TRUE) + mean(data$growth_pct, na.rm = TRUE)) * 100
  pred_data$lower <- (pred_data$lower * sd(data$growth_pct, na.rm = TRUE) + mean(data$growth_pct, na.rm = TRUE)) * 100
  pred_data$upper <- (pred_data$upper * sd(data$growth_pct, na.rm = TRUE) + mean(data$growth_pct, na.rm = TRUE)) * 100

  ## Back-transform scaled variable and interaction to original scale for plotting
  pred_data[[variable]] <- pred_data[[var_s]] * sd(data[[variable]], na.rm = TRUE) + mean(data[[variable]], na.rm = TRUE)
  if (!is.null(interaction)) {
    pred_data[[interaction]] <- pred_data[[int_s]] * sd(data[[interaction]], na.rm = TRUE) + mean(data[[interaction]], na.rm = TRUE)
  }

  if (is.null(interaction)) {
    p <- ggplot(pred_data, aes(x = !!sym(variable), y = fit)) +
      geom_line(linewidth = 2) +
      geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2)
  } else {
    p <- ggplot(pred_data, aes(x = !!sym(variable), y = fit, color = !!sym(interaction), group = !!sym(interaction))) +
      geom_line(linewidth = 2) +
      geom_ribbon(aes(ymin = lower, ymax = upper, fill = !!sym(interaction), group = !!sym(interaction)), alpha = 0.2)
  }

  var_lab  <- variable_name_dict[[variable]]
  int_lab  <- if (!is.null(interaction)) variable_name_dict[[interaction]] else NULL

  p <- p +
    labs(
      x = var_lab,
      y = "Basal area growth (%)"
    ) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    theme_bw()

  if (!is.null(interaction)) {
    p <- p +
      scale_color_gradientn(colors = met.brewer("Archambault", n = 5, type = "continuous"), name = int_lab) +
      scale_fill_gradientn(colors = met.brewer("Archambault", n = 5, type = "continuous"), name = int_lab)
  }

  return(p)
}


plot_multiple_gam_predictions <- function(model, data,
                                          interaction_vars = c("beta_pca", "se_pca", "drought_strength"),
                                          vars_to_plot = NULL,
                                          ci_level = 0.9,
                                          ncol = 4,
                                          align_axes = FALSE,
                                          axis_limits = c(NA, NA)) {

  ## get a list of possible variables to plot from model formula
  model_vars <- all.vars(formula(model))[-1]
  model_vars <- model_vars[!model_vars %in% c("year", "lon", "lat")]
  model_vars <- gsub("_scaled$", "", model_vars)

  ## if a subset of variables is supplied, keep only those
  if (!is.null(vars_to_plot)) {
    model_vars <- intersect(model_vars, vars_to_plot)
  }
  
  ## combine variable list with specified interactions
  vars <- lapply(model_vars, function(v) {
    list(variable = v, interaction = if (v %in% interaction_vars) "ba" else NULL)
  })

  ## we want the interaction panels first if vars_to_plot is not specified
  if (is.null(vars_to_plot)) {
    vars <- c(
      vars[sapply(vars, function(x) !is.null(x$interaction))],
      vars[sapply(vars, function(x) is.null(x$interaction))]
    )
  }
  
  ## we are going to steal the legend element from the first plot
  first_with_interaction <- which(sapply(vars, function(x) !is.null(x$interaction)))[1]
  
  ## make the plots
  plots <- lapply(seq_along(vars), function(i) {
    p <- plot_one_gam_prediction(
      model = model,
      data = data,
      variable = vars[[i]]$variable,
      interaction = vars[[i]]$interaction,
      ci_level = ci_level
    )

    if (!is.na(first_with_interaction) & i != first_with_interaction) {
      p <- p + theme(legend.position = "none")
    }

    if (i != 1 & ((i-1) %% ncol != 0)) {
      p <- p + theme(axis.title.y = element_blank())
    }

    p
  })

  
  if (align_axes) {
    y_limits <- range(sapply(plots, function(p) ggplot_build(p)$layout$panel_params[[1]]$y.range))
    plots <- lapply(plots, function(p) p + coord_cartesian(ylim = y_limits))
    if (!is.na(axis_limits[[1]]) & !is.na(axis_limits[[2]])) {
      plots <- lapply(plots, function(p) p + coord_cartesian(ylim = axis_limits))
    }
  }

  wrap_plots(plots, ncol = ncol, guides = "collect") +
    plot_annotation(tag_levels = "A",
                    theme = theme(plot.tag = element_text(face = "bold", size = 14),
                                  plot.tag.position = c(0.05, 0.95)))  &
    theme(legend.position = "right")
}

