##--------------------------------------------------------------
## Wrappers around plot_growth_predictions and
## plot_mortality_predictions that overlay raw data points.
##
## Each function calls its counterpart from the original source
## files and then prepends a geom_point layer so the raw
## observations appear behind the model line and CI ribbon.
##
## author: jacob levine; jacob.levine@utah.edu
##--------------------------------------------------------------

library(here)
source(here("src", "r", "plot_growth_predictions.r"))
source(here("src", "r", "plot_mortality_predictions.r"))


## ---------------------------------------------------------------
## Growth (GAM) wrappers
## ---------------------------------------------------------------

#' Like plot_one_gam_prediction() but with raw data underneath.
#'
#' @param point_alpha,point_size  appearance of raw data dots
#' @param point_color             fill/colour of dots (NA = inherit interaction colour)
#' @inheritParams plot_one_gam_prediction
plot_one_gam_prediction_wdata <- function(model, data, variable,
                                          interaction = NULL,
                                          n_points = 100,
                                          n_interaction_levels = 3,
                                          ci_level = 0.95,
                                          point_alpha = 0.15,
                                          point_size  = 0.8,
                                          point_color = "grey30") {

  p <- plot_one_gam_prediction(
    model               = model,
    data                = data,
    variable            = variable,
    interaction         = interaction,
    n_points            = n_points,
    n_interaction_levels = n_interaction_levels,
    ci_level            = ci_level
  )

  y_raw <- data$growth_pct * 100

  if (is.null(interaction)) {
    raw_layer <- geom_point(
      data  = data.frame(x = data[[variable]], y = y_raw),
      aes(x = x, y = y),
      alpha  = point_alpha,
      size   = point_size,
      color  = point_color,
      inherit.aes = FALSE
    )
  } else {
    raw_layer <- geom_point(
      data  = data.frame(x   = data[[variable]],
                         y   = y_raw,
                         grp = data[[interaction]]),
      aes(x = x, y = y, color = grp),
      alpha  = point_alpha,
      size   = point_size,
      inherit.aes = FALSE
    )
  }

  ## prepend so points sit behind line and ribbon
  p$layers <- c(list(raw_layer), p$layers)
  p
}


#' Like plot_multiple_gam_predictions() but with raw data underneath
#' each panel.
#'
#' @inheritParams plot_multiple_gam_predictions
#' @param point_alpha,point_size,point_color  passed to each panel
plot_multiple_gam_predictions_wdata <- function(model, data,
                                                interaction_vars = c("beta_pca", "se_pca", "drought_strength"),
                                                vars_to_plot     = NULL,
                                                ci_level         = 0.9,
                                                ncol             = 4,
                                                align_axes       = FALSE,
                                                axis_limits      = c(NA, NA),
                                                point_alpha      = 0.15,
                                                point_size       = 0.8,
                                                point_color      = "grey30") {

  model_vars <- all.vars(formula(model))[-1]
  model_vars <- model_vars[!model_vars %in% c("year", "lon", "lat")]
  model_vars <- gsub("_scaled$", "", model_vars)

  if (!is.null(vars_to_plot)) {
    model_vars <- intersect(model_vars, vars_to_plot)
  }

  vars <- lapply(model_vars, function(v) {
    list(variable = v, interaction = if (v %in% interaction_vars) "ba" else NULL)
  })

  if (is.null(vars_to_plot)) {
    vars <- c(
      vars[sapply(vars, function(x) !is.null(x$interaction))],
      vars[sapply(vars, function(x)  is.null(x$interaction))]
    )
  }

  first_with_interaction <- which(sapply(vars, function(x) !is.null(x$interaction)))[1]

  plots <- lapply(seq_along(vars), function(i) {
    p <- plot_one_gam_prediction_wdata(
      model               = model,
      data                = data,
      variable            = vars[[i]]$variable,
      interaction         = vars[[i]]$interaction,
      ci_level            = ci_level,
      point_alpha         = point_alpha,
      point_size          = point_size,
      point_color         = point_color
    )

    if (!is.na(first_with_interaction) && i != first_with_interaction) {
      p <- p + theme(legend.position = "none")
    }

    if (i != 1 && ((i - 1) %% ncol != 0)) {
      p <- p + theme(axis.title.y = element_blank())
    }

    p
  })

  if (align_axes) {
    y_limits <- range(sapply(plots, function(p) ggplot_build(p)$layout$panel_params[[1]]$y.range))
    plots <- lapply(plots, function(p) p + coord_cartesian(ylim = y_limits))
    if (!is.na(axis_limits[[1]]) && !is.na(axis_limits[[2]])) {
      plots <- lapply(plots, function(p) p + coord_cartesian(ylim = axis_limits))
    }
  }

  wrap_plots(plots, ncol = ncol, guides = "collect") +
    plot_annotation(tag_levels = "A",
                    theme = theme(plot.tag          = element_text(face = "bold", size = 14),
                                  plot.tag.position = c(0.05, 0.95)))
}


## ---------------------------------------------------------------
## Mortality (BEZI) wrappers
## ---------------------------------------------------------------

#' Like plot_one_bezi_prediction() but with raw data underneath.
#'
#' @param point_alpha,point_size,point_color  appearance of raw data dots
#' @inheritParams plot_one_bezi_prediction
plot_one_bezi_prediction_wdata <- function(model, data, variable,
                                           interaction = NULL,
                                           n_points = 100,
                                           n_interaction_levels = 3,
                                           point_alpha = 0.15,
                                           point_size  = 0.8,
                                           point_color = "grey30") {

  p <- plot_one_bezi_prediction(
    model                = model,
    data                 = data,
    variable             = variable,
    interaction          = interaction,
    n_points             = n_points,
    n_interaction_levels = n_interaction_levels
  )

  y_raw <- data$mort_tpa_pct

  if (is.null(interaction)) {
    raw_layer <- geom_point(
      data  = data.frame(x = data[[variable]], y = y_raw),
      aes(x = x, y = y),
      alpha  = point_alpha,
      size   = point_size,
      color  = point_color,
      inherit.aes = FALSE
    )
  } else {
    raw_layer <- geom_point(
      data  = data.frame(x   = data[[variable]],
                         y   = y_raw,
                         grp = data[[interaction]]),
      aes(x = x, y = y, color = grp),
      alpha  = point_alpha,
      size   = point_size,
      inherit.aes = FALSE
    )
  }

  p$layers <- c(list(raw_layer), p$layers)
  p
}


#' Like plot_multiple_bezi_predictions() but with raw data underneath
#' each panel.
#'
#' @inheritParams plot_multiple_bezi_predictions
#' @param point_alpha,point_size,point_color  passed to each panel
plot_multiple_bezi_predictions_wdata <- function(model, data,
                                                 interaction_vars = c("beta_pca", "se_pca", "drought_strength"),
                                                 vars_to_plot     = NULL,
                                                 ci_level         = 0.9,
                                                 ncol             = 4,
                                                 align_axes       = FALSE,
                                                 axis_limits      = c(NA, NA),
                                                 point_alpha      = 0.15,
                                                 point_size       = 0.8,
                                                 point_color      = "grey30") {

  model_vars <- all.vars(formula(model))[-1]
  model_vars <- model_vars[!model_vars %in% c("year", "lon", "lat")]
  model_vars <- gsub("_scaled$", "", model_vars)

  if (!is.null(vars_to_plot)) {
    model_vars <- intersect(model_vars, vars_to_plot)
  }

  vars <- lapply(model_vars, function(v) {
    list(variable = v, interaction = if (v %in% interaction_vars) "ba" else NULL)
  })

  if (is.null(vars_to_plot)) {
    vars <- c(
      vars[sapply(vars, function(x) !is.null(x$interaction))],
      vars[sapply(vars, function(x)  is.null(x$interaction))]
    )
  }

  first_with_interaction <- which(sapply(vars, function(x) !is.null(x$interaction)))[1]

  plots <- lapply(seq_along(vars), function(i) {
    p <- plot_one_bezi_prediction_wdata(
      model                = model,
      data                 = data,
      variable             = vars[[i]]$variable,
      interaction          = vars[[i]]$interaction,
      point_alpha          = point_alpha,
      point_size           = point_size,
      point_color          = point_color
    )

    if (!is.na(first_with_interaction) && i != first_with_interaction) {
      p <- p + theme(legend.position = "none")
    }

    if (i != 1 && ((i - 1) %% ncol != 0)) {
      p <- p + theme(axis.title.y = element_blank())
    }

    p
  })

  if (align_axes) {
    y_limits <- range(sapply(plots, function(p) ggplot_build(p)$layout$panel_params[[1]]$y.range))
    plots <- lapply(plots, function(p) p + coord_cartesian(ylim = y_limits))
    if (!is.na(axis_limits[[1]]) && !is.na(axis_limits[[2]])) {
      plots <- lapply(plots, function(p) p + coord_cartesian(ylim = axis_limits))
    }
  }

  wrap_plots(plots, ncol = ncol, guides = "collect") +
    plot_annotation(tag_levels = "A",
                    theme = theme(plot.tag          = element_text(face = "bold", size = 14),
                                  plot.tag.position = c(0.05, 0.95)))
}
