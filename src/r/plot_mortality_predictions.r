
library(gamlss)
library(here)
library(msm)
library(ggplot2)
library(MetBrewer)
library(patchwork)
library(maps)

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

pred_deltamethod <- function(i, X_mu, X_nu, beta_mu, beta_nu, vcv) {

    eta_mu <- drop(X_mu[i, ] %*% beta_mu)
    eta_nu <- drop(X_nu[i, ] %*% beta_nu)

                                        # Apply inverse links
    mu_linkinv <- BEZI()$mu.linkinv
    nu_linkinv <- BEZI()$nu.linkinv

    mu_hat <- mu_linkinv(eta_mu)
    nu_hat <- nu_linkinv(eta_nu)

    comb_mean <- ~ (1 - (1/(1 + exp(-x2)))) * (1/(1 + exp(-x1)))

    var_eta_mu <- X_mu[i, , drop=FALSE] %*% vcv[names(beta_mu), names(beta_mu)] %*% t(X_mu[i, , drop=FALSE])
    var_eta_nu <- X_nu[i, , drop=FALSE] %*% vcv[names(beta_nu), names(beta_nu)] %*% t(X_nu[i, , drop=FALSE])
    cov_eta    <- X_mu[i, , drop=FALSE] %*% vcv[names(beta_mu), names(beta_nu)] %*% t(X_nu[i, , drop=FALSE])

    V_eta <- matrix(c(var_eta_mu, cov_eta,
                      cov_eta,    var_eta_nu), nrow=2)

    se_mean <- msm::deltamethod(comb_mean, c(eta_mu, eta_nu), V_eta)

                                        # Combined expected mortality
    expected <- (1 - nu_hat) * mu_hat

    return(c(mean_hat = (1 - nu_hat) * mu_hat,
             se_mean  = se_mean))
}


plot_one_bezi_prediction <- function(model, data, variable, interaction = NULL,
                                 n_points = 100, n_interaction_levels = 3) {

    model_vars <- all.vars(formula(model))[-1]

    pred_data <- data.frame(matrix(ncol = length(model_vars), nrow = n_points * max(1, n_interaction_levels)))
    colnames(pred_data) <- model_vars

                                        # Fill other covariates with medians
    vars_to_fix <- setdiff(model_vars, c(variable, interaction))
    for (var in vars_to_fix) {
        if (!is.null(data[[var]])) {
            pred_data[[var]] <- median(data[[var]], na.rm = TRUE)
        }
    }

    scaled_vars <- grep("_scaled", model_vars, value = TRUE)

    for (v in scaled_vars) {
        if (!v %in% names(pred_data)) {
            pred_data[[v]] <- median(data[[v]], na.rm = TRUE)
        }
    }

    var_s <- paste0(variable, "_scaled")
    pred_data[[var_s]] <-
        rep(seq(quantile(data[[var_s]], 0.05, na.rm = TRUE),
                quantile(data[[var_s]], 0.95, na.rm = TRUE),
                length.out = n_points),
            times = ifelse(is.null(interaction), 1, n_interaction_levels))

    if (!is.null(interaction)) {
        int_s <- paste0(interaction, "_scaled")
        interaction_levels <- seq(quantile(data[[int_s]], 0.1, na.rm = TRUE),
                                  quantile(data[[int_s]], 0.9, na.rm = TRUE),
                                  length.out = n_interaction_levels)
        pred_data[[int_s]] <- rep(interaction_levels, each = n_points)
    }
    pred_data$mort_tpa_pct <- 0
    rownames(pred_data) <- NULL

                                        # Design matrices
    X_mu <- model.matrix(delete.response(terms(model, "mu")), data = pred_data)
    X_nu <- model.matrix(delete.response(terms(model, "nu")), data = pred_data)

                                        # Coefficients
    beta_mu <- coef(model, what = "mu")
    beta_nu <- coef(model, what = "nu")

    vcv <- vcov(model)

    for (i in 1:nrow(pred_data)) {

        pred_data[i, c("fit", "se")] <-
            pred_deltamethod(i, X_mu, X_nu, beta_mu, beta_nu, vcv)

    }

    pred_data$lower <- pred_data$fit - 1.96 * pred_data$se
    pred_data$upper <- pred_data$fit + 1.96 * pred_data$se

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
            y = "Expected Annual Mortality (% TPA)"
        ) +
        scale_x_continuous(expand = c(0, 0)) +
        theme_bw()

    if (!is.null(interaction)) {
        p <- p +
            scale_color_gradientn(colors = met.brewer("Archambault", n = 5, type = "continuous"), name = int_lab) +
            scale_fill_gradientn(colors = met.brewer("Archambault", n = 5, type = "continuous"), name = int_lab)
    }

    return(p)

}


plot_multiple_bezi_predictions <- function(model, data,
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
        p <- plot_one_bezi_prediction(
            model = model,
            data = data,
            variable = vars[[i]]$variable,
            interaction = vars[[i]]$interaction
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
        if (!is.na(axis_limits[1]) & !is.na(axis_limits[2])) {
            plots <- lapply(plots, function(p) p + coord_cartesian(ylim = axis_limits))
        }
    }

    wrap_plots(plots, ncol = ncol, guides = "collect") +
        plot_annotation(
            tag_levels = 'A',
            theme = theme(plot.tag = element_text(face = "bold", size = 14),
                          plot.tag.position = c(0.05, 0.95)))  &
        theme(legend.position = "right")
}


plot_spatial_smooth_bezi <- function(model, data, year = 2015, ecoregion = NULL,
                                     lon_range = c(-125, -66), lat_range = c(25, 50),
                                     grid_resolution = 200) {

    # prediction grid
    lon_seq <- seq(lon_range[1], lon_range[2], length.out = grid_resolution)
    lat_seq <- seq(lat_range[1], lat_range[2], length.out = grid_resolution)
    grid <- expand.grid(lon = lon_seq, lat = lat_seq)

    # Extract model variables from formula (mu formula)
    model_vars <- all.vars(formula(model, lhs = NULL))  # variables on RHS

    # Exclude spatial & grouping vars
    exclude_vars <- c("lon", "lat", "year")
    variables_to_mean <- setdiff(model_vars, exclude_vars)

    # Fill grid with mean values of covariates
    for (var in variables_to_mean) {
        if (var %in% names(data)) {
            grid[[var]] <- mean(data[[var]], na.rm = TRUE)
        } else {
            # if variable not in data (e.g. interaction?), set NA or drop
            grid[[var]] <- NA
        }
    }

                                        # Set year and ecoregion
    grid$year_centered <- year

    print(colnames(grid))

    grid$predicted <- predict(model, what = "mu", newdata = grid, type = "response", se.fit = FALSE)
    
                                        # Convert to SpatRaster for plotting
    raster_smooth <- rast(grid[, c("lon", "lat", "predicted")], type = "xyz", crs = "EPSG:4326")

    raster_df <- as.data.frame(raster_smooth, xy = TRUE)


    conus_map <- map_data("state")

                                        # Plot prediction surface
    ggplot() +
        geom_raster(data = raster_df, aes(x = x, y = y, fill = predicted)) +
        scale_fill_viridis_c(name = "Predicted Value", option = "magma", direction = -1) +
        geom_polygon(data = conus_map, aes(x = long, y = lat, group = group),
                     fill = NA, color = "black") +
        coord_fixed(1.3) +
        scale_y_continuous(expand = c(0, 0)) +
        scale_x_continuous(expand = c(0, 0)) +
        labs(title = paste0("Predicted Response Surface (Year = ", year, ")"),
             x = "Longitude",
             y = "Latitude") +
        theme_bw()
}

