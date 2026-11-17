##-------------------------------------------------------------------
## A collection of julia wrapper functions for calling plotting
## utilities written in R
##
## author: jacob levine; email: jacob.levine@utah.edu 
##-------------------------------------------------------------------

module PlottingUtils

## dependencies
using RCall, DataFrames

include(joinpath(@__DIR__, "..", "..", "src", "julia", "utils.jl"))

using .Utils

growth_predictions_path = srcdir("r", "plot_growth_predictions.r")
mortality_predictions_path = srcdir("r", "plot_mortality_predictions.r")
plot_smooths_path = srcdir("r", "plot_spatial_smooths.r")
plot_variograms_path = srcdir("r", "plot_variograms.r")

R"source($growth_predictions_path)"
R"source($mortality_predictions_path)"
R"source($plot_smooths_path)"
R"source($plot_variograms_path)"


## define functions

"""
    plot_multiple_gam_preditions(model::RObject, data::DataFrame;
vars_to_plot = nothing, ci_level = 0.95, ncol = 4, align_axes = false,
save_output = false, output_path = nothing)

Plot the predictions of a growth_model fit using ModelFitting.bam() by calling utility
functions written in R. 
"""
function plot_growth_model_predictions(model::RObject, data::DataFrame;
    interaction_vars = ["beta_pca", "td_pca", "drought_strength"],
    vars_to_plot = nothing, ci_level = 0.95, ncol = 4, align_axes = false,
    save_output = false, output_path = nothing, axis_limits = [missing, missing])

    R"""
    plot_multiple_gam_predictions($model, $data,
        interaction_vars = $interaction_vars,
        vars_to_plot = $(vars_to_plot === nothing ? RObject("NULL") : vars_to_plot),
        ci_level = $ci_level,
        ncol = $ncol,
        align_axes = $align_axes,
        axis_limits = c($axis_limits[1], $axis_limits[2]))
    """

    if save_output && isnothing(output_path)
        error("To save output, you must provide an output_path")
    end

    if save_output
        R"""
            ggsave($output_path, width = ($ncol * 3.5) + 1, height = 4, dpi = 300)
        """
    end
    
end

"""
plot_spatial_smooth(model::RObject, data::DataFrame;
    year = 2015, lon_range = [-125, -66], lat_range = [25, 50], grid_resolution = 200,
save_output = false, output_path = nothing)

Plot the spatial smooth of a growth_model fit using ModelFitting.bam() by calling R.

"""
function plot_spatial_smooth(model::RObject, data::DataFrame;
    year = 2015, lon_range = [-125, -66], lat_range = [25, 50], grid_resolution = 200,
    save_output = false, output_path = nothing)

    R"""
    plot_spatial_smooth($model, $data,
    year = $year, lon_range = $lon_range, lat_range = $lat_range,
    grid_resolution = $grid_resolution)
    """

    if save_output && isnothing(output_path)
        error("To save output, you must provide an output_path")
    end

    if save_output
        R"""
        ggsave($output_path, width = 4.5, height = 4, dpi = 300)
        """
    end
       
end

function plot_residual_variogram(model::RObject,
    save_output::Bool = false, output_path = nothing)

    R"""
    model_frame <- model.frame($model)
    if (inherits($model, "gamlss")) {
    residuals <- residuals($model, type = "simple")
    colnames(model_frame)[colnames(model_frame) == "pb(lon, df = 18)"] <- "lon"
    colnames(model_frame)[colnames(model_frame) == "pb(lat, df = 18)"] <- "lat"
    } else {
    residuals <- residuals($model, type = "pearson")
    }
    model_frame$residuals <- residuals

    print(colnames(model_frame))
    """
    
    R"plot <- residual_variogram(model_frame)"
    
    if save_output && isnothing(output_path)
        error("To save output, you must provide an output_path")
    end

    if save_output
        R"""
        ggsave($output_path, width = 4.5, height = 4, dpi = 300)
        """
    end

    
    
end


function plot_acf(model::RObject,
    save_output::Bool = false, output_path = nothing)

    R"""
    model_frame <- model.frame($model)
    if (inherits($model, "gamlss")) {
    residuals <- residuals($model, type = "simple")
    colnames(model_frame)[colnames(model_frame) == "pb(year_centered)"] <- "year"
    } else {
    residuals <- residuals($model, type = "pearson")
    }
    model_frame$residuals <- residuals

    acf_obj <- acf(model_frame$residuals[order(model_frame$year)],
        plot = FALSE, na.action = na.pass)
    
    acf_df <- with(acf_obj, data.frame(lag = lag[-1], acf = acf[-1]))
    
    ggplot(acf_df, aes(x = lag, y = acf)) +
    geom_col(fill = "steelblue") +
    geom_hline(yintercept = c(-1.96/sqrt(length(model_frame$residuals)),
        1.96/sqrt(length(model_frame$residuals))),
        linetype = "dashed", color = "red") +
    labs(title = "ACF of residuals by year",
        x = "Lag", y = "Autocorrelation") +
    theme_bw()
    """
       
    if save_output && isnothing(output_path)
        error("To save output, you must provide an output_path")
    end

    if save_output
        R"""
        ggsave($output_path, width = 4.5, height = 4, dpi = 300)
        """
    end
    
end


function plot_group_acf(model::RObject;
    group_var::AbstractString = "NA_L2NAME",
    save_output::Bool = false,
    output_path = nothing)

    R"""
    library(ggplot2)

    # get model frame + residuals
    model_frame <- model.frame($model)

    if (inherits($model, "gamlss")) {
      residuals <- residuals($model, type = "simple")
      colnames(model_frame)[colnames(model_frame) == "pb(year_centered)"] <- "year"
    } else {
      residuals <- residuals($model, type = "pearson")
    }
    model_frame$residuals <- residuals

    group_var <- $group_var

    groups <- unique(model_frame[[group_var]])
    acf_list <- list()

    for (g in groups) {
      sub <- model_frame[model_frame[[group_var]] == g, ]
      # order by year
      sub <- sub[order(sub$year), ]
      n <- nrow(sub)
      if (n > 1) {
        ao <- acf(sub$residuals, plot = FALSE, na.action = na.pass)
        df <- data.frame(
          lag = ao$lag[-1],
          acf = ao$acf[-1],
          group = rep(g, length(ao$lag) - 1),
          n = n
        )
        acf_list[[length(acf_list) + 1]] <- df
      }
    }

    acf_df <- do.call(rbind, acf_list)

    p <- ggplot(acf_df, aes(x = lag, y = acf)) +
      geom_col(fill = "steelblue") +
      geom_hline(aes(yintercept = 1.96/sqrt(n)), linetype = "dashed", color = "red") +
      geom_hline(aes(yintercept = -1.96/sqrt(n)), linetype = "dashed", color = "red") +
      facet_wrap(~ group, scales = "free_y") +
      labs(title = paste0("Within-group ACF of residuals (by ", group_var, ")"),
           x = "Lag", y = "Autocorrelation") +
      theme_bw()

    print(p)

    if ($save_output) {
      if (is.null($output_path)) stop("To save output, must provide output_path")
      ggsave($output_path, width = 5, height = 4, dpi = 300)
    }
    """
end



function plot_mortality_model_predictions(model::RObject, data::DataFrame;
    interaction_vars = ["beta_pca", "se_pca", "drought_strength"],
    vars_to_plot = nothing, ci_level = 0.95, ncol = 4, align_axes = false,
    save_output = false, output_path = nothing, axis_limits = [missing, missing])

    R"""
    plot_multiple_bezi_predictions($model, $data,
        interaction_vars = $interaction_vars,
        vars_to_plot = $(vars_to_plot === nothing ? RObject("NULL") : vars_to_plot),
        ci_level = $ci_level,
        ncol = $ncol,
        align_axes = $align_axes,
        axis_limits = c($axis_limits[1], $axis_limits[2]))
    """
    
    if save_output && isnothing(output_path)
        error("To save output, you must provide an output_path")
    end

    if save_output
        R"""
            ggsave($output_path, width = ($ncol * 3.5) + 1, height = 4, dpi = 300)
        """
    end
    
end



function plot_spatial_smooth_bezi(model::RObject, data::DataFrame;
    year = 2015, lon_range = [-125, -66], lat_range = [25, 50], grid_resolution = 200,
    save_output = false, output_path = nothing)

    R"""
    plot_spatial_smooth_bezi($model, $data,
    year = $year, lon_range = $lon_range, lat_range = $lat_range,
    grid_resolution = $grid_resolution)
    """

    if save_output && isnothing(output_path)
        error("To save output, you must provide an output_path")
    end

    if save_output
        R"""
        ggsave($output_path, width = 4.5, height = 4, dpi = 300)
        """
    end
       
end


end
