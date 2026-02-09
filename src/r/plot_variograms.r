##--------------------------------------------------------------
## Functions to generate semi-variograms and assess autocorrelation                                                        
##--------------------------------------------------------------

library(ggplot2)
library(gstat)
library(sf)

## Residual variogram function
residual_variogram <- function(data) {

  ## Remove duplicate spatial locations
  unique_data <- data[!duplicated(data[, c("lon", "lat")]), ]
  complete_data <- unique_data[!is.na(unique_data$residuals), ]

  ## Convert to sf and project to planar CRS
  spatial_data <- st_as_sf(complete_data, coords = c("lon", "lat"), crs = 4326)
  projected_data <- st_transform(spatial_data, crs = 5070)
  coordinates <- st_coordinates(projected_data)

  residuals_df <- data.frame(z = complete_data$residuals, x = coordinates[, "X"], y = coordinates[, "Y"])

  ## Define distance parameters
  min_distance <- 30000         ## meters (adjust if needed)
  max_distance <- 200000        ## meters (e.g. 200 km)
  num_bins <- 20
  bin_width <- (max_distance - min_distance) / num_bins

  ## Define custom breaks to skip the smallest lag class
  distance_breaks <- seq(min_distance, max_distance, by = bin_width)

  ## Compute empirical variogram
  empirical_variogram <- variogram(z ~ 1, data = residuals_df, locations = ~x + y,
                                   cutoff = max_distance, boundaries = distance_breaks)

  ## Robust variogram fitting
  fitted_variogram <- tryCatch(
    fit.variogram(empirical_variogram, vgm("Exp"), fit.method = 2),
    error = function(e) {
      message("Variogram fitting failed: ", e$message)
      return(NULL)
    }
  )

  if (is.null(fitted_variogram)) return(NULL)

  ## Extract estimated range
  estimated_range <- fitted_variogram$range[2]
  message("Estimated range of spatial autocorrelation: ", round(estimated_range / 1000, 1), " km")

  ## Plot
  fitted_values <- variogramLine(fitted_variogram, maxdist = max(empirical_variogram$dist))
  variogram_df <- as.data.frame(empirical_variogram)

  plot <- ggplot(variogram_df, aes(x = dist, y = gamma)) +
    geom_point() +
    geom_line(data = fitted_values, aes(x = dist, y = gamma), color = "blue") +
    geom_vline(xintercept = estimated_range, color = "red", linetype = "dashed") +
    scale_y_continuous(limits = c(0, NA)) +
    scale_x_continuous(limits = c(min_distance, max_distance), expand = c(0, 0)) +
    labs(title = "Residual Variogram",
         x = "Distance (meters)",
         y = "Semivariance") +
    theme_bw()

  return(plot)
}

