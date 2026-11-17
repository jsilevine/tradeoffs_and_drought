##--------------------------------------------------------------
## Functions to visualize the spatial smooth term from a GAM model
##
## author: jacob levine; email: jacob.levine@utah.edu
##--------------------------------------------------------------

library(ggplot2)
library(terra)
library(maps)
library(sf)

plot_spatial_smooth <- function(model, data, year = 2015,
                                     lon_range = c(-125, -66), lat_range = c(25, 50), grid_resolution = 200) {

                                        # Define a grid of lat-lon values covering the CONUS
    lon_seq <- seq(lon_range[1], lon_range[2], length.out = grid_resolution)
    lat_seq <- seq(lat_range[1], lat_range[2], length.out = grid_resolution)
    grid <- expand.grid(lon = lon_seq, lat = lat_seq)

                                        # Extract variables from the model formula
    model_vars <- all.vars(formula(model))

                                        # Exclude variables that are part of the smooth terms (e.g., lon, lat, year, l2_ecoregion)
    exclude_vars <- c("lon", "lat", "year", "l2_ecoregion")
    variables_to_mean <- setdiff(model_vars, exclude_vars)

                                        # Dynamically assign mean values for each variable
    for (var in variables_to_mean) {
        grid[[var]] <- mean(data[[var]], na.rm = TRUE)
    }

                                        # Assign year and ecoregion
    grid$year <- year  # Set the year for prediction

                                        # Predict the smooth term for the grid
    grid$predicted_smooth <- predict(model,
                                     newdata = grid,
                                     type = "terms",
                                     terms = "s(lon,lat)")

                                        # Convert the grid to a SpatRaster
    raster_smooth <- rast(grid[, c("lon", "lat", "predicted_smooth")], type = "xyz", crs = "EPSG:4326")

                                        # Get a map of the CONUS
    states <- maps::map("state", plot = FALSE, fill = TRUE)
    states_sf <- sf::st_as_sf(states)
    states_vect <- terra::vect(states_sf)
    states_vect <- terra::project(states_vect, crs(raster_smooth))

    raster_smooth <- terra::crop(raster_smooth, states_vect)
    raster_smooth <- terra::mask(raster_smooth, states_vect)

    raster_df <- as.data.frame(raster_smooth, xy = TRUE)

    conus_map <- map_data("state")
                                        # Plot the raster over the map
    output <- ggplot() +
        geom_raster(data = raster_df, aes(x = x, y = y, fill = predicted_smooth)) +
        scale_fill_viridis_c(name = "Smooth Value", option = "magma", direction = -1) +
        geom_polygon(data = conus_map, aes(x = long, y = lat, group = group),
                     fill = NA, color = "black", linewidth = 0.3) +
        coord_fixed(1.3) +
        scale_y_continuous(expand = c(0, 0)) +
        scale_x_continuous(expand = c(0, 0)) +
        labs(title = "Lat-Lon Smooth Term",
             x = "Longitude",
             y = "Latitude") +
        theme_bw()

    return(output)
}
