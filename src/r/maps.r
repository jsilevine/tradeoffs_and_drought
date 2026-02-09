library(sf)
library(ggplot2)
library(rnaturalearth)
library(ggspatial)
library(scales)
library(terra)

## create base map of conus
base_map<- function(crs = 5070) {
  conus_map <- ne_states(country = "United States of America", returnclass = "sf")
  conus_map <- conus_map[!(conus_map$postal %in% c("AK", "HI", "PR")), ]
  conus_map <- st_transform(conus_map, crs)

  ggplot() +
    geom_sf(data = conus_map,
            fill = "white", color = "black") +
    scale_y_continuous(expand = c(0,0.1)) +
    scale_x_continuous(expand = c(0.001,0.1)) +
    annotation_scale(location = "bl", width_hint = 0.2) +
    annotation_north_arrow(location = "bl", which_north = "true",
                           pad_x = unit(0.08, "in"), pad_y = unit(0.5, "in"),
                           style = north_arrow_fancy_orienteering) +
    coord_sf(crs = st_crs(crs)) +
    theme_bw(base_size = 12) +
    theme(
      legend.position = "none",
      panel.grid.major = element_line(color = "grey95"),
      panel.background = element_rect(fill = "white", colour = NA),
      axis.title = element_blank(),
      axis.text = element_text(color = "grey30")
    )

}

plot_trend_map <- function(data, var, lims = c(NA, NA), crs = 5070) {

  conus_map <- rnaturalearth::ne_states(country = "United States of America", returnclass = "sf")
  conus_map <- conus_map[!(conus_map$postal %in% c("AK", "HI", "PR")), ]
  conus_aea <- sf::st_transform(conus_map, crs)

  r <- terra::rast(data[, c("lon", "lat", var)], type = "xyz", crs = "EPSG:4326")

  r_proj <- terra::project(r, paste0("EPSG:", crs), method = "bilinear")

  conus_vect <- terra::vect(sf::st_geometry(conus_aea))
  r_masked <- terra::mask(r_proj, conus_vect)

  r_fine <- terra::disagg(r_proj, fact = 3)
  r_masked_fine <- terra::mask(r_fine, conus_vect)
  r_df_fine <- as.data.frame(r_masked_fine, xy = TRUE, na.rm = TRUE)
  colnames(r_df_fine)[3] <- "value"

  r_df <- as.data.frame(r_masked, xy = TRUE, na.rm = TRUE)
  if (ncol(r_df) < 3) stop("Projected raster has no value column after masking.")
  colnames(r_df)[3] <- "value"

  if (all(is.na(lims))) {
    lims_use <- range(r_df$value, na.rm = TRUE)
  } else {
    lims_use <- lims
  }

  p <- base_map() + 
    ggplot2::geom_tile(data = r_df_fine, aes(x = x, y = y, fill = value)) +
    ggplot2::scale_fill_gradient2(
      low = "#B2182B",  high = "#2166AC",
      limits = lims_use,
      oob = scales::squish,
      trans = "asinh",
      name = if (!is.null(name_dict) && !is.null(name_dict[[var]])) name_dict[[var]] else var
    ) +
    ggplot2::labs(title = if (!is.null(name_dict) && !is.null(name_dict[[var]])) name_dict[[var]] else var) +
    theme(legend.position = c(0.15, 0.15))

  return(p)
}

