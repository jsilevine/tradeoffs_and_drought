##---------------------------------------------------------------
## 03_create_grid.r: Create raster grid for analyses
##
## author: jacob levine; jacob.levine@utah.edu
##---------------------------------------------------------------

library(sf)
library(terra)
library(ggplot2)
library(here)
library(tidyterra)

# Load the data
data <- read.csv(here("data", "gm", "gm_full.csv"))

# Remove AK and HI
data <- data[data$STATECD != 2 & data$STATECD != 15,]

# Convert to sf object
sp <- st_as_sf(data[, c("LON", "LAT", "STATECD", "MEASYEAR")],
               coords = c("LON", "LAT"), crs = 4326)

# Bounding box of points
bbox <- st_bbox(sp)

# Create empty raster template with 0.25° resolution
template <- rast(
  xmin = floor(bbox["xmin"]),
  xmax = ceiling(bbox["xmax"]),
  ymin = floor(bbox["ymin"]),
  ymax = ceiling(bbox["ymax"]),
  resolution = c(0.25, 0.25),
  crs = "EPSG:4326"
)

# Visualize the grid
grid_sf <- st_as_sf(as.polygons(template))

# Assign unique grid cell ID to all cells (global indexing)
values(template) <- 1:ncell(template)

# Convert sp to terra vect for raster ops
sp_vect <- vect(sp)

# Assign each plot to a grid cell
data$grid_cell_id <- cellFromXY(template, st_coordinates(sp))

# Rasterize: count number of plots per cell
r_count <- rasterize(sp_vect, template, fun = "count")
plot(r_count, main = "Count of Points per Cell")
writeRaster(r_count, here("data", "template", "plot_count.tif"), overwrite = TRUE)

r_count_proj <- terra::project(r_count, "EPSG:5070")

ggplot() +
  geom_spatraster(data = r_count_proj, na.rm = TRUE) +
  scale_fill_viridis_c(na.value = "transparent", name = "# Obs.") +
  scale_x_continuous(name = "Longitude", expand = c(0,0)) +
  scale_y_continuous(name = "Latitude", expand = c(0,0)) +
  labs(title = "Number of Unique Observations Per Cell") +
  theme_bw()

ggsave(here("figures", "maps", "plot_count_map.png"),
       width = 8, height = 6)

# Binary presence/absence
r_binary <- classify(r_count, cbind(1, Inf, 1))
plot(r_binary, main = "Binary Presence/Absence of Points")
writeRaster(r_binary, here("data", "templates", "plot_bin.tif"), overwrite = TRUE)

# Rasterize min and max measurement year per cell
r_min <- rasterize(sp_vect, template, field = "MEASYEAR", fun = min, na.rm = TRUE)
r_max <- rasterize(sp_vect, template, field = "MEASYEAR", fun = max, na.rm = TRUE)

# Save plots
plot(r_min, main = "Minimum MEASYEAR per Cell")
plot(r_max, main = "Maximum MEASYEAR per Cell")

# Combine raster values into grid_cells dataframe
grid_cells <- as.data.frame(c(r_count, r_min, r_max), xy = TRUE, na.rm = TRUE)
names(grid_cells) <- c("LON", "LAT", "plot_count", "min_year", "max_year")

# Assign matching grid cell ID for each row in grid_cells
grid_cells$grid_cell_id <- cellFromXY(template, grid_cells[, c("LON", "LAT")])

# Save template raster and grid cell summary
writeRaster(template, here("data", "templates", "grid_template.tif"), overwrite = TRUE)
write.csv(grid_cells, here("data", "templates", "grid_cells.csv"), row.names = FALSE)
write.csv(data, here("data", "gm", "gm_full_grid.csv"), row.names = FALSE)

## Save grid to plot mappings:
write.csv(data[,c("PLOT_ID", "grid_cell_id")], here("data", "templates", "plot_to_cell.csv"))

cell_to_plots <- split(data$PLOT_ID, data$grid_cell_id)
saveRDS(cell_to_plots, here("data", "templates", "cell_to_plots.rds"))

