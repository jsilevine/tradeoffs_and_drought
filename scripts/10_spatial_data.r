##---------------------------------------------------------------
## 10_spatial_data.r: Append miscellaneous spatial data to the full dataset
##
## author: jacob levine; jacob.levine@utah.edu
##---------------------------------------------------------------

## Load required libraries
invisible(lapply(c("here", "ggplot2", "dplyr", "sf", "gstat", "sp"), library, character.only = TRUE))

## data average to grid by plots (i.e. tradeoffs fit for each plot, then averaged)
data <- read.csv(here("data", "full_data.csv"))
data_sf <- st_as_sf(data, coords = c("lon", "lat"), crs = 4326)

## read in ecoregion data (level I and II)
l2_ecoregions <- st_read(here("data", "ecoregions", "level_2_ecoregions", "NA_CEC_Eco_Level2.shp"))
l2_ecoregions <- st_transform(l2_ecoregions, crs = 4326)
l2_ecoregions <- st_make_valid(l2_ecoregions)

l1_ecoregions <- st_read(here("data", "ecoregions", "level_1_ecoregions", "NA_CEC_Eco_Level1.shp"))
l1_ecoregions <- st_transform(l1_ecoregions, crs = 4326)
l1_ecoregions <- st_make_valid(l1_ecoregions)

## spatial join to append ecoregion info
data_sf <- st_join(data_sf, l2_ecoregions[,2], join = st_intersects)
data_sf <- st_join(data_sf, l1_ecoregions[,2], join = st_intersects)

## extract coordinates
coords <- st_coordinates(data_sf)
data <- as.data.frame(data_sf)
data$lon <- coords[, 1]
data$lat <- coords[, 2]

## drop geometry column
data$geometry <- NULL
data <- data[,c("lon", "lat", setdiff(names(data), c("lon", "lat")))]

## overwrite full_data.csv with appended ecoregion info
write.csv(data, here("data", "full_data.csv"), row.names = FALSE)


