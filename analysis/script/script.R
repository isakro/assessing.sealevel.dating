library(tidyverse)
library(sf)
library(terra)
library(tmap)



# Read in site and radiocarbon data
sites <- read.csv(here::here('analysis/data/raw_data/sites.csv'))
rcarbon <- read.csv(here::here('analysis/data/raw_data/radiocarbon.csv'))

# Read in site limits
sitesl <- read_sf(here::here('analysis/data/raw_data/site_limits.gpkg'))

# Read in site features
pts <- read_sf(here::here('analysis/data/raw_data/sitefeatures/points.gpkg'))
lines <- read_sf(here::here('analysis/data/raw_data/sitefeatures/lines.gpkg'))
polygons <- read_sf(here::here('analysis/data/raw_data/sitefeatures/polygons.gpkg'))

# Combine into a single object
sitesf <- rbind(pts, lines, polygons)

tm_shape(st_centroid(st_zm(sitesl))) +
  tm_dots()



rcarb_features <- left_join(rcarbon, sitesf, by = c("site_name", "context" = "name", "askeladden_id" = "ask_id"))

View(rcarb_features)

