library(tidyverse)
library(sf)
library(terra)
library(tmap)



# Read in site and radiocarbon data
sites <- read.csv(here::here('analysis/data/raw_data/sites.csv'))
rcarbon <- read.csv(here::here('analysis/data/raw_data/radiocarbon.csv'))

# Read in site limits
sitesl <- st_zm(read_sf(here::here('analysis/data/raw_data/site_limits.gpkg')))

# Read in site features
pts <- read_sf(here::here('analysis/data/raw_data/sitefeatures/points.gpkg'))
lines <- read_sf(here::here('analysis/data/raw_data/sitefeatures/lines.gpkg'))
polygons <- read_sf(here::here('analysis/data/raw_data/sitefeatures/polygons.gpkg'))

# Combine into a single object
sitesf <- st_zm(rbind(pts, lines, polygons))

tm_shape(st_centroid(sitesl)) +
  tm_dots()

# Assign placeholder for NA in Askeladden ID for join
rcarbon <- rcarbon %>%  mutate(askeladden_id = ifelse(is.na(askeladden_id), 0, askeladden_id))

rcarb_features <- left_join(rcarbon, sitesf, by = c("site_name", "context" = "name", "askeladden_id" = "ask_id"))

rcarb_features %>%  filter(is.na(st_dimension(rcarb_features$geom)) == TRUE) %>%
  select(- geom) %>%
 left_join(sitesl, by = c("site_name", "askeladden_id" = "ask_id")) %>% View()


