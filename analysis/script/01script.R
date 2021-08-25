library(tidyverse)
library(sf)
library(tmap)
library(oxcAAR)
quickSetupOxcal()

# Read in site and radiocarbon data
sites <- read.csv(here::here('analysis/data/raw_data/sites.csv'))
rcarbon <- read.csv(here::here('analysis/data/raw_data/radiocarbon.csv'))

# Read in site limits and combine into a single object
sitespol <- read_sf(here::here('analysis/data/raw_data/site_limits/site_limits.gpkg'))
siteslpts <- read_sf(here::here('analysis/data/raw_data/site_limits/site_limits_pts.gpkg'))
sitesl <- st_zm(rbind(sitespol, siteslpts))


# Read in site features and combine
pts <- read_sf(here::here('analysis/data/raw_data/site_features/points.gpkg'))
lines <- read_sf(here::here('analysis/data/raw_data/site_features/lines.gpkg'))
polygons <- read_sf(here::here('analysis/data/raw_data/site_features/polygons.gpkg'))
sitesf <- st_zm(rbind(pts, lines, polygons))

# Assign placeholder for NA in Askeladden ID for joins
rcarbon <- rcarbon %>%  mutate(ask_id = ifelse(ask_id == "", 0, ask_id))
sites <- sites %>% mutate(ask_id = ifelse(ask_id == "", 0, ask_id))
sitesl <- sitesl %>% mutate(ask_id = ifelse(is.na(ask_id), 0, ask_id))
sitesf <- sitesf %>% mutate(ask_id = ifelse(is.na(ask_id), 0, ask_id))

# Combine site data and site limits
site_limits <- st_as_sf(left_join(sites, sitesl,
                                  by = c("name" = "site_name", "ask_id")))

# Combine radiocarbon data and site features
rcarb_features <- left_join(rcarbon, sitesf,
                            by = c("site_name", "context" = "name", "ask_id"))

# Assign radiocarbon dates not specified to any feature/excavation unit to the
# site limit
rcarb_features <- rcarb_features %>%  left_join(sitesl,
                                      by = c("site_name", "ask_id")) %>%
  mutate(geom = ifelse(is.na(st_dimension(geom.x)), geom.y, geom.x),
         .keep = "unused") %>%
  dplyr::select(-name) %>% st_as_sf(crs = st_crs(sitesf))

# Calibrate radiocarbon dates and exclude dates that do not overlap with the
# Stone Age at two sigma.

# Calibrate dates
caldates <- oxcalCalibrate(rcarb_features$c14_bp, rcarb_features$error,
                           rcarb_features$lab_code)

# Retrieve start and end date with 2-sigma confidence
rcarb_features <- rcarb_features %>%
  mutate(sig_2_start_bc = map(caldates, ~ min(.x$sigma_ranges$two_sigma$start)),
         sig_2_end_bc = map(caldates, ~ max(.x$sigma_ranges$two_sigma$end)))

# Exclude date ranges falling outside the Stone Age
rcarb_sa <- rcarb_features %>%  filter(sig_2_end_bc < -1700)

# Use the above to retrieve overview of sites with radiocarbon dates to the
# Stone Age.
sites_sa <- site_limits %>% filter(name %in% unique(rcarb_sa$site_name))


save(caldates, rcarb_sa, sites_sa,
     file = here::here("analysis/data/derived_data/01data.RData"))
