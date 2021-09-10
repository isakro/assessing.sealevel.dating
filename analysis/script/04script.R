library(here)
library(tidyverse)
library(sf)
library(ggridges)
library(tmap)
library(tmaptools)
library(cowplot)
library(terra)
library(oxcAAR)
library(vapour)
library(topoDistance)
quickSetupOxcal()

# Load required functions and prepared data
source(here("analysis/script/03script.R"))
load(here("analysis/data/derived_data/01data.RData"))
load(here("analysis/data/derived_data/02data.RData"))

# Specify path to raster data (currently stored locally)
dtmpath <- "/home/isak/phd/eaa_presentation/dtm1"

# Assign closest isobases to site data (largest = TRUE makes sure that
# sites overlapping isobase polygons are only ascribed the one they overlap the
# most, so as not to create duplicate records).
sites_sa <- st_join(st_make_valid(sites_sa), isopolys,
                    join = st_intersects, largest = TRUE)

# Example site
sitename <- "Langangen Vestgård 1"

sitel <- filter(sites_sa, name == sitename)
siter <- filter(rcarb_sa, site_name == sitename)

# Model dates
datedat <- model_dates(siter)

# Interpolate displacement curve to the site location
sitecurve <- interpolate_curve(years = xvals,
                               isobase1 = sitel$isobase1,
                               isobase2 = sitel$isobase2,
                               target = sitel,
                               dispdat = displacement_curves,
                               isodat = isobases,
                               direction_rel_curve1 = sitel$dir_rel_1)
# Add site name
sitecurve$name <- sitename

# Load correct regional raster
dtm <- load_raster(dtmpath, sitel)

# Create bounding box polygon
location_bbox <- bboxpoly(sitel, 250)

# Use this to clip the dtm to the site area
sitearea <- terra::crop(dtm, location_bbox)

# Simulate sea-level and retrieve distances (takes some time)
output <- sample_shoreline(samps = 1000, sitel, sitecurve, sitearea,
                           posteriorprobs)

# Generate grid with dtm resolution holding number of overlaps for each cell
# (also takes quite some time to execute)
simsea <- sea_overlaps(sitearea, output$seapol)

# Create plots
datplot <- plot_dates(datedat, sitename)
smap <- tmap_grob(shore_plot(simsea, sitel))
distplot <- distance_plot(output$results)

# Arrange plots
plot_grid(datplot, smap, distplot, nrow = 1)
