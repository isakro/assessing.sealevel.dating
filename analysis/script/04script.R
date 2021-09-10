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

# Apply functions (see 03script.R)
output <- apply_functions(sitename, dtmpath, displacement_curves, isobases,
                          nsamp = 1000)

# Create plots
datplot <- plot_dates(output$datedat, sitename)
smap <- tmap_grob(shore_plot(output$simsea, output$sitel))
distplot <- distance_plot(output$results)

# Arrange plots
plot_grid(datplot, smap, distplot, nrow = 1)
