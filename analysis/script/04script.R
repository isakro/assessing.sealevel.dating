library(here)
library(tidyverse)
library(sf)
library(ggridges)
library(ggthemes)
library(ggnewscale)
library(tmap)
library(tmaptools)
library(patchwork)
library(terra)
library(raster)
library(oxcAAR)
library(vapour)
library(topoDistance)
quickSetupOxcal()

library(cowplot)
library(IRanges)

# To install IRanges:
# source("http://bioconductor.org/biocLite.R")
# biocLite("IRanges")


# Load required functions and prepared data
source(here("analysis/script/03script.R"))
load(here("analysis/data/derived_data/01data.RData"))
load(here("analysis/data/derived_data/02data.RData"))

# Load background map
bmap <- st_read(here('analysis/data/raw_data/naturalearth_countries.gpkg'))

# Specify path to raster data (currently stored locally)
dtmpath <- "/home/isak/phd/eaa_presentation/dtm1"

# Assign closest isobases to site data.
sites_sa <- st_join(st_make_valid(sites_sa), isopolys,
                    join = st_intersects, largest = TRUE)

# Both the digital terrain model and the radiocarbon dates required
# manual inspection in case of inconsistencies (i.e. a highway running by the
# site), defining necessary size of the window of analysis, and/or to decide on
# grouping of dates. Each site was therefore first simulated 100 times and
# then rerun at 1000 when these issues were handled.


######### Kvastad A2 #########
sitename <- "Kvastad A2"
date_groups <- c(1, 1, 1, 2, 2)
# output <- apply_functions(sitename, date_groups, dtmpath, displacement_curves,
#                          isobases, nsamp = 1000, loc_bbox = 500)
load(here("analysis/data/derived_data/kvastada2.RData"))

plot_results(sitename, output$sitel, output$datedat, output$sitearea, bmap,
             sites_sa, isobases, output, date_groups, scale_dist = 150)

ggsave(file = here("analysis/figures/kvastad_a2.4.png"), width = 250,
       height = 228, units = "mm")

######### Langangen Vestgård 1 #########
sitename <- "Langangen Vestgård 1"
date_groups <- rep(1, 13)
output <- apply_functions(sitename, date_groups, dtmpath, displacement_curves,
                          isobases, nsamp = 1000, loc_bbox = 250)

save(output,
     file = here::here("analysis/data/derived_data/langangenv1.RData"))

plot_results(sitename, output$sitel, output$datedat, output$sitearea, bmap,
             sites_sa, isobases, output, date_groups, scale_dist = 100)

ggsave(file = here("analysis/figures/langangenv1.png"), width = 250,
       height = 152, units = "mm")

######### Vallermyrene 1a #########
sitename <- "Vallermyrene 1a"
date_groups <- c(1, 1)
output <- apply_functions(sitename, date_groups, dtmpath, displacement_curves,
                          isobases, nsamp = 1000, loc_bbox = 250)

save(output,
     file = here::here("analysis/data/derived_data/vallermyrene1a.RData"))

plot_results(sitename, output$sitel, output$datedat, output$sitearea, bmap,
             sites_sa, isobases, output, date_groups, scale_dist = 100)

ggsave(file = here("analysis/figures/vallermyrene1a.png"), width = 250,
       height = 152, units = "mm")

######### Vallermyrene 1b #########
sitename <- "Vallermyrene 1b"
date_groups <- 1
output <- apply_functions(sitename, date_groups, dtmpath, displacement_curves,
                          isobases, nsamp = 1000, loc_bbox = 250)

save(output,
     file = here::here("analysis/data/derived_data/vallermyrene1b.RData"))

plot_results(sitename, output$sitel, output$datedat, output$sitearea, bmap,
             sites_sa, isobases, output, date_groups, scale_dist = 100)

ggsave(file = here("analysis/figures/vallermyrene1b.png"), width = 250,
       height = 152, units = "mm")

