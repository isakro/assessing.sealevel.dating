library(here)
library(tidyverse)
library(sf)
library(ggridges)
library(ggthemes)
library(ggnewscale)
library(tmap)
library(tmaptools)
library(cowplot)
library(terra)
library(raster)
library(oxcAAR)
library(vapour)
library(topoDistance)
library(IRanges)
quickSetupOxcal()

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

# Assign closest isobases to site data (largest = TRUE makes sure that
# sites overlapping isobase polygons are only ascribed the one they overlap the
# most, so as not to create duplicate records).
sites_sa <- st_join(st_make_valid(sites_sa), isopolys,
                    join = st_intersects, largest = TRUE)

# Example site
sitename <- "Kvastad A2"

# Apply functions (see 03script.R)
output <- apply_functions(sitename, c(1,1,1,2,2), dtmpath, displacement_curves, isobases,
                          nsamp = 1000)

# Create overview plots
datplot <- plot_dates(output$datedat, sitename)
pmap <- site_plot(sitearea, output$sitel, 150)
omap <- overview_plot(bmap, output$sitel, sites_sa, isobases)

# Create simulation plots (one set per group)
datplot1 <- plot_dates(output$datedat, sitename, multigroup = FALSE,
                       group = 1, title = FALSE)
smap1 <- shore_plot(output[[1]]$simsea, output$sitel, 150)
distplot1 <- distance_plot(output[[1]]$results)

datplot2 <- plot_dates(output$datedat, sitename, multigroup = FALSE,
                       group = 2,  title = FALSE)
smap2 <- shore_plot(output[[2]]$simsea, output$sitel, 150)
distplot2 <- distance_plot(output[[2]]$results)

pg1 <- plot_grid(datplot, datplot1, datplot2, nrow = 3)
pg2 <- plot_grid(pmap, smap1, smap2, nrow = 3)
pg3 <- plot_grid(omap, distplot2, distplot2, nrow = 3)
plot_grid(pg1, pg2, pg3, nrow = 1)

(datplot | pmap | omap) /
(datplot1 | smap1 | distplot1) /
(datplot2 | smap2 | distplot2) +
plot_layout(widths = c(1, 1.5, 1))

(p1 | p2 | p3) /
  p4


pg1 <- plot_grid(datplot, pmap, omap, nrow = 1)
pg2 <- plot_grid(datplot1, smap1, distplot1, nrow = 1)
pg3 <- plot_grid(datplot2, smap2, distplot2, nrow = 1)
plot_grid(pg1, pg2, pg3, nrow = 3)

ggsave(file = here("analysis/figures/kvastad_a2.2.png"), width = 450, height = 350, units = "mm")
ggsave(file = here("analysis/figures/kvastad_a2.3.png"), width = 250, height = 230, units = "mm")

save(output,
     file = here("analysis/data/derived_data/kvastada2.RData"))

# Example site
sitename <- "Langangen Vestgård 1"

# Apply functions (see 03script.R) - likely takes a few hours to execute
output <- apply_functions(sitename, dtmpath, displacement_curves, isobases,
                          nsamp = 1000)

# Create plots
datplot <- plot_dates(output$datedat, sitename)
smap <- tmap_grob(shore_plot(output[[1]]$simsea, output$sitel))
distplot <- distance_plot(output[[1]]$results)

# Arrange and call to plot
plot_grid(datplot, smap, distplot, nrow = 1)

save(output,
     file = here::here("analysis/data/derived_data/lv1.RData"))
