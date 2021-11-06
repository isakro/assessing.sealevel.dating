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

# Not currently in use
#library(cowplot)
#library(IRanges)
# To install IRanges:
# source("http://bioconductor.org/biocLite.R")
# biocLite("IRanges")


# Load required functions and prepared data
source(here("analysis/script/04script.R"))
load(here("analysis/data/derived_data/01data.RData"))
load(here("analysis/data/derived_data/02data.RData"))

# Load background map
bmap <- st_read(here("analysis/data/raw_data/naturalearth_countries.gpkg"))

# # Specify path to raster data (currently stored locally)
# dtmpath <- "/home/isak/phd/eaa_presentation/dtm10/dtm10_32.tif"

# Read in raster
dtm <- rast("/home/isak/phd/eaa_presentation/dtm10/dtm10.tif")

# Path to smaller site area rasters
siterpath <- "/home/isak/phd/eaa_presentation/sitearea"


# Assign closest isobases to site data.
sites_sa <- st_join(st_make_valid(sites_sa), isopolys,
                    join = st_intersects, largest = TRUE)

# Prespecified height of plot (to be multiplied by number of phases)
plot_height <- 76

# Both the digital terrain model and the radiocarbon dates required
# manual inspection in case of inconsistencies (i.e. a highway running by the
# site), defining necessary size of the window of analysis, and/or to decide on
# grouping of dates. Each site was therefore first simulated 50 times and
# then rerun at 1000 when these issues were handled.

######### Anvik #########
sitename <- "Anvik"
date_groups <- rep(1, 4)

# output <- apply_functions(sitename, date_groups, dtm, displacement_curves,
#                           isobases, nsamp = 1000, loc_bbox = 400, siterpath)
# save(output,
#      file = here::here("analysis/data/derived_data/anvik.RData"))
load(here("analysis/data/derived_data/anvik.RData"))

sitearea <- rast(file.path(siterpath,
                           paste0(str_replace(sitename, " ", "_"), ".tif")))

plot_results(sitename, output$sitel, output$datedat, sitearea, bmap,
             sites_sa, isobases, output, date_groups, scale_dist = 150,
             s_tdist = 1, s_xpos = 150, s_ypos = 70, s_bheight = 0.75)

ggsave(file = here("analysis/figures/anvik.png"), width = 250,
       height = plot_height + (plot_height * length(unique(date_groups))),
       units = "mm")

######### Bakke #########
sitename <- "Bakke"
date_groups <- 1

# output <- apply_functions(sitename, date_groups, dtm, displacement_curves,
#                           isobases, nsamp = 1000, loc_bbox = 3500, siterpath)
# save(output,
#      file = here::here("analysis/data/derived_data/bakke.RData"))
load(here("analysis/data/derived_data/bakke.RData"))

sitearea <- rast(file.path(siterpath,
                           paste0(str_replace(sitename, " ", "_"), ".tif")))

plot_results(sitename, output$sitel, output$datedat, sitearea, bmap,
             sites_sa, isobases, output, date_groups, scale_dist = 1250,
             s_tdist = 3, s_xpos = 1250, s_ypos = 500,  s_bheight = 1.5)

ggsave(file = here("analysis/figures/bakke.png"), width = 250,
       height = plot_height + (plot_height * length(unique(date_groups))),
       units = "mm")

######### Dørdal #########
sitename <- "Dørdal"
date_groups <- rep(1, 4)

# output <- apply_functions(sitename, date_groups, dtm, displacement_curves,
#                           isobases, nsamp = 1000, loc_bbox = 2000, siterpath)
# save(output,
#      file = here::here("analysis/data/derived_data/dordal.RData"))
load(here("analysis/data/derived_data/dordal.RData"))

sitearea <- rast(file.path(siterpath,
                           paste0(str_replace(sitename, " ", "_"), ".tif")))

plot_results(sitename, output$sitel, output$datedat, sitearea, bmap,
             sites_sa, isobases, output, date_groups, scale_dist = 500,
             s_tdist = 5, s_xpos = 600, s_ypos = 200,  s_bheight = 3)

ggsave(file = here("analysis/figures/dordal.png"), width = 250,
       height = plot_height + (plot_height * length(unique(date_groups))),
       units = "mm")

######### Dybdalshei 1 #########
sitename <- "Dybdalshei 1"
date_groups <- c(1, 1, 1, 1, 2, 2, 3, 3)

# output <- apply_functions(sitename, date_groups, dtm, displacement_curves,
#                           isobases, nsamp = 1000, loc_bbox = 400, siterpath)
# save(output,
#      file = here::here("analysis/data/derived_data/dybdalshei1.RData"))
load(here("analysis/data/derived_data/dybdalshei1.RData"))

sitearea <- rast(file.path(siterpath,
                           paste0(str_replace(sitename, " ", "_"), ".tif")))

plot_results(sitename, output$sitel, output$datedat, sitearea, bmap,
             sites_sa, isobases, output, date_groups, scale_dist = 150,
             s_tdist = 0.5, s_xpos = 135, s_ypos = 65,  s_bheight = 0.25)

ggsave(file = here("analysis/figures/dybdalshei1.png"), width = 250,
       height = plot_height + (plot_height * length(unique(date_groups))),
       units = "mm")

######### Dybdalshei 2 #########
sitename <- "Dybdalshei 2"
date_groups <- 1

# output <- apply_functions(sitename, date_groups, dtm, displacement_curves,
#                            isobases, nsamp = 1000, loc_bbox = 400, siterpath)
# save(output,
#       file = here::here("analysis/data/derived_data/dybdalshei2.RData"))
load(here("analysis/data/derived_data/dybdalshei2.RData"))

sitearea <- rast(file.path(siterpath,
                           paste0(str_replace(sitename, " ", "_"), ".tif")))

plot_results(sitename, output$sitel, output$datedat, sitearea, bmap,
             sites_sa, isobases, output, date_groups, scale_dist = 150,
             s_tdist = 2, s_xpos = 150, s_ypos = 60, s_bheight = 1.5)

ggsave(file = here("analysis/figures/dybdalshei2.png"), width = 250,
       height = plot_height + (plot_height * length(unique(date_groups))),
       units = "mm")

######### Gunnarsrød 4. Edited raster. Two single dates #########
sitename <- "Gunnarsrød 4"
date_groups <- c(2, 1)

# output <- apply_functions(sitename, date_groups, dtm, displacement_curves,
#                           isobases, nsamp = 1000, loc_bbox = 400, siterpath)
# save(output,
#      file = here::here("analysis/data/derived_data/gunnarsrod4.RData"))
load(here("analysis/data/derived_data/gunnarsrod4.RData"))

sitearea <- rast(file.path(siterpath,
                           paste0(str_replace(sitename, " ", "_"), ".tif")))

plot_results(sitename, output$sitel, output$datedat, sitearea, bmap,
             sites_sa, isobases, output, date_groups, scale_dist = 150,
             s_tdist = 9, s_xpos = 150, s_ypos = 70, s_bheight = 5)

ggsave(file = here("analysis/figures/gunnarsrod4.png"), width = 250,
       height = plot_height + (plot_height * length(unique(date_groups))),
       units = "mm")

######### Gunnarsrød 5 Edited raster #########
sitename <- "Gunnarsrød 5"
date_groups <- c(1, 1, 2, 3)

# output <- apply_functions(sitename, date_groups, dtm, displacement_curves,
#                           isobases, nsamp = 1000, loc_bbox = 400, siterpath)
# save(output,
#      file = here::here("analysis/data/derived_data/gunnarsrod5.RData"))
load(here("analysis/data/derived_data/gunnarsrod5.RData"))

sitearea <- rast(file.path(siterpath,
                           paste0(str_replace(sitename, " ", "_"), ".tif")))

plot_results(sitename, output$sitel, output$datedat, sitearea, bmap,
             sites_sa, isobases, output, date_groups, scale_dist = 150,
             s_tdist = 0.5, s_xpos = 135, s_ypos = 65,  s_bheight = 0.25)

ggsave(file = here("analysis/figures/gunnarsrod5.png"), width = 250,
       height = plot_height + (plot_height * length(unique(date_groups))),
       units = "mm")

######### Gunnarsrød 6a #########
sitename <- "Gunnarsrød 6a"
date_groups <- 1

# output <- apply_functions(sitename, date_groups, dtm, displacement_curves,
#                           isobases, nsamp = 1000, loc_bbox = 400, siterpath)
# save(output,
#      file = here::here("analysis/data/derived_data/gunnarsrod6a.RData"))
load(here("analysis/data/derived_data/gunnarsrod6a.RData"))

sitearea <- rast(file.path(siterpath,
                           paste0(str_replace(sitename, " ", "_"), ".tif")))

plot_results(sitename, output$sitel, output$datedat, sitearea, bmap,
             sites_sa, isobases, output, date_groups, scale_dist = 150,
             s_tdist = 1, s_xpos = 135, s_ypos = 65,  s_bheight = 0.5)

ggsave(file = here("analysis/figures/gunnarsrod6a.png"), width = 250,
       height = plot_height + (plot_height * length(unique(date_groups))),
       units = "mm")

######### Gunnarsrød 6b #########
sitename <- "Gunnarsrød 6b"
date_groups <- 1

# output <- apply_functions(sitename, date_groups, dtm, displacement_curves,
#                           isobases, nsamp = 1000, loc_bbox = 400, siterpath)
# save(output,
#      file = here::here("analysis/data/derived_data/gunnarsrod6b.RData"))
load(here("analysis/data/derived_data/gunnarsrod6b.RData"))

sitearea <- rast(file.path(siterpath,
                           paste0(str_replace(sitename, " ", "_"), ".tif")))

plot_results(sitename, output$sitel, output$datedat, sitearea, bmap,
             sites_sa, isobases, output, date_groups, scale_dist = 150,
             s_tdist = 2, s_xpos = 135, s_ypos = 65,  s_bheight = 1)

ggsave(file = here("analysis/figures/gunnarsrod6b.png"), width = 250,
       height = plot_height + (plot_height * length(unique(date_groups))),
       units = "mm")

######### Gunnarsrød 7 #########
sitename <- "Gunnarsrød 7"
date_groups <- 1

# output <- apply_functions(sitename, date_groups, dtm, displacement_curves,
#                           isobases, nsamp = 1000, loc_bbox = 400, siterpath)
# save(output,
#      file = here::here("analysis/data/derived_data/gunnarsrod7.RData"))
load(here("analysis/data/derived_data/gunnarsrod7.RData"))

sitearea <- rast(file.path(siterpath,
                           paste0(str_replace(sitename, " ", "_"), ".tif")))

plot_results(sitename, output$sitel, output$datedat, sitearea, bmap,
             sites_sa, isobases, output, date_groups, scale_dist = 150,
             s_tdist = 1, s_xpos = 135, s_ypos = 65,  s_bheight = 0.5)

ggsave(file = here("analysis/figures/gunnarsrod7.png"), width = 250,
       height = plot_height + (plot_height * length(unique(date_groups))),
       units = "mm")

######### Gunnarsrød 10 #########
sitename <- "Gunnarsrød 10"
date_groups <- 1

# output <- apply_functions(sitename, date_groups, dtm, displacement_curves,
#                           isobases, nsamp = 1000, loc_bbox = 400, siterpath)
# save(output,
#      file = here::here("analysis/data/derived_data/gunnarsrod10.RData"))
load(here("analysis/data/derived_data/gunnarsrod10.RData"))

sitearea <- rast(file.path(siterpath,
                           paste0(str_replace(sitename, " ", "_"), ".tif")))

plot_results(sitename, output$sitel, output$datedat, sitearea, bmap,
             sites_sa, isobases, output, date_groups, scale_dist = 150,
             s_tdist = 3.5, s_xpos = 135, s_ypos = 65,  s_bheight = 2)

ggsave(file = here("analysis/figures/gunnarsrod10.png"), width = 250,
       height = plot_height + (plot_height * length(unique(date_groups))),
       units = "mm")

######### Hegna vest 1 #########
sitename <- "Hegna vest 1"
date_groups <- c(1, 1, 2)

# output <- apply_functions(sitename, date_groups, dtm, displacement_curves,
#                           isobases, nsamp = 1000, loc_bbox = 400, siterpath)
# save(output,
#      file = here::here("analysis/data/derived_data/hegnavest1.RData"))
load(here("analysis/data/derived_data/hegnavest1.RData"))

sitearea <- rast(file.path(siterpath,
                           paste0(str_replace(sitename, " ", "_"), ".tif")))

plot_results(sitename, output$sitel, output$datedat, sitearea, bmap,
             sites_sa, isobases, output, date_groups, scale_dist = 150,
             s_tdist = 0.7, s_xpos = 135, s_ypos = 65,  s_bheight = 0.5)

ggsave(file = here("analysis/figures/hegnavest1.png"), width = 250,
       height = plot_height + (plot_height * length(unique(date_groups))),
       units = "mm")

######### Hegna vest 2 #########
sitename <- "Hegna vest 2"
date_groups <- c(1, 2, 3, 3)

# output <- apply_functions(sitename, date_groups, dtm, displacement_curves,
#                           isobases, nsamp = 1000, loc_bbox = 1200, siterpath)
# save(output,
#      file = here::here("analysis/data/derived_data/hegnavest2.RData"))
load(here("analysis/data/derived_data/hegnavest2.RData"))

sitearea <- rast(file.path(siterpath,
                           paste0(str_replace(sitename, " ", "_"), ".tif")))

plot_results(sitename, output$sitel, output$datedat, sitearea, bmap,
             sites_sa, isobases, output, date_groups, scale_dist = 600,
             s_tdist = 1.25, s_xpos = 500, s_ypos = 250,  s_bheight = 0.8)

ggsave(file = here("analysis/figures/hegnavest2.png"), width = 250,
       height = plot_height + (plot_height * length(unique(date_groups))),
       units = "mm")

######### Hegna vest 3 #########
sitename <- "Hegna vest 3"
date_groups <- 1

# output <- apply_functions(sitename, date_groups, dtm, displacement_curves,
#                           isobases, nsamp = 1000, loc_bbox = 400, siterpath)
# save(output,
#      file = here::here("analysis/data/derived_data/hegnavest3.RData"))
load(here("analysis/data/derived_data/hegnavest3.RData"))

sitearea <- rast(file.path(siterpath,
                           paste0(str_replace(sitename, " ", "_"), ".tif")))

plot_results(sitename, output$sitel, output$datedat, sitearea, bmap,
             sites_sa, isobases, output, date_groups, scale_dist = 150,
             s_tdist = 1, s_xpos = 135, s_ypos = 65,  s_bheight = 0.5)

ggsave(file = here("analysis/figures/hegnavest3.png"), width = 250,
       height = plot_height + (plot_height * length(unique(date_groups))),
       units = "mm")

######### Hesthag C2 #########
sitename <- "Hesthag C2"
date_groups <- c(1, 2)

# output <- apply_functions(sitename, date_groups, dtm, displacement_curves,
#                           isobases, nsamp = 1000, loc_bbox = 400, siterpath)
# save(output,
#      file = here::here("analysis/data/derived_data/hesthagc2.RData"))
load(here("analysis/data/derived_data/hesthagc2.RData"))

sitearea <- rast(file.path(siterpath,
                           paste0(str_replace(sitename, " ", "_"), ".tif")))

plot_results(sitename, output$sitel, output$datedat, sitearea, bmap,
             sites_sa, isobases, output, date_groups, scale_dist = 150,
             s_tdist = 0.5, s_xpos = 135, s_ypos = 65,  s_bheight = 0.25)

ggsave(file = here("analysis/figures/hesthagc2.png"), width = 250,
       height = plot_height + (plot_height * length(unique(date_groups))),
       units = "mm")

######### Hesthag C4 #########
sitename <- "Hesthag C4"
date_groups <- c(1, 2)

# output <- apply_functions(sitename, date_groups, dtm, displacement_curves,
#                           isobases, nsamp = 1000, loc_bbox = 400, siterpath)
# save(output,
#      file = here::here("analysis/data/derived_data/hesthagc4.RData"))
load(here("analysis/data/derived_data/hesthagc4.RData"))

sitearea <- rast(file.path(siterpath,
                           paste0(str_replace(sitename, " ", "_"), ".tif")))

plot_results(sitename, output$sitel, output$datedat, sitearea, bmap,
             sites_sa, isobases, output, date_groups, scale_dist = 150,
             s_tdist = 1, s_xpos = 135, s_ypos = 65,  s_bheight = 0.5)

ggsave(file = here("analysis/figures/hesthagc4.png"), width = 250,
       height = plot_height + (plot_height * length(unique(date_groups))),
       units = "mm")

######### Hovland 1 #########
sitename <- "Hovland 1"
date_groups <- c(1, 1, 1, 2)

# output <- apply_functions(sitename, date_groups, dtm, displacement_curves,
#                           isobases, nsamp = 1000, loc_bbox = 1000, siterpath)
# save(output,
#      file = here::here("analysis/data/derived_data/hovland1.RData"))
load(here("analysis/data/derived_data/hovland1.RData"))

sitearea <- rast(file.path(siterpath,
                           paste0(str_replace(sitename, " ", "_"), ".tif")))

plot_results(sitename, output$sitel, output$datedat, sitearea, bmap,
             sites_sa, isobases, output, date_groups, scale_dist = 250,
             s_tdist = 3, s_xpos = 350, s_ypos = 150,  s_bheight = 2.25)

ggsave(file = here("analysis/figures/hovland1.png"), width = 250,
       height = plot_height + (plot_height * length(unique(date_groups))),
       units = "mm")

######### Hovland 3 #########
sitename <- "Hovland 3"
date_groups <- rep(1, 18)

# output <- apply_functions(sitename, date_groups, dtm, displacement_curves,
#                           isobases, nsamp = 1000, loc_bbox = 400, siterpath)
# save(output,
#      file = here::here("analysis/data/derived_data/hovland3.RData"))
load(here("analysis/data/derived_data/hovland3.RData"))

sitearea <- rast(file.path(siterpath,
                           paste0(str_replace(sitename, " ", "_"), ".tif")))

plot_results(sitename, output$sitel, output$datedat, sitearea, bmap,
             sites_sa, isobases, output, date_groups, scale_dist = 150,
             s_tdist = 1, s_xpos = 135, s_ypos = 65,  s_bheight = 0.5)

ggsave(file = here("analysis/figures/hovland3.png"), width = 250,
       height = plot_height + (plot_height * length(unique(date_groups))),
                               units = "mm")

######### Hovland 4 #########
sitename <- "Hovland 4"
date_groups <- c(rep(1, 4), 2)

# output <- apply_functions(sitename, date_groups, dtm, displacement_curves,
#                           isobases, nsamp = 1000, loc_bbox = 1200, siterpath)
# save(output,
#      file = here::here("analysis/data/derived_data/hovland4.RData"))
load(here("analysis/data/derived_data/hovland4.RData"))

sitearea <- rast(file.path(siterpath,
                           paste0(str_replace(sitename, " ", "_"), ".tif")))

plot_results(sitename, output$sitel, output$datedat, sitearea, bmap,
             sites_sa, isobases, output, date_groups, scale_dist = 400,
             s_tdist = 2.5, s_xpos = 500, s_ypos = 200,  s_bheight = 2)

ggsave(file = here("analysis/figures/hovland4.png"), width = 250,
       height = plot_height + (plot_height * length(unique(date_groups))),
       units = "mm")

######### Hovland 5 #########
sitename <- "Hovland 5"
date_groups <- 1

# output <- apply_functions(sitename, date_groups, dtm, displacement_curves,
#                           isobases, nsamp = 1000, loc_bbox = 1000, siterpath)
# save(output,
#      file = here::here("analysis/data/derived_data/hovland5.RData"))
load(here("analysis/data/derived_data/hovland5.RData"))

sitearea <- rast(file.path(siterpath,
                           paste0(str_replace(sitename, " ", "_"), ".tif")))

plot_results(sitename, output$sitel, output$datedat, sitearea, bmap,
             sites_sa, isobases, output, date_groups, scale_dist = 300,
             s_tdist = 2.5, s_xpos = 500, s_ypos = 200,  s_bheight = 2)

ggsave(file = here("analysis/figures/hovland5.png"), width = 250,
       height = plot_height + (plot_height * length(unique(date_groups))),
       units = "mm")

######### Hydal 4 #########
sitename <- "Hydal 4"
date_groups <- c(1,1)

# output <- apply_functions(sitename, date_groups, dtm, displacement_curves,
#                           isobases, nsamp = 1000, loc_bbox = 400, siterpath)
# save(output,
#      file = here::here("analysis/data/derived_data/hydal4.RData"))
load(here("analysis/data/derived_data/hydal4.RData"))

sitearea <- rast(file.path(siterpath,
                           paste0(str_replace(sitename, " ", "_"), ".tif")))

plot_results(sitename, output$sitel, output$datedat, sitearea, bmap,
             sites_sa, isobases, output, date_groups, scale_dist = 150,
             s_tdist = 1.5, s_xpos = 135, s_ypos = 65,  s_bheight = 1)

ggsave(file = here("analysis/figures/hydal4.png"), width = 250,
       height = plot_height + (plot_height * length(unique(date_groups))),
                               units = "mm")

######### Krøgenes D1 #########
sitename <- "Krøgenes D1"
date_groups <- c(1, rep(2, 4), 3, 4, 5)

output <- apply_functions(sitename, date_groups, dtm, displacement_curves,
                          isobases, nsamp = 1000, loc_bbox = 400, siterpath)
save(output,
     file = here::here("analysis/data/derived_data/krogenesd1.RData"))

sitearea <- rast(file.path(siterpath,
                           paste0(str_replace(sitename, " ", "_"), ".tif")))

plot_results(sitename, output$sitel, output$datedat, sitearea, bmap,
             sites_sa, isobases, output, date_groups, scale_dist = 150)

ggsave(file = here("analysis/figures/krogenesd1.png"), width = 250,
       height = plot_height + (plot_height * length(unique(date_groups))),
       units = "mm")

######### Krøgenes D2 #########
sitename <- "Krøgenes D2"
date_groups <- c(1, 2, 2, 2)

output <- apply_functions(sitename, date_groups, dtm, displacement_curves,
                          isobases, nsamp = 1000, loc_bbox = 400, siterpath)
save(output,
     file = here::here("analysis/data/derived_data/krogenesd2.RData"))

sitearea <- rast(file.path(siterpath,
                           paste0(str_replace(sitename, " ", "_"), ".tif")))

plot_results(sitename, output$sitel, output$datedat, sitearea, bmap,
             sites_sa, isobases, output, date_groups, scale_dist = 150,
             s_tdist = 1.5, s_xpos = 135, s_ypos = 65,  s_bheight = 1)

ggsave(file = here("analysis/figures/krogenesd2.png"), width = 250,
       height = plot_height + (plot_height * length(unique(date_groups))),
       units = "mm")

######### Kvastad A1 #########
sitename <- "Kvastad A1"
date_groups <- 1

output <- apply_functions(sitename, date_groups, dtm, displacement_curves,
                          isobases, nsamp = 1000, loc_bbox = 400, siterpath)

save(output,
     file = here::here("analysis/data/derived_data/kvastada1.RData"))

sitearea <- rast(file.path(siterpath,
                           paste0(str_replace(sitename, " ", "_"), ".tif")))

plot_results(sitename, output$sitel, output$datedat, sitearea, bmap,
             sites_sa, isobases, output, date_groups, scale_dist = 150,
             s_tdist = 0.5, s_xpos = 135, s_ypos = 65,  s_bheight = 0.25)

ggsave(file = here("analysis/figures/kvastada1.png"), width = 250,
       height = plot_height + (plot_height * length(unique(date_groups))),
       units = "mm")

######### Kvastad A2 Both areas collapsed #########
sitename <- "Kvastad A2"
date_groups <- c(1, 1, 1, 2, 2)

output <- apply_functions(sitename, date_groups, dtm, displacement_curves,
                          isobases, nsamp = 1000, loc_bbox = 500, siterpath)
# load(here("analysis/data/derived_data/kvastada2.RData"))

save(output,
     file = here::here("analysis/data/derived_data/kvastada2.RData"))

sitearea <- rast(file.path(siterpath,
                           paste0(str_replace(sitename, " ", "_"), ".tif")))

plot_results(sitename, output$sitel, output$datedat, sitearea, bmap,
             sites_sa, isobases, output, date_groups, scale_dist = 150)

ggsave(file = here("analysis/figures/kvastada2.png"), width = 250,
       height = plot_height + (plot_height * length(unique(date_groups))),
       units = "mm")

######### Kvastad A4 #########
sitename <- "Kvastad A4"
date_groups <- 1

# output <- apply_functions(sitename, date_groups, dtm, displacement_curves,
#                           isobases, nsamp = 1000, loc_bbox = 400, siterpath)
# save(output,
#      file = here::here("analysis/data/derived_data/kvastada4.RData"))
load(here("analysis/data/derived_data/kvastada4.RData"))

sitearea <- rast(file.path(siterpath,
                           paste0(str_replace(sitename, " ", "_"), ".tif")))

plot_results(sitename, output$sitel, output$datedat, sitearea, bmap,
             sites_sa, isobases, output, date_groups, scale_dist = 150,
             s_tdist = 1, s_xpos = 135, s_ypos = 65,  s_bheight = 0.5)

ggsave(file = here("analysis/figures/kvastada4.png"), width = 250,
       height = plot_height + (plot_height * length(unique(date_groups))),
                               units = "mm")

######### Kvastad A9 #########
sitename <- "Kvastad A9"
date_groups <- 1

# output <- apply_functions(sitename, date_groups, dtm, displacement_curves,
#                            isobases, nsamp = 1000, loc_bbox = 400, siterpath)
# save(output,
#      file = here::here("analysis/data/derived_data/kvastada9.RData"))
load(here("analysis/data/derived_data//kvastada9.RData"))

sitearea <- rast(file.path(siterpath,
                           paste0(str_replace(sitename, " ", "_"), ".tif")))

plot_results(sitename, output$sitel, output$datedat, sitearea, bmap,
             sites_sa, isobases, output, date_groups, scale_dist = 150,
             s_tdist = 5, s_xpos = 140, s_ypos = 60,  s_bheight = 3.5)

ggsave(file = here("analysis/figures/kvastada9.png"), width = 250,
       height = plot_height + (plot_height * length(unique(date_groups))),
       units = "mm")

######### Langangen Vestgård 1 #########
sitename <- "Langangen Vestgård 1"
date_groups <- rep(1, 13)

# output <- apply_functions(sitename, date_groups, dtm, displacement_curves,
#                           isobases, nsamp = 1000, loc_bbox = 400, siterpath)
# save(output,
#      file = here::here("analysis/data/derived_data/langangenv1.RData"))
load(here("analysis/data/derived_data/langangenv1.RData"))

sitearea <- rast(file.path(siterpath,
                           paste0(str_replace(sitename, " ", "_"), ".tif")))

plot_results(sitename, output$sitel, output$datedat, sitearea, bmap,
             sites_sa, isobases, output, date_groups, scale_dist = 150)

ggsave(file = here("analysis/figures/langangenv1.png"), width = 250,
       height = plot_height + (plot_height * length(unique(date_groups))),
       units = "mm")

######### Langangen Vestgård 3 #########
sitename <- "Langangen Vestgård 3"
date_groups <- c(2, 2, 1, 2, 2, 2)

# output <- apply_functions(sitename, date_groups, dtm, displacement_curves,
#                            isobases, nsamp = 1000, loc_bbox = 400, siterpath)
# save(output,
#      file = here::here("analysis/data/derived_data/langangenv3.RData"))
load(here("analysis/data/derived_data/langangenv3.RData"))

sitearea <- rast(file.path(siterpath,
                           paste0(str_replace(sitename, " ", "_"), ".tif")))

plot_results(sitename, output$sitel, output$datedat, sitearea, bmap,
             sites_sa, isobases, output, date_groups, scale_dist = 150,
             s_tdist = 1.5, s_xpos = 150, s_ypos = 80,  s_bheight = 1.25)

ggsave(file = here("analysis/figures/langangenv3.png"), width = 250,
       height = plot_height + (plot_height * length(unique(date_groups))),
       units = "mm")

######### Langangen Vestgård 5 #########
sitename <- "Langangen Vestgård 5"
date_groups <- c(1, 1, 2, 2, 2)

output <- apply_functions(sitename, date_groups, dtm, displacement_curves,
                          isobases, nsamp = 1000, loc_bbox = 400, siterpath)
# load(here("analysis/data/derived_data/langangenv5.RData"))

save(output,
     file = here::here("analysis/data/derived_data/langangenv5.RData"))

sitearea <- rast(file.path(siterpath,
                           paste0(str_replace(sitename, " ", "_"), ".tif")))

plot_results(sitename, output$sitel, output$datedat, sitearea, bmap,
             sites_sa, isobases, output, date_groups, scale_dist = 150)

ggsave(file = here("analysis/figures/langangenv5.png"), width = 250,
       height = plot_height + (plot_height * length(unique(date_groups))),
       units = "mm")

######### Langangen Vestgård 6 #########
sitename <- "Langangen Vestgård 6"
date_groups <- c(rep(1, 9), 2)

# output <- apply_functions(sitename, date_groups, dtm, displacement_curves,
#                           isobases, nsamp = 1000, loc_bbox = 400, siterpath)
# save(output,
#      file = here::here("analysis/data/derived_data/langangenv6.RData"))
load(here("analysis/data/derived_data/langangenv6.RData"))

sitearea <- rast(file.path(siterpath,
                           paste0(str_replace(sitename, " ", "_"), ".tif")))

plot_results(sitename, output$sitel, output$datedat, sitearea, bmap,
             sites_sa, isobases, output, date_groups, scale_dist = 150,
             s_tdist = 1, s_xpos = 170, s_ypos = 80,  s_bheight = 0.5)

ggsave(file = here("analysis/figures/langangenv6.png"), width = 250,
       height = plot_height + (plot_height * length(unique(date_groups))),
       units = "mm")

######### Langangen Vestgård 7 Road as site limit #########
sitename <- "Langangen Vestgård 7"
date_groups <- 1

# output <- apply_functions(sitename, date_groups, dtm, displacement_curves,
#                           isobases, nsamp = 1000, loc_bbox = 200, siterpath,
#                           sitelimit = FALSE)
# save(output,
#      file = here::here("analysis/data/derived_data/langangenv7.RData"))
load(here("analysis/data/derived_data/langangenv7.RData"))

sitearea <- rast(file.path(siterpath,
                           paste0(str_replace(sitename, " ", "_"), ".tif")))

plot_results(sitename, output$sitel, output$datedat, sitearea, bmap,
             sites_sa, isobases, output, date_groups, scale_dist = 100,
             s_tdist = 1, s_xpos = 100, s_ypos = 50,  s_bheight = 10)

ggsave(file = here("analysis/figures/langangenv7.png"), width = 250,
       height =  plot_height + (plot_height * length(unique(date_groups))),
       units = "mm")

######### Langemyr #########
sitename <- "Langemyr"
date_groups <- 1

# output <- apply_functions(sitename, date_groups, dtm, displacement_curves,
#                           isobases, nsamp = 1000, loc_bbox = 400, siterpath)
# save(output,
#      file = here::here("analysis/data/derived_data/langemyr.RData"))
load(here("analysis/data/derived_data/langemyr.RData"))

sitearea <- rast(file.path(siterpath,
                           paste0(str_replace(sitename, " ", "_"), ".tif")))

plot_results(sitename, output$sitel, output$datedat, sitearea, bmap,
             sites_sa, isobases, output, date_groups, scale_dist = 150)

ggsave(file = here("analysis/figures/langemyr.png"), width = 250,
       height = plot_height + (plot_height * length(unique(date_groups))),
       units = "mm")

######### Lunaveien #########
sitename <- "Lunaveien"
date_groups <- 1

output <- apply_functions(sitename, date_groups, dtm, displacement_curves,
                          isobases, nsamp = 1000, loc_bbox = 650, siterpath)
save(output,
     file = here::here("analysis/data/derived_data/lunaveien.RData"))
# load(here("analysis/data/derived_data/lunaveien.RData"))

sitearea <- rast(file.path(siterpath,
                           paste0(str_replace(sitename, " ", "_"), ".tif")))

plot_results(sitename, output$sitel, output$datedat, sitearea, bmap,
             scale_dist = 200,  s_tdist = 0.5, s_xpos = 135, s_ypos = 65,
             s_bheight = 0.25)

ggsave(file = here("analysis/figures/lunaveien.png"), width = 250,
       height = plot_height + (plot_height * length(unique(date_groups))),
       units = "mm")

######### Nauen A #########
sitename <- "Nauen A"
date_groups <- c(1, 1, 2)

# output <- apply_functions(sitename, date_groups, dtm, displacement_curves,
#                           isobases, nsamp = 1000, loc_bbox = 1000, siterpath)
# save(output,
#      file = here::here("analysis/data/derived_data/nauena.RData"))
load(here("analysis/data/derived_data/nauena.RData"))

sitearea <- rast(file.path(siterpath,
                           paste0(str_replace(sitename, " ", "_"), ".tif")))

plot_results(sitename, output$sitel, output$datedat, sitearea, bmap,
             sites_sa, isobases, output, date_groups, scale_dist = 400,
             s_tdist = 1, s_xpos = 400, s_ypos = 160,  s_bheight = 0.5)

ggsave(file = here("analysis/figures/nauena.png"), width = 250,
       height = plot_height + (plot_height * length(unique(date_groups))),
       units = "mm")

######### Nauen C #########
sitename <- "Nauen C"
date_groups <- c(1, 1, 1, 2)

# output <- apply_functions(sitename, date_groups, dtm, displacement_curves,
#                           isobases, nsamp = 1000, loc_bbox = 1000, siterpath)
# save(output,
#      file = here::here("analysis/data/derived_data/nauenc.RData"))
load(here("analysis/data/derived_data/nauenc.RData"))

sitearea <- rast(file.path(siterpath,
                           paste0(str_replace(sitename, " ", "_"), ".tif")))

plot_results(sitename, output$sitel, output$datedat, sitearea, bmap,
             sites_sa, isobases, output, date_groups, scale_dist = 400,
             s_tdist = 1.5, s_xpos = 400, s_ypos = 160,  s_bheight = 1)

ggsave(file = here("analysis/figures/nauenc.png"), width = 250,
       height = plot_height + (plot_height * length(unique(date_groups))),
       units = "mm")

######### Nordby 1 #########
sitename <- "Nordby 1"
date_groups <- rep(1, 5)

# output <- apply_functions(sitename, date_groups, dtm, displacement_curves,
#                           isobases, nsamp = 1000, loc_bbox = 400, siterpath)
#
# save(output,
#      file = here::here("analysis/data/derived_data/nordby1.RData"))
load(here("analysis/data/derived_data/nordby1.RData"))

sitearea <- rast(file.path(siterpath,
                           paste0(str_replace(sitename, " ", "_"), ".tif")))

plot_results(sitename, output$sitel, output$datedat, sitearea, bmap,
             sites_sa, isobases, output, date_groups, scale_dist = 150,
             s_tdist = 0.5, s_xpos = 110, s_ypos = 50,  s_bheight = 0.25)

ggsave(file = here("analysis/figures/nordby1.png"), width = 250,
       height = plot_height + (plot_height * length(unique(date_groups))),
       units = "mm")

######### Nordby 52 #########
sitename <- "Nordby 52"
date_groups <- c(rep(1, 7), 2)

# output <- apply_functions(sitename, date_groups, dtm, displacement_curves,
#                           isobases, nsamp = 1000, loc_bbox = 400, siterpath)
# save(output,
#      file = here::here("analysis/data/derived_data/nordby52.RData"))
load(here("analysis/data/derived_data/nordby52.RData"))

sitearea <- rast(file.path(siterpath,
                           paste0(str_replace(sitename, " ", "_"), ".tif")))

plot_results(sitename, output$sitel, output$datedat, sitearea, bmap,
             sites_sa, isobases, output, date_groups, scale_dist = 150,
             s_tdist = 2, s_xpos = 150, s_ypos = 80,  s_bheight = 1)

ggsave(file = here("analysis/figures/nordby52.png"), width = 250,
       height = plot_height + (plot_height * length(unique(date_groups))),
       units = "mm")

######### Pauler 1 #########
sitename <- "Pauler 1"
date_groups <- 1

# output <- apply_functions(sitename, date_groups, dtm, displacement_curves,
#                           isobases, nsamp = 1000, loc_bbox = 2500, siterpath)
# save(output,
#      file = here::here("analysis/data/derived_data/pauler1.RData"))
load(here("analysis/data/derived_data/pauler1.RData"))


sitearea <- rast(file.path(siterpath,
                           paste0(str_replace(sitename, " ", "_"), ".tif")))

plot_results(sitename, output$sitel, output$datedat, sitearea, bmap,
             sites_sa, isobases, output, date_groups, scale_dist = 1000,
             s_tdist = 2, s_xpos = 110, s_ypos = 50,  s_bheight = 1)

ggsave(file = here("analysis/figures/pauler1.png"), width = 250,
       height = plot_height + (plot_height * length(unique(date_groups))),
       units = "mm")

######### Pauler 2 #########
sitename <- "Pauler 2"
date_groups <- c(1, 1, 2)

output <- apply_functions(sitename, date_groups, dtm, displacement_curves,
                          isobases, nsamp = 1000, loc_bbox = 2500, siterpath)

save(output,
     file = here::here("analysis/data/derived_data/pauler2.RData"))

sitearea <- rast(file.path(siterpath,
                           paste0(str_replace(sitename, " ", "_"), ".tif")))

plot_results(sitename, output$sitel, output$datedat, sitearea, bmap,
             sites_sa, isobases, output, date_groups, scale_dist = 1000,
             s_tdist = 2, s_xpos = 110, s_ypos = 50,  s_bheight = 1)

ggsave(file = here("analysis/figures/pauler2.png"), width = 250,
       height = plot_height + (plot_height * length(unique(date_groups))),
       units = "mm")

######### Pjonkerød R1 #########
sitename <- "Pjonkerød R1"
date_groups <- c(2, 1)

output <- apply_functions(sitename, date_groups, dtm, displacement_curves,
                          isobases, nsamp = 1000, loc_bbox = 400, siterpath)

save(output,
     file = here::here("analysis/data/derived_data/pjonkerodr1.RData"))

sitearea <- rast(file.path(siterpath,
                           paste0(str_replace(sitename, " ", "_"), ".tif")))

plot_results(sitename, output$sitel, output$datedat, sitearea, bmap,
             sites_sa, isobases, output, date_groups, scale_dist = 150,
             s_tdist = 1, s_xpos = 110, s_ypos = 50,  s_bheight = 0.5)

ggsave(file = here("analysis/figures/pjonkerodr1.png"), width = 250,
       height = 152, units = "mm")

######### Prestemoen 1 #########
sitename <- "Prestemoen 1"
date_groups <- c(1, 1, 1)

output <- apply_functions(sitename, date_groups, dtm, displacement_curves,
                          isobases, nsamp = 1000, loc_bbox = 400, siterpath)

save(output,
     file = here::here("analysis/data/derived_data/prestemoen1.RData"))

sitearea <- rast(file.path(siterpath,
                           paste0(str_replace(sitename, " ", "_"), ".tif")))

plot_results(sitename, output$sitel, output$datedat, sitearea, bmap,
             sites_sa, isobases, output, date_groups, scale_dist = 150,
             s_tdist = 1, s_xpos = 110, s_ypos = 50,  s_bheight = 0.5)

ggsave(file = here("analysis/figures/prestemoen1.png"), width = 250,
       height = 152, units = "mm")

######### Ragnhildrød #########
sitename <- "Ragnhildrød"
date_groups <- c(1, 2)

output <- apply_functions(sitename, date_groups, dtm, displacement_curves,
                          isobases, nsamp = 1000, loc_bbox = 400, siterpath)

save(output,
     file = here::here("analysis/data/derived_data/ragnhildrod.RData"))

sitearea <- rast(file.path(siterpath,
                           paste0(str_replace(sitename, " ", "_"), ".tif")))

plot_results(sitename, output$sitel, output$datedat, sitearea, bmap,
             sites_sa, isobases, output, date_groups, scale_dist = 150,
             s_tdist = 1, s_xpos = 110, s_ypos = 50,  s_bheight = 0.5)

ggsave(file = here("analysis/figures/ragnhildrod.png"), width = 250,
       height = 152, units = "mm")

######### Rødbøl 54 #########
sitename <- "Rødbøl 54"
date_groups <- 1

output <- apply_functions(sitename, date_groups, dtm, displacement_curves,
                          isobases, nsamp = 1000, loc_bbox = 400, siterpath)

save(output,
     file = here::here("analysis/data/derived_data/rodbol54.RData"))

sitearea <- rast(file.path(siterpath,
                           paste0(str_replace(sitename, " ", "_"), ".tif")))

plot_results(sitename, output$sitel, output$datedat, sitearea, bmap,
             sites_sa, isobases, output, date_groups, scale_dist = 150,
             s_tdist = 1, s_xpos = 110, s_ypos = 50,  s_bheight = 0.5)

ggsave(file = here("analysis/figures/rodbol54.png"), width = 250,
       height = 152, units = "mm")

######### Rognlien #########
sitename <- "Rognlien"
date_groups <- c(1, 1, 2)

output <- apply_functions(sitename, date_groups, dtm, displacement_curves,
                          isobases, nsamp = 100, loc_bbox = 400, siterpath)

save(output,
     file = here::here("analysis/data/derived_data/rognlien.RData"))

sitearea <- rast(file.path(siterpath,
                           paste0(str_replace(sitename, " ", "_"), ".tif")))

plot_results(sitename, output$sitel, output$datedat, sitearea, bmap,
             sites_sa, isobases, output, date_groups, scale_dist = 150,
             s_tdist = 2, s_xpos = 150, s_ypos = 60, s_bheight = 1.5)

ggsave(file = here("analysis/figures/rognlien.png"), width = 250,
       height = 152, units = "mm")

######### Sagene B1 #########
sitename <- "Sagene B1"
date_groups <- c(1, 2)

output <- apply_functions(sitename, date_groups, dtm, displacement_curves,
                          isobases, nsamp = 100, loc_bbox = 400, siterpath)

save(output,
     file = here::here("analysis/data/derived_data/sageneb1.RData"))

sitearea <- rast(file.path(siterpath,
                           paste0(str_replace(sitename, " ", "_"), ".tif")))

plot_results(sitename, output$sitel, output$datedat, sitearea, bmap,
             sites_sa, isobases, output, date_groups, scale_dist = 150,
             s_tdist = 2, s_xpos = 150, s_ypos = 60, s_bheight = 1.5)

ggsave(file = here("analysis/figures/sageneb1.png"), width = 250,
       height = 152, units = "mm")

######### Sagene B2 #########
sitename <- "Sagene B2"
date_groups <- 1

output <- apply_functions(sitename, date_groups, dtm, displacement_curves,
                          isobases, nsamp = 100, loc_bbox = 400, siterpath)

save(output,
     file = here::here("analysis/data/derived_data/sageneb2.RData"))

sitearea <- rast(file.path(siterpath,
                           paste0(str_replace(sitename, " ", "_"), ".tif")))

plot_results(sitename, output$sitel, output$datedat, sitearea, bmap,
             sites_sa, isobases, output, date_groups, scale_dist = 150,
             s_tdist = 2, s_xpos = 150, s_ypos = 60, s_bheight = 1.5)

ggsave(file = here("analysis/figures/sageneb2.png"), width = 250,
       height = 152, units = "mm")

######### Solum 3 Will need a lot of processing time due to distance #########
sitename <- "Solum 3"
date_groups <- 1

output <- apply_functions(sitename, date_groups, dtm, displacement_curves,
                          isobases, nsamp = 1000, loc_bbox = 2500, siterpath)
save(output,
     file = here::here("analysis/data/derived_data/solum3.RData"))

sitearea <- rast(file.path(siterpath,
                           paste0(str_replace(sitename, " ", "_"), ".tif")))

plot_results(sitename, output$sitel, output$datedat, sitearea, bmap,
             sites_sa, isobases, output, date_groups, scale_dist = 1000,
             s_tdist = 15, s_xpos = 800, s_ypos = 300,  s_bheight = 8)

ggsave(file = here("analysis/figures/solum3.png"), width = 250,
       height = 152, units = "mm")

######### Stokke/Polland 1 #########
sitename <- "Stokke/Polland 1"
date_groups <- c(1, 2, 3, 3)

output <- apply_functions(sitename, date_groups, dtm, displacement_curves,
                          isobases, nsamp = 100, loc_bbox = 400, siterpath)

save(output,
     file = here::here("analysis/data/derived_data/stokkep1.RData"))

sitearea <- rast(file.path(siterpath,
                           paste0(str_replace(sitename, " ", "_"), ".tif")))

plot_results(sitename, output$sitel, output$datedat, sitearea, bmap,
             sites_sa, isobases, output, date_groups, scale_dist = 150,
             s_tdist = 2, s_xpos = 150, s_ypos = 60, s_bheight = 1.5)

ggsave(file = here("analysis/figures/stokkep1.png"), width = 250,
       height = 152, units = "mm")

######### Stokke/Polland 5 #########
sitename <- "Stokke/Polland 5"
date_groups <- c(1, 1, 1)

output <- apply_functions(sitename, date_groups, dtm, displacement_curves,
                          isobases, nsamp = 100, loc_bbox = 400, siterpath)

save(output,
     file = here::here("analysis/data/derived_data/stokkep5.RData"))

sitearea <- rast(file.path(siterpath,
                           paste0(str_replace(sitename, " ", "_"), ".tif")))

plot_results(sitename, output$sitel, output$datedat, sitearea, bmap,
             sites_sa, isobases, output, date_groups, scale_dist = 150,
             s_tdist = 2, s_xpos = 150, s_ypos = 60, s_bheight = 1.5)

ggsave(file = here("analysis/figures/stokkep5.png"), width = 250,
       height = 152, units = "mm")

######### Stokke/Polland 8 #########
sitename <- "Stokke/Polland 8"
date_groups <- 1

output <- apply_functions(sitename, date_groups, dtm, displacement_curves,
                          isobases, nsamp = 100, loc_bbox = 400, siterpath)

save(output,
     file = here::here("analysis/data/derived_data/stokkep8.RData"))

sitearea <- rast(file.path(siterpath,
                           paste0(str_replace(sitename, " ", "_"), ".tif")))

plot_results(sitename, output$sitel, output$datedat, sitearea, bmap,
             sites_sa, isobases, output, date_groups, scale_dist = 150,
             s_tdist = 2, s_xpos = 150, s_ypos = 60, s_bheight = 1.5)

ggsave(file = here("analysis/figures/stokkep8.png"), width = 250,
       height = 152, units = "mm")

######### Torstvet #########
sitename <- "Torstvet"
date_groups <- c(1, 1)

output <- apply_functions(sitename, date_groups, dtm, displacement_curves,
                          isobases, nsamp = 100, loc_bbox = 400, siterpath)

save(output,
     file = here::here("analysis/data/derived_data/torstvet.RData"))

sitearea <- rast(file.path(siterpath,
                           paste0(str_replace(sitename, " ", "_"), ".tif")))

plot_results(sitename, output$sitel, output$datedat, sitearea, bmap,
             sites_sa, isobases, output, date_groups, scale_dist = 150,
             s_tdist = 2, s_xpos = 150, s_ypos = 60, s_bheight = 1.5)

ggsave(file = here("analysis/figures/torstvet.png"), width = 250,
       height = 152, units = "mm")

######### Trondalen #########
sitename <- "Trondalen"
date_groups <- c(1, 1, 1)

output <- apply_functions(sitename, date_groups, dtm, displacement_curves,
                          isobases, nsamp = 100, loc_bbox = 400, siterpath)

save(output,
     file = here::here("analysis/data/derived_data/trondalen.RData"))

sitearea <- rast(file.path(siterpath,
                           paste0(str_replace(sitename, " ", "_"), ".tif")))

plot_results(sitename, output$sitel, output$datedat, sitearea, bmap,
             sites_sa, isobases, output, date_groups, scale_dist = 150,
             s_tdist = 2, s_xpos = 150, s_ypos = 60, s_bheight = 1.5)

ggsave(file = here("analysis/figures/trondalen.png"), width = 250,
       height = 152, units = "mm")

######### Tverdal #########
sitename <- "Tverdal"
date_groups <- rep(1, 4)

output <- apply_functions(sitename, date_groups, dtm, displacement_curves,
                          isobases, nsamp = 1000, loc_bbox = 400, siterpath)
save(output,
     file = here::here("analysis/data/derived_data/tverdal.RData"))
# load(here("analysis/data/derived_data/tverdal.RData"))

sitearea <- rast(file.path(siterpath,
                           paste0(str_replace(sitename, " ", "_"), ".tif")))

plot_results(sitename, output$sitel, output$datedat, sitearea, bmap,
             sites_sa, isobases, output, date_groups, scale_dist = 150,
             s_tdist = 0.5, s_xpos = 175, s_ypos = 85,  s_bheight = 0.25)

ggsave(file = here("analysis/figures/tverdal.png"), width = 250,
       height = 152, units = "mm")

######### Vallermyrene 1a #########
sitename <- "Vallermyrene 1a"
date_groups <- c(1, 1)

output <- apply_functions(sitename, date_groups, dtm, displacement_curves,
                           isobases, nsamp = 1000, loc_bbox = 400, siterpath)
# load(here("analysis/data/derived_data/vallermyrene1a.RData"))

save(output,
     file = here::here("analysis/data/derived_data/vallermyrene1a.RData"))

sitearea <- rast(file.path(siterpath,
                           paste0(str_replace(sitename, " ", "_"), ".tif")))

plot_results(sitename, output$sitel, output$datedat, sitearea, bmap,
             sites_sa, isobases, output, date_groups, scale_dist = 150)

ggsave(file = here("analysis/figures/vallermyrene1a.png"), width = 250,
       height = 152, units = "mm")

######### Vallermyrene 1b #########
sitename <- "Vallermyrene 1b"
date_groups <- 1

output <- apply_functions(sitename, date_groups, dtm, displacement_curves,
                           isobases, nsamp = 1000, loc_bbox = 400, siterpath)
# load(here("analysis/data/derived_data/vallermyrene1b.RData"))

save(output,
     file = here::here("analysis/data/derived_data/vallermyrene1b.RData"))

sitearea <- rast(file.path(siterpath,
                           paste0(str_replace(sitename, " ", "_"), ".tif")))

plot_results(sitename, output$sitel, output$datedat, sitearea, bmap,
             sites_sa, isobases, output, date_groups, scale_dist = 150)

ggsave(file = here("analysis/figures/vallermyrene1b.png"), width = 250,
       height = 152, units = "mm")

######### Vallermyrene 2 #########
sitename <- "Vallermyrene 2"
date_groups <- 1

# Use dated feature instead of site limit
output <- apply_functions(sitename, date_groups, dtm, displacement_curves,
                          isobases, nsamp = 1000, loc_bbox = 200, siterpath,
                          sitelimit = FALSE)
save(output,
     file = here::here("analysis/data/derived_data/vallermyrene2.RData"))

sitearea <- rast(file.path(siterpath,
                           paste0(str_replace(sitename, " ", "_"), ".tif")))

plot_results(sitename, output$sitel, output$datedat, sitearea, bmap,
             sites_sa, isobases, output, date_groups, scale_dist = 100)

ggsave(file = here("analysis/figures/vallermyrene2.png"), width = 250,
       height = 152, units = "mm")

######### Vallermyrene 4a #########
sitename <- "Vallermyrene 4a"
date_groups <- c(1, 1)

output <- apply_functions(sitename, date_groups, dtm, displacement_curves,
                           isobases, nsamp = 1000, loc_bbox = 400, siterpath)
# load(here("analysis/data/derived_data/vallermyrene4a.RData"))

save(output,
     file = here::here("analysis/data/derived_data/vallermyrene4a.RData"))

sitearea <- rast(file.path(siterpath,
                           paste0(str_replace(sitename, " ", "_"), ".tif")))

plot_results(sitename, output$sitel, output$datedat, sitearea, bmap,
             sites_sa, isobases, output, date_groups, scale_dist = 100)

ggsave(file = here("analysis/figures/vallermyrene4a.png"), width = 250,
       height = 152, units = "mm")

######### Vallermyrene 4b #########
sitename <- "Vallermyrene 4b"
date_groups <- c(1, 1)

output <- apply_functions(sitename, date_groups, dtm, displacement_curves,
                          isobases, nsamp = 1000, loc_bbox = 400, siterpath)
save(output,
     file = here::here("analysis/data/derived_data/vallermyrene4b.RData"))

sitearea <- rast(file.path(siterpath,
                           paste0(str_replace(sitename, " ", "_"), ".tif")))

plot_results(sitename, output$sitel, output$datedat, sitearea, bmap,
             sites_sa, isobases, output, date_groups, scale_dist = 150)

ggsave(file = here("analysis/figures/vallermyrene4b.png"), width = 250,
       height = 152, units = "mm")
