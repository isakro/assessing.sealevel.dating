library(here)
library(tidyverse)
library(sf)
library(ggridges)
library(ggthemes)
library(ggnewscale)
library(patchwork)
library(terra)
library(raster)
library(oxcAAR)
library(topoDistance)
library(IRanges)
quickSetupOxcal()

# To install IRanges:
# source("http://bioconductor.org/biocLite.R")
# biocLite("IRanges")

# If that fails, try:
# install.packages("BiocManager")
# BiocManager::install("IRanges")

# For reproducibility
set.seed(123)

# Load required functions and prepared data
source(here("analysis/script/04analysis_functions.R"))
load(here("analysis/data/derived_data/01data.RData"))
load(here("analysis/data/derived_data/02data.RData"))

# Load background map
bmap <- st_read(here("analysis/data/raw_data/naturalearth_countries.gpkg"))

# # Specify path to raster data (currently stored locally)
dtm <- rast("/home/isak/phd/eaa_presentation/dtm10/dtm10_32.tif")

# This tries to the load the edited DTM. Due to the file sizes involved,
# this is distributed as tiles which have to be merged by first
# running 00dtm_prep.R
dtm <- rast(here("analysis/data/derived_data/dtm/dtm10.tif"))

# Path to directory to hold smaller site area rasters
siterpath <- here("analysis/data/derived_data/sitearea_temp")

# Site to directory to hold simulated sea level raster
searasterpath <- here("analysis/data/derived_data/sealevel_rasters")

# Prespecified height of plot (to be multiplied by number of phases)
plot_height <- 76

# Specify number of simulation runs based on convergence below
# (1000 seems adequate as only minute variation is observable beyond this point)
nsamp = 1000

# Identifying number of simulation runs to perform per site by evaluating
# when mean distance values converges for a chosen test site with complex
# topographic surroundings and uncertain date range.
sitename <- "Hovland 5"
date_groups <- group_dates(rcarb_sa, sitename)

###### Uncomment to rerun simulation
# sims <- 5000
# output <- apply_functions(sitename, date_groups, dtm, displacement_curves,
#                             isobases, nsamp = sims, loc_bbox = 1000,
#                             siterpath)
# simresults <- data.frame(output[[1]]$results, "simn" = 1:sims)
#
# simdat <- simresults %>%
#   mutate(cmean_h = cummean(hordist),
#          cmean_t = cummean(topodist),
#          cmean_v = cummean(vertdist))

# Load simulation results
load(here("analysis/data/derived_data/05simdat.RData"))

simplt1 <- ggplot(simdat) +
  geom_line(aes(x = simn, y = cmean_h)) +
  labs(x = "Number of simulations",
       y = "Mean horisontal distance between site and shoreline") +
  theme_classic()

simplt2 <- ggplot(simdat) +
  geom_line(aes(x = simn, y = cmean_t)) +
  labs(x = "Number of simulations",
     y = "Mean topographic distance between site and shoreline") +
  theme_classic()

simplt3 <- ggplot(simdat) +
  geom_line(aes(x = simn, y = cmean_v)) +
  labs(x = "Number of simulations",
       y = "Mean vertical distance between site and shoreline") +
  theme_classic()

simplt1 + simplt2 + simplt3 +
plot_annotation("Test-run on Hovland 5 to inform number of simulation runs")

ggsave(file = here("analysis/figures/sample_size.png"), width = 300,
       height = 120,
       units = "mm")

# Some manual inspection was necessary in case of inconsistencies in the DTM
# (i.e. a highway running by the site) and for defining necessary size of the
# window of analysis. Each site was therefore first simulated 50 times and then
# rerun at 1000 when these issues were handled and parameters specified.

######### Adal vestre 1 #########
sitename <- "Adal vestre 1"
date_groups <- group_dates(rcarb_sa, sitename)

output <- apply_functions(sitename, date_groups, dtm, displacement_curves,
                          isobases, nsamp = nsamp, loc_bbox = 400, siterpath,
                          searasterpath)

save(output,
      file = here::here("analysis/data/derived_data/adalvestre1.RData"))

sitearea <- rast(file.path(siterpath,
                           paste0(str_replace(sitename, " ", "_"), ".tif")))

plot_results(sitename, output$sitel, output$datedat, sitearea, bmap,
             sites_sa, isobases, output, date_groups, scale_dist = 150,
             s_tdist = 0.5, s_xpos = 150, s_ypos = 70, s_bheight = 0.25,
             searasterpath = searasterpath)

ggsave(file = here("analysis/figures/adalvestre1.png"), width = 250,
       height = plot_height + (plot_height * length(unique(date_groups))),
       units = "mm")

######### Adal vestre 2 #########
sitename <- "Adal vestre 2"
date_groups <- group_dates(rcarb_sa, sitename)

output <- apply_functions(sitename, date_groups, dtm, displacement_curves,
                          isobases, nsamp = nsamp, loc_bbox = 400, siterpath,
                          searasterpath)
save(output,
     file = here::here("analysis/data/derived_data/adalvestre2.RData"))

sitearea <- rast(file.path(siterpath,
                           paste0(str_replace(sitename, " ", "_"), ".tif")))

plot_results(sitename, output$sitel, output$datedat, sitearea, bmap,
             sites_sa, isobases, output, date_groups, scale_dist = 150,
             s_tdist = 1, s_xpos = 150, s_ypos = 70, s_bheight = 0.5,
             searasterpath = searasterpath)

ggsave(file = here("analysis/figures/adalvestre2.png"), width = 250,
       height = plot_height + (plot_height * length(unique(date_groups))),
       units = "mm")

######### Alveberget 8 #########
sitename <- "Alveberget 8"
date_groups <- group_dates(rcarb_sa, sitename)

output <- apply_functions(sitename, date_groups, dtm, displacement_curves,
                          isobases, nsamp = nsamp, loc_bbox = 400, siterpath,
                          searasterpath)
save(output,
     file = here::here("analysis/data/derived_data/alveberget8.RData"))

sitearea <- rast(file.path(siterpath,
                           paste0(str_replace(sitename, " ", "_"), ".tif")))

plot_results(sitename, output$sitel, output$datedat, sitearea, bmap,
             sites_sa, isobases, output, date_groups, scale_dist = 150,
             s_tdist = 1, s_xpos = 150, s_ypos = 70, s_bheight = 0.75,
             searasterpath = searasterpath)

ggsave(file = here("analysis/figures/alveberget8.png"), width = 250,
       height = plot_height + (plot_height * length(unique(date_groups))),
       units = "mm")

######### Anvik #########
sitename <- "Anvik"
date_groups <- group_dates(rcarb_sa, sitename)

output <- apply_functions(sitename, date_groups, dtm, displacement_curves,
                          isobases, nsamp = nsamp, loc_bbox = 400, siterpath,
                          searasterpath)
save(output,
     file = here::here("analysis/data/derived_data/anvik.RData"))

sitearea <- rast(file.path(siterpath,
                           paste0(str_replace(sitename, " ", "_"), ".tif")))

plot_results(sitename, output$sitel, output$datedat, sitearea, bmap,
             sites_sa, isobases, output, date_groups, scale_dist = 150,
             s_tdist = 1, s_xpos = 150, s_ypos = 70, s_bheight = 0.75,
             searasterpath = searasterpath)

ggsave(file = here("analysis/figures/anvik.png"), width = 250,
       height = plot_height + (plot_height * length(unique(date_groups))),
       units = "mm")

######### Auve #########
sitename <- "Auve"
date_groups <- group_dates(rcarb_sa, sitename)

output <- apply_functions(sitename, date_groups, dtm, displacement_curves,
                          isobases, nsamp = nsamp, loc_bbox = 400, siterpath,
                          searasterpath)
save(output,
     file = here::here("analysis/data/derived_data/auve.RData"))

sitearea <- rast(file.path(siterpath,
                           paste0(str_replace(sitename, " ", "_"), ".tif")))

plot_results(sitename, output$sitel, output$datedat, sitearea, bmap,
             sites_sa, isobases, output, date_groups, scale_dist = 150,
             s_tdist = 0.35, s_xpos = 150, s_ypos = 70, s_bheight = 0.25,
             searasterpath = searasterpath)

ggsave(file = here("analysis/figures/auve.png"), width = 250,
       height = plot_height + (plot_height * length(unique(date_groups))),
       units = "mm")

######### Bakke #########
sitename <- "Bakke"
date_groups <- group_dates(rcarb_sa, sitename)

output <- apply_functions(sitename, date_groups, dtm, displacement_curves,
                          isobases, nsamp = nsamp, loc_bbox = 3500, siterpath,
                          searasterpath)
save(output,
     file = here::here("analysis/data/derived_data/bakke.RData"))

sitearea <- rast(file.path(siterpath,
                           paste0(str_replace(sitename, " ", "_"), ".tif")))

plot_results(sitename, output$sitel, output$datedat, sitearea, bmap,
             sites_sa, isobases, output, date_groups, scale_dist = 1250,
             s_tdist = 3, s_xpos = 1250, s_ypos = 500,  s_bheight = 1.5,
             searasterpath = searasterpath)

ggsave(file = here("analysis/figures/bakke.png"), width = 250,
       height = plot_height + (plot_height * length(unique(date_groups))),
       units = "mm")

######### Brunstad 24 #########
sitename <- "Brunstad 24"
date_groups <- group_dates(rcarb_sa, sitename)

output <- apply_functions(sitename, date_groups, dtm, displacement_curves,
                          isobases, nsamp = nsamp, loc_bbox = 400, siterpath,
                          searasterpath)
save(output,
     file = here::here("analysis/data/derived_data/brunstad24.RData"))

sitearea <- rast(file.path(siterpath,
                           paste0(str_replace(sitename, " ", "_"), ".tif")))

plot_results(sitename, output$sitel, output$datedat, sitearea, bmap,
             sites_sa, isobases, output, date_groups,  scale_dist = 150,
             s_tdist = 0.5, s_xpos = 150, s_ypos = 70, s_bheight = 0.25,
             searasterpath = searasterpath)

ggsave(file = here("analysis/figures/brunstad24.png"), width = 250,
       height = plot_height + (plot_height * length(unique(date_groups))),
       units = "mm")

######### Brunstad 25 #########
sitename <- "Brunstad 25"
date_groups <- group_dates(rcarb_sa, sitename)

output <- apply_functions(sitename, date_groups, dtm, displacement_curves,
                          isobases, nsamp = nsamp, loc_bbox = 400, siterpath,
                          searasterpath)
save(output,
     file = here::here("analysis/data/derived_data/brunstad25.RData"))

sitearea <- rast(file.path(siterpath,
                           paste0(str_replace(sitename, " ", "_"), ".tif")))

plot_results(sitename, output$sitel, output$datedat, sitearea, bmap,
             sites_sa, isobases, output, date_groups,  scale_dist = 150,
             s_tdist = 0.5, s_xpos = 150, s_ypos = 70, s_bheight = 0.25,
             searasterpath = searasterpath)

ggsave(file = here("analysis/figures/brunstad25.png"), width = 250,
       height = plot_height + (plot_height * length(unique(date_groups))),
       units = "mm")

######### Dørdal #########
sitename <- "Dørdal"
date_groups <- group_dates(rcarb_sa, sitename)

output <- apply_functions(sitename, date_groups, dtm, displacement_curves,
                          isobases, nsamp = nsamp, loc_bbox = 2000, siterpath,
                          searasterpath)
save(output,
     file = here::here("analysis/data/derived_data/dordal.RData"))

sitearea <- rast(file.path(siterpath,
                           paste0(str_replace(sitename, " ", "_"), ".tif")))

plot_results(sitename, output$sitel, output$datedat, sitearea, bmap,
             sites_sa, isobases, output, date_groups, scale_dist = 500,
             s_tdist = 5, s_xpos = 600, s_ypos = 200,  s_bheight = 3,
             searasterpath = searasterpath)

ggsave(file = here("analysis/figures/dordal.png"), width = 250,
       height = plot_height + (plot_height * length(unique(date_groups))),
       units = "mm")

######### Dybdalshei 1 #########
sitename <- "Dybdalshei 1"
date_groups <- group_dates(rcarb_sa, sitename)

output <- apply_functions(sitename, date_groups, dtm, displacement_curves,
                          isobases, nsamp = nsamp, loc_bbox = 400, siterpath,
                          searasterpath)
save(output,
     file = here::here("analysis/data/derived_data/dybdalshei1.RData"))

sitearea <- rast(file.path(siterpath,
                           paste0(str_replace(sitename, " ", "_"), ".tif")))

plot_results(sitename, output$sitel, output$datedat, sitearea, bmap,
             sites_sa, isobases, output, date_groups, scale_dist = 150,
             s_tdist = 0.5, s_xpos = 135, s_ypos = 65,  s_bheight = 0.25,
             searasterpath = searasterpath)

ggsave(file = here("analysis/figures/dybdalshei1.png"), width = 250,
       height = plot_height + (plot_height * length(unique(date_groups))),
       units = "mm")

######### Dybdalshei 2 #########
sitename <- "Dybdalshei 2"
date_groups <- group_dates(rcarb_sa, sitename)

output <- apply_functions(sitename, date_groups, dtm, displacement_curves,
                          isobases, nsamp = nsamp, loc_bbox = 400, siterpath,
                          searasterpath)
save(output,
      file = here::here("analysis/data/derived_data/dybdalshei2.RData"))

sitearea <- rast(file.path(siterpath,
                           paste0(str_replace(sitename, " ", "_"), ".tif")))

plot_results(sitename, output$sitel, output$datedat, sitearea, bmap,
             sites_sa, isobases, output, date_groups, scale_dist = 150,
             s_tdist = 2, s_xpos = 150, s_ypos = 60, s_bheight = 1.5,
             searasterpath = searasterpath)

ggsave(file = here("analysis/figures/dybdalshei2.png"), width = 250,
       height = plot_height + (plot_height * length(unique(date_groups))),
       units = "mm")

######### Gunnarsrød 4 #########
sitename <- "Gunnarsrød 4"
# Oldest interpreted as unrelated in original report, manually seperated
date_groups <- c(2, 1)

output <- apply_functions(sitename, date_groups, dtm, displacement_curves,
                          isobases, nsamp = nsamp, loc_bbox = 400, siterpath,
                          searasterpath)
save(output,
     file = here::here("analysis/data/derived_data/gunnarsrod4.RData"))

sitearea <- rast(file.path(siterpath,
                           paste0(str_replace(sitename, " ", "_"), ".tif")))

plot_results(sitename, output$sitel, output$datedat, sitearea, bmap,
             sites_sa, isobases, output, date_groups, scale_dist = 150,
             s_tdist = 9, s_xpos = 150, s_ypos = 70, s_bheight = 5,
             searasterpath = searasterpath)

ggsave(file = here("analysis/figures/gunnarsrod4.png"), width = 250,
       height = plot_height + (plot_height * length(unique(date_groups))),
       units = "mm")

######### Gunnarsrød 5 #########
sitename <- "Gunnarsrød 5"
date_groups <- group_dates(rcarb_sa, sitename)

output <- apply_functions(sitename, date_groups, dtm, displacement_curves,
                          isobases, nsamp = nsamp, loc_bbox = 400, siterpath,
                          searasterpath)
save(output,
     file = here::here("analysis/data/derived_data/gunnarsrod5.RData"))

sitearea <- rast(file.path(siterpath,
                           paste0(str_replace(sitename, " ", "_"), ".tif")))

plot_results(sitename, output$sitel, output$datedat, sitearea, bmap,
             sites_sa, isobases, output, date_groups, scale_dist = 150,
             s_tdist = 0.5, s_xpos = 135, s_ypos = 65,  s_bheight = 0.25,
             searasterpath = searasterpath)

ggsave(file = here("analysis/figures/gunnarsrod5.png"), width = 250,
       height = plot_height + (plot_height * length(unique(date_groups))),
       units = "mm")

######### Gunnarsrød 6a #########
sitename <- "Gunnarsrød 6a"
date_groups <- group_dates(rcarb_sa, sitename)

output <- apply_functions(sitename, date_groups, dtm, displacement_curves,
                          isobases, nsamp = nsamp, loc_bbox = 400, siterpath,
                          searasterpath)
save(output,
     file = here::here("analysis/data/derived_data/gunnarsrod6a.RData"))

sitearea <- rast(file.path(siterpath,
                           paste0(str_replace(sitename, " ", "_"), ".tif")))

plot_results(sitename, output$sitel, output$datedat, sitearea, bmap,
             sites_sa, isobases, output, date_groups, scale_dist = 150,
             s_tdist = 1, s_xpos = 135, s_ypos = 65,  s_bheight = 0.5,
             searasterpath = searasterpath)

ggsave(file = here("analysis/figures/gunnarsrod6a.png"), width = 250,
       height = plot_height + (plot_height * length(unique(date_groups))),
       units = "mm")

######### Gunnarsrød 6b #########
sitename <- "Gunnarsrød 6b"
date_groups <- 1

output <- apply_functions(sitename, date_groups, dtm, displacement_curves,
                          isobases, nsamp = nsamp, loc_bbox = 400, siterpath,
                          searasterpath)
save(output,
     file = here::here("analysis/data/derived_data/gunnarsrod6b.RData"))

sitearea <- rast(file.path(siterpath,
                           paste0(str_replace(sitename, " ", "_"), ".tif")))

plot_results(sitename, output$sitel, output$datedat, sitearea, bmap,
             sites_sa, isobases, output, date_groups, scale_dist = 150,
             s_tdist = 2, s_xpos = 135, s_ypos = 65,  s_bheight = 1,
             searasterpath = searasterpath)

ggsave(file = here("analysis/figures/gunnarsrod6b.png"), width = 250,
       height = plot_height + (plot_height * length(unique(date_groups))),
       units = "mm")

######### Gunnarsrød 7 #########
sitename <- "Gunnarsrød 7"
date_groups <- group_dates(rcarb_sa, sitename)

output <- apply_functions(sitename, date_groups, dtm, displacement_curves,
                          isobases, nsamp = nsamp, loc_bbox = 400, siterpath,
                          searasterpath)
save(output,
     file = here::here("analysis/data/derived_data/gunnarsrod7.RData"))

sitearea <- rast(file.path(siterpath,
                           paste0(str_replace(sitename, " ", "_"), ".tif")))

plot_results(sitename, output$sitel, output$datedat, sitearea, bmap,
             sites_sa, isobases, output, date_groups, scale_dist = 150,
             s_tdist = 1, s_xpos = 135, s_ypos = 65,  s_bheight = 0.5,
             searasterpath = searasterpath)

ggsave(file = here("analysis/figures/gunnarsrod7.png"), width = 250,
       height = plot_height + (plot_height * length(unique(date_groups))),
       units = "mm")

######### Gunnarsrød 10 #########
sitename <- "Gunnarsrød 10"
date_groups <- group_dates(rcarb_sa, sitename)

output <- apply_functions(sitename, date_groups, dtm, displacement_curves,
                          isobases, nsamp = nsamp, loc_bbox = 400, siterpath,
                          searasterpath)
save(output,
     file = here::here("analysis/data/derived_data/gunnarsrod10.RData"))

sitearea <- rast(file.path(siterpath,
                           paste0(str_replace(sitename, " ", "_"), ".tif")))

plot_results(sitename, output$sitel, output$datedat, sitearea, bmap,
             sites_sa, isobases, output, date_groups, scale_dist = 150,
             s_tdist = 3.5, s_xpos = 135, s_ypos = 65,  s_bheight = 2,
             searasterpath = searasterpath)

ggsave(file = here("analysis/figures/gunnarsrod10.png"), width = 250,
       height = plot_height + (plot_height * length(unique(date_groups))),
       units = "mm")

######### Hegna vest 1 #########
sitename <- "Hegna vest 1"
date_groups <- group_dates(rcarb_sa, sitename)

output <- apply_functions(sitename, date_groups, dtm, displacement_curves,
                          isobases, nsamp = nsamp, loc_bbox = 400, siterpath,
                          searasterpath)
save(output,
     file = here::here("analysis/data/derived_data/hegnavest1.RData"))

sitearea <- rast(file.path(siterpath,
                           paste0(str_replace(sitename, " ", "_"), ".tif")))

plot_results(sitename, output$sitel, output$datedat, sitearea, bmap,
             sites_sa, isobases, output, date_groups, scale_dist = 150,
             s_tdist = 0.7, s_xpos = 135, s_ypos = 65,  s_bheight = 0.5,
             searasterpath = searasterpath)

ggsave(file = here("analysis/figures/hegnavest1.png"), width = 250,
       height = plot_height + (plot_height * length(unique(date_groups))),
       units = "mm")

######### Hegna vest 2 #########
sitename <- "Hegna vest 2"
date_groups <- group_dates(rcarb_sa, sitename)

output <- apply_functions(sitename, date_groups, dtm, displacement_curves,
                          isobases, nsamp = nsamp, loc_bbox = 1200, siterpath,
                          searasterpath)
save(output,
     file = here::here("analysis/data/derived_data/hegnavest2.RData"))

sitearea <- rast(file.path(siterpath,
                           paste0(str_replace(sitename, " ", "_"), ".tif")))

plot_results(sitename, output$sitel, output$datedat, sitearea, bmap,
             sites_sa, isobases, output, date_groups, scale_dist = 600,
             s_tdist = 1.25, s_xpos = 500, s_ypos = 250,  s_bheight = 0.8,
             searasterpath = searasterpath)

ggsave(file = here("analysis/figures/hegnavest2.png"), width = 250,
       height = plot_height + (plot_height * length(unique(date_groups))),
       units = "mm")

######### Hegna vest 3 #########
sitename <- "Hegna vest 3"
date_groups <- group_dates(rcarb_sa, sitename)

output <- apply_functions(sitename, date_groups, dtm, displacement_curves,
                          isobases, nsamp = nsamp, loc_bbox = 400, siterpath,
                          searasterpath)
save(output,
     file = here::here("analysis/data/derived_data/hegnavest3.RData"))

sitearea <- rast(file.path(siterpath,
                           paste0(str_replace(sitename, " ", "_"), ".tif")))

plot_results(sitename, output$sitel, output$datedat, sitearea, bmap,
             sites_sa, isobases, output, date_groups, scale_dist = 150,
             s_tdist = 1, s_xpos = 135, s_ypos = 65,  s_bheight = 0.5,
             searasterpath = searasterpath)

ggsave(file = here("analysis/figures/hegnavest3.png"), width = 250,
       height = plot_height + (plot_height * length(unique(date_groups))),
       units = "mm")

######### Hesthag C2 #########
sitename <- "Hesthag C2"
date_groups <- group_dates(rcarb_sa, sitename)

output <- apply_functions(sitename, date_groups, dtm, displacement_curves,
                          isobases, nsamp = nsamp, loc_bbox = 400, siterpath,
                          searasterpath)
save(output,
     file = here::here("analysis/data/derived_data/hesthagc2.RData"))

sitearea <- rast(file.path(siterpath,
                           paste0(str_replace(sitename, " ", "_"), ".tif")))

plot_results(sitename, output$sitel, output$datedat, sitearea, bmap,
             sites_sa, isobases, output, date_groups, scale_dist = 150,
             s_tdist = 0.5, s_xpos = 135, s_ypos = 65,  s_bheight = 0.25,
             searasterpath = searasterpath)

ggsave(file = here("analysis/figures/hesthagc2.png"), width = 250,
       height = plot_height + (plot_height * length(unique(date_groups))),
       units = "mm")

######### Hesthag C4 #########
sitename <- "Hesthag C4"
date_groups <- group_dates(rcarb_sa, sitename)

output <- apply_functions(sitename, date_groups, dtm, displacement_curves,
                          isobases, nsamp = nsamp, loc_bbox = 400, siterpath,
                          searasterpath)
save(output,
     file = here::here("analysis/data/derived_data/hesthagc4.RData"))

sitearea <- rast(file.path(siterpath,
                           paste0(str_replace(sitename, " ", "_"), ".tif")))

plot_results(sitename, output$sitel, output$datedat, sitearea, bmap,
             sites_sa, isobases, output, date_groups, scale_dist = 150,
             s_tdist = 1, s_xpos = 135, s_ypos = 65,  s_bheight = 0.5,
             searasterpath = searasterpath)

ggsave(file = here("analysis/figures/hesthagc4.png"), width = 250,
       height = plot_height + (plot_height * length(unique(date_groups))),
       units = "mm")

######### Hovland 1 #########
sitename <- "Hovland 1"
date_groups <- group_dates(rcarb_sa, sitename)

output <- apply_functions(sitename, date_groups, dtm, displacement_curves,
                          isobases, nsamp = nsamp, loc_bbox = 1000, siterpath,
                          searasterpath)
save(output,
     file = here::here("analysis/data/derived_data/hovland1.RData"))

sitearea <- rast(file.path(siterpath,
                           paste0(str_replace(sitename, " ", "_"), ".tif")))

plot_results(sitename, output$sitel, output$datedat, sitearea, bmap,
             sites_sa, isobases, output, date_groups, scale_dist = 250,
             s_tdist = 3, s_xpos = 350, s_ypos = 150,  s_bheight = 2.25,
             searasterpath = searasterpath)

ggsave(file = here("analysis/figures/hovland1.png"), width = 250,
       height = plot_height + (plot_height * length(unique(date_groups))),
       units = "mm")

######### Hovland 3 #########
sitename <- "Hovland 3"
date_groups <- group_dates(rcarb_sa, sitename)

output <- apply_functions(sitename, date_groups, dtm, displacement_curves,
                          isobases, nsamp = 10, loc_bbox = 400, siterpath,
                          searasterpath)
save(output,
     file = here::here("analysis/data/derived_data/hovland3.RData"))

sitearea <- rast(file.path(siterpath,
                           paste0(str_replace(sitename, " ", "_"), ".tif")))

plot_results(sitename, output$sitel, output$datedat, sitearea, bmap,
             sites_sa, isobases, output, date_groups, scale_dist = 150,
             s_tdist = 1, s_xpos = 135, s_ypos = 65,  s_bheight = 0.5,
             searasterpath = searasterpath)

ggsave(file = here("analysis/figures/hovland3.png"), width = 250,
       height = plot_height + (plot_height * length(unique(date_groups))),
                               units = "mm")

######### Hovland 4 #########
sitename <- "Hovland 4"
date_groups <- group_dates(rcarb_sa, sitename)

output <- apply_functions(sitename, date_groups, dtm, displacement_curves,
                          isobases, nsamp = nsamp, loc_bbox = 1200, siterpath,
                          searasterpath)
save(output,
     file = here::here("analysis/data/derived_data/hovland4.RData"))

sitearea <- rast(file.path(siterpath,
                           paste0(str_replace(sitename, " ", "_"), ".tif")))

plot_results(sitename, output$sitel, output$datedat, sitearea, bmap,
             sites_sa, isobases, output, date_groups, scale_dist = 400,
             s_tdist = 2.5, s_xpos = 500, s_ypos = 200,  s_bheight = 2,
             searasterpath = searasterpath)

ggsave(file = here("analysis/figures/hovland4.png"), width = 250,
       height = plot_height + (plot_height * length(unique(date_groups))),
       units = "mm")

######### Hovland 5 #########
sitename <- "Hovland 5"
date_groups <- group_dates(rcarb_sa, sitename)

output <- apply_functions(sitename, date_groups, dtm, displacement_curves,
                          isobases, nsamp = nsamp, loc_bbox = 1000, siterpath,
                          searasterpath)
save(output,
     file = here::here("analysis/data/derived_data/hovland5.RData"))

sitearea <- rast(file.path(siterpath,
                           paste0(str_replace(sitename, " ", "_"), ".tif")))

plot_results(sitename, output$sitel, output$datedat, sitearea, bmap,
             sites_sa, isobases, output, date_groups, scale_dist = 300,
             s_tdist = 2.5, s_xpos = 500, s_ypos = 200,  s_bheight = 2,
             searasterpath = searasterpath)

ggsave(file = here("analysis/figures/hovland5.png"), width = 250,
       height = plot_height + (plot_height * length(unique(date_groups))),
       units = "mm")

######### Hydal 4 #########
sitename <- "Hydal 4"
date_groups <- group_dates(rcarb_sa, sitename)

output <- apply_functions(sitename, date_groups, dtm, displacement_curves,
                          isobases, nsamp = nsamp, loc_bbox = 400, siterpath,
                          searasterpath)
save(output,
     file = here::here("analysis/data/derived_data/hydal4.RData"))

sitearea <- rast(file.path(siterpath,
                           paste0(str_replace(sitename, " ", "_"), ".tif")))

plot_results(sitename, output$sitel, output$datedat, sitearea, bmap,
             sites_sa, isobases, output, date_groups, scale_dist = 150,
             s_tdist = 1.5, s_xpos = 135, s_ypos = 65,  s_bheight = 1,
             searasterpath = searasterpath)

ggsave(file = here("analysis/figures/hydal4.png"), width = 250,
       height = plot_height + (plot_height * length(unique(date_groups))),
                               units = "mm")

######### Krøgenes D1 #########
sitename <- "Krøgenes D1"
date_groups <- group_dates(rcarb_sa, sitename)

output <- apply_functions(sitename, date_groups, dtm, displacement_curves,
                          isobases, nsamp = nsamp, loc_bbox = 400, siterpath,
                          searasterpath)
save(output,
     file = here::here("analysis/data/derived_data/krogenesd1.RData"))

sitearea <- rast(file.path(siterpath,
                           paste0(str_replace(sitename, " ", "_"), ".tif")))

plot_results(sitename, output$sitel, output$datedat, sitearea, bmap,
             sites_sa, isobases, output, date_groups, scale_dist = 150,
             searasterpath = searasterpath)

ggsave(file = here("analysis/figures/krogenesd1.png"), width = 250,
       height = plot_height + (plot_height * length(unique(date_groups))),
       units = "mm")

######### Krøgenes D2 #########
sitename <- "Krøgenes D2"
date_groups <- group_dates(rcarb_sa, sitename)

output <- apply_functions(sitename, date_groups, dtm, displacement_curves,
                          isobases, nsamp = nsamp, loc_bbox = 400, siterpath,
                          searasterpath)
save(output,
     file = here::here("analysis/data/derived_data/krogenesd2.RData"))

sitearea <- rast(file.path(siterpath,
                           paste0(str_replace(sitename, " ", "_"), ".tif")))

plot_results(sitename, output$sitel, output$datedat, sitearea, bmap,
             sites_sa, isobases, output, date_groups, scale_dist = 150,
             s_tdist = 1.5, s_xpos = 135, s_ypos = 65,  s_bheight = 1,
             searasterpath = searasterpath)

ggsave(file = here("analysis/figures/krogenesd2.png"), width = 250,
       height = plot_height + (plot_height * length(unique(date_groups))),
       units = "mm")

######### Kvastad A1 #########
sitename <- "Kvastad A1"
date_groups <- group_dates(rcarb_sa, sitename)

output <- apply_functions(sitename, date_groups, dtm, displacement_curves,
                          isobases, nsamp = nsamp, loc_bbox = 400, siterpath,
                          searasterpath)
save(output,
     file = here::here("analysis/data/derived_data/kvastada1.RData"))

sitearea <- rast(file.path(siterpath,
                           paste0(str_replace(sitename, " ", "_"), ".tif")))

plot_results(sitename, output$sitel, output$datedat, sitearea, bmap,
             sites_sa, isobases, output, date_groups, scale_dist = 150,
             s_tdist = 0.5, s_xpos = 135, s_ypos = 65,  s_bheight = 0.25,
             searasterpath = searasterpath)

ggsave(file = here("analysis/figures/kvastada1.png"), width = 250,
       height = plot_height + (plot_height * length(unique(date_groups))),
       units = "mm")

######### Kvastad A2 Both areas collapsed #########
sitename <- "Kvastad A2"
date_groups <- group_dates(rcarb_sa, sitename)

output <- apply_functions(sitename, date_groups, dtm, displacement_curves,
                          isobases, nsamp = nsamp, loc_bbox = 500, siterpath,
                          searasterpath)
save(output,
     file = here::here("analysis/data/derived_data/kvastada2.RData"))

sitearea <- rast(file.path(siterpath,
                           paste0(str_replace(sitename, " ", "_"), ".tif")))

plot_results(sitename, output$sitel, output$datedat, sitearea, bmap,
             sites_sa, isobases, output, date_groups, scale_dist = 150,
             searasterpath = searasterpath)

ggsave(file = here("analysis/figures/kvastada2.png"), width = 250,
       height = plot_height + (plot_height * length(unique(date_groups))),
       units = "mm")

######### Kvastad A4 #########
sitename <- "Kvastad A4"
date_groups <- group_dates(rcarb_sa, sitename)

output <- apply_functions(sitename, date_groups, dtm, displacement_curves,
                          isobases, nsamp = nsamp, loc_bbox = 400, siterpath,
                          searasterpath)
save(output,
     file = here::here("analysis/data/derived_data/kvastada4.RData"))

sitearea <- rast(file.path(siterpath,
                           paste0(str_replace(sitename, " ", "_"), ".tif")))

plot_results(sitename, output$sitel, output$datedat, sitearea, bmap,
             sites_sa, isobases, output, date_groups, scale_dist = 150,
             s_tdist = 1, s_xpos = 135, s_ypos = 65,  s_bheight = 0.5,
             searasterpath = searasterpath)

ggsave(file = here("analysis/figures/kvastada4.png"), width = 250,
       height = plot_height + (plot_height * length(unique(date_groups))),
                               units = "mm")

######### Kvastad A9 #########
sitename <- "Kvastad A9"
date_groups <- group_dates(rcarb_sa, sitename)

output <- apply_functions(sitename, date_groups, dtm, displacement_curves,
                           isobases, nsamp = nsamp, loc_bbox = 400, siterpath,
                          searasterpath)
save(output,
     file = here::here("analysis/data/derived_data/kvastada9.RData"))

sitearea <- rast(file.path(siterpath,
                           paste0(str_replace(sitename, " ", "_"), ".tif")))

plot_results(sitename, output$sitel, output$datedat, sitearea, bmap,
             sites_sa, isobases, output, date_groups, scale_dist = 150,
             s_tdist = 5, s_xpos = 140, s_ypos = 60,  s_bheight = 3.5,
             searasterpath = searasterpath)

ggsave(file = here("analysis/figures/kvastada9.png"), width = 250,
       height = plot_height + (plot_height * length(unique(date_groups))),
       units = "mm")

######### Langangen Vestgård 1 #########
sitename <- "Langangen Vestgård 1"
date_groups <- group_dates(rcarb_sa, sitename)

output <- apply_functions(sitename, date_groups, dtm, displacement_curves,
                          isobases, nsamp = nsamp, loc_bbox = 400, siterpath,
                          searasterpath)
save(output,
     file = here::here("analysis/data/derived_data/langangenv1.RData"))

sitearea <- rast(file.path(siterpath,
                           paste0(str_replace(sitename, " ", "_"), ".tif")))

plot_results(sitename, output$sitel, output$datedat, sitearea, bmap,
             sites_sa, isobases, output, date_groups, scale_dist = 150,
             searasterpath = searasterpath)

ggsave(file = here("analysis/figures/langangenv1.png"), width = 250,
       height = plot_height + (plot_height * length(unique(date_groups))),
       units = "mm")

######### Langangen Vestgård 3 #########
sitename <- "Langangen Vestgård 3"
date_groups <- group_dates(rcarb_sa, sitename)

output <- apply_functions(sitename, date_groups, dtm, displacement_curves,
                          isobases, nsamp = nsamp, loc_bbox = 400, siterpath,
                          searasterpath)
save(output,
     file = here::here("analysis/data/derived_data/langangenv3.RData"))

sitearea <- rast(file.path(siterpath,
                           paste0(str_replace(sitename, " ", "_"), ".tif")))

plot_results(sitename, output$sitel, output$datedat, sitearea, bmap,
             sites_sa, isobases, output, date_groups, scale_dist = 150,
             s_tdist = 1.5, s_xpos = 150, s_ypos = 80,  s_bheight = 1.25,
             searasterpath = searasterpath)

ggsave(file = here("analysis/figures/langangenv3.png"), width = 250,
       height = plot_height + (plot_height * length(unique(date_groups))),
       units = "mm")

######### Langangen Vestgård 5 #########
sitename <- "Langangen Vestgård 5"
date_groups <- group_dates(rcarb_sa, sitename)

output <- apply_functions(sitename, date_groups, dtm, displacement_curves,
                          isobases, nsamp = nsamp, loc_bbox = 400, siterpath,
                          searasterpath)
save(output,
     file = here::here("analysis/data/derived_data/langangenv5.RData"))

sitearea <- rast(file.path(siterpath,
                           paste0(str_replace(sitename, " ", "_"), ".tif")))

plot_results(sitename, output$sitel, output$datedat, sitearea, bmap,
             sites_sa, isobases, output, date_groups, scale_dist = 150,
             searasterpath = searasterpath)

ggsave(file = here("analysis/figures/langangenv5.png"), width = 250,
       height = plot_height + (plot_height * length(unique(date_groups))),
       units = "mm")

######### Langangen Vestgård 6 #########
sitename <- "Langangen Vestgård 6"
date_groups <- group_dates(rcarb_sa, sitename)

output <- apply_functions(sitename, date_groups, dtm, displacement_curves,
                          isobases, nsamp = nsamp, loc_bbox = 400, siterpath,
                          searasterpath)
save(output,
     file = here::here("analysis/data/derived_data/langangenv6.RData"))

sitearea <- rast(file.path(siterpath,
                           paste0(str_replace(sitename, " ", "_"), ".tif")))

plot_results(sitename, output$sitel, output$datedat, sitearea, bmap,
             sites_sa, isobases, output, date_groups, scale_dist = 150,
             s_tdist = 1, s_xpos = 170, s_ypos = 80,  s_bheight = 0.5,
             searasterpath = searasterpath)

ggsave(file = here("analysis/figures/langangenv6.png"), width = 250,
       height = plot_height + (plot_height * length(unique(date_groups))),
       units = "mm")

######### Langangen Vestgård 7 Road as site limit #########
sitename <- "Langangen Vestgård 7"
date_groups <- group_dates(rcarb_sa, sitename)

# Use dated feature instead of site limit
output <- apply_functions(sitename, date_groups, dtm, displacement_curves,
                          isobases, nsamp = nsamp, loc_bbox = 200, siterpath,
                          sitelimit = FALSE, searasterpath = searasterpath)
save(output,
     file = here::here("analysis/data/derived_data/langangenv7.RData"))

sitearea <- rast(file.path(siterpath,
                           paste0(str_replace(sitename, " ", "_"), ".tif")))

plot_results(sitename, output$sitel, output$datedat, sitearea, bmap,
             sites_sa, isobases, output, date_groups, scale_dist = 100,
             s_tdist = 1, s_xpos = 100, s_ypos = 50,  s_bheight = 10,
             searasterpath = searasterpath)

ggsave(file = here("analysis/figures/langangenv7.png"), width = 250,
       height =  plot_height + (plot_height * length(unique(date_groups))),
       units = "mm")

######### Langemyr #########
sitename <- "Langemyr"
date_groups <- group_dates(rcarb_sa, sitename)

output <- apply_functions(sitename, date_groups, dtm, displacement_curves,
                          isobases, nsamp = nsamp, loc_bbox = 400, siterpath,
                          searasterpath)
save(output,
     file = here::here("analysis/data/derived_data/langemyr.RData"))

sitearea <- rast(file.path(siterpath,
                           paste0(str_replace(sitename, " ", "_"), ".tif")))

plot_results(sitename, output$sitel, output$datedat, sitearea, bmap,
             sites_sa, isobases, output, date_groups, scale_dist = 150,
             searasterpath = searasterpath)

ggsave(file = here("analysis/figures/langemyr.png"), width = 250,
       height = plot_height + (plot_height * length(unique(date_groups))),
       units = "mm")

######### Lunaveien #########
sitename <- "Lunaveien"
date_groups <- group_dates(rcarb_sa, sitename)

output <- apply_functions(sitename, date_groups, dtm, displacement_curves,
                          isobases, nsamp = nsamp, loc_bbox = 3000, siterpath,
                          searasterpath)
save(output,
     file = here::here("analysis/data/derived_data/lunaveien.RData"))

sitearea <- rast(file.path(siterpath,
                           paste0(str_replace(sitename, " ", "_"), ".tif")))

plot_results(sitename, output$sitel, output$datedat, sitearea, bmap,
             sites_sa, isobases, output, date_groups,  scale_dist = 1250,
             s_tdist = 9, s_xpos = 1250, s_ypos = 500,  s_bheight = 4,
             searasterpath = searasterpath)

ggsave(file = here("analysis/figures/lunaveien.png"), width = 250,
       height = plot_height + (plot_height * length(unique(date_groups))),
       units = "mm")

######### Løvås 1 #########
sitename <- "Løvås 1"
date_groups <- group_dates(rcarb_sa, sitename)

output <- apply_functions(sitename, date_groups, dtm, displacement_curves,
                          isobases, nsamp = nsamp, loc_bbox = 400, siterpath,
                          searasterpath)
save(output,
     file = here::here("analysis/data/derived_data/lovas1.RData"))

sitearea <- rast(file.path(siterpath,
                           paste0(str_replace(sitename, " ", "_"), ".tif")))

plot_results(sitename, output$sitel, output$datedat, sitearea, bmap,
             sites_sa, isobases, output, date_groups, scale_dist = 150,
             s_tdist = 1, s_xpos = 150, s_ypos = 70, s_bheight = 0.5,
             searasterpath = searasterpath)

ggsave(file = here("analysis/figures/lovas1.png"), width = 250,
       height = plot_height + (plot_height * length(unique(date_groups))),
       units = "mm")

######### Løvås 2 #########
sitename <- "Løvås 2"
date_groups <- group_dates(rcarb_sa, sitename)

output <- apply_functions(sitename, date_groups, dtm, displacement_curves,
                          isobases, nsamp = nsamp, loc_bbox = 400, siterpath,
                          searasterpath)
save(output,
     file = here::here("analysis/data/derived_data/lovas2.RData"))

sitearea <- rast(file.path(siterpath,
                           paste0(str_replace(sitename, " ", "_"), ".tif")))

plot_results(sitename, output$sitel, output$datedat, sitearea, bmap,
             sites_sa, isobases, output, date_groups, scale_dist = 150,
             s_tdist = 1, s_xpos = 150, s_ypos = 70, s_bheight = 0.5,
             searasterpath = searasterpath)

ggsave(file = here("analysis/figures/lovas2.png"), width = 250,
       height = plot_height + (plot_height * length(unique(date_groups))),
       units = "mm")

######### Løvås 3 #########
sitename <- "Løvås 3"
date_groups <- group_dates(rcarb_sa, sitename)

output <- apply_functions(sitename, date_groups, dtm, displacement_curves,
                          isobases, nsamp = nsamp, loc_bbox = 400, siterpath,
                          searasterpath)
save(output,
     file = here::here("analysis/data/derived_data/lovas3.RData"))

sitearea <- rast(file.path(siterpath,
                           paste0(str_replace(sitename, " ", "_"), ".tif")))

plot_results(sitename, output$sitel, output$datedat, sitearea, bmap,
             sites_sa, isobases, output, date_groups, scale_dist = 150,
             s_tdist = 0.5, s_xpos = 150, s_ypos = 70, s_bheight = 0.25,
             searasterpath = searasterpath)

ggsave(file = here("analysis/figures/lovas3.png"), width = 250,
       height = plot_height + (plot_height * length(unique(date_groups))),
       units = "mm")

######### Nauen A #########
sitename <- "Nauen A"
date_groups <- group_dates(rcarb_sa, sitename)

output <- apply_functions(sitename, date_groups, dtm, displacement_curves,
                          isobases, nsamp = nsamp, loc_bbox = 1000, siterpath,
                          searasterpath)
save(output,
     file = here::here("analysis/data/derived_data/nauena.RData"))

sitearea <- rast(file.path(siterpath,
                           paste0(str_replace(sitename, " ", "_"), ".tif")))

plot_results(sitename, output$sitel, output$datedat, sitearea, bmap,
             sites_sa, isobases, output, date_groups, scale_dist = 400,
             s_tdist = 1, s_xpos = 400, s_ypos = 160,  s_bheight = 0.5,
             searasterpath = searasterpath)

ggsave(file = here("analysis/figures/nauena.png"), width = 250,
       height = plot_height + (plot_height * length(unique(date_groups))),
       units = "mm")

######### Nauen C #########
sitename <- "Nauen C"
date_groups <- group_dates(rcarb_sa, sitename)

output <- apply_functions(sitename, date_groups, dtm, displacement_curves,
                          isobases, nsamp = nsamp, loc_bbox = 1000, siterpath,
                          searasterpath)
save(output,
     file = here::here("analysis/data/derived_data/nauenc.RData"))

sitearea <- rast(file.path(siterpath,
                           paste0(str_replace(sitename, " ", "_"), ".tif")))

plot_results(sitename, output$sitel, output$datedat, sitearea, bmap,
             sites_sa, isobases, output, date_groups, scale_dist = 400,
             s_tdist = 1.5, s_xpos = 400, s_ypos = 160,  s_bheight = 1,
             searasterpath = searasterpath)

ggsave(file = here("analysis/figures/nauenc.png"), width = 250,
       height = plot_height + (plot_height * length(unique(date_groups))),
       units = "mm")

######### Nordby 1 #########
sitename <- "Nordby 1"
date_groups <- group_dates(rcarb_sa, sitename)

output <- apply_functions(sitename, date_groups, dtm, displacement_curves,
                          isobases, nsamp = nsamp, loc_bbox = 400, siterpath,
                          searasterpath)
save(output,
     file = here::here("analysis/data/derived_data/nordby1.RData"))

sitearea <- rast(file.path(siterpath,
                           paste0(str_replace(sitename, " ", "_"), ".tif")))

plot_results(sitename, output$sitel, output$datedat, sitearea, bmap,
             sites_sa, isobases, output, date_groups, scale_dist = 150,
             s_tdist = 0.5, s_xpos = 130, s_ypos = 70,  s_bheight = 0.25,
             searasterpath = searasterpath)

ggsave(file = here("analysis/figures/nordby1.png"), width = 250,
       height = plot_height + (plot_height * length(unique(date_groups))),
       units = "mm")

######### Nordby 52 #########
sitename <- "Nordby 52"
date_groups <- group_dates(rcarb_sa, sitename)

output <- apply_functions(sitename, date_groups, dtm, displacement_curves,
                          isobases, nsamp = nsamp, loc_bbox = 400, siterpath,
                          searasterpath)
save(output,
     file = here::here("analysis/data/derived_data/nordby52.RData"))

sitearea <- rast(file.path(siterpath,
                           paste0(str_replace(sitename, " ", "_"), ".tif")))

plot_results(sitename, output$sitel, output$datedat, sitearea, bmap,
             sites_sa, isobases, output, date_groups, scale_dist = 150,
             s_tdist = 2, s_xpos = 150, s_ypos = 80,  s_bheight = 1,
             searasterpath = searasterpath)

ggsave(file = here("analysis/figures/nordby52.png"), width = 250,
       height = plot_height + (plot_height * length(unique(date_groups))),
       units = "mm")

######### Pauler 1 #########
sitename <- "Pauler 1"
date_groups <- group_dates(rcarb_sa, sitename)

output <- apply_functions(sitename, date_groups, dtm, displacement_curves,
                          isobases, nsamp = nsamp, loc_bbox = 2500, siterpath,
                          searasterpath)
save(output,
     file = here::here("analysis/data/derived_data/pauler1.RData"))


sitearea <- rast(file.path(siterpath,
                           paste0(str_replace(sitename, " ", "_"), ".tif")))

plot_results(sitename, output$sitel, output$datedat, sitearea, bmap,
             sites_sa, isobases, output, date_groups, scale_dist = 1000,
             s_tdist = 2, s_xpos = 1200, s_ypos = 500,  s_bheight = 1.5,
             searasterpath = searasterpath)

ggsave(file = here("analysis/figures/pauler1.png"), width = 250,
       height = plot_height + (plot_height * length(unique(date_groups))),
       units = "mm")

######### Pauler 2 #########
sitename <- "Pauler 2"
date_groups <- group_dates(rcarb_sa, sitename)

output <- apply_functions(sitename, date_groups, dtm, displacement_curves,
                          isobases, nsamp = nsamp, loc_bbox = 2500, siterpath,
                          searasterpath)
save(output,
     file = here::here("analysis/data/derived_data/pauler2.RData"))

sitearea <- rast(file.path(siterpath,
                           paste0(str_replace(sitename, " ", "_"), ".tif")))

plot_results(sitename, output$sitel, output$datedat, sitearea, bmap,
             sites_sa, isobases, output, date_groups, scale_dist = 1000,
             s_tdist = 10, s_xpos = 1200, s_ypos = 500,  s_bheight = 5,
             searasterpath = searasterpath)

ggsave(file = here("analysis/figures/pauler2.png"), width = 250,
       height = plot_height + (plot_height * length(unique(date_groups))),
       units = "mm")

######### Pjonkerød R1 #########
sitename <- "Pjonkerød R1"
date_groups <- group_dates(rcarb_sa, sitename)

output <- apply_functions(sitename, date_groups, dtm, displacement_curves,
                          isobases, nsamp = nsamp, loc_bbox = 400, siterpath,
                          searasterpath)
save(output,
     file = here::here("analysis/data/derived_data/pjonkerodr1.RData"))

sitearea <- rast(file.path(siterpath,
                           paste0(str_replace(sitename, " ", "_"), ".tif")))

plot_results(sitename, output$sitel, output$datedat, sitearea, bmap,
             sites_sa, isobases, output, date_groups, scale_dist = 150,
             s_tdist = 1, s_xpos = 150, s_ypos = 80,  s_bheight = 0.5,
             searasterpath = searasterpath)

ggsave(file = here("analysis/figures/pjonkerodr1.png"), width = 250,
       height = plot_height + (plot_height * length(unique(date_groups))),
       units = "mm")

######### Prestemoen 1 #########
sitename <- "Prestemoen 1"
date_groups <- group_dates(rcarb_sa, sitename)

output <- apply_functions(sitename, date_groups, dtm, displacement_curves,
                          isobases, nsamp = nsamp, loc_bbox = 400, siterpath,
                          searasterpath)
save(output,
     file = here::here("analysis/data/derived_data/prestemoen1.RData"))

sitearea <- rast(file.path(siterpath,
                           paste0(str_replace(sitename, " ", "_"), ".tif")))

plot_results(sitename, output$sitel, output$datedat, sitearea, bmap,
             sites_sa, isobases, output, date_groups, scale_dist = 150,
             s_tdist = 1, s_xpos = 110, s_ypos = 50,  s_bheight = 0.5,
             searasterpath = searasterpath)

ggsave(file = here("analysis/figures/prestemoen1.png"), width = 250,
       height = plot_height + (plot_height * length(unique(date_groups))),
       units = "mm")

######### Ragnhildrød #########
sitename <- "Ragnhildrød"
date_groups <- group_dates(rcarb_sa, sitename)

output <- apply_functions(sitename, date_groups, dtm, displacement_curves,
                          isobases, nsamp = nsamp, loc_bbox = 4000, siterpath,
                          searasterpath)
save(output,
     file = here::here("analysis/data/derived_data/ragnhildrod.RData"))

sitearea <- rast(file.path(siterpath,
                           paste0(str_replace(sitename, " ", "_"), ".tif")))

plot_results(sitename, output$sitel, output$datedat, sitearea, bmap,
             sites_sa, isobases, output, date_groups, scale_dist = 1500,
             s_tdist = 5, s_xpos = 1500, s_ypos = 700,  s_bheight = 3,
             searasterpath = searasterpath)

ggsave(file = here("analysis/figures/ragnhildrod.png"), width = 250,
       height = plot_height + (plot_height * length(unique(date_groups))),
       units = "mm")

######### Rødbøl 54 #########
sitename <- "Rødbøl 54"
date_groups <- group_dates(rcarb_sa, sitename)

output <- apply_functions(sitename, date_groups, dtm, displacement_curves,
                          isobases, nsamp = nsamp, loc_bbox = 400, siterpath,
                          searasterpath)
save(output,
     file = here::here("analysis/data/derived_data/rodbol54.RData"))

sitearea <- rast(file.path(siterpath,
                           paste0(str_replace(sitename, " ", "_"), ".tif")))

plot_results(sitename, output$sitel, output$datedat, sitearea, bmap,
             sites_sa, isobases, output, date_groups, scale_dist = 150,
             s_tdist = 0.5, s_xpos = 150, s_ypos = 50,  s_bheight = 0.25,
             searasterpath = searasterpath)

ggsave(file = here("analysis/figures/rodbol54.png"), width = 250,
       height = plot_height + (plot_height * length(unique(date_groups))),
       units = "mm")

######### Rognlien #########
sitename <- "Rognlien"
date_groups <- group_dates(rcarb_sa, sitename)

output <- apply_functions(sitename, date_groups, dtm, displacement_curves,
                          isobases, nsamp = nsamp, loc_bbox = 400, siterpath,
                          searasterpath)
save(output,
     file = here::here("analysis/data/derived_data/rognlien.RData"))

sitearea <- rast(file.path(siterpath,
                           paste0(str_replace(sitename, " ", "_"), ".tif")))

plot_results(sitename, output$sitel, output$datedat, sitearea, bmap,
             sites_sa, isobases, output, date_groups, scale_dist = 150,
             s_tdist = 2, s_xpos = 150, s_ypos = 70, s_bheight = 1.5,
             searasterpath = searasterpath)

ggsave(file = here("analysis/figures/rognlien.png"), width = 250,
       height = plot_height + (plot_height * length(unique(date_groups))),
       units = "mm")

######### Sagene B1 #########
sitename <- "Sagene B1"
date_groups <- group_dates(rcarb_sa, sitename)

output <- apply_functions(sitename, date_groups, dtm, displacement_curves,
                          isobases, nsamp = nsamp, loc_bbox = 1000, siterpath,
                          searasterpath)
save(output,
     file = here::here("analysis/data/derived_data/sageneb1.RData"))

sitearea <- rast(file.path(siterpath,
                           paste0(str_replace(sitename, " ", "_"), ".tif")))

plot_results(sitename, output$sitel, output$datedat, sitearea, bmap,
             sites_sa, isobases, output, date_groups, scale_dist = 500,
             s_tdist = 1.5, s_xpos = 400, s_ypos = 150, s_bheight = 1,
             searasterpath = searasterpath)

ggsave(file = here("analysis/figures/sageneb1.png"), width = 250,
       height = plot_height + (plot_height * length(unique(date_groups))),
       units = "mm")

######### Sagene B2 #########
sitename <- "Sagene B2"
date_groups <- group_dates(rcarb_sa, sitename)

output <- apply_functions(sitename, date_groups, dtm, displacement_curves,
                          isobases, nsamp = nsamp, loc_bbox = 1000, siterpath,
                          searasterpath)
save(output,
     file = here::here("analysis/data/derived_data/sageneb2.RData"))

sitearea <- rast(file.path(siterpath,
                           paste0(str_replace(sitename, " ", "_"), ".tif")))

plot_results(sitename, output$sitel, output$datedat, sitearea, bmap,
             sites_sa, isobases, output, date_groups, scale_dist = 500,
             s_tdist = 1, s_xpos = 300, s_ypos = 120, s_bheight = 0.5,
             searasterpath = searasterpath)

ggsave(file = here("analysis/figures/sageneb2.png"), width = 250,
       height = plot_height + (plot_height * length(unique(date_groups))),
       units = "mm")

######### Solum 3 #########
sitename <- "Solum 3"
date_groups <- group_dates(rcarb_sa, sitename)

output <- apply_functions(sitename, date_groups, dtm, displacement_curves,
                          isobases, nsamp = nsamp, loc_bbox = 3500, siterpath,
                          searasterpath)
save(output,
     file = here::here("analysis/data/derived_data/solum3.RData"))

sitearea <- rast(file.path(siterpath,
                           paste0(str_replace(sitename, " ", "_"), ".tif")))

plot_results(sitename, output$sitel, output$datedat, sitearea, bmap,
             sites_sa, isobases, output, date_groups, scale_dist = 1000,
             s_tdist = 15, s_xpos = 1200, s_ypos = 500,  s_bheight = 8,
             searasterpath = searasterpath)

ggsave(file = here("analysis/figures/solum3.png"), width = 250,
       height = plot_height + (plot_height * length(unique(date_groups))),
       units = "mm")

######### Stokke/Polland 1 #########
sitename <- "Stokke/Polland 1"

# Typological indicators indicate seperate Late Mesolithic and a Early Neolithic
# visits. Date groups were manually set to reflect this.
date_groups <- c(1, 2, 2, 2)

output <- apply_functions(sitename, date_groups, dtm, displacement_curves,
                          isobases, nsamp = nsamp, loc_bbox = 400, siterpath,
                          searasterpath)
save(output,
     file = here::here("analysis/data/derived_data/stokkep1.RData"))

sitearea <- rast(file.path(siterpath,
                           paste0(str_replace(str_replace(sitename, " ", "_"),
                                              "/", "_"), ".tif")))

plot_results(sitename, output$sitel, output$datedat, sitearea, bmap,
             sites_sa, isobases, output, date_groups, scale_dist = 150,
             s_tdist = 0.5, s_xpos = 150, s_ypos = 60, s_bheight = 0.25,
             searasterpath = searasterpath)

ggsave(file = here("analysis/figures/stokkep1.png"), width = 250,
       height = plot_height + (plot_height * length(unique(date_groups))),
       units = "mm")

######### Stokke/Polland 5 #########
sitename <- "Stokke/Polland 5"
date_groups <- group_dates(rcarb_sa, sitename)

output <- apply_functions(sitename, date_groups, dtm, displacement_curves,
                          isobases, nsamp = nsamp, loc_bbox = 400, siterpath,
                          searasterpath)
save(output,
     file = here::here("analysis/data/derived_data/stokkep5.RData"))

sitearea <- rast(file.path(siterpath,
                           paste0(str_replace(str_replace(sitename, " ", "_"),
                                                        "/", "_"), ".tif")))

plot_results(sitename, output$sitel, output$datedat, sitearea, bmap,
             sites_sa, isobases, output, date_groups, scale_dist = 150,
             s_tdist = 1, s_xpos = 150, s_ypos = 60, s_bheight = 1,
             searasterpath = searasterpath)

ggsave(file = here("analysis/figures/stokkep5.png"), width = 250,
       height = plot_height + (plot_height * length(unique(date_groups))),
       units = "mm")

######### Stokke/Polland 8 #########
sitename <- "Stokke/Polland 8"
date_groups <- group_dates(rcarb_sa, sitename)

output <- apply_functions(sitename, date_groups, dtm, displacement_curves,
                          isobases, nsamp = nsamp, loc_bbox = 400, siterpath,
                          searasterpath)
save(output,
     file = here::here("analysis/data/derived_data/stokkep8.RData"))

sitearea <- rast(file.path(siterpath,
                           paste0(str_replace(str_replace(sitename, " ", "_"),
                                                        "/", "_"), ".tif")))

plot_results(sitename, output$sitel, output$datedat, sitearea, bmap,
             sites_sa, isobases, output, date_groups, scale_dist = 150,
             s_tdist = 0.5, s_xpos = 150, s_ypos = 60, s_bheight = 0.25,
             searasterpath = searasterpath)

ggsave(file = here("analysis/figures/stokkep8.png"), width = 250,
       height = plot_height + (plot_height * length(unique(date_groups))),
       units = "mm")

######### Torstvet #########
sitename <- "Torstvet"
date_groups <- group_dates(rcarb_sa, sitename)

output <- apply_functions(sitename, date_groups, dtm, displacement_curves,
                          isobases, nsamp = nsamp, loc_bbox = 1000, siterpath,
                          searasterpath)
save(output,
     file = here::here("analysis/data/derived_data/torstvet.RData"))

sitearea <- rast(file.path(siterpath,
                           paste0(str_replace(sitename, " ", "_"), ".tif")))

plot_results(sitename, output$sitel, output$datedat, sitearea, bmap,
             sites_sa, isobases, output, date_groups, scale_dist = 150,
             s_tdist = 1, s_xpos = 150, s_ypos = 60, s_bheight = 1,
             searasterpath = searasterpath)

ggsave(file = here("analysis/figures/torstvet.png"), width = 250,
       height = plot_height + (plot_height * length(unique(date_groups))),
       units = "mm")

######### Trondalen #########
sitename <- "Trondalen"
date_groups <- group_dates(rcarb_sa, sitename)

output <- apply_functions(sitename, date_groups, dtm, displacement_curves,
                          isobases, nsamp = nsamp, loc_bbox = 400, siterpath,
                          searasterpath)
save(output,
     file = here::here("analysis/data/derived_data/trondalen.RData"))

sitearea <- rast(file.path(siterpath,
                           paste0(str_replace(sitename, " ", "_"), ".tif")))

plot_results(sitename, output$sitel, output$datedat, sitearea, bmap,
             sites_sa, isobases, output, date_groups, scale_dist = 150,
             s_tdist = 1, s_xpos = 150, s_ypos = 60, s_bheight = 1,
             searasterpath = searasterpath)

ggsave(file = here("analysis/figures/trondalen.png"), width = 250,
       height = plot_height + (plot_height * length(unique(date_groups))),
       units = "mm")

######### Tverdal #########
sitename <- "Tverdal"
date_groups <- group_dates(rcarb_sa, sitename)

output <- apply_functions(sitename, date_groups, dtm, displacement_curves,
                          isobases, nsamp = nsamp, loc_bbox = 400, siterpath,
                          searasterpath)
save(output,
     file = here::here("analysis/data/derived_data/tverdal.RData"))

sitearea <- rast(file.path(siterpath,
                           paste0(str_replace(sitename, " ", "_"), ".tif")))

plot_results(sitename, output$sitel, output$datedat, sitearea, bmap,
             sites_sa, isobases, output, date_groups, scale_dist = 150,
             s_tdist = 1, s_xpos = 175, s_ypos = 85,  s_bheight = 1,
             searasterpath = searasterpath)

ggsave(file = here("analysis/figures/tverdal.png"), width = 250,
       height = plot_height + (plot_height * length(unique(date_groups))),
       units = "mm")

######### Vallermyrene 1a #########
sitename <- "Vallermyrene 1a"
date_groups <- group_dates(rcarb_sa, sitename)

output <- apply_functions(sitename, date_groups, dtm, displacement_curves,
                           isobases, nsamp = nsamp, loc_bbox = 400, siterpath,
                          searasterpath)
save(output,
     file = here::here("analysis/data/derived_data/vallermyrene1a.RData"))

sitearea <- rast(file.path(siterpath,
                           paste0(str_replace(sitename, " ", "_"), ".tif")))

plot_results(sitename, output$sitel, output$datedat, sitearea, bmap,
             sites_sa, isobases, output, date_groups, scale_dist = 150,
             s_tdist = 3, s_xpos = 175, s_ypos = 85,  s_bheight = 1.5,
             searasterpath = searasterpath)

ggsave(file = here("analysis/figures/vallermyrene1a.png"), width = 250,
       height = plot_height + (plot_height * length(unique(date_groups))),
       units = "mm")

######### Vallermyrene 1b #########
sitename <- "Vallermyrene 1b"
date_groups <- group_dates(rcarb_sa, sitename)

output <- apply_functions(sitename, date_groups, dtm, displacement_curves,
                          isobases, nsamp = nsamp, loc_bbox = 400, siterpath,
                          searasterpath)
save(output,
     file = here::here("analysis/data/derived_data/vallermyrene1b.RData"))

sitearea <- rast(file.path(siterpath,
                           paste0(str_replace(sitename, " ", "_"), ".tif")))

plot_results(sitename, output$sitel, output$datedat, sitearea, bmap,
             sites_sa, isobases, output, date_groups, scale_dist = 150,
             s_tdist = 2, s_xpos = 175, s_ypos = 85,  s_bheight = 1,
             searasterpath = searasterpath)

ggsave(file = here("analysis/figures/vallermyrene1b.png"), width = 250,
       height = plot_height + (plot_height * length(unique(date_groups))),
       units = "mm")

######### Vallermyrene 2 #########
sitename <- "Vallermyrene 2"
date_groups <- group_dates(rcarb_sa, sitename)

# Use dated feature instead of site limit
output <- apply_functions(sitename, date_groups, dtm, displacement_curves,
                          isobases, nsamp = nsamp, loc_bbox = 200, siterpath,
                          sitelimit = FALSE, searasterpath = searasterpath)
save(output,
     file = here::here("analysis/data/derived_data/vallermyrene2.RData"))

sitearea <- rast(file.path(siterpath,
                           paste0(str_replace(sitename, " ", "_"), ".tif")))

plot_results(sitename, output$sitel, output$datedat, sitearea, bmap,
             sites_sa, isobases, output, date_groups, scale_dist = 100,
             s_tdist = 20, s_xpos = 85, s_ypos = 40,  s_bheight = 10,
             searasterpath = searasterpath)

ggsave(file = here("analysis/figures/vallermyrene2.png"), width = 250,
       height = plot_height + (plot_height * length(unique(date_groups))),
       units = "mm")

######### Vallermyrene 4a #########
sitename <- "Vallermyrene 4a"
date_groups <- group_dates(rcarb_sa, sitename)

output <- apply_functions(sitename, date_groups, dtm, displacement_curves,
                           isobases, nsamp = nsamp, loc_bbox = 400, siterpath,
                          searasterpath)
save(output,
     file = here::here("analysis/data/derived_data/vallermyrene4a.RData"))

sitearea <- rast(file.path(siterpath,
                           paste0(str_replace(sitename, " ", "_"), ".tif")))

plot_results(sitename, output$sitel, output$datedat, sitearea, bmap,
             sites_sa, isobases, output, date_groups, scale_dist = 100,
             s_tdist = 1, s_xpos = 175, s_ypos = 85,  s_bheight = 0.5,
             searasterpath = searasterpath)

ggsave(file = here("analysis/figures/vallermyrene4a.png"), width = 250,
       height = plot_height + (plot_height * length(unique(date_groups))),
       units = "mm")

######### Vallermyrene 4b #########
sitename <- "Vallermyrene 4b"
date_groups <- group_dates(rcarb_sa, sitename)

output <- apply_functions(sitename, date_groups, dtm, displacement_curves,
                          isobases, nsamp = nsamp, loc_bbox = 400, siterpath,
                          searasterpath)
save(output,
     file = here::here("analysis/data/derived_data/vallermyrene4b.RData"))

sitearea <- rast(file.path(siterpath,
                           paste0(str_replace(sitename, " ", "_"), ".tif")))

plot_results(sitename, output$sitel, output$datedat, sitearea, bmap,
             sites_sa, isobases, output, date_groups, scale_dist = 150,
             s_tdist = 1.5, s_xpos = 175, s_ypos = 85,  s_bheight = 1,
             searasterpath = searasterpath)

ggsave(file = here("analysis/figures/vallermyrene4b.png"), width = 250,
       height = plot_height + (plot_height * length(unique(date_groups))),
       units = "mm")

# Delete all site area rasters stored in sitearea_temp
do.call(file.remove, list(list.files(siterpath, full.names = TRUE)))
