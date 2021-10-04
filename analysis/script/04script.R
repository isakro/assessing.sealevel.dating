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
# Path to smaller site area rasters
siterpath <- "/home/isak/phd/eaa_presentation/sitearea"


# Assign closest isobases to site data.
sites_sa <- st_join(st_make_valid(sites_sa), isopolys,
                    join = st_intersects, largest = TRUE)

# Both the digital terrain model and the radiocarbon dates required
# manual inspection in case of inconsistencies (i.e. a highway running by the
# site), defining necessary size of the window of analysis, and/or to decide on
# grouping of dates. Each site was therefore first simulated 50 times and
# then rerun at 1000 when these issues were handled.



######### Vallermyrene 1a #########
sitename <- "Vallermyrene 1a"
date_groups <- c(1, 1)

# output <- apply_functions(sitename, date_groups, dtmpath, displacement_curves,
#                           isobases, nsamp = 1000, loc_bbox = 250)
load(here("analysis/data/derived_data/vallermyrene1a.RData"))

sitearea <- rast(file.path(siterpath,
                           paste0(str_replace(sitename, " ", "_"), ".tif")))

plot_results(sitename, output$sitel, output$datedat, sitearea, bmap,
             sites_sa, isobases, output, date_groups, scale_dist = 100)

ggsave(file = here("analysis/figures/vallermyrene1a.png"), width = 250,
       height = 152, units = "mm")

######### Vallermyrene 1b #########
sitename <- "Vallermyrene 1b"
date_groups <- 1

# output <- apply_functions(sitename, date_groups, dtmpath, displacement_curves,
#                           isobases, nsamp = 1000, loc_bbox = 250)
load(here("analysis/data/derived_data/vallermyrene1b.RData"))

sitearea <- rast(file.path(siterpath,
          paste0(str_replace(sitename, " ", "_"), ".tif")))

plot_results(sitename, output$sitel, output$datedat, sitearea, bmap,
             sites_sa, isobases, output, date_groups, scale_dist = 100)

ggsave(file = here("analysis/figures/vallermyrene1b.png"), width = 250,
       height = 152, units = "mm")

######### Vallermyrene 4a #########
sitename <- "Vallermyrene 4a"
date_groups <- c(1, 1)

# output <- apply_functions(sitename, date_groups, dtmpath, displacement_curves,
#                           isobases, nsamp = 1000, loc_bbox = 250, siterpath)
load(here("analysis/data/derived_data/vallermyrene4a.RData"))

sitearea <- rast(file.path(siterpath,
                           paste0(str_replace(sitename, " ", "_"), ".tif")))

plot_results(sitename, output$sitel, output$datedat, sitearea, bmap,
             sites_sa, isobases, output, date_groups, scale_dist = 100)

ggsave(file = here("analysis/figures/vallermyrene4a.png"), width = 250,
       height = 152, units = "mm")
######### Tverdal #########
sitename <- "Tverdal"
date_groups <- rep(1, 4)

# output <- apply_functions(sitename, date_groups, dtmpath, displacement_curves,
#                           isobases, nsamp = 1000, loc_bbox = 250, siterpath)
load(here("analysis/data/derived_data/tverdal.RData"))

sitearea <- rast(file.path(siterpath,
           paste0(str_replace(sitename, " ", "_"), ".tif")))

plot_results(sitename, output$sitel, output$datedat, sitearea, bmap,
             sites_sa, isobases, output, date_groups, scale_dist = 100,
             s_tdist = 0.5, s_xpos = 175, s_ypos = 85,  s_bheight = 0.25)

ggsave(file = here("analysis/figures/tverdal.png"), width = 250,
       height = 152, units = "mm")

######### Kvastad A1 #########
sitename <- "Kvastad A1"
date_groups <- 1

load(here("analysis/data/derived_data/kvastada1.RData"))
# output <- apply_functions(sitename, date_groups, dtmpath, displacement_curves,
#                           isobases, nsamp = 1000, loc_bbox = 250, siterpath)

sitearea <- rast(file.path(siterpath,
                           paste0(str_replace(sitename, " ", "_"), ".tif")))

plot_results(sitename, output$sitel, output$datedat, sitearea, bmap,
             sites_sa, isobases, output, date_groups, scale_dist = 100,
             s_tdist = 0.5, s_xpos = 135, s_ypos = 65,  s_bheight = 0.25)

ggsave(file = here("analysis/figures/kvastada1.png"), width = 250,
       height = 152, units = "mm")

######### Kvastad A2 #########
sitename <- "Kvastad A2"
date_groups <- c(1, 1, 1, 2, 2)

# output <- apply_functions(sitename, date_groups, dtmpath, displacement_curves,
#                          isobases, nsamp = 1000, loc_bbox = 500, siterpath)
load(here("analysis/data/derived_data/kvastada2.RData"))

sitearea <- rast(file.path(siterpath,
                           paste0(str_replace(sitename, " ", "_"), ".tif")))

plot_results(sitename, output$sitel, output$datedat, sitearea, bmap,
             sites_sa, isobases, output, date_groups, scale_dist = 150)

ggsave(file = here("analysis/figures/kvastad_a2.png"), width = 250,
       height = 228, units = "mm")
######### Kvastad A4 #########
sitename <- "Kvastad A4"
date_groups <- 1

load(here("analysis/data/derived_data/kvastada4.RData"))
# output <- apply_functions(sitename, date_groups, dtmpath, displacement_curves,
#                           isobases, nsamp = 1000, loc_bbox = 250, siterpath)

sitearea <- rast(file.path(siterpath,
                           paste0(str_replace(sitename, " ", "_"), ".tif")))

plot_results(sitename, output$sitel, output$datedat, sitearea, bmap,
             sites_sa, isobases, output, date_groups, scale_dist = 100,
             s_tdist = 0.5, s_xpos = 135, s_ypos = 65,  s_bheight = 0.25)

ggsave(file = here("analysis/figures/kvastada4.png"), width = 250,
       height = 152, units = "mm")

######### Kvastad A9 #########
sitename <- "Kvastad A9"
date_groups <- 1

output <- apply_functions(sitename, date_groups, dtmpath, displacement_curves,
                          isobases, nsamp = 1000, loc_bbox = 250, siterpath)
save(output,
     file = here::here("analysis/data/derived_data/kvastada9.RData"))

sitearea <- rast(file.path(siterpath,
                           paste0(str_replace(sitename, " ", "_"), ".tif")))

plot_results(sitename, output$sitel, output$datedat, sitearea, bmap,
             sites_sa, isobases, output, date_groups, scale_dist = 100,
             s_tdist = 3.5, s_xpos = 100, s_ypos = 40,  s_bheight = 2.5)

ggsave(file = here("analysis/figures/kvastada9.png"), width = 250,
       height = 152, units = "mm")


######### Langangen Vestgård 1 #########
sitename <- "Langangen Vestgård 1"
date_groups <- rep(1, 13)

# output <- apply_functions(sitename, date_groups, dtmpath, displacement_curves,
#                           isobases, nsamp = 1000, loc_bbox = 250, siterpath)
load(here("analysis/data/derived_data/langangenv1.RData"))

sitearea <- rast(file.path(siterpath,
                           paste0(str_replace(sitename, " ", "_"), ".tif")))

plot_results(sitename, output$sitel, output$datedat, sitearea, bmap,
             sites_sa, isobases, output, date_groups, scale_dist = 100)

ggsave(file = here("analysis/figures/langangenv1.png"), width = 250,
       height = 152, units = "mm")
######### Langangen Vestgård 3 #########
sitename <- "Langangen Vestgård 3"
date_groups <- rep(1, 6)

# output <- apply_functions(sitename, date_groups, dtmpath, displacement_curves,
#                           isobases, nsamp = 1000, loc_bbox = 250, siterpath)
load(here("analysis/data/derived_data/langangenv3.RData"))

sitearea <- rast(file.path(siterpath,
                           paste0(str_replace(sitename, " ", "_"), ".tif")))

plot_results(sitename, output$sitel, output$datedat, sitearea, bmap,
             sites_sa, isobases, output, date_groups, scale_dist = 100,
             s_tdist = 1, s_xpos = 110, s_ypos = 50,  s_bheight = 0.5)

ggsave(file = here("analysis/figures/langangenv3.png"), width = 250,
       height = 152, units = "mm")
######### Langangen Vestgård 5 #########
sitename <- "Langangen Vestgård 5"
date_groups <- c(1, 1, 2, 2, 2)

# output <- apply_functions(sitename, date_groups, dtmpath, displacement_curves,
#                           isobases, nsamp = 1000, loc_bbox = 250, siterpath)
load(here("analysis/data/derived_data/langangenv5.RData"))

sitearea <- rast(file.path(siterpath,
                           paste0(str_replace(sitename, " ", "_"), ".tif")))

plot_results(sitename, output$sitel, output$datedat, sitearea, bmap,
             sites_sa, isobases, output, date_groups, scale_dist = 100)

ggsave(file = here("analysis/figures/langangenv5.png"), width = 250,
       height = 228, units = "mm")

######### Langangen Vestgård 6 #########
sitename <- "Langangen Vestgård 6"
date_groups <- rep(1, 9)

# output <- apply_functions(sitename, date_groups, dtmpath, displacement_curves,
#                           isobases, nsamp = 1000, loc_bbox = 250, siterpath,
#                           rcarbcor_true = TRUE)
load(here("analysis/data/derived_data/langangenv6.RData"))

sitearea <- rast(file.path(siterpath,
                           paste0(str_replace(sitename, " ", "_"), ".tif")))

plot_results(sitename, output$sitel, output$datedat, sitearea, bmap,
             sites_sa, isobases, output, date_groups, scale_dist = 100,
             s_tdist = 0.5, s_xpos = 110, s_ypos = 50,  s_bheight = 0.3)

ggsave(file = here("analysis/figures/langangenv6.png"), width = 250,
       height = 152, units = "mm")
######### Rognlien #########
sitename <- "Rognlien"
date_groups <- 1

# output <- apply_functions(sitename, date_groups, dtmpath, displacement_curves,
#                           isobases, nsamp = 100, loc_bbox = 250, siterpath,
#                           rcarbcor_true = TRUE)
load(here("analysis/data/derived_data/rognlien.RData"))

sitearea <- rast(file.path(siterpath,
                           paste0(str_replace(sitename, " ", "_"), ".tif")))

plot_results(sitename, output$sitel, output$datedat, sitearea, bmap,
             sites_sa, isobases, output, date_groups, scale_dist = 100,
             s_tdist = 1, s_xpos = 110, s_ypos = 50,  s_bheight = 0.5)

ggsave(file = here("analysis/figures/rognlien.png"), width = 250,
       height = 152, units = "mm")

######### Dybdalshei 2 #########
sitename <- "Dybdalshei 2"
date_groups <- 1

# output <- apply_functions(sitename, date_groups, dtmpath, displacement_curves,
#                           isobases, nsamp = 1000, loc_bbox = 250, siterpath)
load(here("analysis/data/derived_data/dybdalshei2.RData"))

sitearea <- rast(file.path(siterpath,
                           paste0(str_replace(sitename, " ", "_"), ".tif")))

plot_results(sitename, output$sitel, output$datedat, sitearea, bmap,
             sites_sa, isobases, output, date_groups, scale_dist = 100)

ggsave(file = here("analysis/figures/dybdalshei2.png"), width = 250,
       height = 152, units = "mm")

######### Dybdalshei 1 #########
sitename <- "Dybdalshei 1"
date_groups <- c(1, 1, 1, 1, 2, 2, 3, 3)

# output <- apply_functions(sitename, date_groups, dtmpath, displacement_curves,
#                           isobases, nsamp = 1000, loc_bbox = 250, siterpath)
load(here("analysis/data/derived_data/dybdalshei1.RData"))

sitearea <- rast(file.path(siterpath,
                           paste0(str_replace(sitename, " ", "_"), ".tif")))

plot_results(sitename, output$sitel, output$datedat, sitearea, bmap,
             sites_sa, isobases, output, date_groups, scale_dist = 100,
             s_tdist = 0.5, s_xpos = 135, s_ypos = 65,  s_bheight = 0.25)

ggsave(file = here("analysis/figures/dybdalshei1.png"), width = 250,
       height = 304, units = "mm")

######### Gunnarsrød 4. Need to edit raster #########
sitename <- "Gunnarsrød 4"
date_groups <- c(1, 1)

output <- apply_functions(sitename, date_groups, dtmpath, displacement_curves,
                          isobases, nsamp = 1000, loc_bbox = 250, siterpath)
save(output,
     file = here::here("analysis/data/derived_data/gunnarsrod4.RData"))

sitearea <- rast(file.path(siterpath,
                           paste0(str_replace(sitename, " ", "_"), ".tif")))

plot_results(sitename, output$sitel, output$datedat, sitearea, bmap,
             sites_sa, isobases, output, date_groups, scale_dist = 100,
             s_tdist = 0.5, s_xpos = 135, s_ypos = 65,  s_bheight = 0.25)

ggsave(file = here("analysis/figures/gunnarsrod4.png"), width = 250,
       height = 152, units = "mm")

######### Gunnarsrød 5 Need to edit raster #########
sitename <- "Gunnarsrød 5"
date_groups <- c(1, 1, 2)

output <- apply_functions(sitename, date_groups, dtmpath, displacement_curves,
                          isobases, nsamp = 1000, loc_bbox = 400, siterpath,
                          rcarbcor_true = TRUE)
save(output,
     file = here::here("analysis/data/derived_data/gunnarsrod5.RData"))

sitearea <- rast(file.path(siterpath,
                           paste0(str_replace(sitename, " ", "_"), ".tif")))

plot_results(sitename, output$sitel, output$datedat, sitearea, bmap,
             sites_sa, isobases, output, date_groups, scale_dist = 100,
             s_tdist = 0.5, s_xpos = 135, s_ypos = 65,  s_bheight = 0.25)

ggsave(file = here("analysis/figures/gunnarsrod5.png"), width = 250,
       height = 152, units = "mm")

######### Langangen Vestgård 7 Road as site limit #########
sitename <- "Langangen Vestgård 7"
date_groups <- 1

output <- apply_functions(sitename, date_groups, dtmpath, displacement_curves,
                          isobases, nsamp = 1000, loc_bbox = 250, siterpath)
save(output,
     file = here::here("analysis/data/derived_data/langangenv7.RData"))

sitearea <- rast(file.path(siterpath,
                           paste0(str_replace(sitename, " ", "_"), ".tif")))

plot_results(sitename, output$sitel, output$datedat, sitearea, bmap,
             sites_sa, isobases, output, date_groups, scale_dist = 100,
             s_tdist = 1, s_xpos = 110, s_ypos = 50,  s_bheight = 0.5)

ggsave(file = here("analysis/figures/langangenv7.png"), width = 250,
       height = 152, units = "mm")

######### Hesthag C2 Two individual dates #########
sitename <- "Hesthag C2"
date_groups <- c(1, 2)

output <- apply_functions(sitename, date_groups, dtmpath, displacement_curves,
                          isobases, nsamp = 1000, loc_bbox = 250, siterpath)
save(output,
     file = here::here("analysis/data/derived_data/hesthagc2.RData"))

sitearea <- rast(file.path(siterpath,
                           paste0(str_replace(sitename, " ", "_"), ".tif")))

plot_results(sitename, output$sitel, output$datedat, sitearea, bmap,
             sites_sa, isobases, output, date_groups, scale_dist = 100)

ggsave(file = here("analysis/figures/hesthagc2.png"), width = 250,
       height = 152, units = "mm")

######### Hegna vest 1 Need to edit raster #########
sitename <- "Hegna vest 1"
date_groups <- c(1, 1)

output <- apply_functions(sitename, date_groups, dtmpath, displacement_curves,
                          isobases, nsamp = 1000, loc_bbox = 250, siterpath,
                          rcarbcor_true = TRUE)
save(output,
     file = here::here("analysis/data/derived_data/hegnavest1.RData"))

sitearea <- rast(file.path(siterpath,
                           paste0(str_replace(sitename, " ", "_"), ".tif")))

plot_results(sitename, output$sitel, output$datedat, sitearea, bmap,
             sites_sa, isobases, output, date_groups, scale_dist = 100)

ggsave(file = here("analysis/figures/hegnavest1.png"), width = 250,
       height = 152, units = "mm")

######### Hegna vest 3 Need to edit raster #########
sitename <- "Hegna vest 3"
date_groups <- 1

output <- apply_functions(sitename, date_groups, dtmpath, displacement_curves,
                          isobases, nsamp = 1000, loc_bbox = 250, siterpath)
save(output,
     file = here::here("analysis/data/derived_data/hegnavest3.RData"))

sitearea <- rast(file.path(siterpath,
                           paste0(str_replace(sitename, " ", "_"), ".tif")))

plot_results(sitename, output$sitel, output$datedat, sitearea, bmap,
             sites_sa, isobases, output, date_groups, scale_dist = 100)

ggsave(file = here("analysis/figures/hegnavest3.png"), width = 250,
       height = 152, units = "mm")

######### Hesthag C4 Need to edit raster #########
sitename <- "Stokke/Polland 8"
date_groups <- 1

output <- apply_functions(sitename, date_groups, dtmpath, displacement_curves,
                          isobases, nsamp = 1000, loc_bbox = 250, siterpath,
                          rcarbcor_true = TRUE)
save(output,
     file = here::here("analysis/data/derived_data/stokkepolland8.RData"))

sitearea <- rast(file.path(siterpath,
                           paste0(str_replace(sitename, " ", "_"), ".tif")))

plot_results(sitename, output$sitel, output$datedat, sitearea, bmap,
             sites_sa, isobases, output, date_groups, scale_dist = 100)

ggsave(file = here("analysis/figures/stokkepolland8.png"), width = 250,
       height = 152, units = "mm")

######### Vallermyrene 4b Need to edit raster #########
sitename <- "Vallermyrene 4b"
date_groups <- c(1, 1)

output <- apply_functions(sitename, date_groups, dtmpath, displacement_curves,
                          isobases, nsamp = 1000, loc_bbox = 250, siterpath)
save(output,
     file = here::here("analysis/data/derived_data/vallermyrene4b.RData"))

sitearea <- rast(file.path(siterpath,
                           paste0(str_replace(sitename, " ", "_"), ".tif")))

plot_results(sitename, output$sitel, output$datedat, sitearea, bmap,
             sites_sa, isobases, output, date_groups, scale_dist = 100)

ggsave(file = here("analysis/figures/vallermyrene4b.png"), width = 250,
       height = 152, units = "mm")
######### Solum 3 #########
sitename <- "Solum 3"
date_groups <- 1

output <- apply_functions(sitename, date_groups, dtmpath, displacement_curves,
                          isobases, nsamp = 1000, loc_bbox = 250, siterpath,
                          rcarbcor_true = TRUE)
save(output,
     file = here::here("analysis/data/derived_data/solum3.RData"))

sitearea <- rast(file.path(siterpath,
                           paste0(str_replace(sitename, " ", "_"), ".tif")))

plot_results(sitename, output$sitel, output$datedat, sitearea, bmap,
             sites_sa, isobases, output, date_groups, scale_dist = 100)

ggsave(file = here("analysis/figures/solum3.png"), width = 250,
       height = 152, units = "mm")


######### Langemyr Need to solve raster loading #########
sitename <- "Langemyr"
date_groups <- 1

output <- apply_functions(sitename, date_groups, dtmpath, displacement_curves,
                          isobases, nsamp = 1000, loc_bbox = 250, siterpath)
save(output,
     file = here::here("analysis/data/derived_data/langemyr.RData"))

sitearea <- rast(file.path(siterpath,
                           paste0(str_replace(sitename, " ", "_"), ".tif")))

plot_results(sitename, output$sitel, output$datedat, sitearea, bmap,
             sites_sa, isobases, output, date_groups, scale_dist = 100)

ggsave(file = here("analysis/figures/langemyr.png"), width = 250,
       height = 152, units = "mm")

######### Hovland 5 Need to solve raster loading #########
sitename <- "Hovland 5"
date_groups <- 1

output <- apply_functions(sitename, date_groups, dtmpath, displacement_curves,
                          isobases, nsamp = 1000, loc_bbox = 250, siterpath)
save(output,
     file = here::here("analysis/data/derived_data/hovland5.RData"))

sitearea <- rast(file.path(siterpath,
                           paste0(str_replace(sitename, " ", "_"), ".tif")))

plot_results(sitename, output$sitel, output$datedat, sitearea, bmap,
             sites_sa, isobases, output, date_groups, scale_dist = 100)

ggsave(file = here("analysis/figures/hovland5.png"), width = 250,
       height = 152, units = "mm")

######### Vallermyrene 2 Perhaps not use the site limit, but rather dated feature #########
sitename <- "Vallermyrene 2"
date_groups <- 1

output <- apply_functions(sitename, date_groups, dtmpath, displacement_curves,
                          isobases, nsamp = 1000, loc_bbox = 250, siterpath)
save(output,
     file = here::here("analysis/data/derived_data/vallermyrene1.RData"))

sitearea <- rast(file.path(siterpath,
                           paste0(str_replace(sitename, " ", "_"), ".tif")))

plot_results(sitename, output$sitel, output$datedat, sitearea, bmap,
             sites_sa, isobases, output, date_groups, scale_dist = 100)

ggsave(file = here("analysis/figures/vallermyrene1.png"), width = 250,
       height = 152, units = "mm")

