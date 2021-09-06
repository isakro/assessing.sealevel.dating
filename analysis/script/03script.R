library(tidyverse)
library(sf)
library(tmap)
library(tmaptools)
library(cowplot)
library(terra)
library(oxcAAR)
library(topoDistance)
library(vapour)
quickSetupOxcal()

source(here::here("analysis/script/02script.R"))
load(here::here("analysis/data/derived_data/01data.RData"))
load(here::here("analysis/data/derived_data/02data.RData"))


dtmpath <- "/home/isak/phd/eaa_presentation/dtm1/data/dtm1_33_120_109.tif"

rinfo <- vapour_raster_info(dtmpath)




buf <- sites_sa %>%
  filter(st_is_empty(.) == FALSE) %>%
  st_buffer(dist = 1000)

# hoydedata.no only accepts 10 or less polygons per shape.

tm_shape(st_union(buf, by_feature = FALSE)) +
  tm_polygons()

buffer <- st_union(buf, by_feature = FALSE)
st_write(buf[1,], here::here("analysis/buffer1.shp"))
nrow(buf)


