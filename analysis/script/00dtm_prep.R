library(raster)

# As the raster data is very large, only the edited version is distributed with
# the research compendium as a series of raster tiles. These have to be merged
# for the rest of the scripts to run.

# The unedited DTM is not available here, but the polygon data used to edit the
# raster are, and the script run to do the editing is available in 03dtm_edit.R,
# making it possible to re-run the editing by retrieving raster data for the
# study area from hoydedata.no

# Get path to DTM tiles
dtmtiles <- list.files(here::here("analysis/data/raw_data/tiled_dtm"),
                       pattern= ".*.tif$", full.names = TRUE)

# Load the tiles and merge
tiles <- sapply(dtmtiles, raster::raster)
names(tiles)[1:2] <- c('x', 'y')
dtm <- do.call(raster::merge, tiles)

# Save merged raster to the derived data folder
dtmpath <-  here::here("analysis/data/derived_data/dtm")
edrast <- file.path(dtmpath, "dtm10.tif")
writeRaster(dtm, edrast,
            overwrite = TRUE)
