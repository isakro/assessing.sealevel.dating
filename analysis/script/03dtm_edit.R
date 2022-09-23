library(terra)
library(sf)
library(rgrass7)

# Script to remove roads and other disturbances in the DTM, and interpolate
# new values. Polygon data with which to remove these features has been
# created manually based on an inspection of satellite data and the DTM.

# This path has to be set to the local installation of GRASS
# grasspath <- "/usr/lib/grass78"

try(!is.na(grasspath)) # Uncomment and det path to local GRASS installation
                       #  above or this will throw and error.

# Due to the file size (c. 1.6 GB), the DTM is stored in 400 tiles that are
# first merged and then saved at this destination
dtmpath <-  here::here("analysis/data/derived_data")

# Get path to DTM tiles
dtmtiles <- list.files(here::here("analysis/data/raw_data/tiled_dtm10"),
                       pattern= ".*.tif$", full.names = TRUE)

# Load the tiles and merge
tiles <- sapply(dtmtiles, raster::raster)
names(tiles)[1:2] <- c('x', 'y')
dtm <- rast(do.call(raster::merge, tiles))

# Read in vector data defining problem areas to be edited
clip <- st_read(here::here("analysis/data/raw_data/clipping_poly.gpkg"))

# Make the vector feature a raster
rclip <- terra::rasterize(vect(clip), dtm)

# Make the raster NA at the polygons
dmask <- mask(dtm, rclip, inverse = TRUE)

edrast <- file.path(dtmpath, "dtm10.tif")
writeRaster(dmask, edrast,
            overwrite = TRUE)

# Initiate GRASS
loc <- initGRASS(gisBase = grasspath,
                 mapset = "PERMANENT", override = TRUE)

# Set projection
execGRASS("g.proj", flags = "c",
          proj4 = as.character(crs(clip, proj = TRUE)))

# Load raster into GRASS
execGRASS("r.in.gdal", flags = c("overwrite"),
          parameters = list(input = edrast, output = "tmprast"))

# Set computational region to the raster
execGRASS("g.region", parameters = list(raster = "tmprast") )

# Perform regularized spline interpolation with tension
# using default settings for the tension parameter.
execGRASS("r.fillnulls", flags = "overwrite",
          parameters = list(input = "tmprast",
                            output = "intrast",
                            method = "rst"))

# An error occurs when loading the raster to R unless this is specified
use_sp()

# Load raster back into R
intrast <- readRAST("intrast")

# Save the raster
writeRaster(raster::raster(intrast), edrast,
            overwrite = TRUE)
