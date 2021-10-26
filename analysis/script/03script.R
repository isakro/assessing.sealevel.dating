library(terra)
library(sf)
library(rgrass7)

# Script to remove roads and other disturbances in the DTM, and interpolate
# new values. Polygon data with which to remove these features has been
# created manually based on an inspection of satellite data and the DTM.

# This path has to be set to local installation of GRASS
grasspath <- "/usr/lib/grass78"
# Path to raster data. Stored locally for now
rpath <- "/home/isak/phd/eaa_presentation/dtm10/"

# Read in data
dtm <- rast(file.path(rpath, "dtm10_32.tif"))
clip <- st_read(here::here("analysis/data/raw_data/clipping_poly.gpkg"))

# Make the vector feature a raster
rclip <- rasterize(vect(clip), dtm)

# Make the raster NA at the polygons
dmask <- mask(dtm, rclip, inverse = TRUE)

edrast <- file.path(rpath, "dtm10.tif")
writeRaster(dmask, edrast,
            overwrite = TRUE)

# Initiate GRASS
loc <- initGRASS(gisBase = grasspath,
                 mapset = "PERMANENT", override = TRUE)

# Set projection
execGRASS("g.proj", flags = "c",
          proj4 = crs(dtm, proj = TRUE))

# Load raster into GRASS
execGRASS("r.in.gdal", flags = c("overwrite"),
          parameters = list(input = edrast, output = "tmprast"))

# Set computational region to the raster
execGRASS("g.region", parameters = list(raster = "tmprast") )

# Perform interpolation using regularized spline interpolation with tension
# with default settings for the tension parameter.
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
