library(terra)
library(sf)
library(rgrass7)

# Script to remove roads and other disturbances in the DTM, and interpolate
# new values. Polygon data with which to remove these features has been
# created manually based on an inspection of satellite data and the DTM.


# --- Only the edited DTM is distributed with this research compendium ---
#                         (see 00dtm_prep.R)
# To run, this script requires a local installation of GRASS GIS as well as
# the manual download of unedited DTM from hoydedata.no


##  This path has to be set to the local installation of GRASS
# grasspath <- "/usr/lib/grass78"

try(!is.na(grasspath)) # Uncomment and set path to local GRASS installation
                       #  above or this will throw and error.

# Path to locally stored unedited version of the DTM from the
# Norwegian Mapping Authority (freely available from hoydedata.no)
# dtm <- rast("/home/isak/spatial_data/dtm10.tif")

try(!is.na(dtm)) # Uncomment and read in locally stored unedited DTM
                 # above or this will throw and error.

# Read in vector data defining problem areas to be edited
clip <- st_read(here::here("analysis/data/raw_data/clipping_poly.gpkg"))

# Make the vector feature a raster
rclip <- terra::rasterize(vect(clip), dtm)

# Make the raster NA at the polygons
dmask <- mask(dtm, rclip, inverse = TRUE)

dtmpath <-  here::here("analysis/data/derived_data")
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
