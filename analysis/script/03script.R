library(terra)
library(sf)

# Script to remove roads and other disturbances in the DTM, and interpolate
# new values. Polygon data with which to remove these features has been
# created manually based on an inspection of satellite data and the DTM.

