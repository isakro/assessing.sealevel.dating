library(tidyverse)
library(here)
library(sf)
library(tmap)
library(cowplot)
library(era)

# Load radiocarbon and site data from the first script
load(here("analysis/data/derived_data/01data.RData"))

# Read in displacement curves
arendal_raw <- read.csv(
  here("analysis/data/raw_data/displacement_curves/bjornebu.csv"))
tvedestrand_raw <- read.csv(
  here("analysis/data/raw_data/displacement_curves/hanto.csv"))
larvik_raw <- read.csv(
  here("analysis/data/raw_data/displacement_curves/gunnarsrod.csv"))
skoppum_raw <- read.csv(
  here("analysis/data/raw_data/displacement_curves/skoppum.csv"))

# Read in coast polgyon
coast <- read_sf(here("analysis/data/raw_data/coast.gpkg"))

# Read in isobase centre points
centpts <- read_sf(here("analysis/data/raw_data/isobase_centrepts.gpkg"))

#### Setting up displacement curves ####

# Define function to fit displacement curve to a set interval of years.
# (This ensures that curves consist of vectors of the same length, necessary
# for operations below to work).
curve_interval <- function(years, curve){
  upperelev <- approx(curve[,"upperyr"],
                      curve[,"upperelev"],
                      xout = years)[["y"]]
  lowerelev <- approx(curve[,"loweryr"],
                      curve[,"lowerelev"],
                      xout = years)[["y"]]
  values <- cbind.data.frame(years, upperelev, lowerelev)
  return(values)
}

# Define sequence of years for which to define displacement curves (cal BP)
xvals <- seq(0, 12500, 1)

# Define curves on a per year interval, using the sequence of years defined
# above
larvik <- curve_interval(xvals, larvik_raw)
larvik$name <- "Larvik"
tvedestrand <- curve_interval(xvals, tvedestrand_raw)
tvedestrand$name <- "Tvedestrand"
arendal <- curve_interval(xvals, arendal_raw)
arendal$name <- "Arendal"
horten <-  curve_interval(xvals, skoppum_raw)
horten$name <- "Horten"

# # Transform to BCE/CE using Joe Roe's era package and make BCE negative
# # to corresponds to oxcAAR later on. [This gave some issue with rounding
# errors]
# displacement_curves <- rbind(larvik, tvedestrand, arendal, horten)
# displacement_curves$years <- yr(displacement_curves$years, "cal BP")
# displacement_curves$years <- yr_transform(displacement_curves$years, "BCE") * -1
# xvals <- as.integer(yr_transform(yr(xvals, "cal BP"), "BCE") * -1)

displacement_curves <- rbind(larvik, tvedestrand, arendal, horten)
displacement_curves$years <- (displacement_curves$years - 1950) * -1
xvals <- (xvals - 1950) * -1

# Specify arbitrarily long distance of the line to represent the isobase at
# each centre-point.
isodist <- 9000000

# Specify direction of gradient
deg <- 327

# Loop over points and create isobase
for (i in 1:nrow(centpts)){
  # Find x and y coords
  x <- st_coordinates(centpts[i, ])[1]
  y <- st_coordinates(centpts[i, ])[2]

  # Find coords at the specified distance from the point at deg degree angle,
  # adding isodist
  xx <- x + isodist * (cos(deg * pi / 180))
  yy <- y + isodist * (sin(deg * pi / 180))

  # Find coords at the specified distance from the point at deg degree angle,
  # subtracting isodist
  xx2 <- x - isodist * (cos(deg * pi / 180))
  yy2 <- y - isodist * (sin(deg * pi / 180))

  # Create points from the identified coordinates
  pts <- st_sfc(st_multipoint(rbind(st_coordinates(centpts[i,])[1,],
                                    c(xx, yy), c(xx2, yy2)))) %>%
    st_set_crs(32632) %>% # WGS 84 / UTM 32N
    st_sf()

  # Combine these into a line
  assign(paste0("isobase", i), st_cast(pts, to = 'LINESTRING'))
}

# Combine the results and add the name of the centre-points
isobases <- rbind(isobase1, isobase2, isobase3, isobase4)
isobases$name <- centpts$name

# Find bounding box for centre-points and increase this for plotting
centpts_bb <- st_bbox(centpts)
centpts_bb[1:2] <- centpts_bb[1:2] - 20000
centpts_bb[3:4] <- centpts_bb[3:4] + 20000

# Split study area into a series of polygons to assign correct displacement
# curve to sites to be analysed
# First find study area outline and make this a polygon
outline <- st_as_sf(st_as_sfc(centpts_bb))
# Then split this by the isobases
isopolys <- st_collection_extract(lwgeom::st_split(outline, isobases))
# Then assign isobasenames to the polygons and the direction (in case both are
# to the north-east).
isopolys$isobase1 <- c("Arendal", "Tvedestrand", "Larvik", "Horten", "Horten")
isopolys$isobase2 <- c("Tvedestrand", "Arendal", "Tvedestrand", "Larvik",
                       "Larvik")
isopolys$dir_rel_1 <- c("sw", "sw", "sw", "sw", "nw")

# Save required data
save(displacement_curves, isobases, isopolys, xvals,
     file = here("analysis/data/derived_data/02data.RData"))
