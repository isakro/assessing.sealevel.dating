library(tidyverse)
library(sf)
library(tmap)
library(tmaptools)
library(terra)
library(oxcAAR)
library(topoDistance)
quickSetupOxcal()

# Check github page for installation of gdalio:
# https://github.com/hypertidy/gdalio

source(here::here("analysis/script/02script.R"))
load(here::here("analysis/data/derived_data/01data.RData"))
load(here::here("analysis/data/derived_data/02data.RData"))

dtm <- rast("/home/isak/phd/eaa_presentation/dtm1/dtm1_33_120_109.tif")

start_time <- Sys.time()

# Example site
sitename <- "Langangen Vestgård 1"

sitel <- filter(sites_sa, name == sitename)
siter <- filter(rcarb_sa, site_name == sitename)


# Function for creating a bounding box polygon around a feature,
# the size of which can be adjusted. Used for plotting purposes
# and ease of raster handling.
bboxpoly <- function(feature, xy_adjustment) {
  feature_bbox <- st_bbox(feature)
  feature_bbox[1:2] <- feature_bbox[1:2] - xy_adjustment
  feature_bbox[3:4] <-feature_bbox[3:4] + xy_adjustment
  feature_bboxpoly <- st_as_sf(st_as_sfc(feature_bbox))
  return(feature_bboxpoly)
}

# Create bounding box with 10 meter clearing around site limit
location_bbox <- bboxpoly(sitel, 10)

# Plot site with features
tm_shape(location_bbox, unit = "m") +
  tm_borders() +
  tm_shape(sitel) +
  tm_borders() +
  tm_shape(siter) +
  tm_fill(col = "black") +
  tm_scale_bar(breaks = c(0, 10, 20), text.size = 0.8) +
  tm_layout(legend.outside = TRUE,
            main.title = sitename,
            main.title.size = 1)

# Calibrating and modelling radiocarbon dates
dates <- R_Date(siter$lab_code, siter$c14_bp, siter$error)
oxcal_code <- paste0('
  Plot()
   {
    Sequence()
    {
     Boundary("Start");
     Phase()
     {
     ',
     dates,
     '
     Sum("Sum");
     };
     Boundary("End");
    };
  };')

oxcal_exec <- executeOxcalScript(oxcal_code)
oxcal_read <- readOxcalOutput(oxcal_exec)
oxcal_data <- parseOxcalOutput(oxcal_read, only.R_Date = FALSE)
oxcal_data_dates <- parseOxcalOutput(oxcal_read, only.R_Date = TRUE)

plot(oxcal_data_dates)

posterior <- data.frame(
  probabilities = oxcal_data$Sum$posterior_probabilities$probabilities,
  date_grid = oxcal_data$Sum$posterior_probabilities$dates,
  site = siter$site_name[1])


# Idenfity closest displacement curves and direction
closestiso <- closest_isobases(feature = sitel, isobases = isobases)
sitel <- st_as_sf(cbind(sitel, closestiso))

sitecurve <- interpolate_curve(years = xvals,
                               isobase1 = sitel$isobase1,
                               isobase2 = sitel$isobase2,
                               target = sitel,
                               dispdat = displacement_curves,
                               isodat = isobases,
                               direction_rel_curve1 = sitel$dir_to_1)

sitecurve$name <- sitename
dispcurves <- displacement_curves %>%
  filter(name %in% c(sitel$isobase1, sitel$isobase2)) %>%
  rbind(sitecurve)

ggplot(dispcurves, aes(x = years, col = name)) +
  geom_line(aes(y =  upperelev)) +
  geom_line(aes(y = lowerelev)) +
  ylab("Meters above present sea level") +
  xlab("Cal. BCE/CE") +
  scale_colour_manual(values = c("red", "blue", "seagreen2")) +
  theme_bw() +
  theme(legend.title = element_blank(), legend.position = "bottom",
        legend.direction = "horizontal")

location_bbox <- bboxpoly(sitel, 250)
locationarea <- terra::crop(dtm, location_bbox)

samps <- 1000
output <- list(length = samps)
seapolygons <- list(length = samps)
topopaths <- list(length = samps)

for(i in 1:samps){
  # Set up sampling frame from sum of posterior densities
  samplingframe <- list(x = posterior[["date_grid"]],
                        y = posterior[["probabilities"]])
  # Draw single sample date
  sdate <- with(samplingframe, sample(x, size = 1, prob = y))

  # Find sealevel-range at sampled date
  uppers <- approx(sitecurve[,"years"], sitecurve[,"upperelev"],
          xout = sdate)[["y"]]
  lowers <- approx(sitecurve[,"years"], sitecurve[,"lowerelev"],
         xout = sdate)[["y"]]

  # Set up empty
  results <- data.frame(matrix(ncol = 3, nrow = samps))
  names(results) <- c("vertdist", "hordist", "topodist")
  spolygons <- list(length = samps)
  tpaths <- list(length = samps)

  # Loop over and sample sea levels at 5 centimetre intervals
  for(j in 1:samps){
    sealevel <- sample(seq(lowers, uppers, 0.05), 1)

    # Create polygon representing the elevation of the sea level
    rclmat <- matrix(c(sealevel, Inf, NA, 0, sealevel, 1),
                     ncol = 3, byrow = TRUE)
    classarea <- terra::classify(locationarea, rclmat)
    seapoly <- st_as_sf(terra::as.polygons(classarea))
    tm_shape(seapoly) + tm_fill() + tm_shape(sitel) + tm_borders()

    spolygons[[j]] <- seapoly

    # Find the points nearest eachother on site- and sea-polygons
    distline <- st_nearest_points(sitel, seapoly)
    # Retrieve the two points
    distpts <- st_cast(distline, "POINT")

    # If these are at the same coordinate (i.e. the polygons are overlapping),
    # return all values as 0
    if((st_coordinates(distpts[1]) == st_coordinates(distpts[2]))[1]){
      results$vertdist[j] <- 0
      results$hordist[j] <- 0
      results$topodist[j] <- 0

      # If not find vertical, horisontal and topographical distance
    } else{
      # Find the lowest elevation on the site polygon and subtract the elevation
      # of the sealevel
      siteelev <- extract(locationarea, vect(sitel), fun = min)[2]
      results$vertdist[j] <- siteelev - sealevel

      # Find the distance of the line between the closest points
      results$hordist[j] <- st_length(distline)

      # Use the
      coords <- rbind(st_coordinates(distpts[1]),
                      st_coordinates(distpts[2]))
      rownames(coords) <- c("loc", "sea")

      # Find topographic distance (retrieve path)
      topodist <- topoDist(raster::raster(locationarea), coords,
                                      paths = TRUE)
      # Store distance
      results$topodist[j] <- topodist[[1]]["loc", "sea"]
      # Store path geometry
      tpaths[[j]] <- topodist[[2]]
    }
  }
  # Return results
  output[[i]] <- results
  seapolygons[[i]] <- do.call(rbind, spolygons)
  topopaths[[i]] <- do.call(rbind, tpaths)
}

output <- do.call(rbind, results)
seapol <- do.call(rbind, seapolygons)
topop <- do.call(rbind, topopaths)

# Save output here in case of crash in the following data handling and
# visualisation.
save(output, seapol, topop, locationarea,
     file = here::here("analysis/data/derived_data/eaa.RData"))

end_time <- Sys.time()
end_time - start_time

grid <- st_make_grid(locationarea,
                   cellsize = c(res(locationarea)[1], res(locationarea)[2]),
                   n = c(ncol(locationarea), nrow(locationarea)))

grid_intersects <- st_intersects(tst, seapol)
overlaps <- sapply(grid_intersects, length)
overlapgrid <- data.frame(grid, overlaps) %>%
  st_as_sf()

tm_shape(overlapgrid) + tm_fill(col = "overlaps") +
  tm_shape(sitel) + tm_borders(lwd = 1.5) +
  tm_legend(legend.outside = TRUE)

plot(tst, colour = overlaps)

tm_shape(overlapgrid, unit = "m") + tm_fill(col = "overlaps", title = "") +
  tm_shape(sitel) + tm_borders(col = "black", lwd = 1) +
  tm_legend(legend.outside = TRUE) + tm_scale_bar(text.size = 0.8)

# tm_shape(locationarea, unit = "m") +
#   tm_raster(palette = ,
#             title = "m a.s.l.",
#             legend.show = TRUE, style = "cont")


# Trying a solution for "glimpsing" rasters. Inspired by @joshualerickson:
# https://twitter.com/joshualerickson/status/1422662520519593986/photo/1

# library(vapour)
# library(gdalio)

# dtmpath <- "/home/isak/phd/eaa_presentation/dtm1/data/dtm1_33_120_109.tif"

# ri <- vapour::vapour_raster_info(dtmpath)
# aoi <- mapedit::drawFeatures()
# aoi <- aoi %>% st_transform(crs = ri$projection)
#
# aoi_bbox <- function(aoi){
#   aoi_bbox <- st_bbox(aoi)
#   extent <- c(aoi_bbox[[1]], aoi_bbox[[3]], aoi_bbox[[2]], aoi_bbox[[4]])
# }
#
# aoi_ext <- aoi_bbox(aoi)
#
# grid0 <- gdalio_set_default_grid(list(extent = aoi_ext,
#                                       dimension = c(800, 800),
#                                       projection = ri$projection))
# gdalio_set_default_grid(grid0)
#
# gdalio_raster <- function(dsn, ...){
#   v - gdalio_data(dsn, ...)
#   g <- gdalio_get_default_grid()
#   r <- raster::raster(raster::extent(g$extent), nrows = g$dimension[2],
#                       ncols = g$dimension[1], crs = g$projection)
#   if(length(v) > 1){
#     r <- raster::brick(replicate(length(v), r, simplify = FALSE))
#   }
#   raster::setValues(r, do.call(cbind, v))
# }
#
# cpg <- gdalio_raster(f)
# raster::plot(cpg, col = hcl.colours(26))
