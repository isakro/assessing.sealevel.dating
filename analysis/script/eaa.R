library(tidyverse)
library(sf)
library(tmap)
library(tmaptools)
library(cowplot)
library(terra)
library(oxcAAR)
library(topoDistance)
quickSetupOxcal()

# Check github page for installation of gdalio:
# https://github.com/hypertidy/gdalio

source(here::here("analysis/script/02script.R"))
load(here::here("analysis/data/derived_data/01data.RData"))
load(here::here("analysis/data/derived_data/02data.RData"))

dtm1 <- rast("/home/isak/phd/eaa_presentation/dtm1/dtm1_33_120_109.tif") # LV1
dtm2 <- rast("/home/isak/phd/eaa_presentation/dtm1/dtm1_33_119_108.tif")
dtm3 <- rast("/home/isak/phd/eaa_presentation/dtm1/dtm_33_115_116_105.tif")
dtm1[dtm1 <= 0] <- NA

# Trondal, Sandvigen 1, Langangen Vestgård 1, Hesthag C2,
# Stokke/Polland 1 (partly highway), Stokke/Polland 5 (highway)
# Vallermyrene 1a

countries <- st_read(here::here('analysis/data/raw_data/naturalearth_countries.gpkg'))

# assign isobases to sites (using st_centroid so that sites crossing isobase
# only get a single occurence)
sites_sa <- st_join(st_make_valid(sites_sa), isopolys, join = st_intersects, largest = TRUE)
site_pts <- st_centroid(sites_sa)

# Create bounding box around the sites, both for extent indicator on the
# overview map and to create basemap for the site overview.
sitbbox2 <- st_bbox(site_pts)
sitbbox2[1:2] <- sitbbox2[1:2] - 15000
sitbbox2[3:4] <- sitbbox2[3:4] + 15000
boundingpoly2 <- st_as_sf(st_as_sfc(sitbbox2))

# same as above, using the smaller bounding box
bound_reproj2 <- st_transform(boundingpoly2, st_crs(countries))
mapcountries2 <- countries %>%
  filter(st_intersects(., bound_reproj2, sparse = FALSE))
count_reproj2 <- st_transform(mapcountries2, st_crs(site_pts))

# Assemble site map. anc defined first is used for placing the scalebar,
# which required some tweaking to get right.
anc <- as.numeric(c(sitbbox2$ymin, sitbbox2$xmax))
sitm1 <- ggplot() +
  geom_sf(data = count_reproj2, fill = "grey", colour = NA) +
  geom_sf(data = site_pts, fill = "red",
          size = 2.25, shape = 21,
          colour = "black") +
  labs(title = substitute(paste("Excavated Stone Age sites (n = ", rws, ")"),
                          list(rws = nrow(site_pts)))) +
  ggsn::scalebar(data = site_pts, dist = 20, dist_unit = "km",
                 transform = FALSE, st.size = 4, height = 0.02,
                 border.size = 0.1, st.dist = 0.03,
                 anchor = c(x = anc[2] - 15500, y = anc[1]) + 8000) +
  coord_sf(xlim = c(sitbbox2[1], sitbbox2[3]),
           ylim = c(sitbbox2[2], sitbbox2[4]),
           expand = FALSE) +
  theme_bw() + theme(axis.title=element_blank(),
                     axis.text.y = element_blank(),
                     axis.text.x = element_blank(),
                     rect = element_rect(),
                     axis.ticks = element_blank(),
                     panel.grid.major = element_blank()) +
  guides(fill = guide_legend(ncol = 2,
                             title.position = "top", title.hjust = 0.5))

sites_rcarb <- site_pts %>% filter(name %in% unique(rcarb_sa$site_name))

sitm2 <- ggplot() +
  geom_sf(data = count_reproj2, fill = "grey", colour = NA) +
  geom_sf(data = sites_rcarb, fill = "red",
          size = 2.25, shape = 21,
          colour = "black") +
  labs(title = substitute(paste("Sites with "^"14",
                                "C-date (n = ", rws, ")"),
                          list(rws = nrow(sites_rcarb)))) +
  ggsn::scalebar(data = site_pts, dist = 20, dist_unit = "km",
                 transform = FALSE, st.size = 4, height = 0.02,
                 border.size = 0.1, st.dist = 0.03,
                 anchor = c(x = anc[2] - 15500, y = anc[1]) + 8000) +
  coord_sf(xlim = c(sitbbox2[1], sitbbox2[3]),
           ylim = c(sitbbox2[2], sitbbox2[4]),
           expand = FALSE) +
  theme_bw() + theme(axis.title=element_blank(),
                     axis.text.y = element_blank(),
                     axis.text.x = element_blank(),
                     rect = element_rect(),
                     axis.ticks = element_blank(),
                     panel.grid.major = element_blank()) +
  guides(fill = guide_legend(ncol = 2,
                             title.position = "top", title.hjust = 0.5))

sites_sa_pts <- site_pts %>% filter(name %in% unique(rcarb_sa$site_name))

sitm3 <- ggplot() +
  geom_sf(data = count_reproj2, fill = "grey", colour = NA) +
  geom_sf(data = sites_sa_pts, fill = "red",
          size = 2.25, shape = 21,
          colour = "black") +
  labs(title = substitute(paste("Sites with "^"14",
                                "C-date to the Stone Age (n = ", rws, ")"),
                          list(rws = nrow(sites_sa_pts)))) +
  ggsn::scalebar(data = site_pts, dist = 20, dist_unit = "km",
                 transform = FALSE, st.size = 4, height = 0.02,
                 border.size = 0.1, st.dist = 0.03,
                 anchor = c(x = anc[2] - 15500, y = anc[1]) + 8000) +
  coord_sf(xlim = c(sitbbox2[1], sitbbox2[3]),
           ylim = c(sitbbox2[2], sitbbox2[4]),
           expand = FALSE) +
  theme_bw() + theme(axis.title=element_blank(),
                     axis.text.y = element_blank(),
                     axis.text.x = element_blank(),
                     rect = element_rect(),
                     axis.ticks = element_blank(),
                     panel.grid.major = element_blank()) +
  guides(fill = guide_legend(ncol = 2,
                             title.position = "top", title.hjust = 0.5))

plot_grid(sitm1, sitm3, nrow = 1)


isomap <- ggplot() +
  geom_sf(data = count_reproj2, fill = "grey", colour = NA) +
  geom_sf(data = sites_sa_pts, fill = "red",
          size = 2.25, shape = 21,
          colour = "black") +
  geom_sf(data = isobases, aes(colour = name)) +
  scale_colour_manual(values = c("Arendal" = "black",
                        "Larvik" = "darkgreen",
                        "Tvedestrand" = "blue",
                        "Halden" = 'red')) +
  ggsn::scalebar(data = site_pts, dist = 20, dist_unit = "km",
                 transform = FALSE, st.size = 4, height = 0.02,
                 border.size = 0.1, st.dist = 0.03,
                 anchor = c(x = anc[2] - 15500, y = anc[1]) + 8000) +
  coord_sf(xlim = c(sitbbox2[1], sitbbox2[3]),
           ylim = c(sitbbox2[2], sitbbox2[4]),
           expand = FALSE) +
  theme_bw() + theme(axis.title=element_blank(),
                     axis.text.y = element_blank(),
                     axis.text.x = element_blank(),
                     rect = element_rect(),
                     axis.ticks = element_blank(),
                     panel.grid.major = element_blank(),
                     legend.position = "none")


plot_grid(isomap, isocurves_plot)

start_time <- Sys.time()

# Example site
sitename <- "Vallermyrene 1a"

sitel <- filter(sites_sa, name == sitename)
siter <- filter(rcarb_sa, site_name == sitename)


# Kvastad A2
siter <- siter[1:5,]

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
sitmap <- tmap_grob(tm_shape(location_bbox, unit = "m") +
                      tm_borders() +
                      tm_shape(sitel) +
                      tm_borders() +
                      tm_shape(siter) +
                      tm_fill(col = "black") +
                      tm_scale_bar(breaks = c(0, 10, 20), text.size = 0.8) +
                      tm_layout(legend.outside = TRUE))

# siter[class(st_geometry(siter) == "sfc_POINT"]

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

posteriorprobs <- data.frame(
  probabilities = oxcal_data$Sum$posterior_probabilities$probabilities,
  date_grid = oxcal_data$Sum$posterior_probabilities$dates,
  site = siter$site_name[1])

prior <- get_raw_probabilities(oxcal_data)
names(prior) <- paste0(get_name(oxcal_data))
prior <- prior[!(names(prior) %in% c("Sum", "End", "Start",
                                     "character(0)"))]
prior <-  bind_rows(prior, .id = "name")
prior$class <- "aprior"

posterior <- get_posterior_probabilities(oxcal_data)
names(posterior) <- get_name(oxcal_data)
posterior <- bind_rows(posterior[names(posterior) != "character(0)"],
                       .id = "name")
posterior$class <- "bposterior"
posterior <- mutate(posterior, class = ifelse(name == "Sum", "Posterior sum", class),
                    name = ifelse(name == "Sum", "Posterior sum", name))

datedat <- rbind(prior, posterior)

sitdat <- datedat %>% filter(!(name %in% c("Start", "End"))) %>%
  arrange(class) %>%
  ggplot() +
  ggridges::geom_ridgeline(aes(x = dates, y = name, height = probabilities * 100, fill = class, alpha = class)) +
  scale_fill_manual(values = c("Posterior sum" = "grey35", "aprior" = "white", "bposterior" = "grey")) +
  scale_alpha_manual(values = c("Posterior sum" = 1, "aprior" = 0.1, "bposterior" = 0.7)) +
  labs(title = sitename, y = "", x = "Cal. BCE") + theme_bw() + theme(legend.position = "none")

plot_grid(sitdat, sitmap)

# Interpolate displacement curve
sitecurve <- interpolate_curve(years = xvals,
                               isobase1 = sitel$isobase1,
                               isobase2 = sitel$isobase2,
                               target = sitel,
                               dispdat = displacement_curves,
                               isodat = isobases,
                               direction_rel_curve1 = sitel$dir_rel_1)

sitecurve$name <- sitename
dispcurves <- displacement_curves %>%
  filter(name %in% c(sitel$isobase1, sitel$isobase2)) %>%
  rbind(sitecurve)

ggplot(dispcurves, aes(x = years, col = name)) +
  geom_line(aes(y =  upperelev)) +
  geom_line(aes(y = lowerelev)) +
  ggridges::geom_ridgeline(data = filter(datedat, name == "Posterior sum"), aes(x = dates, y = min(siteelev[1:nrow(siteelev),]), height = probabilities * 150))
  ylab("Meters above present sea level") +
  xlab("Cal. BCE/CE") +
  scale_colour_manual(values = c("red", "blue", "seagreen2")) +
  theme_bw() +
  theme(legend.title = element_blank(), legend.position = "bottom",
        legend.direction = "horizontal")


samplingframe <- list(x = posteriorprobs[["date_grid"]],
                      y = posteriorprobs[["probabilities"]])

sdate <- round(with(samplingframe, sample(x, size = 1, prob = y)))
library(ggridges)
postsum <- filter(datedat, name == "Posterior sum")
postplot <- ggplot(postsum) +
  ggridges::geom_ridgeline(aes(x = dates, y = name, height = probabilities * 150), fill = "grey") +
  #geom_segment(aes(x = sdate, xend = sdate, y = as.numeric(name), yend = as.numeric(name) + 0.9), colour = "red") +
  labs(y = "", x = "Cal. BCE")  +
  scale_y_discrete(expand = c(0.01, 0)) + theme_bw() + theme(legend.position = "none",
                                                             panel.grid.major = element_blank(),
                                                             panel.grid.minor = element_blank(),
                                                             panel.border = element_blank(),
                                                             axis.line = element_line("black"),
                                                             axis.title.y = element_blank(),
                                                             axis.text.y = element_blank(),
                                                             axis.ticks.y = element_blank())

curveplot <- ggplot(sitecurve, aes(x = years)) +
  geom_line(aes(y =  upperelev)) +
  geom_line(aes(y = lowerelev)) +
  ylab("Meters above present sea level") +
  xlab("Cal. BCE/CE") +
  theme_bw() +
  theme(legend.title = element_blank(), legend.position = "bottom",
        legend.direction = "horizontal")


plot_grid(postplot, curveplot)

samp1plot <- postplot +
  geom_vline(xintercept = sdate, colour = "red") +
  scale_x_continuous(breaks = sort(c(-7000, -6800, sdate, -6600, -6400)))

plot_grid(samp1plot, curveplot)

uppers <- round(approx(sitecurve[,"years"], sitecurve[,"upperelev"],
                       xout = sdate)[["y"]], 1)
lowers <- round(approx(sitecurve[,"years"], sitecurve[,"lowerelev"],
                       xout = sdate)[["y"]], 1)

curvesamp1 <- curveplot + geom_segment(aes(x = sdate, xend = sdate, y = uppers,
                                           yend = -Inf), col = "red", linetype = "dashed") +
  geom_segment(aes(x = -Inf, xend = sdate, y = uppers,
                   yend = uppers), col = "red") +
  geom_segment(aes(x = -Inf, xend = sdate, y = lowers,
                   yend = lowers), col = "red") +
  scale_x_continuous(breaks = sort(c(-8000, sdate, -4000, 0)))


plot_grid(samp1plot, curvesamp1)

masls <- data.frame(x = 0:150, y = 0:150)
ggplot(masls, aes(x, y)) +
  geom_vline(xintercept = uppers, colour = "red") +
  geom_text(aes(x = uppers + 5, label = as.character(uppers), y = 10)) +
  geom_text(aes(x = lowers - 5, label = as.character(lowers), y = 10)) +
  geom_vline(xintercept = lowers, colour = "red") +
  scale_x_continuous(limits = c(0, 150),
                     name = "Meters above present sea level") +
  theme_bw()  +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line.x = element_line("black"),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())


location_bbox <- bboxpoly(sitel, 250)
sitearea <- terra::crop(dtm1, location_bbox)

sample_shoreline <- function(samps, sitel, sitecurve, sitearea, posteriorprobs){

  results <- data.frame(matrix(ncol = 4, nrow = samps))
  names(results) <- c("vertdist", "hordist", "topodist", "year")
  seapolygons <- list(length = samps)
  topopaths <- list(length = samps)

  for(i in 1:samps){
    # Set up sampling frame from sum of posterior densities
    samplingframe <- list(x = posteriorprobs[["date_grid"]],
                        y = posteriorprobs[["probabilities"]])
    # Draw single sample date
    sdate <- with(samplingframe, sample(x, size = 1, prob = y))
    results$year[i] <- sdate

    # Find sealevel-range at sampled date
    uppers <- approx(sitecurve[,"years"], sitecurve[,"upperelev"],
            xout = sdate)[["y"]]
    lowers <- approx(sitecurve[,"years"], sitecurve[,"lowerelev"],
           xout = sdate)[["y"]]

    sealevel <- sample(seq(lowers, uppers, 0.05), 1)

    # Create polygon representing the elevation of the sea level
    rclmat <- matrix(c(sealevel, Inf, NA, 0, sealevel, 1),
                       ncol = 3, byrow = TRUE)
    classarea <- terra::classify(sitearea, rclmat)
    seapoly <- st_combine(st_as_sf(terra::as.polygons(classarea)))

    seapolygons[[i]] <- seapoly

    # Retrieve lowest elevation on polygon
    siteelev <- extract(sitearea, vect(sitel), fun = min)[2]
    # Find difference to sea level elevation
    results$vertdist[i] <-  min(siteelev[1:nrow(siteelev),]) - sealevel

    # If these are at the same coordinate (i.e. the polygons are overlapping),
    # return all values as 0
    if(st_overlaps(seapoly, sitel, sparse = FALSE)){
      results$hordist[i] <- 0
      results$topodist[i] <- 0

      # Else...if the entire site polygon is in the sea
    } else if(st_contains(seapoly, sitel, sparse = FALSE)){

      # Find the nearest points on site- and seapolygons, making the sea a
      # multiline feature as it contains the site
      distline <- st_nearest_points(sitel, st_cast(seapoly,"MULTILINESTRING"))

      # Retrieve the two points
      distpts <- st_cast(distline, "POINT")

      # Make the value negative to indicate that the site is located below sea level
      results$hordist[i] <- st_length(distline) * -1

      # Find topographic distance
      coords <- rbind(st_coordinates(distpts[1]),
                      st_coordinates(distpts[2]))
      rownames(coords) <- c("loc", "sea")

      # Find topographic distance (retrieve path)
      topodist <- topoDist(raster::raster(sitearea), coords,
                           paths = TRUE)
      # Store distance (negative to indicate below sea-level)
      results$topodist[i] <- topodist[[1]]["loc", "sea"] * -1
      # Store path geometry
      topopaths[[i]] <- st_as_sf(topodist[[2]])

      # Else, that is not overlapping slightly or completely, i.e. the site is
      # located some distance from the sea, perform above steps bu without
      # negative values
    } else {
      # Find the nearest points on the site- and seapolygons,
      distline <- st_nearest_points(sitel, seapoly)

      # Retrieve the two points
      distpts <- st_cast(distline, "POINT")

      # Retrieve horisontal distance
      results$hordist[i] <- st_length(distline)

      # Find topographic distance
      coords <- rbind(st_coordinates(distpts[1]),
                      st_coordinates(distpts[2]))
      rownames(coords) <- c("loc", "sea")

      # Find topographic distance (retrieve path)
      topodist <- topoDist(raster::raster(sitearea), coords,
                           paths = TRUE)
      # Store distance
      results$topodist[i] <- topodist[[1]]["loc", "sea"]
      # Store path geometry
      topopaths[[i]] <- st_as_sf(topodist[[2]])
    }
  }
  # If there are no topopaths at all, do not try to return these
  if(any(!(summary(topopaths)[,2]) %in% c("-none-", "sf"))){
    return(list(
      results = results,
      seapol = st_as_sf(st_sfc(do.call(rbind, seapolygons)),
                        crs = st_crs(sitel)),
      topop = st_as_sf(st_sfc(do.call(rbind, topopaths)),
                       crs = st_crs(sitel))
    ))
  } else
    return(list(
      results = results,
      seapol = st_as_sf(st_sfc(do.call(rbind, seapolygons)),
                        crs = st_crs(sitel))))
}

# Function returns a list holding numerical results, polygons representing sea,
# and
output <- sample_shoreline(samps = 1000, sitel = sitel, sitecurve = sitecurve,
                           sitearea = sitearea, posteriorprobs = posteriorprobs)

# Save output here in case of crash in the following data handling and
# visualisation.
save(output,
     file = here::here("analysis/data/derived_data/lv1.RData"))

seapol <- output$seapol
topopath <- output$topop

end_time <- Sys.time()
end_time - start_time

grid <- st_make_grid(sitearea,
                   cellsize = c(res(sitearea)[1], res(sitearea)[2]),
                   n = c(ncol(sitearea), nrow(sitearea)))

grid_intersects <- st_intersects(grid, seapol)
overlaps <- sapply(grid_intersects, length)
overlapgrid <- data.frame(grid, overlaps) %>%
  st_as_sf()

shoremaplv1 <- tmap_grob(tm_shape(overlapgrid, unit = "m") + tm_fill(col = "overlaps", title = "", palette = "Greys", alpha = 0.5, legend.show = FALSE) +
  tm_shape(sitel) + tm_borders(col = "black", lwd = 1) +
  #tm_shape(topopath) + tm_lines() +
  tm_legend(legend.outside = TRUE) + tm_scale_bar(text.size = 0.8))

distplotlv1 <- ggplot() + geom_boxplot(aes(x = "Vertical distance", y = output$results$vertdist)) +
  geom_boxplot(aes(x = "Horisontal distance", y = output$results$hordist)) +
  geom_boxplot(aes(x = "Topographic distance", y = output$results$topodist)) +
  labs(y = "Distance from shoreline (m)", x = "") +
  theme_bw()

lv1plot <- plot_grid(lv1dat, shoremaplv1, distplotlv1, nrow = 1)

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
