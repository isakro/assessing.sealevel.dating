library(tidyverse)
library(sf)
library(ggridges)
library(tmap)
library(tmaptools)
library(cowplot)
library(terra)
library(oxcAAR)
library(topoDistance)
quickSetupOxcal()

# This script is an absolute mess and will not run. The main components are
# available in the other scripts and defintion of relevant functions mainly in
# 03script.R

source(here::here("analysis/script/02script.R"))
load(here::here("analysis/data/derived_data/01data.RData"))
load(here::here("analysis/data/derived_data/02data.RData"))
load(here::here("analysis/data/derived_data/sp1.RData"))

dtm1 <- rast("/home/isak/phd/eaa_presentation/dtm1/dtm1_33_120_109.tif") # LV1
dtm2 <- rast("/home/isak/phd/eaa_presentation/dtm1/dtm1_33_119_108.tif")
dtm3 <- rast("/home/isak/phd/eaa_presentation/dtm1/dtm_33_115_116_105.tif")
dtm4 <- rast("/home/isak/phd/eaa_presentation/dtm1/vm1/dtm1_33_119_109.tif")
dtm5 <- rast("/home/isak/phd/eaa_presentation/dtm1/trondal/trondal.tif")
dtm1[dtm1 <= 0] <- NA

# Tverdal, Sandvigen 1, Langangen Vestgård 1, Hesthag C2,
# Stokke/Polland 1 (partly highway), Stokke/Polland 5 (highway)
# Vallermyrene 1a

countries <- st_read(here::here('analysis/data/raw_data/naturalearth_countries.gpkg'))

# assign isobases to sites
sites_sa <- st_join(st_make_valid(sites_sa), isopolys, join = st_intersects, largest = TRUE)
site_pts <- st_centroid(sites_sa)


model_rcarb <- function(data){

  # Format individual dates for R_Date in Oxcal
  dates <- R_Date(data[["lab_code"]], data[["c14_bp"]], data[["error"]])

  # Manually compile rest of OxCal code
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

  # Run OxCal code - wrapped this in some functions to  not print to
  # rmarkdown
  invisible(capture.output(oxcal_exec <- executeOxcalScript(oxcal_code)))
  oxcal_read <- readOxcalOutput(oxcal_exec)
  # Note that parseOxcalOutput() defaults to only returning individual dates,
  # specified with only.R_Date
  oxcal_data <- parseOxcalOutput(oxcal_read, only.R_Date = FALSE)

  # Assemble data frame to be returned
  posterior <- data.frame(
    probabilities = oxcal_data$Sum$posterior_probabilities$probabilities,
    date_grid = oxcal_data$Sum$posterior_probabilities$dates,
    site = data[["site_name"]][1])

  # Return data
  return(posterior)
}

# Apply modelling function to the radiocarbon dates, grouped by site name
posteriors <- rcarb_sa %>%
  group_by(site_name) %>%
  group_map(~ model_rcarb(.), .keep = TRUE)


posteriors_t <- rcarb_sa %>%
  filter(rcarb_cor == "t") %>%
  group_by(site_name) %>%
  group_map(~ model_rcarb(.), .keep = TRUE)

rcarb_data <- do.call(rbind, posteriors)
rcarb_data_t <- do.call(rbind, posteriors_t)

ggplot(rcarb_data_t) +
  geom_ridgeline(aes(x = date_grid, y = site, height = probabilities * 100,
                     fill = site,
                     colour = site)) +
  xlab("Cal. BCE") +
  theme_classic() +
  theme(axis.title.y = element_blank()) +
  guides(fill =  guide_legend(title.position = "top", title.hjust = 0.5),
         colour = FALSE)

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
sitename <- "Kvastad A2"

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

sitmap1 <- ggplot() +
  geom_sf(data = count_reproj2, fill = "grey", colour = NA) +
  geom_sf(data = isobases, aes(colour = name)) +
  geom_sf(data = st_centroid(sitel), size = 3, shape = 21, fill = "red",
          colour = "black") +
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

# Plot site with features
sitmap <- tmap_grob(tm_shape(sitearea, unit = "m") +
                      tm_raster(title = "Elevation (m)", palette = terrain.colors(252)) +
                      tm_shape(sitel) +
                      tm_fill(col = "red") +
                      tm_shape(siter) +
                      tm_dots(col = "black") +
                      tm_scale_bar(breaks = c(0, 100, 200), text.size = 0.8))

# siter[class(st_geometry(siter) == "sfc_POINT"]

tst <- model_dates(siter)

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

plot_dates(tst, sitename)

sitdat <- datedat %>% filter(!(name %in% c("Start", "End"))) %>%
  arrange(class) %>%
  ggplot() +
  ggridges::geom_ridgeline(aes(x = dates, y = name, height = probabilities * 100, fill = class, alpha = class)) +
  scale_fill_manual(values = c("Posterior sum" = "grey35", "aprior" = "white", "bposterior" = "grey")) +
  scale_alpha_manual(values = c("Posterior sum" = 1, "aprior" = 0.1, "bposterior" = 0.7)) +
  labs(title = sitename, y = "", x = "cal BCE") + theme_bw() + theme(legend.position = "none")

plot_grid(sitdat, sitmap, sitmap1, nrow = 1)

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
  ggridges::geom_ridgeline(data = filter(datedat, name == "Posterior sum"), aes(x = dates, y = min(siteelev[1:nrow(siteelev),]), height = probabilities * 100))
  ylab("Meters above present sea level") +
  xlab("Cal. BCE/CE") +
  scale_colour_manual(values = c("red", "blue", "seagreen2")) +
  theme_bw() +
  theme(legend.title = element_blank(), legend.position = "bottom",
        legend.direction = "horizontal")


samplingframe <- list(x = posteriorprobs[["date_grid"]],
                      y = posteriorprobs[["probabilities"]])

#sdate <- round(with(samplingframe, sample(x, size = 1, prob = y)))
sdate <- -3016
postsum <- filter(datedat, name == "Posterior sum")
postplot <- ggplot(postsum) +
  ggridges::geom_ridgeline(aes(x = dates, y = name, height = probabilities * 100), fill = "grey") +
  #geom_segment(aes(x = sdate, xend = sdate, y = as.numeric(name), yend = as.numeric(name) + 0.9), colour = "red") +
  labs(y = "", x = "cal BCE")  +
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
  xlab("cal BCE/CE") +
  theme_bw() +
  theme(legend.title = element_blank(), legend.position = "bottom",
        legend.direction = "horizontal")


plot_grid(postplot, curveplot)

samp1plot <- postplot +
  geom_vline(xintercept = sdate, colour = "red") +
  scale_x_continuous(breaks = sort(c(-3750, -3250, sdate, -2750)))

plot_grid(samp1plot, curveplot)

uppers <- round(approx(sitecurve[,"years"], sitecurve[,"upperelev"],
                       xout = sdate)[["y"]], 1)
lowers <- round(approx(sitecurve[,"years"], sitecurve[,"lowerelev"],
                       xout = sdate)[["y"]], 1)

curvesamp1 <- curveplot + geom_segment(aes(x = sdate, xend = sdate, y = uppers,
                                           yend = -Inf), col = "red", linetype = "dashed") +
  geom_segment(aes(x = -Inf, xend = sdate, y = uppers,
                   yend = uppers), col = "red", size = 0.25) +
  geom_segment(aes(x = -Inf, xend = sdate, y = lowers,
                   yend = lowers), col = "red", size = 0.25) +
  annotate("text", x = -8000, y = uppers + 1, label = as.character(uppers)) +
  annotate("text", x = -8000, y = lowers - 1, label = as.character(lowers)) +
  scale_x_continuous(breaks = sort(c(-8000, sdate, -4000, 0)))


plot_grid(samp1plot, curvesamp1)

sealevel <- sample(seq(lowers, uppers, 0.05), 1)

masls <- data.frame(x = 0:50, y = 0:50)
elevsamp1 <- ggplot(masls, aes(x, y)) +
  annotate("text", x = uppers + 0.5, y = 10, label = as.character(uppers)) +
  annotate("text", x = lowers - 0.5, y = 10, label = as.character(lowers)) +
  geom_vline(xintercept = uppers, colour = "red") +
  geom_vline(xintercept = lowers, colour = "red") +
  geom_vline(xintercept = sealevel, colour = "red", linetype = "dashed") +
  scale_x_continuous(limits = c(5, 20), breaks = c(0, 5, 10, sealevel, 20),
                     name = "Meters above present sea level") +
  theme_bw()  +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line.x = element_line("black"),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())

rclmat <- matrix(c(sealevel, Inf, NA, 0, sealevel, 1),
                 ncol = 3, byrow = TRUE)
classarea <- terra::classify(sitearea, rclmat)
seapoly <- st_combine(st_as_sf(terra::as.polygons(classarea)))


locbbox <- st_bbox(sitel)
locbbox[1:2] <- locbbox[1:2] - 100
locbbox[3:4] <- locbbox[3:4] + 100

anc <- as.numeric(c(locbbox$ymin, locbbox$xmax))


seaplot1 <- ggplot() +
  geom_sf(data = seapoly, fill = "grey", colour = NA) +
  geom_sf(data = sitel, fill = NA, colour = "red") +
  labs(title =  paste0("Sealevel at ", sealevel, " m above present")) +
  coord_sf(xlim = c(locbbox[1], locbbox[3]),
           ylim = c(locbbox[2], locbbox[4]), expand = FALSE) +
  ggsn::scalebar(data = sitel, dist = 20, dist_unit = "m",
                 transform = FALSE, st.size = 4, height = 0.1,
                 border.size = 0.1, st.dist = 0.2,
                 anchor = c(x = anc[2] - 10, y = anc[1] + 10)) +
  theme_bw() + theme(axis.title=element_blank(),
                     axis.text.y = element_blank(),
                     axis.text.x = element_blank(),
                     rect = element_rect(),
                     axis.ticks = element_blank(),
                     panel.grid.major = element_blank(),
                     legend.position = "none")

elevsamp1

seaplot1

# Vertical distance
siteelev <- extract(sitearea, vect(sitel), fun = min)[2]
vertdist <-  min(siteelev[1:nrow(siteelev),]) - sealevel

# Horisontal distance
distline <- st_nearest_points(sitel, st_cast(seapoly,"MULTILINESTRING"))
distpts <- st_cast(distline, "POINT")
hordist <- st_length(distline)

# Topographic distance
coords <- rbind(st_coordinates(distpts[1]),
                st_coordinates(distpts[2]))
rownames(coords) <- c("loc", "sea")
topodist <- topoDist(raster::raster(sitearea), coords,
                     paths = TRUE)
topdist <- topodist[[1]]["loc", "sea"]


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
      results$vertdist[i] <- 0

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
     file = here::here("analysis/data/derived_data/vm1a.RData"))

seapol <- output$seapol
topopath <- output$topop

end_time <- Sys.time()
end_time - start_time

# Define function to count the number of times the sea is simulated as present
# for all raster cells
sea_overlaps <- function(sitearea, seapolygons){
  grid <- st_make_grid(sitearea,
                       cellsize = c(res(sitearea)[1], res(sitearea)[2]),
                       n = c(ncol(sitearea), nrow(sitearea)))

  grid_intersects <- st_intersects(grid, seapol)
  overlaps <- sapply(grid_intersects, length)
  overlapgrid <- data.frame(grid, overlaps) %>%
    st_as_sf()

}

shoremaplv1 <- tmap_grob(tm_shape(overlapgrid, unit = "m") + tm_fill(col = "overlaps", title = "", legend.show = FALSE, palette = "Greys", alpha = 0.5) +
  tm_shape(sitel) + tm_borders(col = "black", lwd = 1) +
  tm_legend(legend.outside = TRUE) + tm_scale_bar(text.size = 0.8))

distplotlv1 <- ggplot() + geom_boxplot(aes(x = "Vertical distance", y = output$results$vertdist)) +
  geom_boxplot(aes(x = "Horisontal distance", y = output$results$hordist)) +
  geom_boxplot(aes(x = "Topographic distance", y = output$results$topodist)) +
  labs(y = "Distance from shoreline (m)", x = "") +
  theme_bw()

tdal <- plot_grid(shoremaplv1, distplotlv1)

vm1aplot <- plot_grid(sitdat, shoremaplv1, distplotlv1, nrow = 1)

# Model comparisons drawing on:
# https://groups.google.com/g/oxcal/c/YRGxmMKWIIo/m/lOUbzBsP7_oJ

# See also discussion at:
# https://groups.google.com/g/oxcal/c/mDVhh8FNnGY/m/6ugk93lyFrEJ
# https://groups.google.com/g/oxcal/c/g5AsmhUBat8/m/RB_zye5fZosJ
# And Bronk Ramsey 2015 Bayesian Approaches to the Building of
# Archaeological Chronologies



# Initial model groups all dates using the Boundary function
sitedates <- siter
dates <- R_Date(sitedates$lab_code, sitedates$c14_bp, sitedates$error)
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
     };
     Boundary("End");
    };
  };')

oxcal_exec <- executeOxcalScript(oxcal_code)
oxcal_read <- readOxcalOutput(oxcal_exec)
oxcal_data <- parseOxcalOutput(oxcal_read, only.R_Date = TRUE)
indices <- parseFullOxcalOutput(oxcal_read)
amodel <- indices$model$modelAgreement


# Second model group all dates overlapping with two-sigma date ranges
ir <- IRanges(as.numeric(siter$sig_2_start_bc), as.numeric(siter$sig_2_end_bc))
siter$group <- subjectHits(findOverlaps(ir, reduce(ir)))


phase_model <- vector()
for(i in 1:length(unique(siter$group))){
  # All dates in present group i
  dats <- filter(siter, group == i)
  txt <- paste0('
         Boundary("', i, '");
         Phase("', i, '")
         {
         ',

         R_Date(dats$lab_code, dats$c14_bp, dats$error),
         '
         };
         '
         )
  phase_model[i] <- txt
}

twosig_model <- paste0('
  Plot()
   {
    Sequence()
    {', paste(phase_model, collapse = " "),
    ' Boundary("End");
    };
  };')

oxcal_exec2 <- executeOxcalScript(twosig_model)
oxcal_read2 <- readOxcalOutput(oxcal_exec2)
indices2 <- parseFullOxcalOutput(oxcal_read2)
amodel2 <- indices2$model$modelAgreement

# Third model combines all dates with overlapping three sigma date ranges
ir <- IRanges(as.numeric(siter$sig_3_start_bc), as.numeric(siter$sig_3_end_bc))
siter$group <- subjectHits(findOverlaps(ir, reduce(ir)))

siter <- siter %>% mutate(group = ifelse(lab_code %in% g1, 1, 2))

g1 <- c("Ua-52878", "Ua-52879", "Ua-52880")
g2 <- c("Ua-52926", "Ua-52925")

phase_model <- vector()
for(i in 1:length(unique(siter$group))){
  # All dates in present group i
  dats <- filter(siter, group == i)
  txt <- paste0('
         Boundary("Start ', i, '");
         Phase("', i, '")
         {
         ',

                R_Date(dats$lab_code, dats$c14_bp, dats$error),
                '
                KDE_Plot("KDE_Plot_', i,'");
         };
         Boundary("End ', i, '");
         '
  )
  phase_model[i] <- txt
}

threesig_model <- paste0('
  Plot()
   {
    Sequence()
    {', paste(phase_model, collapse = " "),
    '};
    };')

oxcal_exec3 <- executeOxcalScript(threesig_model)
oxcal_read3 <- readOxcalOutput(oxcal_exec3)
indices3 <- parseFullOxcalOutput(oxcal_read3)
amodel4 <- indices3$model$modelAgreement

# Find F-model for each model
fmodel <- (amodel/100)^sqrt(nrow(siter))
fmodel2 <-(amodel2/100)^sqrt(nrow(siter))
fmodel3 <-(amodel3/100)^sqrt(nrow(siter))

fmodels <- c(fmodel, fmodel2, fmodel3)

# Choose whichever has the highest value
if(which(fmodels == max(fmodels)) == 1){
  siter$group <- 1
} else if(which(fmodels == max(fmodels)) == 2){
  ir <- IRanges(as.numeric(siter$sig_2_start_bc), as.numeric(siter$sig_2_end_bc))
  siter$group <- subjectHits(findOverlaps(ir, reduce(ir)))
} else{
  ir <- IRanges(as.numeric(siter$sig_3_start_bc), as.numeric(siter$sig_3_end_bc))
  siter$group <- subjectHits(findOverlaps(ir, reduce(ir)))
}

# Rerun final model, this time summing each group
phase_model <- vector()
for(i in 1:length(unique(siter$group))){
  # All dates in present group i
  dats <- filter(siter, group == i)
  txt <- paste0('
         Boundary("', i, '");
         Phase("', i, '")
         {
         ',
          R_Date(dats$lab_code, dats$c14_bp, dats$error),
        if(nrow(dats) > 1){
         paste0('
         Sum("Sum ', i, '");
          };
         ')
        } else{
         paste0('
         };
         ')
        }
  )
  phase_model[i] <- txt
}

final_model <- paste0('
  Plot()
   {
    Sequence()
    {', paste(phase_model, collapse = " "),
                       ' Boundary("End");
    };
  };')

oxcal_exec <- executeOxcalScript(final_model)
oxcal_read <- readOxcalOutput(oxcal_exec)
oxcal_data <- parseOxcalOutput(oxcal_read3, only.R_Date = FALSE)

# Retrieve the raw probabilities for individual 14C-dates
prior <- get_raw_probabilities(oxcal_data)
names(prior) <- paste0(get_name(oxcal_data))
prior <- prior[!(names(prior) %in% c(unique(siter$group),
                                    names(prior)[grep("Start", names(prior))],
                                    names(prior)[grep("End", names(prior))],
                                    names(prior)[grep("KDE", names(prior))],
                                     "character(0)"))]
prior <-  bind_rows(prior, .id = "name")
prior <- left_join(prior, siter, by = c("name" = "lab_code"))
prior$class <- "aprior"

# Retrieve the posterior probabilities for individual 14C-dates and the
# posterior sum.
posterior <- get_posterior_probabilities(oxcal_data)
names(posterior) <- get_name(oxcal_data)
posterior <- posterior[!(names(posterior) %in%
                           c(names(posterior)[grep("Kernel", names(posterior))],
                            names(posterior)[grep("Scale", names(posterior))],
                            names(posterior)[grep("Start", names(posterior))],
                            names(posterior)[grep("End", names(posterior))]))]


posterior <- bind_rows(posterior[!(names(posterior) %in% c(unique(siter$group),
                       "character(0)"))], .id = "name")
posterior <- left_join(posterior, siter, by = c("name" = "lab_code"))
posterior$class <- "bposterior"

# Could not for the life of me get this to work:
# tst <- posterior %>% mutate(group = ifelse(grepl("Sum", name),
#        strsplit(name)[[1]][2], group))
# So I gave up and went for the for-loop
kdes <- unique(posterior[grep("KDE", posterior$name),]$name)
for(i in 1:length(kdes)){
  posterior[posterior$name == kdes[i], "group"] <- as.numeric(
                                      strsplit(kdes[i], "_")[[1]][3])
  posterior[posterior$name == kdes[i], "class"] <- "kde"
}


datedat <- rbind(prior, posterior)
#forcats::fct_relevel

datedat %>%
  arrange(class) %>%
  ggplot() +
  ggridges::geom_ridgeline(aes(x = dates, y = fct_reorder(name, group, .desc = TRUE),
                               height = probabilities * 50,
                               fill = as.factor(group),
                               alpha = fct_reorder(class, group)))+
  scale_alpha_manual(values = c("aprior" = 0.1,
                                "bposterior" = 0.5,
                                "kde" = 1)) +
  theme_bw()


phase_model <- vector()
for(i in 1:length(unique(siter$group))){
  # All dates in present group i
  dats <- filter(siter, group == i)
  txt <- paste0('
         Boundary("Start ', i, '");
         Phase("', i, '")
         {
         ',

                R_Date(dats$lab_code, dats$c14_bp, dats$error),
                '
                Sum("Sum ', i,'");
         };
         Boundary("End ', i, '");
         '
  )
  phase_model[i] <- txt
}

threesig_model <- paste0('
  Plot()
   {
    Sequence()
    {', paste(phase_model, collapse = " "),
                         '};
    };')

oxcal_exec3 <- executeOxcalScript(threesig_model)
oxcal_read3 <- readOxcalOutput(oxcal_exec3)
indices3 <- parseFullOxcalOutput(oxcal_read3)

oxcal_data <- parseOxcalOutput(oxcal_read3, only.R_Date = FALSE)

# Retrieve the raw probabilities for individual 14C-dates
prior <- get_raw_probabilities(oxcal_data)
names(prior) <- paste0(get_name(oxcal_data))
prior <- prior[!(names(prior) %in% c(unique(siter$group),
                                     names(prior)[grep("Start", names(prior))],
                                     names(prior)[grep("End", names(prior))],
                                     names(prior)[grep("Sum", names(prior))],
                                     "character(0)"))]
prior <-  bind_rows(prior, .id = "name")
prior <- left_join(prior, siter, by = c("name" = "lab_code"))
prior$class <- "aprior"

# Retrieve the posterior probabilities for individual 14C-dates and the
# posterior sum.
posterior <- get_posterior_probabilities(oxcal_data)
names(posterior) <- get_name(oxcal_data)
posterior <- posterior[!(names(posterior) %in%
                           c(names(posterior)[grep("Kernel", names(posterior))],
                             names(posterior)[grep("Scale", names(posterior))],
                             names(posterior)[grep("Start", names(posterior))],
                             names(posterior)[grep("End", names(posterior))]))]
posterior <- bind_rows(posterior[!(names(posterior) %in% c(unique(siter$group),
                                                           "character(0)"))], .id = "name")
posterior <- left_join(posterior, siter, by = c("name" = "lab_code"))
posterior$class <- "bposterior"


sums <- unique(posterior[grep("Sum", posterior$name),]$name)
for(i in 1:length(sums)){
  posterior[posterior$name == sums[i], "group"] <- as.numeric(
    strsplit(kdes[i], "_")[[1]][3])
  posterior[posterior$name == sums[i], "class"] <- "sum"
}


datedat <- rbind(prior, posterior)
#forcats::fct_relevel

datedat %>%
  arrange(class) %>%
  ggplot() +
  ggridges::geom_ridgeline(aes(x = dates, y = fct_reorder(name, group, .desc = TRUE),
                               height = probabilities * 50,
                               fill = as.factor(group),
                               alpha = fct_reorder(class, group)))+
  scale_alpha_manual(values = c("aprior" = 0.1,
                                "bposterior" = 0.5,
                                "sum" = 1)) +
  theme_bw()


# Assemble new model, this time summing the modelled dates


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
