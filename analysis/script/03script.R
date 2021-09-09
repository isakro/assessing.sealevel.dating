library(here)
library(tidyverse)
library(sf)
library(ggridges)
library(tmap)
library(tmaptools)
library(cowplot)
library(terra)
library(oxcAAR)
library(vapour)
library(topoDistance)
quickSetupOxcal()

source(here("analysis/script/02script.R"))
load(here("analysis/data/derived_data/01data.RData"))
load(here("analysis/data/derived_data/02data.RData"))

# Function to identify and load the correct raster
load_raster <- function(dtmfolder, sitelimit) {

  # Number of raster files in total
  nrasts <- length(list.files(dtmfolder))

  # Empty data frame to hold raster data
  rasterexts <- data.frame(matrix(ncol = 5, nrow = nrasts))
  names(rasterexts) <- c("xmin", "xmax", "ymin", "ymax", "path")

  # Loop over and retrieve raster info using the package vapour() without
  # actually loading the raster
  for (i in 1:nrasts){
    # Path to the current raster
    path <- here(list.files(dtmfolder, full.names = TRUE)[i])

    # Retrieve raster extent
    rextent <- vapour_raster_info(path)$extent

    # Populate row in the data frame
    rasterexts[i,] <- c(rextent[1], rextent[2], rextent[3], rextent[4], path)
  }

  sitebbox <- st_bbox(sitelimit)

  # Identify correct raster
  craster <- rasterexts %>% filter(xmin < sitebbox$xmin &
                          xmax > sitebbox$xmax &
                          ymin < sitebbox$ymin &
                          ymax > sitebbox$ymax)

  # Return the correct raster
  return(rast(craster$path))
}


# Define function to model a series of dates using the Boundary
# function from oxcal, sum the posterior density estimate, and return
# all the data for plotting and further analysis (including raw probabilites).
model_dates <- function(sitedates){
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
     Sum("Sum");
     };
     Boundary("End");
    };
  };')

  # Execute the OxCal script and retrieve the results
  oxcal_exec <- executeOxcalScript(oxcal_code)
  oxcal_read <- readOxcalOutput(oxcal_exec)
  oxcal_data <- parseOxcalOutput(oxcal_read, only.R_Date = FALSE)

  # Retrieve the raw probabilities for individual 14C-dates
  prior <- get_raw_probabilities(oxcal_data)
  names(prior) <- paste0(get_name(oxcal_data))
  prior <- prior[!(names(prior) %in% c("Sum", "End", "Start",
                                       "character(0)"))]
  prior <-  bind_rows(prior, .id = "name")
  prior$class <- "aprior"

  # Retrieve the posterior probabilities for individual 14C-dates and the
  # posterior sum.
  posterior <- get_posterior_probabilities(oxcal_data)
  names(posterior) <- get_name(oxcal_data)
  posterior <- bind_rows(posterior[names(posterior) != "character(0)"],
                         .id = "name")
  posterior$class <- "bposterior"
  posterior <- mutate(posterior, class =
                      ifelse(name == "Sum", "Posterior sum", class),
                      name = ifelse(name == "Sum", "Posterior sum", name))

  # Return results
  return(rbind(prior, posterior))
}

# Define function to plot modelled and unmodelled 14C dates
plot_dates <- function(datedata, sitename){
 datedata %>%
  filter(!(name %in% c("Start", "End"))) %>%
  arrange(class) %>%
  ggplot() +
  ggridges::geom_ridgeline(aes(x = dates, y = name,
                               height = probabilities * 100,
                               fill = class, alpha = class)) +
  scale_fill_manual(values = c("Posterior sum" = "grey35",
                               "aprior" = "white", "bposterior" = "grey")) +
  scale_alpha_manual(values = c("Posterior sum" = 1, "aprior" = 0.1,
                                "bposterior" = 0.7)) +
  labs(title = sitename, y = "", x = "cal BCE")
  + theme_bw() + theme(legend.position = "none")
}

# Define function to simulate shoreline and measure distances given a single
# radiocarbon probability density function, site limit, and adjusted
# displacement curve
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

      # Make the value negative to indicate that the site is located below
      # sea level
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

# Define function to count the number of times the sea is simulated as present
# for all raster cells
sea_overlaps <- function(sitearea, seapolygons){
  grid <- st_make_grid(sitearea,
                       cellsize = c(res(sitearea)[1], res(sitearea)[2]),
                       n = c(ncol(sitearea), nrow(sitearea)))

  grid_intersects <- st_intersects(grid, seapolygons)
  overlaps <- sapply(grid_intersects, length)
  overlapgrid <- data.frame(grid, overlaps) %>%
    st_as_sf()
  return(overlapgrid)
}

