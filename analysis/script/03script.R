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

# Define function to compare three models of the radiocarbon dates, compare
# F_model values, implement the one with the highest value, and

model_phases <- function(sitedates, manual_groups = NULL){

  # Initial model groups all dates using the Boundary function
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
  ir <- IRanges(as.numeric(sitedates$sig_3_start_bc),
                as.numeric(sitedates$sig_3_end_bc))
  sitedates$group <- subjectHits(findOverlaps(ir, reduce(ir)))

  # Skip model comparison if there is only a single group at three sigma
  if(length(unique(sitedates$group)) > 1){

    # Second model combines all dates with overlapping three sigma date ranges
    ir <- IRanges(as.numeric(sitedates$sig_3_start_bc),
                  as.numeric(sitedates$sig_3_end_bc))
    sitedates$group <- subjectHits(findOverlaps(ir, reduce(ir)))

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


    oxcal_exec2 <- executeOxcalScript(threesig_model)
    oxcal_read2 <- readOxcalOutput(oxcal_exec2)
    indices2 <- parseFullOxcalOutput(oxcal_read2)
    amodel2 <- indices2$model$modelAgreement
  }

  if(!is.null(manual_groups)){
    sitedates$groups <- manual_groups

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

    manual_model <- paste0('
    Plot()
    {
    Sequence()
    {', paste(phase_model, collapse = " "),
                             '};
    };')

    oxcal_exec3 <- executeOxcalScript(threesig_model)
    oxcal_read3 <- readOxcalOutput(oxcal_exec3)
    indices3 <- parseFullOxcalOutput(oxcal_read3)
    amodel3 <- indices3$model$modelAgreement

  }

    # Find F-model for each model
    fmodel <- (amodel/100)^sqrt(nrow(sitedates))
    fmodel2 <-(amodel2/100)^sqrt(nrow(sitedates))
    fmodel3 <-(amodel3/100)^sqrt(nrow(sitedates))

    fmodels <- c(fmodel, fmodel2, fmodel3)

    # Find pseudo bayes factor for all model comparisons
    bfactor <- function(x, y){(x / y) ^ sqrt(nrow(sitedates))}
    bfactors <- outer(fmodels, fmodels, FUN = bfactor)
    rownames(bfactors) <- c("Model 1", "Model 2", "Model 3")

    # Check which model outcompetes the others (no Bayes factor below 1),
    # and return the name
    model <- names(which(!(apply(bfactors, 1, function(x) any(x < 1)))))


    # Choose whichever has the highest value
    if(model == "Model 1"){
      sitedates$group <- 1
    } else if(model == "Model 2"){
      sitedates$group <- subjectHits(findOverlaps(ir, reduce(ir)))
    } else if(model == "Model 3"){
      sitedates$group <- manual_groups
    }
  } else {
    sitedates$group <- 1
  }


  # Rerun final model, this time summing each group of more than one date
  phase_model <- vector()
  for(i in 1:length(unique(sitedates$group))){
    # All dates in present group i
    dats <- filter(sitedates, group == i)
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

  # Add start and end of OxCal code
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
  oxcal_data <- parseOxcalOutput(oxcal_read, only.R_Date = FALSE)

  # Retrieve the raw probabilities for individual 14C-dates
  prior <- get_raw_probabilities(oxcal_data)
  names(prior) <- paste0(get_name(oxcal_data))
  prior <- prior[!(names(prior) %in% c(unique(sitedates$group),
                                       names(prior)[grep("Sum", names(prior))],
                                       "End", "Start",
                                       "character(0)"))]
  prior <-  bind_rows(prior, .id = "name")
  prior <- left_join(prior, sitedates, by = c("name" = "lab_code"))
  prior$class <- "aprior"

  # Retrieve the posterior probabilities for individual 14C-dates and the
  # posterior sums.
  posterior <- get_posterior_probabilities(oxcal_data)
  names(posterior) <- get_name(oxcal_data)
  posterior <- bind_rows(posterior[!(names(posterior) %in%
                    c(unique(sitedates$group), "character(0)"))], .id = "name")
  posterior <- left_join(posterior, sitedates, by = c("name" = "lab_code"))
  posterior$class <- "bposterior"

  # Could not for the life of me get something like this to work:
  # posterior <- posterior %>% mutate(group = ifelse(grepl("Sum", name),
  #        strsplit(name)[[1]][2], group))
  # So I gave up and went for a for-loop
  sums <- unique(posterior[grep("Sum", posterior$name),]$name)
  for(i in 1:length(sums)){
    posterior[posterior$name == sums[i], "group"] <- as.numeric(
      strsplit(sums[i], " ")[[1]][2])
    posterior[posterior$name == sums[i], "class"] <- "sum"
  }

  return(rbind(prior, posterior))

}

# Define function to interpolate displacement curve for a given location,
# based on distance to curves on two provided isobases on the
# north-east to south-west axis perpendicular to the isobases.
interpolate_curve <- function(years, isobase1, isobase2, target, dispdat,
                              isodat, direction_rel_curve1){

  # Distance between "northern" and "southern" isobase
  dist <- st_distance(filter(isodat, name == isobase1),
                      filter(isodat, name == isobase2))

  curve1 <- filter(dispdat, name == isobase1)
  curve2 <- filter(dispdat, name == isobase2)

  # Difference in displacement per meter between the the curves,
  # upper confidence limit
  prm_u <- (dplyr::select(curve1, upperelev) -
              dplyr::select(curve2, upperelev))/ as.numeric(dist)
  # Difference in difference per meter, lower confidence limit
  prm_l <- (dplyr::select(curve1, lowerelev) -
              dplyr::select(curve2, lowerelev))/ as.numeric(dist)

  # Distance to target isobase from isobase of curve 1
  distfrom1 <- st_distance(filter(isodat, name == isobase1),
                           target)

  # Find and return values for the target isobase
  uppervals <- prm_u * as.numeric(distfrom1)
  lowervals <- prm_l * as.numeric(distfrom1)

  # If the direction relative to curve1 is southwest the values are subtracted,
  # if not the values are added.
  if (direction_rel_curve1 == "sw"){
    upperelev <- dplyr::select(curve1, upperelev) - uppervals
    lowerelev <- dplyr::select(curve1, lowerelev) - lowervals
  } else {
    upperelev <- dplyr::select(curve1, upperelev) + uppervals
    lowerelev <- dplyr::select(curve1, lowerelev) + lowervals
  }
  values <- cbind.data.frame(years, upperelev, lowerelev)
  return(values)
}

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
    samplingframe <- list(x = posteriorprobs[["dates"]],
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

# Function that combines all of the above
apply_functions <- function(sitename, dtmpath, displacement_curves, isobases,
                            nsamp = 1000){

  sitel <- filter(sites_sa, name == sitename)
  siter <- filter(rcarb_sa, site_name == sitename)

  # Model dates
  datedat <- model_phases(sitedates = siter)

  # Interpolate displacement curve to the site location
  sitecurve <- interpolate_curve(years = xvals,
                                 isobase1 = sitel$isobase1,
                                 isobase2 = sitel$isobase2,
                                 target = sitel,
                                 dispdat = displacement_curves,
                                 isodat = isobases,
                                 direction_rel_curve1 = sitel$dir_rel_1)
  # Add site name
  sitecurve$name <- sitename

  # Load correct regional raster
  dtm <- load_raster(dtmpath, sitel)

  # Create bounding box polygon
  location_bbox <- bboxpoly(sitel, 250)

  # Use this to clip the dtm to the site area
  sitearea <- terra::crop(dtm, location_bbox)

  # Retrieve the posterior density estimate
  posteriorprobs <- filter(datedat, datedat$class == "Sum")

  # Simulate sea-level and retrieve distances (takes some time)
  output <- sample_shoreline(samps = nsamp, sitel, sitecurve, sitearea,
                             posteriorprobs)

  # Generate grid with dtm resolution holding number of overlaps for each cell
  # (also takes quite some time to execute)
  simsea <- sea_overlaps(sitearea, output$seapol)

  output$simsea <- simsea
  output$datedat <- datedat
  output$sitel <- sitel
  return(output)
}

# Define function to plot site relative to simulated sea-levels
shore_plot <- function(overlapgrid, sitelimit) {
  tm_shape(overlapgrid, unit = "m") +
    tm_fill(col = "overlaps", title = "",
            legend.show = FALSE, palette = "Greys", alpha = 0.5) +
    tm_shape(sitelimit) +
    tm_borders(col = "black", lwd = 1) +
    tm_scale_bar(text.size = 0.8)

}

# Define function to plot boxplots of distance from site to shoreline across
# simulation runs
distance_plot <- function(data) {
  ggplot(data) +
    geom_boxplot(aes(x = "Vertical distance", y = vertdist)) +
    geom_boxplot(aes(x = "Horisontal distance", y = hordist)) +
    geom_boxplot(aes(x = "Topographic distance", y = topodist)) +
    labs(y = "Distance from shoreline (m)", x = "") +
    theme_bw()
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
    labs(title = sitename, y = "", x = "cal BCE") +
  theme_bw() + theme(legend.position = "none")
}
