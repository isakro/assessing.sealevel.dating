# Define function to model a series of dates using the Boundary
# function from oxcal, sum the posterior density estimate, and return
# all the data for plotting and further analysis (including raw probabilities).
model_dates <- function(sitedates, date_groups){

  # Assign group to dates
  sitedates$group <- date_groups

    # as.data.frame(table(date_groups))$Freq # number of dates per group
    # Assemble OxCal code for modelling dates, summing for each group of dates
    phase_model <- vector()
    for(i in 1:length(unique(sitedates$group))){

    dats <- filter(sitedates, group == i)
    txt <- paste0('
      Sequence()
           {
           Boundary("', i, '");
           Phase("', i, '")
           {
           ',
                  R_Date(dats$lab_code, dats$c14_bp, dats$error),
           if(nrow(dats) > 1){
                    paste0('
           Sum("Sum ', i, '");
            };
            Boundary("End ', i, '");
           };
           ')
                  } else{
                    paste0('
           };
           Boundary("End ', i, '");
            };
           ')
                  }
    )
    phase_model[i] <- txt
  }

  # Add start and end of OxCal code
  # This can be reviewed by calling cat(model)
  model <- paste0('
    Phase()
     {', paste(phase_model, collapse = " "),
                  '};')

    # Execture the code and retrieve results
    oxcal_exec <- executeOxcalScript(model)
    oxcal_read <- readOxcalOutput(oxcal_exec)
    oxcal_data <- parseOxcalOutput(oxcal_read, only.R_Date = FALSE)

    # Retrieve the raw probabilities for individual 14C-dates
    prior <- get_raw_probabilities(oxcal_data)
    names(prior) <- paste0(get_name(oxcal_data))
    prior <- prior[!(names(prior) %in% c(unique(sitedates$group),
                                        names(prior)[grep("Sum", names(prior))],
                                        names(prior)[grep("End", names(prior))],
                                        "End", "Start",
                                        "character(0)"))]
    prior <-  bind_rows(prior, .id = "name")
    prior <- left_join(prior, sitedates, by = c("name" = "lab_code"))
    prior$class <- "aprior"
    prior$order <- 1

    # Retrieve the posterior probabilities for individual 14C-dates and the
    # posterior sums.
    posterior <- get_posterior_probabilities(oxcal_data)
    names(posterior) <- get_name(oxcal_data)
    posterior <- bind_rows(posterior[!(names(posterior) %in%
                                        c(unique(sitedates$group),
                                        names(prior)[grep("End", names(prior))],
                                       "Start", "character(0)"))],
                                        .id = "name")
    posterior <- left_join(posterior, sitedates, by = c("name" = "lab_code"))
    posterior$class <- "bposterior"
    posterior$order <- 1

    # Add group number to sums and whether or not it is based on 14-dates
    # believed to correlate with the lithic inventory.
    # Could not for the life of me get something like this to work:
    # posterior <- posterior %>% mutate(group = ifelse(grepl("Sum", name),
    #        strsplit(name)[[1]][2], group))
    # So I gave up and went for a for-loop
    sums <- unique(posterior[grep("Sum", posterior$name),]$name)
    if(length(sums) > 0){
      for(i in 1:length(sums)){
        group <- as.numeric(strsplit(sums[i], " ")[[1]][2])
        # Assign group number
        posterior[posterior$name == sums[i], "group"] <- group
        # Assign the class to be "sum"
        posterior[posterior$name == sums[i], "class"] <- "sum"
        # Retrieve the rcarb_cor value of 14-dates for associated with the sum
        posterior[posterior$name == sums[i], "rcarb_cor"] <-
                           unique(posterior[!is.na(posterior$rcarb_cor) &
                           posterior$group == group,"rcarb_cor"])
        # Add order = 2 for the sum to be plotted after the rest of the group in
        # plot_dates() defined below
        posterior[posterior$name == sums[i], "order"] <- 2
      }
    }


    # Return results
    return(rbind(prior, posterior))

    # If there is only a singe date per group: Calibrate and return
  # } else {
  #   caldat <- oxcalCalibrate(sitedates$c14_bp,
  #                            sitedates$error, sitedates$lab_code)
  #   dat <- get_raw_probabilities(caldat)
  #   names(dat) <- get_name(caldat)
  #   dat_df <- bind_rows(dat, .id = "name")
  #   caldate <- left_join(dat_df, sitedates, by = c("name" = "lab_code"))
  #   return(caldate)
  # }
}

# Define utility function to compare three models of radiocarbon dates
# associated with a site. Either as all dates grouped, only dates overlapping
# at three sigma and alternatively a manually defined model
model_phases <- function(sitedates, date_groups = NULL){

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
  threesgroup <- subjectHits(findOverlaps(ir, reduce(ir)))

  # Skip model comparison if there is only a single group at three sigma
  if(length(unique(threesgroup)) > 1){

    # Second model combines all dates with overlapping three sigma date ranges
    sitedates$group <- threesgroup

    phase_model <- vector()
    for(i in 1:length(unique(sitedates$group))){
      # All dates in present group i
      dats <- filter(sitedates, group == i)
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

  if(!is.null(date_groups)){
    sitedates$group <- date_groups

    phase_model <- vector()
    for(i in 1:length(unique(sitedates$group))){
      # All dates in present group i
      dats <- filter(sitedates, group == i)
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

    oxcal_exec3 <- executeOxcalScript(manual_model)
    oxcal_read3 <- readOxcalOutput(oxcal_exec3)
    indices3 <- parseFullOxcalOutput(oxcal_read3)
    amodel3 <- indices3$model$modelAgreement

   }

    # Find F-model for each model
   fmodels <- c("Model 1" = (amodel/100)^sqrt(nrow(sitedates)))

   if(length(unique(threesgroup)) > 1){
      fmodel2 <-(amodel2/100)^sqrt(nrow(sitedates))
      fmodels <- c(fmodels, "Model 2" = fmodel2)
   }

   if(!is.null(date_groups)){
      fmodel3 <-(amodel3/100)^sqrt(nrow(sitedates))
      fmodels <- c(fmodels, "Model 3" = fmodel3)
   }

    # Find pseudo bayes factor for all model comparisons
   bfactor <- function(x, y){(x / y) ^ sqrt(nrow(sitedates))}
   bfactors <- outer(fmodels, fmodels, FUN = bfactor)
   rownames(bfactors) <- names(fmodels)

    # Check which model out-competes the others (no Bayes factor below 1)
   model <- names(which(!(apply(bfactors, 1, function(x) any(x < 1)))))


    # Assign model groupings to 14c-dates
   if(model == "Model 1"){
      sitedates$group <- 1
   } else if(model == "Model 2"){
      sitedates$group <- threesgroup
   } else if(model == "Model 3"){
      sitedates$group <- date_groups
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
                    c(unique(sitedates$group), "End", "Start",
                      "character(0)"))], .id = "name")
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

  # Return result and model number
  return(list(model, rbind( prior, posterior)))
}

# Group dates based on whether or not they overlap at 99.7% prob
group_dates <- function(data, sitename){
  siter <-  filter(data, site_name == sitename)
  ir <- IRanges(as.numeric(siter$sig_3_start_bc),
                  as.numeric(siter$sig_3_end_bc))
  date_groups <- subjectHits(findOverlaps(ir, reduce(ir)))
  return(date_groups)
}

# Function to interpolate displacement curve using IDW
interpolate_curve <- function(years, target, dispdat, isodat){

  dists <- as.data.frame(st_distance(target, isobases))
  names(dists) <- isobases$name

  values <- data.frame(matrix(ncol = 3, nrow = length(years)))
  names(values) <- c("years", "lowerelev", "upperelev")

  # In the case that a site is on the isobase of a
  # displacement curve (e.g. Alveberget 8), simply return that displacement
  # curve
  if(any(as.numeric(dists) == 0)){
    values <-
      dispdat[dispdat$name == names(dists)[which(as.numeric(dists) == 0)],]

  } else { for(i in 1:length(years)){
    for(j in 1:ncol(dists)){
      le <- dispdat[which(dispdat$name == names(dists)[j] &
                            dispdat$years == years[i]),
                    "lowerelev"]

      ue <- dispdat[which(dispdat$name == names(dists)[j] &
                            dispdat$years == years[i]),
                    "upperelev"]

      dists[2, j] <- le
      dists[3, j] <- ue
    }
    distdat <- as.data.frame(t(dists))
    names(distdat) <- c("distance", "lower", "upper")

      # No sites are older than the lowest limit of any displacement curve
      # so in case of NA values, simply assign NA
    if(any(is.na(distdat))){
      lowerval <- upperval <- NA
    } else {
      # Inverse distance weighting
      lowerval <- sum(apply(distdat, 1,
                            function(x) x["lower"] * x["distance"]^-2)) /
        sum(apply(distdat, 1, function(x) x["distance"] ^-2))
      upperval <- sum(apply(distdat, 1,
                            function(x) x["upper"] * x["distance"]^-2)) /
        sum(apply(distdat, 1, function(x) x["distance"] ^-2))

    }
    values[i,] <- c(years[i], lowerval, upperval)
  }
  }
    return(values)
}

# Function to identify and load the correct raster if it is split up into
#  multiple tiles (not in use).
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
sample_shoreline <- function(samps, sitel, sitecurve, sitearea, posteriorprobs,
                             paths = FALSE){

  results <- data.frame(matrix(ncol = 6, nrow = samps))
  names(results) <- c("vertdist", "hordist", "topodist",
                      "year", "rcarb_cor", "sitename")
  results$rcarb_cor <- unique(na.omit(posteriorprobs$rcarb_cor))
  results$sitename <- unique(sitel$name)
  seapolygons <- list(length = samps)

  # If topographic paths are to be retrieved
  if(paths == TRUE){
    topopaths <- list()
  }

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

    # If these are at the same coordinate (i.e. the polygons are overlapping),
    # return all values as 0
    if(st_overlaps(seapoly, sitel, sparse = FALSE)){
      results$hordist[i] <- 0
      results$topodist[i] <- 0
      results$vertdist[i] <- 0

      # Else...if the entire site polygon is in the sea
    } else if(st_contains(seapoly, sitel, sparse = FALSE)){

      # Retrieve highest elevation on polygon
      siteelev <- extract(sitearea, vect(sitel), fun = max)[2]
      # Find difference to sea level elevation
      results$vertdist[i] <-  max(siteelev[1:nrow(siteelev),]) - sealevel

      # Find the nearest points on site- and seapolygons, making the sea a
      # multiline feature as it contains the site
      distline <- st_nearest_points(sitel, st_cast(seapoly,"MULTILINESTRING"))

      # Retrieve the two points
      distpts <- st_cast(distline, "POINT")

      # Find topographic distance
      coords <- rbind(st_coordinates(distpts[1]),
                      st_coordinates(distpts[2]))
      rownames(coords) <- c("loc", "sea")

      # Find topographic distance (retrieve path)
      topodist <- topoDist(raster::raster(sitearea), coords,
                           paths = TRUE)

      if(paths == TRUE){
        topopaths[i] <- topodist[[2]]
      }


      # Store distance (negative to indicate below sea-level)
      results$topodist[i] <- topodist[[1]]["loc", "sea"] * -1

      # If statement required as the code will fail if the line is of
      # length 0
      if(results$topodist[i] != 0){

        # Retrieve path
        tpath <- st_as_sf(topodist[[2]])
        # topopaths[[i]] <- tpath

        # Find end-points
        endptns <- st_line_sample(tpath, sample = c(0, 1))

        # Cast as line
        hordistl <- st_cast(endptns, "LINESTRING")

        # Find horisontal distance
        results$hordist[i] <- st_length(hordistl) * -1

      } else{
        results$hordist[i] <- 0
      }

      # Else, not overlapping slightly or completely, i.e. the site is
      # located some distance from the sea, perform above steps but without
      # negative values
    } else {

      # Retrieve lowest elevation on polygon
      siteelev <- extract(sitearea, vect(sitel), fun = min)[2]
      # Find difference to sea level elevation
      results$vertdist[i] <-  min(siteelev[1:nrow(siteelev),]) - sealevel

      # Find the nearest points on the site- and seapolygons,
      distline <- st_nearest_points(sitel, seapoly)

      # Retrieve the two points
      distpts <- st_cast(distline, "POINT")

      # Find topographic distance
      coords <- rbind(st_coordinates(distpts[1]),
                      st_coordinates(distpts[2]))
      rownames(coords) <- c("loc", "sea")

      # Find topographic distance (retrieve path)
      topodist <- topoDist(raster::raster(sitearea), coords,
                           paths = TRUE)

      if(paths == TRUE){
        topopaths[i] <- topodist[[2]]
      }

      # Store distance
      results$topodist[i] <- topodist[[1]]["loc", "sea"]

      if(results$topodist[i] != 0){

        # Retrieve path
        tpath <- st_as_sf(topodist[[2]])
        # topopaths[[i]] <- tpath

        # Find end-points
        endptns <- st_line_sample(tpath, sample = c(0, 1))

        # Cast as line
        hordistl <- st_cast(endptns, "LINESTRING")

        # Find horisontal distance
        results$hordist[i] <- st_length(hordistl)
      } else{
        results$hordist[i] <- 0
      }
  }
  }

  if(paths == TRUE){
    return(list(
      results = results,
      seapol = st_as_sf(st_sfc(do.call(rbind, seapolygons)),
                        crs = st_crs(sitel)),
      topopaths = topopaths))
  } else {
    return(list(
      results = results,
      seapol = st_as_sf(st_sfc(do.call(rbind, seapolygons)),
                        crs = st_crs(sitel))))
  }
  # # If there are no topopaths at all, do not try to return these
  # if(any(!(summary(topopaths)[,2]) %in% c("-none-", "sf"))){
  #   return(list(
  #     results = results,
  #     seapol = st_as_sf(st_sfc(do.call(rbind, seapolygons)),
  #                       crs = st_crs(sitel)),
  #     topop = st_as_sf(st_sfc(do.call(rbind, topopaths)),
  #                      crs = st_crs(sitel))
  #   ))
  # } else
  #   return(list(
  #     results = results,
  #     seapol = st_as_sf(st_sfc(do.call(rbind, seapolygons)),
  #                       crs = st_crs(sitel))))
}

# Function to only simulate the sea without finding distances
simulate_sea <- function(nsim, sitel, sitecurve, sitearea, posteriorprobs){

  searasters <- list(length = nsim)

  for(i in 1:nsim){
    # Set up sampling frame from sum of posterior densities
    samplingframe <- list(x = posteriorprobs[["dates"]],
                          y = posteriorprobs[["probabilities"]])
    # Draw single sample date
    sdate <- with(samplingframe, sample(x, size = 1, prob = y))

    # Find sealevel-range at sampled date
    uppers <- approx(sitecurve[,"years"], sitecurve[,"upperelev"],
                     xout = sdate)[["y"]]
    lowers <- approx(sitecurve[,"years"], sitecurve[,"lowerelev"],
                     xout = sdate)[["y"]]

    sealevel <- sample(seq(lowers, uppers, 0.05), 1)

    # Create polygon representing the elevation of the sea level
    rclmat <- matrix(c(sealevel, Inf, NA, 0, sealevel, 1),
                     ncol = 3, byrow = TRUE)
    searasters[[i]] <- raster(terra::classify(sitearea, rclmat))
  }
  rs <- stack(searasters)
  searaster <- calc(rs, sum, na.rm = TRUE)
  return(searaster)
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
apply_functions <- function(sitename, date_groups, dtm, displacement_curves,
                            isobases, nsamp = 1000, loc_bbox, siterpath,
                            rcarbcor_true = FALSE, sitelimit = TRUE){

  # Retrieve site features
  siter <- filter(rcarb_sa, site_name == sitename)

  # If the site limit is to be used
  if(sitelimit == TRUE){
    sitel <- filter(sites_sa, name == sitename)
    siteu <- st_union(sitel)
    sitel <- st_as_sf(cbind(siteu, st_drop_geometry(sitel[1,])))
  } else {
    # If not create a convex hull around the dated features
    sitel <- st_convex_hull(siter)
    # And assign the values from the limit to the convex hull
    sitel <- left_join(sitel, st_drop_geometry(filter(sites_sa,
             name == sitename)), by = c("site_name" = "name", "ask_id"))
    names(sitel)[which(names(sitel) == "site_name")] <- "name"
    sitel <- st_as_sf(sitel)
  }

  # If this option is set to true, only include radiocarbon dates which were
  # evaluated by the original excavators to correlate with typological
  # indicators of the assemblages. (Not in use)
  if(rcarbcor_true == TRUE){
    siter <- filter(siter, rcarb_cor == "t")
  }

  # Model dates
  datedat <- model_dates(sitedates = siter, date_groups = date_groups)

  # Interpolate displacement curve to the site location
  sitecurve <- interpolate_curve(years = xvals,
                                 target = sitel,
                                 dispdat = displacement_curves,
                                 isodat = isobases)
  # Add site name
  sitecurve$name <- sitename

  # Create bounding box polygon
  location_bbox <- bboxpoly(sitel, loc_bbox)

  # Use this to clip the dtm to the site area
  sitearea <- terra::crop(dtm, location_bbox)

  # Define empty list to hold results
  output <- list()
  # Retrieve the posterior density estimate for each date group
  for(i in 1:length(unique(date_groups))){

    # If there is only a single date from the site
    if(nrow(siter) == 1){
      posteriorprobs <- datedat
    # If the length of the group is not longer than one, there is no sum
    # associated with the date
    } else if(length(date_groups[date_groups == i]) <= 1){
      posteriorprobs <- filter(datedat, group == i & class == "bposterior")
    # If the group length is longer than one, use the posterior sum
    } else{
      posteriorprobs <- filter(datedat, group == i & class == "sum")
    }

    # Simulate sea-level and retrieve distances (takes time)
    output[[i]] <- sample_shoreline(samps = nsamp, sitel, sitecurve, sitearea,
                               posteriorprobs)

    # Generate grid with dtm resolution holding number of overlaps for each cell
    # (also takes quite some time to execute)
    simsea <- sea_overlaps(sitearea, output[[i]]$seapol)

    output[[i]]$simsea <- simsea
  }

  output$datedat <- datedat
  output$sitel <- sitel
  output$sitecurve <- sitecurve

  # Store raster (replacing spaces and dashes in the file name)
  writeRaster(sitearea, file.path(siterpath,
                        paste0(str_replace(str_replace(sitename, " ", "_"),
                                           "/", "_"), ".tif")),
              overwrite = TRUE)

  # Return results
  return(output)
}

# Define function to plot site in present day landscape. The various
# "s" parameters are to adjust the scale bar, which needed some manual
# adjustment for nearly every plot
site_plot <- function(locationraster, sitelimit, dist, date_groups,
                      s_tdist, s_xpos, s_ypos, s_bheight, s_tsize) {

  # Make present day sea-level 0 (for hillshade)
  locationraster[locationraster <= 0] <- 0

  # Retrieve bounding box for scale bar placement
  bboxraster <- st_bbox(locationraster)
  anc <- as.numeric(c(bboxraster$ymin, bboxraster$xmax))

  # Create hillshade
  # For some reason the terra version returned pathchy hillshade
  slope <- raster::terrain(raster(locationraster), 'slope')
  aspect <- raster::terrain(raster(locationraster), 'aspect')
  hill <- raster::hillShade(slope, aspect, 40, 270)

  # Now make sea-level NA
  locationraster[locationraster <= 0] <- NA

  # Make the raster amenable for plotting with ggplot (as.data.frame with terra
  # returned some strange errors)
  raster_df <- raster::as.data.frame(raster(locationraster), xy = TRUE)
  names(raster_df) <- c("x", "y", "value")

  # Code partly taken from
  # https://gist.github.com/dirkseidensticker/ce98c6adfe16d5e4590e95c587ea0432
  ggplot() +
    geom_raster(data = raster::as.data.frame(hill, xy = TRUE),
                aes(x = x, y = y, fill = layer)) +
    scale_fill_gradient(low = "black", high = "grey40", na.value = NA) +
    new_scale_fill() +
    geom_raster(data = raster_df, aes(x = x, y = y, fill = value),
                alpha = 0.3) +
    scale_fill_gradient(low = "grey", high = "white", na.value = NA) +
    new_scale_fill() +
    geom_raster(data = raster_df, aes(x = x, y = y, fill = value)) +
    scale_fill_gradient(low = NA, high = NA, na.value = "white") +
    geom_sf(data = sitelimit, fill = NA,
            colour = "red") +
    ggsn::scalebar(data = sitelimit, dist = dist, dist_unit = "m",
                   transform = FALSE, st.size = s_tsize, height = s_bheight,
                   border.size = 0.1, st.dist = s_tdist,
                   anchor = c(x = anc[2] - s_xpos, y = anc[1]) + s_ypos) +
    coord_sf(expand = FALSE) +
    theme_bw() + theme(axis.title=element_blank(),
                       axis.text.y = element_blank(),
                       axis.text.x = element_blank(),
                       rect = element_rect(),
                       axis.ticks = element_blank(),
                       panel.grid.major = element_blank(),
                       legend.position="none")
}

overview_plot <- function(background_map, sitelimit, sites, isobases,
                          os_tdist, os_xpos, os_ypos, os_bheight, os_tsize) {
  bboxsites <- st_bbox(sites)
  bboxsites[1] <- bboxsites[1] - 15000
  bboxsites[3] <- bboxsites[3] + 15000
  bboxsites[2] <- bboxsites[2] - 5000
  bboxsites[4] <- bboxsites[4] + 5000
  bboxsitespoly <- st_as_sf(st_as_sfc(bboxsites))

  bound_reproj <- st_transform(bboxsitespoly, st_crs(background_map))
  bmap2 <- background_map %>%
    filter(st_intersects(., bound_reproj, sparse = FALSE))
  bmap_reproj <- st_transform(bmap2, st_crs(sites))

  anc <- as.numeric(c(bboxsites$ymin, bboxsites$xmax))

  ggplot() +
    geom_sf(data = bmap_reproj, fill = "grey", colour = NA) +
    geom_sf(data = isobases, aes(colour = name)) +
    geom_sf(data = st_centroid(sitelimit), size = 1.5, shape = 21,
            colour = "black", fill = "red") +
    ggsn::scalebar(data = sites, dist = 20, dist_unit = "km",
                   transform = FALSE, st.size = os_tsize, height = os_bheight,
                   border.size = 0.1, st.dist = os_tdist,
                   anchor = c(x = anc[2] - os_xpos, y = anc[1]) + os_ypos) +
    coord_sf(xlim = c(bboxsites[1], bboxsites[3]),
             ylim = c(bboxsites[2], bboxsites[4]),
             expand = FALSE) +
    scale_colour_manual(values = c("Arendal" = "black",
                                   "Larvik" = "darkgreen",
                                   "Tvedestrand" = "blue",
                                   "Horten" = "darkorange")) +
    theme_bw() + theme(axis.title=element_blank(),
                       axis.text.y = element_blank(),
                       axis.text.x = element_blank(),
                       rect = element_rect(),
                       axis.ticks = element_blank(),
                       panel.grid.major = element_blank(),
                       legend.position = "none")
}

# Define function to plot site relative to simulated sea-levels
shore_plot <- function(locationraster, overlapgrid, sitelimit, dist,
                       date_groups, s_tdist, s_xpos, s_ypos, s_bheight,
                       s_tsize) {

  bboxgrid <- st_bbox(overlapgrid)
  anc <- as.numeric(c(bboxgrid$ymin, bboxgrid$xmax))


  # Make 0 overlaps NA to remove colour using na.values in call to plot
  overlapgrid <- mutate(overlapgrid, overlaps = ifelse(overlaps == 0,
                                                       NA, overlaps))

  # Make present day sea-level 0 (for hillshade)
  locationraster[locationraster <= 0] <- 0

  # Retrieve bounding box for scale bar placement
  bboxraster <- st_bbox(locationraster)
  anc <- as.numeric(c(bboxraster$ymin, bboxraster$xmax))

  # Create hillshade
  # For some reason the terra version returned pathchy hillshade
  slope <- raster::terrain(raster(locationraster), 'slope')
  aspect <- raster::terrain(raster(locationraster), 'aspect')
  hill <- raster::hillShade(slope, aspect, 40, 270)

  # Now make sea-level NA
  locationraster[locationraster <= 0] <- NA

  # Make the raster amenable for plotting with ggplot (as.data.frame with terra
  # returned some strange errors)
  raster_df <- raster::as.data.frame(raster(locationraster), xy = TRUE)
  names(raster_df) <- c("x", "y", "value")

  ggplot() +
    geom_raster(data = raster::as.data.frame(hill, xy = TRUE),
                aes(x = x, y = y, fill = layer)) +
    scale_fill_gradient(low = "black", high = "grey40", na.value = NA) +
    new_scale_fill() +
    geom_raster(data = raster_df, aes(x = x, y = y, fill = value),
                alpha = 0.3) +
    scale_fill_gradient(low = "grey", high = "white", na.value = NA) +
    new_scale_fill() +
    geom_raster(data = raster_df, aes(x = x, y = y, fill = value)) +
    scale_fill_gradient(low = NA, high = NA, na.value = "white") +
    new_scale_fill() +
    geom_sf(data = overlapgrid, aes(alpha = overlaps), col = NA,
            fill = "#dbe3f3") + ##B6D0E2 ##bfe6ff
    geom_sf(data = sitelimit, colour = "red", fill = NA) +
    scale_alpha_continuous(range = c(0.01, 1), na.value = 0) +
    ggsn::scalebar(data = sitelimit, dist = dist, dist_unit = "m",
                   transform = FALSE, st.size = s_tsize, height = s_bheight,
                   border.size = 0.1, st.dist = s_tdist,
                   anchor = c(x = anc[2] - s_xpos, y = anc[1]) + s_ypos) +
    coord_sf(expand = FALSE) +
    theme_bw() + theme(axis.title=element_blank(),
                       axis.text.y = element_blank(),
                       axis.text.x = element_blank(),
                       rect = element_rect(),
                       axis.ticks = element_blank(),
                       panel.grid.major = element_blank(),
                       legend.position = "none")
}

# Define function to plot violin plot of distance from site to shoreline across
# simulation runs
distance_plot <- function(data) {
  pivot_longer(data, c(vertdist, hordist, topodist)) %>%
  ggplot(aes(name, value, fill = name)) +
    geom_violin(alpha = 0.5, scale='width') +
    # ggdist::stat_halfeye(alpha = 0.7, adjust = 0.5, width = 0.6,
    #                   justification = -0.2, .width = 0, point_colour = NA) +
    # geom_boxplot(width = 0.12, outlier.size = 0.5) +
    labs(y = "Distance from shoreline (m)", x = "") +
    scale_x_discrete(labels = c('Horizontal', 'Topographic', 'Vertical')) +
    # viridis::scale_fill_viridis(discrete = TRUE) +
    scale_fill_manual(values = c("hordist" = "#00ba38",
                                  "topodist" = "#fad510",
                                  "vertdist" = "#046c9a")) +
    theme_bw() +
    theme(legend.position = "none")
}

# Define function to plot modelled and unmodelled 14C-dates
plot_dates <- function(datedata, sitename, multigroup = TRUE, groupn = NA,
                       title = TRUE){

  # If only one of the groups of dates is being plotted
  if(!is.na(groupn)){
    datedata <- filter(datedata, group == groupn)
  }

  if(length(unique(datedata$name)) > 1){

    plot <-
      datedata %>%
      filter(!grepl("Start|End", name)) %>%
      ggplot() +
      ggridges::geom_ridgeline(aes(x = dates,
      # 3 nested ordering levels: median age, group and sum/individual dates
       y = fct_reorder(fct_reorder(fct_reorder(
          name, sig_3_median_bc, .desc = TRUE),
          order, .desc = TRUE),
          group, .desc = TRUE),
       height = probabilities * 50,
       fill = as.factor(group),
       colour = as.factor(rcarb_cor),
       alpha = class))+
      scale_alpha_manual(values = c("aprior" = 0.05,
                                    "bposterior" = 0.3, #0.4
                                    "sum" = 1)) +
      scale_colour_manual(breaks = c("f", "t", NA),
                          values = c("red", "black", "yellow")) +
      theme_bw() +
      labs(y = "", x = "") +
      theme(legend.position = "none")
  } else{
    plot <- datedata %>%
      filter(!grepl("Start|End", name)) %>%
      ggplot() +
      ggridges::geom_ridgeline(aes(x = dates, y = name,
                                   height = probabilities * 50,
                                   fill = as.factor(group),
                                   alpha = class,
                                   colour = as.factor(rcarb_cor))) +
      scale_colour_manual(breaks = c("f", "t", NA),
                         values = c("red", "black", "yellow")) +
      scale_alpha_manual(values = c("aprior" = 0.05,
                                    "bposterior" = 0.3, #0.4
                                    "sum" = 1)) +
      theme_bw() +
      labs(y = "", x = "") +
      theme(legend.position = "none")
  }

  if(multigroup == TRUE){
    plot <- plot + scale_fill_colorblind()
  } else {
    plot <- plot + scale_fill_manual(
      values = as.character(ggthemes_data$colorblind[groupn, 2]))
  }

  if(title == TRUE){
    plot <- plot + ggtitle(sitename)
  }

  return(plot)
}

# Define function to plot results
plot_results <- function(sitename, sitelimit, datedata, sitearea,
                         background_map, sites, isobases, simulation_output,
                         date_groups, scale_dist, s_tdist = 0.5, s_xpos = 175,
                         s_ypos = 85,  s_bheight = 0.3, s_tsize = 3,
                         os_tdist = 0.03, os_xpos = 20000, os_ypos = 8000,
                         os_bheight = 0.02, os_tsize = 3, adjust_widths = NA){

  # Create list to hold all plots
  plots <- list()

  # Then create plots that are always needed
  plots$dplot <- plot_dates(datedata, sitename, title = FALSE)
  plots$pmap <- site_plot(sitearea, sitelimit, scale_dist, date_groups,
                          s_tdist, s_xpos, s_ypos, s_bheight, s_tsize)
  plots$omap <- overview_plot(background_map, sitelimit, sites, isobases,
                              os_tdist, os_xpos, os_ypos, os_bheight, os_tsize)

  # Then loop over and create plots per group of dates
  for(i in 1:length(unique(date_groups))){

      # Radiocarbon dates
      plots[[paste0("datplot_", i)]] <- plot_dates(
        datedata, sitename, multigroup = FALSE, group = i, title = FALSE)

      # Map of simulated sea-levels
      plots[[paste0("seaplot_", i)]] <- shore_plot(sitearea,
        simulation_output[[i]]$simsea, sitelimit, scale_dist, date_groups,
        s_tdist, s_xpos, s_ypos, s_bheight, s_tsize)

      # Boxplot of distances from site to sea
      plots[[paste0("distplot_", i)]] <- distance_plot(
        simulation_output[[i]]$results)
  }

  # If the number of date groups is more than one the plots are assembled
  #  a bit differently (mainly so that each individual group of dates get their
  # own plot for inspection)
  if(length(unique(date_groups)) > 1){
    # For some reason the widths are sometimes completely off, so I
    # added an option to specify these.
    if(is.na(adjust_widths)){
      wrap_plots(plots,
      nrow = length(unique(date_groups)) + 1, ncol = 3) +
      plot_annotation(title = sitename) &
      theme(plot.title = element_text(hjust = 0.5))
    } else {
      wrap_plots(plots,
        nrow = length(unique(date_groups)) + 1, ncol = 3,
        widths = adjust_widths) +
        plot_annotation(title = sitename) &
        theme(plot.title = element_text(hjust = 0.5))
    }
  } else{
      ( plots$pmap | plots$omap )
      ( plots$pmap | plots$omap ) /
      ( plots$datplot_1 | plots$seaplot_1 | plots$distplot_1 ) +
      plot_annotation(title = sitename) &
      theme(plot.title = element_text(hjust = 0.5))
    }
}

# Shoreline date, using exponential function for elevation offset
shoreline_date <- function(sitename, elev = dtm,
                           disp_curves = displacement_curves,
                           sites = sites_sa,
                           iso = isobases,
                           reso = 0.1, exponential = TRUE,
                           modelfit, siteelev = "mean",
                           specified_elev = NA,
                           sitelimit = TRUE,
                           features = NA){

  # If the site limit is to be used
  if(sitelimit == TRUE){
    sitel <- filter(sites, name == sitename)
    siteu <- st_union(sitel)
    sitel <- st_as_sf(cbind(siteu, st_drop_geometry(sitel[1,])))
  } else {
    # If site features and not site limit is to be used
    siter <- filter(features, site_name == sitename)

    # Create a convex hull around the dated features
    sitel <- st_convex_hull(siter)
    # And assign the values from the limit to the convex hull
    sitel <- left_join(sitel, st_drop_geometry(filter(sites_sa,
                  name == sitename)), by = c("site_name" = "name", "ask_id"))
    names(sitel)[which(names(sitel) == "site_name")] <- "name"
    sitel <- st_as_sf(sitel)
  }

  sitecurve <- interpolate_curve(years = xvals,
                                 target = sitel,
                                 dispdat = disp_curves,
                                 isodat = iso)

  if(!(is.na(specified_elev))){
    siteelev <- specified_elev
  } else{
    if(siteelev == "mean") {
      siteelev <- terra::extract(elev, vect(sitel), fun = mean)[2]
    } else if(siteelev == "min"){
      siteelev <- terra::extract(elev, vect(sitel), fun = min)[2]
    }
  }

  if(exponential){
  inc <- seq(0, as.numeric(siteelev), reso)
  expdat <- data.frame(
    offset = inc,
    px = pexp(inc, rate = modelfit$estimate)) %>%
    mutate(probs = px - lag(px, default =  dplyr::first(px))) %>%
    tail(-1) %>%
    filter(px < 0.99999) # Probability cut-off
  }else{
  inc <- seq(0.001, as.numeric(siteelev), reso)
  expdat <- data.frame(
    offset = inc,
    px = pgamma(inc, shape = modelfit$estimate[1],
                     rate = modelfit$estimate[2])) %>%
    mutate(probs = px - lag(px, default =  dplyr::first(px))) %>%
    tail(-1) %>%
    filter(px < 0.99999) # Probability cut-off
  }

  dategrid <- data.frame(
    years = seq(-10000, 2000, 1),
    probability = 0)

  for(i in 1:nrow(expdat)){
    adjusted_elev <- as.numeric(siteelev - expdat$offset[i])
    if(!(adjusted_elev > 0)) {
      adjusted_elev <- reso # Just to make sure sea-level is never
                            #lower than present (unnecessary due to inc?)
    }
    # Find lower date, subtracting offset (defaults to 0)
    lowerd <- round(approx(sitecurve[,"lowerelev"],
                           xvals, xout = adjusted_elev)[['y']])

    # Find upper date, subtracting offset (defaults to 0)
    upperd <- round(approx(sitecurve[,"upperelev"],
                           xvals, xout = adjusted_elev)[['y']])

    # Find youngest and oldest date
    earliest <- min(c(lowerd, upperd))
    latest <- max(c(lowerd, upperd))

    # Add probability to each year in range
    if(!is.na(latest) && !is.na(earliest)){

      year_range <- seq(earliest, latest, 1)
      prob <- 1/length(year_range)*expdat$probs[i]

      dategrid[dategrid$years %in% year_range, "probability"] <-
        dategrid[dategrid$years %in% year_range, "probability"] + prob
    }
  }
  dategrid$site_name <- sitename

  return(dategrid)
}

# Function to find age difference between shoreline and radiocarbondates
age_diff <- function(shorelinedata, radiocarbondata, nsamp = 100000){

  shorelinedata <- filter(shorelinedata,
                          site_name %in% unique(radiocarbondata$site_name))

  prob95 <- radiocarbondata %>% group_by(site_name) %>%
    summarise(
      year_min = min(hdrcde::hdr(den = list("x" = dates, "y" = probabilities),
                                 prob = 95)$hdr),
      year_max = max(hdrcde::hdr(den = list("x" = dates, "y" = probabilities),
                                 prob = 95)$hdr)
    )

  agedif <- setNames(data.frame(matrix(ncol = 6, nrow = 0)),
                     c("age_diff", "site_name", "year_min",
                       "year_max", "phase", "cross"))

  for(i in 1:length(unique(shorelinedata[,"site_name"]))){
    samplingframe <- list(x = shorelinedata[shorelinedata$site_name ==
                                        unique(shorelinedata[,"site_name"])[i],
                                            "years"],
                          y = shorelinedata[shorelinedata$site_name ==
                                        unique(shorelinedata[,"site_name"])[i],
                                            "probability"])
    sdate <- with(samplingframe, sample(x, size = nsamp, prob = y,
                                        replace = TRUE))

    samplingframe <- list(x = radiocarbondata[radiocarbondata$site_name ==
                                        unique(shorelinedata[,"site_name"])[i],
                                              "dates"]$dates,
                          y = radiocarbondata[radiocarbondata$site_name ==
                                        unique(shorelinedata[,"site_name"])[i],
                                              "probabilities"]$probabilities)

    rdate <- with(samplingframe, sample(x, size = nsamp, prob = y,
                                        replace = TRUE))
    age_diff <- rdate - sdate
    agedif <- rbind(agedif, data.frame(age_diff,
                    site_name =  unique(shorelinedata$site_name)[i],
                    year_min = prob95[prob95$site_name ==
                                        unique(shorelinedata$site_name)[i],
                                      "year_min"][[1]],
                    year_max = prob95[prob95$site_name ==
                                        unique(shorelinedata$site_name)[i],
                                      "year_max"][[1]],
                    phase = ifelse(prob95[prob95$site_name ==
                                            unique(shorelinedata$site_name)[i],
                                          "year_min"][[1]] < -4000, "meso",
                    ifelse(prob95[prob95$site_name ==
                                    unique(shorelinedata$site_name)[i],
                                  "year_min"][[1]] < -2500, "pre_ln", "ln")),
                    cross = ifelse(between(0, min(hdrcde::hdr(age_diff,
                                                              prob = 95)$hdr),
                                           max(hdrcde::hdr(age_diff,
                                       prob = 95)$hdr)), "true", "false")))
  }

  agedif$phase <- factor(agedif$phase, levels=c("meso", "pre_ln", "ln"))
  return(agedif)
}
