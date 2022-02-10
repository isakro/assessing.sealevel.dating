library(ggplot2)
library(tidyverse)
library(here)
library(terra)
library(sf)

# Load required functions and prepared data
source(here("analysis/script/04script.R"))
load(here("analysis/data/derived_data/01data.RData"))
load(here("analysis/data/derived_data/02data.RData"))
load(here("analysis/data/derived_data/06data.RData"))

# Read in raster
dtm <- rast("/home/isak/phd/eaa_presentation/dtm10/dtm10.tif")

# List all files except those starting with a number.
# That is, all data files resulting from analysis in 05script.R
datfiles <- grep("^[0-9]", list.files(here("analysis/data/derived_data")),
                 invert = TRUE, value = TRUE)

# Create empty list to hold dat data
posterior_list <- list()

# Loop over, load and assign them to the list
for(i in 1:length(datfiles)){
  # Load date data
  load(file.path(here("analysis/data/derived_data", datfiles[i])))
  posterior_list[[i]] <- output[!(names(output) %in%
                                  c("", "sitel", "sitecurve"))]$datedat %>%
    filter(rcarb_cor == "t") %>%
    filter(group == min(group)) %>%
    # Make sure site name is included on all rows (sums have NA)
    mutate(site_name = unique(na.omit(site_name)))
}

# Unpack list of lists
rdates <- posterior_list %>% bind_rows()
rdates <- rdates %>%  group_by(site_name) %>%
  # If there is a sum, exclude the rest
  filter(class == 'sum' | all(class != 'sum')) %>%
  ungroup() %>%
  filter(class != "aprior")


# Assign closest isobases to site data and exclude sites not analysed
sites_sa <- st_join(st_make_valid(sites_sa), isopolys,
                    join = st_intersects, largest = TRUE) %>%
  filter(!(name %in% c("Dybdalshei 2", "Lunaveien", "Frebergsvik C")))

shorelinedates <- list()
for(i in 1:nrow(sites_sa)){
  print(sites_sa$name[i])
  shorelinedates[[i]] <- shoreline_date_exp_avg(sitename = sites_sa$name[i],
                                            dtm = dtm,
                                            displacement_curves = displacement_curves,
                                            features = rcarb_sa,
                                            sites = sites_sa,
                                            isobases = isobases,
                                            sitelimit = TRUE,
                                            expratio =  expfit$estimate)
}

shorelinedates <- list()
for(i in 1:nrow(sites_sa)){
  print(sites_sa$name[i])
  shorelinedates[[i]] <- shoreline_date_exp(sitename = sites_sa$name[i],
                                            elev = dtm,
                                            disp_curves = displacement_curves,
                                            sites = sites_sa,
                                            iso = isobases,
                                            expratio = expfit$estimate,
                                            siteelev = "mean",
                                            specified_elev = NA)
}

sdates <- bind_rows(shorelinedates)

# Radiocarbon dates corresponding to site inventory and older than 2500 BCE
corsites <- rdates %>%
  # Excluded, see supplmenetary
  filter(!(site_name %in% c("Dybdalshei 2", "Lunaveien"))) %>%
  group_by(site_name) %>%
  # Exclude Late Neolithic sites
  filter(n() != sum(is.na(name)) & min(dates) < -2500)

# Retrieve shoreline date data for the same sites
corshore <- sdates %>%  filter(site_name %in% unique(corsites$site_name))

# Find 95 % probability range for shoreline dates and median shoreline date
# for ordering in the plot
hdr <- corshore %>%  group_by(site_name) %>%
  filter(cumsum(replace_na(probability, 0)) < 0.95) %>%
  summarise(
    year_min = min(year, na.rm = TRUE),
    year_max = max(year, na.rm = TRUE),
    year_median = median(year, na.rm = TRUE))

# Call to plot
shrplt <- ggplot(data = hdr,
                 aes(x = year_median, y = reorder(site_name, -year_median))) +
  geom_segment(data = hdr, aes(x = year_min, xend = year_max,
                               yend = site_name), col = "red", size = 1) +
  ggridges::geom_ridgeline(data = corshore,
                           aes(x = year, y = site_name,
                               height = probability*200),
                           colour = "grey", fill = "grey") +
  ggridges::geom_ridgeline(data = corsites,
                           aes(x = dates, y = site_name,
                               height = probabilities*50),
                           colour = "black", fill = NA) +
  labs(x = "cal BCE/CE", y = "", title = paste("\U03BB =",
                                  as.numeric(round(expfit$estimate, 3)))) +
  theme_bw()

ggsave(file = here("analysis/figures/shoredate.png"), shrplt,
       width = 200, height = 230, units = "mm")


# Example site for development
sitename <- "Alveberget 8" #
sitel <- filter(sites_sa, name == sitename)

sitecurve <- interpolate_curve(years = xvals,
                               target = sitel,
                               dispdat = displacement_curves,
                               isodat = isobases)
sitecurve$name <- sitename

siteelev <- terra::extract(dtm, vect(sitel), fun = min)[2]

# Exponential decay: y = a(1-b)x
# b =  decay factor (identified ratio)
# a = original amount before decay (1 here)
# x = time (distance here)

inc <- seq(0, 70, 0.001)

expdat <- data.frame(
  offset = inc,
  px = pexp(inc, rate = expfit$estimate)) %>%
  mutate(probs = px - lag(px, default =  dplyr::first(px))) %>%
  tail(-1)

dategrid <- data.frame(
  years = seq(-10000, 2000, 1),
  probability = 0)

for(i in 1:nrow(expdat)){
  negative_offset <- as.numeric(siteelev - expdat$offset[i])
  if(!(negative_offset > 0)) {
    negative_offset <- 0.01
  }
  # Find lower date, subtracting offset (defaults to 0)
  lowerd <- round(approx(sitecurve[,"lowerelev"],
                          xvals, xout = negative_offset)[['y']])

  # Find upper date, subtracting offset (defaults to 0)
  upperd <- round(approx(sitecurve[,"upperelev"],
                          xvals, xout =  negative_offset)[['y']])

  # Find youngest and oldest date
  earliest <- min(c(lowerd, upperd))
  latest <- max(c(lowerd, upperd))

  print(c(earliest, latest))

  # Add probability to each year in range
  if(!is.na(latest) && !is.na(earliest)){

    year_range <- seq(earliest, latest, 1)
    prob <- 1/length(year_range)*expdat$probs[i]

    dategrid[dategrid$years %in% year_range, "probability"] <-
      dategrid[dategrid$years %in% year_range, "probability"] + prob
  }
}

dategrid %>%
ggplot() +
  ggridges::geom_ridgeline(aes(x = years, y = 0, height = probability),
                           colour = "black", fill = "grey") +
  theme_bw()



cd <- ggplot() +
  geom_hline(yintercept = as.numeric(siteelev), linetype = "dashed", col = "red") +
  geom_line(data = sitecurve,
            aes(x = years, y = upperelev)) +
  geom_line(data = sitecurve,
            aes(x = years, y = lowerelev)) +
  ylab("Meters above present sea-level") +
  xlab("cal BCE/CE") +
  theme_bw() +
  theme(legend.title = element_blank(), legend.position = "bottom",
        legend.direction = "horizontal")

cd + dt

tst <- shoreline_date("Krøgenes D2", expratio = expfit$estimate)

tst %>% ggplot() +
  ggridges::geom_ridgeline(aes(x = years, y = 0, height = probability),
                           colour = "black", fill = "grey")

dates <- data.frame(matrix(ncol = 2))
names(dates) <- c("year", "probability")

for(i in 1:length(offsets)){
  negative_offset <- as.numeric(siteelev - offsets[i])
  if(!(negative_offset > 0)) {
    negative_offset <- 0.01
  }
  positive_offset <- as.numeric(siteelev + offsets[i])

  # Find lower date, subtracting offset (defaults to 0)
  lowerd1 <- round(approx(sitecurve[,"lowerelev"],
                          xvals, xout = negative_offset)[['y']])

  # Find upper date, subtracting offset (defaults to 0)
  upperd1 <- round(approx(sitecurve[,"upperelev"],
                          xvals, xout =  negative_offset)[['y']])

  lowerd2 <- round(approx(sitecurve[,"lowerelev"],
                          xvals, xout = positive_offset)[['y']])

  upperd2 <- round(approx(sitecurve[,"upperelev"],
                          xvals, xout =  positive_offset)[['y']])

  # Find youngest and oldest date
  earliest <- min(c(lowerd1, upperd1))
  latest <- max(c(lowerd1, upperd1))

  if(!is.na(earliest)){
    yrs <- seq(earliest, latest, 1)
    yrs <- yrs[!(yrs %in% dates$year)]
  } else if(is.na(earliest && latest)) {
    yrs <- NA
  } else {
    seq(max(dates$year, na.rm = TRUE), latest, 1)
    yrs <- yrs[!(yrs %in% dates$year)]
  }

  if(length(yrs) > 0){
    dates <- rbind(dates, cbind(year = yrs, probability = as.numeric(probs[i])))
  }
}



dat1 <- dates %>%
  dplyr::mutate(probability = probability/sum(probability)) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(combined = mean(c(earliest_date, latest_date)))

dat1 <- dates %>%
  mutate(probability = probability/sum(probability),
                combined = (earliest_date + latest_date)/2)

dat1 %>%
  # filter(probability >= 0.05) %>%
  filter(cumsum(probability) < 0.95) %>%
ggplot(aes(y = 0)) +
  ggridges::geom_ridgeline(aes(x = combined, height = probability),
                           colour = "black", fill = "grey") +
  ggridges::geom_ridgeline(aes(x = earliest_date, height = probability),
                           colour = "blue", fill = "grey", alpha = 0.3) +
  ggridges::geom_ridgeline(aes(x = latest_date, height = probability),
                           colour = "red", fill = "grey", alpha = 0.3) +
  labs(y = "Density", x = "BCE")

dat <- dates %>%
  mutate(probability = probability/sum(probability, na.rm = TRUE))
  ggplot() +
  ggridges::geom_ridgeline(aes(x = year, y = 0, height = probability),
                           colour = "black", fill = "grey")


is.na(dates) <- do.call(cbind,lapply(dates, is.infinite))
yrs <- seq(-10500, 1950, 1)
agegrid <- data.frame(year = yrs,
                      earliest_prob = numeric(length(yrs)),
                      latest_prob = numeric(length(yrs)))

erange <- seq(min(dates$earliest_date, na.rm = TRUE),
              max(dates$earliest_date, na.rm = TRUE), 1)
efirst <- dates[match(unique(dates$earliest_date), dates$earliest_date),]

for (i in 1:length(erange)) {
  if(erange[i] %in% efirst$earliest_date) {
    yr <- efirst$earliest_date[which(efirst$earliest_date == erange[i])]
    agegrid[which(agegrid$year == yr),
            "earliest_prob"] <- efirst[which(efirst$earliest_date == yr),"probability"]
    j <- 1
  } else {
    yr <- efirst$earliest_date[which(efirst$earliest_date == erange[i - j])]
    agegrid[which(agegrid$year == yr + j),
            "earliest_prob"] <- efirst[which(efirst$earliest_date == yr),"probability"]
    j <- j + 1
  }
}

lrange <- seq(min(dates$latest_date, na.rm = TRUE),
              max(dates$latest_date, na.rm = TRUE), 1)
lfirst <- dates[match(unique(dates$latest_date), dates$latest_date),]

j <- 1
for (i in 1:length(lrange)) {
  if(lrange[i] %in% lfirst$latest_date) {
    yr <- lfirst$latest_date[which(lfirst$latest_date == lrange[i])]
    agegrid[which(agegrid$year == yr),
            "latest_prob"] <- lfirst[which(lfirst$latest_date == yr),"probability"]
    j <- 1
  } else {
    yr <- lfirst$latest_date[which(lfirst$latest_date == lrange[i - j])]
    agegrid[which(agegrid$year == yr + j),
            "latest_prob"] <- lfirst[which(lfirst$latest_date == yr),"probability"]
    j <- j +1
  }
}


tst <- agegrid %>%
  mutate(probability = earliest_prob + latest_prob) %>%
  mutate(probability = probability/sum(probability))

tst %>%
  ggplot(aes(x= year, y = 0)) +
  ggridges::geom_ridgeline(aes(height = probability), fill = "grey", col = "blue")


agegrid %>%
  # filter(probability >= 0.05) %>%
  ggplot(aes(y = 0)) +
  ggridges::geom_ridgeline(aes(x = year, height = earliest_prob * 50),
                           fill = NA, col = "blue") +
  ggridges::geom_ridgeline(aes(x = year, height = latest_prob * 50),
                           fill = NA, col = "red")



dat2 <- dates %>%
  group_by(r = row_number()) %>%
  mutate(date_range = list(earliest_date:latest_date)) %>%
  ungroup %>% dplyr::select(-r) %>%
  # Unnest the list column to get the desired "long" data frame
  unnest(date_range) %>%
  mutate(probability = probability/sum(probability)) %>%
  dplyr::select(-earliest_date, -latest_date)

dat2 %>%
  # filter(probability >= 0.05) %>%
  filter(cumsum(probability) < 0.95) %>%
  ggplot(aes(y = 0)) +
  ggridges::geom_ridgeline(aes(x = date_range, height = probability * 50),
                           fill = "grey", col = "grey")
