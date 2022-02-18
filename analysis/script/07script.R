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


# Exclude sites not analysed
sites_sa <- sites_sa %>%
  filter(!(name %in% c("Dybdalshei 2", "Lunaveien", "Frebergsvik C",
                       "Bratsberg", "Lar\xf8nningen")))

# Specify for which sites to not shoreline date the site limit
# (see supplementary material for more on this)
feature_sites <- c("Langangen Vestgård 7", "Vallermyrene 2")

shorelinedates <- list()
for(i in 1:nrow(sites_sa)){
  print(c(i, sites_sa$name[i]))
  if(sites_sa$name[i] %in% feature_sites){
    shorelinedates[[i]] <- shoreline_date(sitename = sites_sa$name[i],
                                          elev = dtm,
                                          disp_curves = displacement_curves,
                                          sites = sites_sa,
                                          iso = isobases,
                                          expratio = expfit$estimate,
                                          siteelev = "min",
                                          reso = 0.001,
                                          specified_elev = NA,
                                          sitelimit = FALSE,
                                          features = rcarb_sa)
  } else {
  shorelinedates[[i]] <- shoreline_date(sitename = sites_sa$name[i],
                                        elev = dtm,
                                        disp_curves = displacement_curves,
                                        sites = sites_sa,
                                        iso = isobases,
                                        expratio = expfit$estimate,
                                        siteelev = "min",
                                        reso = 0.001,
                                        specified_elev = NA)
  }
}
l
sdates <- bind_rows(shorelinedates) %>% group_by(site_name) %>%
  filter(cumsum(replace_na(probability, 0)) < 0.95)

# Radiocarbon dates corresponding to site inventory and older than 2500 BCE
corsites <- rdates %>%
  # Excluded, see supplementary
  filter(!(site_name %in% c("Dybdalshei 2", "Lunaveien"))) %>%
  group_by(site_name) %>%
  # Exclude Late Neolithic sites
  filter(n() != sum(is.na(name)) & min(dates) < -2500)

# Retrieve shoreline date data for the same sites
corshore <- sdates %>%  filter(site_name %in% unique(corsites$site_name))

# Find 95 % probability range for shoreline dates and median shoreline date
# for ordering in the plot
hdr <- corshore %>%  group_by(site_name) %>%
  filter(cumsum(replace_na(probability, 0)) < 0.95 &
           probability != 0) %>%
  summarise(
    year_min = min(years, na.rm = TRUE),
    year_max = max(years, na.rm = TRUE),
    year_median = median(years, na.rm = TRUE))

# Call to plot
shrplt <- ggplot(data = hdr,
                 aes(x = year_median, y = reorder(site_name, -year_median))) +
  geom_segment(data = hdr, aes(x = year_min, xend = year_max,
                               yend = site_name), col = "red", size = 0.75) +
  ggridges::geom_ridgeline(data = corshore,
                           aes(x = years, y = site_name,
                               height = probability * 200),
                           colour = "grey", fill = "grey") +
  ggridges::geom_ridgeline(data = corsites,
                           aes(x = dates, y = site_name,
                               height = probabilities * 50),
                           colour = "black", fill = NA) +
  labs(x = "cal BCE/CE", y = "", title = paste("\U03BB =",
                                  as.numeric(round(expfit$estimate, 3)))) +
  theme_bw()

ggsave(file = here("analysis/figures/shoredate.png"), shrplt,
       width = 200, height = 230, units = "mm")


# Example site for development
sitename <- "Dybdalshei 1" #
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

inc <- seq(0, 80, 0.1)

expdat <- data.frame(
  offset = inc,
  px = pexp(inc, rate = expfit$estimate)) %>%
  mutate(probs = px - lag(px, default =  dplyr::first(px))) %>%
  tail(-1) %>%
  filter(px < 0.99999) # Probability cut-off

dategrid <- data.frame(
  years = seq(-10000, 2000, 1),
  probability = 0)

for(i in 1:nrow(expdat)){
  adjusted_elev <- as.numeric(siteelev - expdat$offset[i])
  if(!(adjusted_elev > 0)) {
    adjusted_elev <- 0.01
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

  print(c(earliest, latest))

  # Add probability to each year in range
  if(!is.na(latest) && !is.na(earliest)){

    year_range <- seq(earliest, latest, 1)
    prob <- 1/length(year_range)*expdat$probs[i]

    dategrid[dategrid$years %in% year_range, "probability"] <-
      dategrid[dategrid$years %in% year_range, "probability"] + prob
  }
}

system.time(dategrid <- shoreline_date("Dybdalshei 1", expratio = expfit$estimate, res = 0.001))

dategrid %>%
ggplot() +
  ggridges::geom_ridgeline(aes(x = years, y = 0, height = cumsum(probability)),
                           colour = "black", fill = "grey") +
  theme_bw()

dategrid %>%
  ggplot() +
  ggridges::geom_ridgeline(aes(x = years, y = 0, height = probability),
                           colour = "black", fill = "grey") +
  theme_bw()
