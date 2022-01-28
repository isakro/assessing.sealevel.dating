library(ggplot2)
library(dplyr)
library(magrittr)
library(stringr)
library(here)
library(terra)
library(sf)

source(here("analysis/script/04script.R"))
load(here("analysis/data/derived_data/01data.RData"))
load(here("analysis/data/derived_data/02data.RData"))
load(here("analysis/data/derived_data/06data.RData"))

# Read in all site limits
sitespol <- read_sf(
  here('analysis/data/raw_data/site_limits/site_limits.gpkg'))

# Read in raster
dtm <- rast("/home/isak/phd/eaa_presentation/dtm10/dtm10.tif")

# Retrieve the Pauler sites
brunlanes <- sitespol %>%
  filter(str_detect(site_name, 'Pauler|Bakke|Sky 1')) %>%
  st_zm() %>%
  rename("name" = "site_name",
         "limit" = "name")

brunlanes <- st_join(st_make_valid(brunlanes), isopolys,
                    join = st_intersects, largest = TRUE)

# Use elevations provided by Jaksland
brunlanes$elevation <- c(108, 98, 127, 124, 114, 110.5, 98, 95, 108)

sdatesp <- list()
for(i in 1:nrow(brunlanes)){
  print(pauler$name[i])
  sdatesp[[i]] <- shoreline_date_exp(sitename = brunlanes$name[i],
                                     dtm = dtm,
                                     displacement_curves = displacement_curves,
                                     sites = pauler,
                                     isobases = isobases,
                                     expratio =  expfit2$estimate,
                                     siteelev = "min",
                                     specified_elev = brunlanes$elevation[i])
}

bdates <- bind_rows(sdatesp)

# Find 95 % probability range for shoreline dates and median shoreline date
# for ordering in the plot
hdr <- bdates %>%  group_by(site_name) %>%
  filter(cumsum(probability) < 0.95) %>%
  summarise(
    comb_min = min(combined, na.rm = TRUE),
    comb_max = max(combined, na.rm = TRUE),
    comb_median = median(combined, na.rm = TRUE))

# Lasse Jaksland's dates of the Brunlanes sites. He provided a maximum age,
# and argued that the relative uncertainty would be around 50 years.
jaksland <- data.frame(site_name = sort(unique(bdates$site_name)),
                       earliest_date = c(-8850, -9200, -9150, -9000,
                                         -8950, -8975, -8850, -8800, -8950)) %>%
  mutate(latest_date1 = earliest_date + 50,
         latest_date2 = earliest_date + 200)

# Call to plot
ggplot(data = hdr,
                 aes(x = comb_median, y = reorder(site_name, -comb_min))) +
  geom_segment(data = hdr, aes(x = comb_min, xend = comb_max,
                               yend = site_name), col = "red", size = 2) +
  ggridges::geom_ridgeline(data = bdates,
                           aes(x = combined, y = site_name,
                               height = probability*50),
                           colour = "grey", fill = "grey") +
  geom_segment(data = hdr, aes(x = comb_min, xend = comb_max,
                               yend = site_name), col = "red", size = 1) +
  geom_linerange(data = jaksland, aes(xmin = earliest_date, xmax = latest_date1,
                                      y = site_name), size = 1.75,
                 position = position_dodge(width = 0.2, preserve = 'single'),
                 inherit.aes = FALSE) +
  geom_linerange(data = jaksland, aes(xmin = earliest_date, xmax = latest_date2,
                                      y = site_name), size = 0.5,
                 position = position_dodge(width = 0.2, preserve = 'single'),
                 inherit.aes = FALSE) +
  # geom_segment(data = jaksland, aes(x = earliest_date, xend = latest_date,
  #                                   y = site_name, yend = site_name),
  #              col = "black", size = 5) +
  labs(y = "", x = "") +
  theme_bw()


plot(st_centroid(st_zm(st_geometry(brunlanes))))

