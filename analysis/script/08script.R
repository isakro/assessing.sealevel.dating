library(ggplot2)
library(ggrepel)
library(patchwork)
library(ggnewscale)
library(dplyr)
library(magrittr)
library(stringr)
library(here)
library(raster)
library(terra)
library(sf)
library(topoDistance)

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
  print(brunlanes$name[i])
  sdatesp[[i]] <- shoreline_date_exp(sitename = brunlanes$name[i],
                                     dtm = dtm,
                                     displacement_curves = displacement_curves,
                                     sites = brunlanes,
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
    comb_median = median(combined, na.rm = TRUE)) %>%
  mutate(siteno = ifelse(site_name == "Pauler 1", "p1", "other"))


# Lasse Jaksland's dates of the Brunlanes sites. He provided a maximum age,
# and argued that the relative uncertainty would be around 50 years.
jaksland <- data.frame(site_name = sort(unique(bdates$site_name)),
                       earliest_date = c(-8850, -9200, -9150, -9000,
                                         -8950, -8975, -8850, -8800, -8950)) %>%
  mutate(latest_date1 = earliest_date + 50,
         latest_date2 = earliest_date + 200)

# Call to plot
bdates <- ggplot(data = hdr, aes(x = comb_median, y = reorder(site_name, -comb_min))) +
  geom_segment(data = hdr, aes(x = comb_min, xend = comb_max,
                               yend = site_name), col = "red", size = 2) +
  ggridges::geom_ridgeline(data = bdates,
                           aes(x = combined, y = site_name,
                               height = probability*50),
                           colour = "grey", fill = "grey") +
  geom_segment(data = hdr, aes(x = comb_min, xend = comb_max,
                               yend = site_name, colour = as.factor(siteno)),
                               size = 2) +
  geom_linerange(data = jaksland, aes(xmin = earliest_date, xmax = latest_date1,
                                      y = site_name), size = 1.75,
                 position = position_dodge(width = 0.3, preserve = 'single'),
                 inherit.aes = FALSE) +
  geom_linerange(data = jaksland, aes(xmin = earliest_date, xmax = latest_date2,
                                      y = site_name), size = 0.5,
                 position = position_dodge(width = 0.3, preserve = 'single'),
                 inherit.aes = FALSE) +
  # geom_segment(data = jaksland, aes(x = earliest_date, xend = latest_date,
  #                                   y = site_name, yend = site_name),
  #              col = "black", size = 5) +
  labs(y = "", x = "BCE") +
  theme_bw() +
  scale_color_manual(values = c("p1" = "gold", "other" = "red")) +
  theme(legend.position = "None")

pauler <- filter(brunlanes, str_detect(name, 'Pauler|Sky'))

sitename <- "Pauler 1"
sitel <- filter(pauler, name == sitename)

sitecurve <- interpolate_curve(years = xvals,
                               isobase1 = sitel$isobase1,
                               isobase2 = sitel$isobase2,
                               target = sitel,
                               dispdat = displacement_curves,
                               isodat = isobases,
                               direction_rel_curve1 = sitel$dir_rel_1)
# Add site name
sitecurve$name <- sitename

# Create bounding box polygon
location_bbox <- bboxpoly(sitel, 1000)

# Use this to clip the dtm to the site area
sitearea <- terra::crop(dtm, location_bbox)

# Retrieve the posterior density estimate for each date group
samplingframe <- filter(bdates, site_name == sitename) %>%
  rename("dates" = "combined", "probabilities" = "probability")
samplingframe$rcarb_cor = "t"

# Simulate sea-level and retrieve distances (uncomment to re-run)
# output <- sample_shoreline(1000, sitel, sitecurve, sitearea,
#                               posteriorprobs = samplingframe)
# save(output,
#      file = here("analysis/data/derived_data/08data.RData"))

load(here("analysis/data/derived_data/08data.RData"))

# Generate grid with dtm resolution holding number of overlaps for each cell
# (also takes quite some time to execute)
simsea <- sea_overlaps(sitearea, output$seapol)

# Code below repurposed from shore_plot() in 04script.R
bboxgrid <- st_bbox(simsea)
anc <- as.numeric(c(bboxgrid$ymin, bboxgrid$xmax))


# Make 0 overlaps NA to remove colour using na.values in call to plot
simsea <- mutate(simsea, overlaps = ifelse(simsea == 0,  NA, simsea))

# Make present day sea-level 0 (for hillshade)
sitearea[sitearea <= 0] <- 0

# Retrieve bounding box for scale bar placement
bboxraster <- st_bbox(sitearea)
anc <- as.numeric(c(bboxraster$ymin, bboxraster$xmax))

# Create hillshade
# For some reason the terra version returned pathchy hillshade
slope <- raster::terrain(raster(sitearea), 'slope')
aspect <- raster::terrain(raster(sitearea), 'aspect')
hill <- raster::hillShade(slope, aspect, 40, 270)

# Now make sea-level NA
sitearea[sitearea <= 0] <- NA

# Make the raster amenable for plotting with ggplot (as.data.frame with terra
# returned some strange errors)
raster_df <- raster::as.data.frame(raster(sitearea), xy = TRUE)
names(raster_df) <- c("x", "y", "value")

pauler1 <- pauler %>%  filter(name == "Pauler 1")

# For labels
site_coords <- as.data.frame(st_coordinates(st_centroid(pauler)))
site_coords$name <- pauler$name

bmap <- ggplot() +
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
  geom_sf(data = simsea, aes(alpha = overlaps), col = NA,
          fill = "#dbe3f3") + ##B6D0E2 ##bfe6ff #"#dbe3f3"
  geom_sf(data = st_centroid(pauler), size = 2, shape = 21, fill = "red") +
  # geom_label_repel(data = site_coords, aes(x = X, y = Y, label = name),
  #                  force = 50, min.segment.length = 0) +
  geom_sf(data = st_centroid(pauler1), size = 2, shape = 21, fill = "gold") +
  scale_alpha_continuous(range = c(0.01, 1), na.value = 0) +
  ggsn::scalebar(data = pauler, dist = 200, dist_unit = "m",
                 transform = FALSE, st.size = 3, height = 0.06,
                 border.size = 0.1, st.dist = 0.07,
                 anchor = c(x = anc[2] - 175, y = anc[1]) + 85) +
  coord_sf(expand = FALSE) +
  theme_bw() + theme(axis.title=element_blank(),
                     axis.text.y = element_blank(),
                     axis.text.x = element_blank(),
                     rect = element_rect(),
                     axis.ticks = element_blank(),
                     panel.grid.major = element_blank(),
                     legend.position = "none")

bdates + bmap

ggsave(file = here("analysis/figures/brunlanes.png"), width = 250,
       height = 120,
       units = "mm")

