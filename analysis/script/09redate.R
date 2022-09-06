library(ggplot2)
library(ggrepel)
library(patchwork)
library(ggnewscale)
library(dplyr)
library(magrittr)
library(stringr)
library(tidyr)
library(here)
library(raster)
library(terra)
library(sf)
library(topoDistance)
library(cowplot)

source(here("analysis/script/04functions.R"))
load(here("analysis/data/derived_data/01data.RData"))
load(here("analysis/data/derived_data/02data.RData"))
load(here("analysis/data/derived_data/06data.RData"))

# Read in all site limits
sitespol <- read_sf(
  here('analysis/data/raw_data/site_limits/site_limits.gpkg'))

# Read in previous shoreline dates and reported elevations
prevdates <- read.csv((here(
  "analysis/data/raw_data/previous_shoreline_dates.csv"))) %>%
  mutate(start = start * -1,
         end = end * -1) %>%
  dplyr::rename("site_name" = "name")

# Read in raster
dtm <- rast("/home/isak/phd/eaa_presentation/dtm10/dtm10.tif")

# Retrieve the Pauler sites
brunlanes <- sitespol %>%
  filter(str_detect(site_name, 'Pauler|Bakke|Sky 1')) %>%
  st_zm() %>%
  dplyr::rename("name" = "site_name",
         "limit" = "name")

brun_rep <- prevdates %>% filter(str_detect(site_name, 'Pauler|Bakke|Sky 1'))

brunlanes <- full_join(brunlanes, brun_rep, by = c("name" = "site_name"))

sdatesp <- list()
for(i in 1:nrow(brunlanes)){
  print(brunlanes$name[i])
  sdatesp[[i]] <- shoreline_date(sitename = brunlanes$name[i],
                                        elev = dtm,
                                        disp_curves = displacement_curves,
                                        sites = brunlanes,
                                        iso = isobases,
                                        expratio = expfit3$estimate,
                                        siteelev = "min",
                                        reso = 0.001,
                                 specified_elev = brunlanes$min_elev[i])
}

bdates <- as_tibble(bind_rows(sdatesp))

# Find 95 % probability range for shoreline dates and median shoreline date
# for ordering in the plot
hdrdat <- data.frame()
for(i in 1:length(unique(bdates$site_name))){
  dat <- hdrcde::hdr(den = list("x" = bdates[bdates$site_name ==
                                         unique(bdates$site_name)[i],
                                       "years"][[1]],
                        "y" = bdates[bdates$site_name ==
                                         unique(bdates$site_name)[i],
                                       "probability"][[1]]), prob = 95)
  segdat <- data.frame(t(dat$hdr))
  groupdat <- data.frame(matrix(nrow = length(segdat[seq(1, nrow(segdat), 2),]),
                                ncol = 5))
  names(groupdat) <- c("site_name", "start", "end", "group", "year_median")

  groupdat$site_name <-  unique(bdates$site_name)[i]
  groupdat$start <- segdat[seq(1, nrow(segdat), 2), "X95."]
  groupdat$end <- segdat[seq(2, nrow(segdat), 2), "X95."]
  groupdat$group <- seq(1:length(dat$hdr[c(TRUE, FALSE)]))
  groupdat$year_median <- dat$mode

  hdrdat <- rbind(hdrdat, groupdat)
}

bdates <- bdates %>%  group_by(site_name) %>%
  filter(cumsum(replace_na(probability, 0)) < 0.99999  &
           probability != 0)

# Lasse Jaksland's dates of the Brunlanes sites. He provided a maximum age,
# and argued that the relative uncertainty would be around 50 years.
jaksland <- brunlanes %>%
  mutate(end2 = start + 200)

# Call to plot
dateplot <- ggplot(data = hdrdat, aes(x = year_median,
                                      y = reorder(site_name, -year_median))) +
  geom_segment(data = hdrdat, aes(x = start, xend = end,
                               yend = site_name), col = NA) +
  ggridges::geom_ridgeline(data = bdates,
                           aes(x = years, y = site_name,
                               height = probability*50),
                           colour = "grey", fill = "grey") +
  geom_linerange(data = hdrdat, aes(xmin = start, xmax = end,
                                      y = site_name), size = 0.5, col = "black",
                 position = position_dodge(width = 0.12, preserve = 'single'),
                 inherit.aes = FALSE) +
  geom_linerange(data = jaksland, aes(xmin = start, xmax = end,
                                      y = name), size = 2, col = "red",
                 position = position_dodge(width = 0.5, preserve = 'single'),
                 inherit.aes = FALSE) +
  geom_linerange(data = jaksland, aes(xmin = start, xmax = end2,
                                      y = name), size = 0.75, col = "red",
                 position = position_dodge(width = 0.5, preserve = 'single'),
                 inherit.aes = FALSE) +
  labs(y = "", x = "BCE", title = paste("\U03BB =",
                          as.numeric(round(expfit3$estimate, 3)))) +
  xlim(c(-10000, -7000)) +
  theme_bw()

pauler <- filter(brunlanes, str_detect(name, 'Pauler|Sky'))

sitename <- "Pauler 1"
sitel <- filter(pauler, name == sitename)

sitecurve <- interpolate_curve(years = xvals,
                               target = sitel,
                               dispdat = displacement_curves,
                               isodat = isobases)

# Add site name
sitecurve$name <- sitename

# Create bounding box polygon
location_bbox <- bboxpoly(sitel, 1000)

# Use this to clip the dtm to the site area
sitearea <- terra::crop(dtm, location_bbox)

samplingframe <- dplyr::filter(bdates, site_name == sitename) %>%
  rename("dates" = "years", "probabilities" = "probability") %>%
  na.omit()
samplingframe$rcarb_cor = "t"

# Simulate sea-level and retrieve distances (uncomment to re-run)
output <- sample_shoreline(1000, sitel, sitecurve, sitearea,
                              posteriorprobs = samplingframe)
# Generate grid with dtm resolution holding number of overlaps for each cell
# (also takes quite some time to execute)
simsea <- sea_overlaps(sitearea, output$seapol)

# Create overview map for inset
imap <- st_read(here('analysis/data/raw_data/naturalearth_countries.gpkg'))

bboxsites <- st_bbox(sites_sa)
bboxsites[1] <- bboxsites[1] - 15000
bboxsites[3] <- bboxsites[3] + 15000
bboxsites[2] <- bboxsites[2] - 5000
bboxsites[4] <- bboxsites[4] + 5000
bboxsitespoly <- st_as_sf(st_as_sfc(bboxsites))

bound_reproj <- st_transform(bboxsitespoly, st_crs(imap))
imap2 <- imap %>%
  dplyr::filter(st_intersects(., bound_reproj, sparse = FALSE))
imap_reproj <- st_transform(imap2, st_crs(sites_sa))

overview <- ggplot() +
  geom_sf(data = imap_reproj, fill = "grey", colour = NA) +
  geom_sf(data = isobases, aes(colour = name)) +
  geom_sf(data = st_centroid(location_bbox), fill = "red",
          size = 2.25,  shape = 21, colour = "black") +
  coord_sf(xlim = c(bboxsites[1], bboxsites[3]),
           ylim = c(bboxsites[2], bboxsites[4]),
           expand = FALSE) +
  scale_colour_manual(values = c("Arendal" = "black",
                                 "Larvik" = "darkgreen",
                                 "Tvedestrand" = "blue",
                                 "Horten" = "darkorange")) +
  theme_nothing() + theme(axis.title=element_blank(),
                     axis.text.y = element_blank(),
                     axis.text.x = element_blank(),
                     rect = element_rect(),
                     axis.ticks = element_blank(),
                     panel.grid.major = element_blank(),
                     legend.position = "none",
                     plot.background = element_rect(fill = "white"),
                     panel.border = element_rect(fill=NA)) +
  scale_x_continuous(expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0)) +
  labs(x = NULL, y = NULL)




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
  geom_sf(data = st_centroid(pauler), size = 1, shape = 21, fill = NA) +
  # geom_label_repel(data = site_coords, aes(x = X, y = Y, label = name),
  #                  force = 50, min.segment.length = 0) +
  geom_sf(data = st_centroid(pauler1), size = 1, shape = 21, fill = "red") +
  scale_alpha_continuous(range = c(0.01, 1), na.value = 0) +
  ggsn::scalebar(data = pauler, dist = 200, dist_unit = "m",
                 transform = FALSE, st.size = 3, height = 0.06,
                 border.size = 0.1, st.dist = 0.09,
                 anchor = c(x = anc[2] - 250, y = anc[1]) + 100) +
  coord_sf(expand = FALSE) +
  theme_bw() + theme(axis.title=element_blank(),
                     axis.text.y = element_blank(),
                     axis.text.x = element_blank(),
                     rect = element_rect(),
                     axis.ticks = element_blank(),
                     panel.grid.major = element_blank(),
                     legend.position = "none")

inset_map <- ggdraw() +
  draw_plot(bmap) +
  draw_plot(overview, x = 0.088, y = 0.64, width = 0.35, height = 0.35)

dateplot + bmap

ggsave(file = here("analysis/figures/brunlanes.png"),
       width = 200, height = 150, units = "mm")

psdates <- prevdates %>%  filter(!str_detect(site_name, 'Pauler|Bakke|Sky 1'))

sites_sl <- site_limits %>% filter(radiocarbon == "f" | radiocarbon == "t" &
                                     rcarb_cor == "f" | name == "Viulsrød 2"
                                   | name == "Anvik")

sites_sl <- sites_sl %>%
  filter(name %in% psdates$site_name)

sites_sl <- inner_join(sites_sl,
                       dplyr::select(psdates, site_name, min_elev, max_elev),
                       by = c("name" = "site_name"))

sitdates <- list()
for(i in 1:nrow(sites_sl)){
  print(c(i, sites_sl$name[i]))

  sitdates[[i]] <- shoreline_date(sitename = sites_sl$name[i],
                                 elev = dtm,
                                 disp_curves = displacement_curves,
                                 sites = sites_sl,
                                 iso = isobases,
                                 expratio = expfit$estimate,
                                 reso = 0.001,
                                 specified_elev = (sites_sl$min_elev[i]+
                                                 sites_sl$max_elev[i])/2)
}

bdates2 <- bind_rows(sitdates) %>% group_by(site_name) %>%
  filter(probability != 0)

# save(bdates, bdates2,
#      file = here("analysis/data/derived_data/09data.RData"))
#
# load(here("analysis/data/derived_data/09data.RData"))

# Find 95 % HDR for shoreline dates and median shoreline date
# for ordering in the plot
hdrdat <- data.frame()
for(i in 1:length(unique(bdates2$site_name))){
  dat <-hdrcde::hdr(den = list("x" = bdates2[bdates2$site_name ==
                                                unique(bdates2$site_name)[i],
                                              "years"][[1]],
                               "y" = bdates2[bdates2$site_name ==
                                                unique(bdates2$site_name)[i],
                                              "probability"][[1]]), prob = 95)
  segdat <- data.frame(t(dat$hdr))
  groupdat <- data.frame(matrix(nrow = length(segdat[seq(1, nrow(segdat), 2),]),
                                ncol = 5))
  names(groupdat) <- c("site_name", "start", "end", "group", "year_median")

  groupdat$site_name <-  unique(bdates2$site_name)[i]
  groupdat$start <- segdat[seq(1, nrow(segdat), 2), "X95."]
  groupdat$end <- segdat[seq(2, nrow(segdat), 2), "X95."]
  groupdat$group <- seq(1:length(dat$hdr[c(TRUE, FALSE)]))
  groupdat$year_median <- dat$mode

  hdrdat <- rbind(hdrdat, groupdat)
}


hdrdat <- hdrdat[order(hdrdat$year_median, decreasing = TRUE),]
hdrdat$site_name <- factor(hdrdat$site_name, levels = unique(hdrdat$site_name))

psdates <- psdates[order(match(psdates$site_name, hdrdat$site_name)),]
psdates$site_name <- factor(psdates$site_name,
                            levels = unique(psdates$site_name))

psdates1 <- psdates %>%  filter(start != end)
psdates2 <- psdates %>%  filter(start == end)

# Call to plot
redateplt <- ggplot(data = hdrdat, aes(x = year_median, y = site_name)) +
  geom_segment(data = hdrdat, aes(x = start, xend = end,
                      yend = site_name), col = NA) +
  # ggridges::geom_ridgeline(data = bdates2,
  #                          aes(x = years, y = site_name,
  #                              height = probability*50),
  #                          colour = "grey", fill = "grey") +
  geom_linerange(data = hdrdat, aes(xmin = start, xmax = end,
                                    y = site_name), col = "grey", size = 2.5) +
# position = position_dodge(width = 0.4, preserve = 'single'),
      # inherit.aes = FALSE
  geom_segment(data = psdates1, aes(x = start, xend = end, yend = site_name),
               size = 1, col = "red", alpha = 1) +
  geom_point(data = psdates2, aes(x = start), col = "red",
             size = 1, alpha = 1) +
  labs(x = "BCE/CE", y = "", title = paste("\U03BB =",
                    as.numeric(round(expfit$estimate, 3)))) +
  scale_x_continuous(breaks = c(seq(-10000, -4000, 2000), -2000, 0, 2000) ) +
  theme_bw()

ggsave(file = here("analysis/figures/redate.png"), redateplt,
       width = 180, height = 250, units = "mm")
