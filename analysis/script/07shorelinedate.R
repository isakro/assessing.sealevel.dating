library(tidyverse)
library(patchwork)
library(here)
library(terra)
library(sf)
library(hdrcde)

# For reproducibility
set.seed(321)

# Load required functions and prepared data
source(here("analysis/script/04functions.R"))
load(here("analysis/data/derived_data/01data.RData"))
load(here("analysis/data/derived_data/02data.RData"))
load(here("analysis/data/derived_data/06data.RData"))

# This tries to the load the edited DTM. Due to the file sizes involved,
# this is distributed as tiles which have to be merged by first
# running 00dtm_prep.R
dtm <- rast(here("analysis/data/derived_data/dtm10.tif"))

# Example site (effectively stepping through the code from shoreline_date in
# 04.functions.R)
sitename <- "Hegna vest 1"
sitel <- filter(sites_sa, name == sitename)

sitecurve <- interpolate_curve(years = xvals,
                               target = sitel,
                               dispdat = displacement_curves,
                               isodat = isobases)
sitecurve$name <- sitename
siteelev <- terra::extract(dtm, vect(sitel), fun = mean)[2]

inc <- seq(0, as.numeric(siteelev), 0.001)

expdat <- data.frame(
  offset = inc,
  # Cumulative exponential distribution function
  px = pexp(inc, rate = expfit$estimate)) %>%
  # Probability for each elevation increment is given by subtracting the
  # previous cumulative probability
  mutate(probs = px - lag(px, default =  dplyr::first(px))) %>%
  tail(-1) %>%
  filter(px < 0.99999) # Probability cut-off

dategrid <- data.frame(
  years = seq(-10000, 2000, 1),
  probability = 0)

for(i in 1:nrow(expdat)){
  adjusted_elev <- as.numeric(siteelev - expdat$offset[i])
  if(!(adjusted_elev > 0)) { # Do not allow sea-level to go below present
    adjusted_elev <- 0.01
  }
  # Find lower date
  lowerd <- round(approx(sitecurve[,"lowerelev"],
                         xvals, xout = adjusted_elev)[['y']])

  # Find upper date
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

# Find 95% HDR
dathdr <- hdrcde::hdr(den = list("x" = dategrid$years,
                       "y" = dategrid$probability, prob = 95))
segdat <- data.frame(t(dathdr$hdr))
groupdat <- data.frame(matrix(nrow = 4,
                              ncol = 4))
names(groupdat) <- c("start", "end", "group", "year_median")
groupdat$start <- segdat[seq(1, nrow(segdat), 2), "X95."]
groupdat$end <- segdat[seq(2, nrow(segdat), 2), "X95."]
groupdat$group <- seq(1:length(segdat$hdr[c(TRUE, FALSE)]))
groupdat$year_median <- segdat$mode


plt <- dategrid %>%
  filter(probability != 0) %>%
ggplot() +
  geom_ribbon(data= sitecurve,
              aes(x = years, ymin = lowerelev, ymax = upperelev), fill = "red",
              alpha = 0.2) +
  geom_line(data = sitecurve, aes(x = years, y = upperelev, col = "red")) +
  geom_line(data = sitecurve, aes(x = years, y = lowerelev, col = "red")) +
  ggridges::geom_ridgeline(aes(x = years, y = 1, height = probability * 18000),
                              colour = NA, fill = "darkgrey", alpha = 0.7) +
  geom_segment(data = groupdat, aes(x = start, xend = end,
                                  y = 0, yend = 0), col = "black") +
  labs(y = "Meters above present sea-level", x = "Shoreline date (BCE)",
       title = paste("\U03BB =", as.numeric(round(expfit$estimate, 3)))) +
  scale_x_continuous(breaks = (seq(-9000, -5000, 1000)), expand = c(0, 0),
                     limits = c(-9500, -4500)) +
  coord_cartesian(ylim = c(0, 80)) +
  theme_bw() +
  theme(legend.position = "None")

# , breaks = c(0, 20, 40, as.numeric(siteelev), 80),
# labels = c("0", "20", "40", expression(alpha), "80")

#  For plotting purposes to close the geom_polygon on the y-axis
expdatt <- rbind(c(0, 0, 0), expdat)

# Code taken from oxcAAR to plot density to y-axis
x_extent <- ggplot_build(plt)$layout$panel_scales_x[[1]]$range$range
expdatt$probs_scaled <- expdatt$probs / max(expdatt$probs) *
  diff(x_extent)/6 + x_extent[1]

dplt <- plt + geom_polygon(data = expdatt, aes(x =  probs_scaled,
                                       y = as.numeric(siteelev) - offset),
                   fill = "#046c9a", alpha = 0.6) +
  geom_hline(yintercept = as.numeric(siteelev), linetype = "dashed",
             col = "#046c9a", alpha = 0.6)

ggsave(file = here("analysis/figures/sdate.png"), dplt,
       width = 100, height = 100, units = "mm")
