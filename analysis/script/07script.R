library(ggplot2)
library(dplyr)
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
  shorelinedates[[i]] <- shoreline_date_exp(sitename = sites_sa$name[i],
                                            dtm = dtm,
                                            displacement_curves = displacement_curves,
                                            features = rcarb_sa,
                                            sites = sites_sa,
                                            isobases = isobases,
                                            sitelimit = TRUE,
                                            expratio =  expfit$estimate)
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
  filter(cumsum(probability) < 0.95) %>%
  summarise(
    comb_min = min(combined, na.rm = TRUE),
    comb_max = max(combined, na.rm = TRUE),
    comb_median = median(combined, na.rm = TRUE))

# Call to plot
shrplt <- ggplot(data = hdr,
                 aes(x = comb_median, y = reorder(site_name, -comb_median))) +
  geom_segment(data = hdr, aes(x = comb_min, xend = comb_max,
                               yend = site_name), col = "red", size = 1) +
  ggridges::geom_ridgeline(data = corshore,
                           aes(x = combined, y = site_name,
                               height = probability*50),
                           colour = "grey", fill = "grey") +
  ggridges::geom_ridgeline(data = corsites,
                           aes(x = dates, y = site_name,
                               height = probabilities*50),
                           colour = "black", fill = NA) +
  labs(x = "cal BCE/CE", y = "", title = paste("\U03BB =",
                                  as.numeric(round(expfit$estimate, 3)))) +
  theme_bw()

ggsave(file = here("analysis/figures/shoredate.png"), shrplt,
       width = 200, height = 200, units = "mm")








  # Call to plot
  ggplot(aes(x = combined, y = site_name)) +
  ggridges::geom_ridgeline(aes(x = dates, height = probability), colour = "grey", fill = "grey") +
  ggridges::geom_ridgeline(aes(x = dates, height = probabilities), colour = "red", fill = "white") +
  labs(x = "cal BCE/CE", y = "") +
  scale_x_continuous(breaks = seq(-9000, 1000, by = 2000)) +
  theme_bw() +
  theme(legend.position = "None")

datedata %>%  ggplot(aes(x = median_date, y = reorder(site_name, -median_date))) +
  ggridges::geom_ridgeline(aes(x = dates, height = probability), colour = "black", fill = "white")






shorelinedates <- list()
for(i in 1:nrow(sites_sa)){
  print(sites_sa$name[i])
  sdates <- shoreline_date(sitename = sites_sa$name[i],
                 dtm = dtm,
                 displacement_curves = displacement_curves,
                 features = rcarb_sa,
                 sites = sites_sa,
                 isobases = isobases,
                 sitelimit = TRUE,
                 xvals = xvals,
                 offset = c(16, 4))
  shorelinedates[[i]] <- c("site_name" = sites_sa$name[i], sdates)
}

sdates <- shorelinedates %>% bind_rows() %>%
  mutate(date_range = as.numeric(latest_date) - as.numeric(earliest_date)) %>%
  rowwise() %>%
  mutate(median_date = median(c(as.numeric(earliest_date), as.numeric(latest_date))))

datedata <- full_join(rdates, sdates, by = "site_name")

shrplt <- datedata %>%
  # Excluded, see supplmenetary
  filter(!(site_name %in% c("Dybdalshei 2", "Lunaveien"))) %>%
  group_by(site_name) %>%
  # Exclude Late Neolithic sites
  filter(n() != sum(is.na(name)) & min(dates) < -2500) %>%
  # Call to plot
ggplot(aes(x = median_date, y = reorder(site_name, -median_date))) +
  geom_segment(aes(x = as.numeric(earliest_date), xend = as.numeric(latest_date),
                   yend = site_name), size = 0.5) +
  ggridges::geom_ridgeline(aes(x = dates, height = probabilities * 100), colour = "red", fill = "white") +
  labs(x = "cal BCE/CE", y = "") +
  scale_x_continuous(breaks = seq(-9000, 1000, by = 2000)) +
  theme_bw() +
  theme(legend.position = "None")

ggsave(file = here("analysis/figures/shoredate.png"), shrplt,
       width = 200, height = 200, units = "mm")



# Exponential decay: y = a(1-b)x
# b =  decay factor
# a = original amount before decay
# x = time
tst <- 1:19
prob = 0.14

plot(tst, (1 * (1 - fitexp$estimate))^tst)


offset = c(0, 19)

sitename <- "Dybdalshei 1"

# Retrieve site features
siter <- filter(features, site_name == sitename)

# If the site limit is to be used
if(sitelimit == TRUE){
  sitel <- filter(sites, name == sitename)
  siteu <- st_union(sitel)
  sitel <- st_as_sf(cbind(siteu, st_drop_geometry(sitel[1,])))
} else {
  # If not create a convex hull around the dated features
  sitel <- st_convex_hull(siter)
  # And assign the values from the limit to the convex hull
  sitel <- left_join(sitel, st_drop_geometry(filter(sites_sa,
                                                    name == sitename)), by = c("site_name" = "name", "ask_id"))
}

sitecurve <- interpolate_curve(years = xvals,
                               isobase1 = sitel$isobase1,
                               isobase2 = sitel$isobase2,
                               target = sitel,
                               dispdat = displacement_curves,
                               isodat = isobases,
                               direction_rel_curve1 = sitel$dir_rel_1)
sitecurve$name <- sitename

siteelev <- terra::extract(dtm, vect(sitel), fun = min)[2]

probs <- c()
prob <- 1
i <- 1
while(prob > 0.00001) {
  prob <- (1 * (1 - (expfit$estimate/10)))^i
  probs[i] <- as.numeric(prob)
  i <- i + 1
}


offsets <- seq(1, length(probs), 1)/10
dates <- data.frame(matrix(ncol = 3, nrow = length(probs)))
names(dates) <- c("earliest_date", "latest_date", "probability")
dates$probability <- probs

for(i in 1:length(offsets)){
  negative_offset <- as.numeric(siteelev - offsets[i])
  if(negative_offset < 1) {
    negative_offset <- 1
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
  earliest <- min(c(lowerd1, upperd1), na.rm = TRUE)
  latest <- max(c(lowerd1, upperd1), na.rm = TRUE)

  dates[i, 1:2] <- cbind(earliest, latest)
}

dates$site_name <- sitename

dates %>%
  # filter(cumsum(probability) < 0.95) %>%
  filter(probability >  0.05) %>%
  ggplot(aes(y = site_name)) +
  # ggridges::geom_ridgeline(aes(x = date, height = probability * 100),
  #                          colour = "black", fill = NA) +
  ggridges::geom_ridgeline(aes(x = earliest_date, height = probability),
                           colour = "blue", fill = NA) +
  ggridges::geom_ridgeline(aes(x = latest_date, height = probability),
                           colour = "red", fill = NA)

dates <- dates %>%
  mutate(combined = (earliest_date + latest_date)/2,
         probability = probability/sum(probability))


dates %>%
  # filter(probability >= 0.05) %>%
  filter(cumsum(probability) < 0.95) %>%
ggplot(aes(y = probability)) +
  ggridges::geom_ridgeline(aes(x = combined, height = probability * 100),
                           colour = "black", fill = "grey") +
  ggridges::geom_ridgeline(aes(x = earliest_date, height = probability * 100),
                           colour = "blue", fill = "blue", alpha = 0.3) +
  ggridges::geom_ridgeline(aes(x = latest_date, height = probability * 100),
                           colour = "red", fill = "red", alpha = 0.3)


nrow(dates)
dates %>%
  filter(cumsum(probability) < 0.95) %>%  nrow()

# Find lower date, subtracting offset (defaults to 0)
lowerd1 <- round(approx(sitecurve[,"lowerelev"],
                        xvals, xout = negative_offset)[['y']])

# Find upper date, subtracting offset (defaults to 0)
upperd1 <- round(approx(sitecurve[,"upperelev"],
                        xvals, xout = negative_offset)[['y']])

# Find lower date, adding upper offset (defaults to 16)
lowerd2 <- round(approx(sitecurve[,"lowerelev"],
                        xvals, xout = positive_offset)[['y']])

# Find upper date, adding upper offset (defaults to 16)
upperd2 <- round(approx(sitecurve[,"upperelev"],
                        xvals, xout = positive_offset)[['y']])

# Find youngest and oldest date
latest <- max(c(lowerd1, upperd1, lowerd2, upperd2), na.rm = TRUE)
earliest <- min(c(lowerd1, upperd1, lowerd2, upperd2), na.rm = TRUE)

shoredates <- c("latest_date" = latest, "earliest_date" = earliest)
return(shoredates)


tst <- shorelinedates[1:2]


tst <-  shoreline_date_exp(sitename = "Dybdalshei 1",
                           dtm = dtm,
                           displacement_curves = displacement_curves,
                           features = rcarb_sa,
                           sites = sites_sa,
                           isobases = isobases,
                           sitelimit = TRUE,
                           expratio =  expfit$estimate)



sdates <- tst %>% bind_rows() %>%
  group_by(site_name) %>%
  mutate(median_date = median(c(as.numeric(earliest_date), as.numeric(latest_date))))





ggsave(file = here("analysis/figures/shoredate.png"), shrplt,
       width = 200, height = 200, units = "mm")



