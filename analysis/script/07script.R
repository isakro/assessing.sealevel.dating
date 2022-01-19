library(ggplot2)
library(dplyr)
library(here)
library(terra)
library(sf)

# Load required functions and prepared data
source(here("analysis/script/04script.R"))
load(here("analysis/data/derived_data/01data.RData"))
load(here("analysis/data/derived_data/02data.RData"))

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
# b =  decay proportion (decay factor)
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
  prob <- (1 * (1 - (fitexp$estimate/100)))^i
  probs[i] <- as.numeric(prob)
  i <- i + 1
}


offsets <- seq(1, length(probs), 1)/100
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

dates <- dates %>%
  mutate(probability = probability/sum(probability),
         date = coalesce(latest_date, earliest_date, ))



dates %>%
  filter(cumsum(probability) < 0.95) %>%
ggplot(aes(y = site_name)) +
  ggridges::geom_ridgeline(aes(x = date, height = probability * 100),
                           colour = "black", fill = NA) +
  ggridges::geom_ridgeline(aes(x = earliest_date, height = probability * 100),
                           colour = "blue", fill = NA) +
  ggridges::geom_ridgeline(aes(x = latest_date, height = probability * 100),
                           colour = "red", fill = NA)


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


