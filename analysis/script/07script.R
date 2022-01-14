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
