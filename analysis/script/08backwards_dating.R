library(tidyverse)
library(patchwork)
library(here)
library(terra)
library(sf)
library(hdrcde)

# For reproducibility
set.seed(1)

# Load required functions and prepared data
source(here("analysis/script/04analysis_functions.R"))
load(here("analysis/data/derived_data/01data.RData"))
load(here("analysis/data/derived_data/02data.RData"))
load(here("analysis/data/derived_data/06data.RData"))

# This tries to the load the edited DTM. Due to the file sizes involved,
# this is distributed as tiles which have to be merged by first
# running 00dtm_prep.R
dtm <- rast(here("analysis/data/derived_data/dtm10.tif"))
dtm <- rast("/home/isak/phd/eaa_presentation/dtm10/dtm10.tif")

# List all files except those starting with a number.
# That is, all data files resulting from analysis in 05script.R
datfiles <- grep("^[0-9]", list.files(here("analysis/data/derived_data")),
                 invert = TRUE, value = TRUE)
# Exclude the temp directory
datfiles <- datfiles[ !grepl("sitearea_temp", datfiles) ]

# Create empty list to hold dat data
posterior_list <- list()
# Loop over, load and assign them to the list
for(i in 1:length(datfiles)){
  # Load date data
  load(file.path(here("analysis/data/derived_data", datfiles[i])))
  posterior_list[[i]] <- output[!(names(output) %in%
                                  c("", "sitel", "sitecurve"))]$datedat %>%
    filter(rcarb_cor == "t") %>%
    # filter(group == min(group)) %>%
    # Make sure site name is included on all rows (sums have NA)
    mutate(site_name = unique(na.omit(site_name)))
}

# Unpack list of lists
rdates <- posterior_list %>% bind_rows()
rdates <- rdates %>%  group_by(site_name, group) %>%
  # If there is a sum, exclude the rest
  filter(class == 'sum' | all(class != 'sum')) %>%
  ungroup() %>%
  filter(class != "aprior") %>%
  filter(!(site_name %in% c("Dybdalshei 2", "Lunaveien", "Frebergsvik C")))

# Exclude sites not analysed
# (see supplementary material for more on this)
sites_sa <- sites_sa %>%
  filter(!(name %in% c("Dybdalshei 2", "Lunaveien", "Frebergsvik C")))

# Specify for which sites to not shoreline date the site limit but rather
# use features (see supplementary material)
feature_sites <- c("Langangen Vestgård 7", "Vallermyrene 2")

# Load shoreline date results (uncomment below to rerun)
load(here("analysis/data/derived_data/08data.RData"))

shorelinedates <- list()
for(i in 1:nrow(sites_sa)){
  print(paste(i, sites_sa$name[i]))
  if(sites_sa$name[i] %in% feature_sites){
    shorelinedates[[i]] <- shoreline_date(sitename = sites_sa$name[i],
                                          elev = dtm,
                                          disp_curves = displacement_curves,
                                          sites = sites_sa,
                                          iso = isobases,
                                          exponential = FALSE,
                                          modelfit = gammafit,
                                          # expratio = expfit$estimate,
                                          siteelev = "mean",
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
                                        exponential = FALSE,
                                        modelfit = gammafit,
                                        # expratio = expfit$estimate,
                                        siteelev = "mean",
                                        reso = 0.001,
                                        specified_elev = NA)
  }
}
# save(shorelinedates,
#      file = here("analysis/data/derived_data/08data.RData"))

sdates <- bind_rows(shorelinedates) %>% group_by(site_name) %>%
  filter(probability != 0)

# Radiocarbon dates corresponding to site inventory.
# First group
corsites <- rdates %>%
  # Excluded, see supplementary
  filter(!(site_name %in% c("Dybdalshei 2", "Lunaveien"))) %>%
  group_by(site_name) %>%
  # Exclude Late Neolithic sites
  filter(n() != sum(is.na(name)) & min(dates) < -1500 & group == min(group))

# Regroup after min(group), some sites have non-correlating first 14C-dates
corgroupsn <- rdates %>%
  # Excluded, see supplementary
  filter(!(site_name %in% c("Dybdalshei 2", "Lunaveien"))) %>%
  group_by(site_name) %>%
  filter(group != min(group)) %>%
  mutate(group2 = ifelse(group == min(group), 2,
                         ifelse(group == min(group) + 1 , 3, 4))) %>%
  filter(n() != sum(is.na(name)))

corgroupsn %>% group_by(site_name) %>% distinct(group)

# Struggled a bit with aes() for geom_ridgeline, so instead
# of aes(col = group2), all three remaining groups are separated here
corgroupsn2 <- corgroupsn %>%  filter(group2 == 2)
corgroupsn3 <- corgroupsn %>%  filter(group2 == 3)
corgroupsn4 <- corgroupsn %>%  filter(group2 == 4)

# Retrieve shoreline date data for the same sites
corshore <- sdates %>%  filter(site_name %in% unique(corsites$site_name))

# Find 95% HDR for shoreline dates and median shoreline date
# for ordering in the plot
hdrdat <- data.frame()
for(i in 1:length(unique(corshore$site_name))){
  dat <-hdrcde::hdr(den = list("x" = corshore[corshore$site_name ==
                                         unique(corshore$site_name)[i],
                                       "years"][[1]],
                        "y" = corshore[corshore$site_name ==
                                         unique(corshore$site_name)[i],
                                       "probability"][[1]]), prob = 95)
  segdat <- data.frame(t(dat$hdr))
  groupdat <- data.frame(matrix(nrow = length(segdat[seq(1, nrow(segdat), 2),]),
                                ncol = 5))
  names(groupdat) <- c("site_name", "start", "end", "group", "year_median")

  groupdat$site_name <-  unique(corshore$site_name)[i]
  groupdat$start <- segdat[seq(1, nrow(segdat), 2), "X95."]
  groupdat$end <- segdat[seq(2, nrow(segdat), 2), "X95."]
  groupdat$group <- seq(1:length(dat$hdr[c(TRUE, FALSE)]))
  groupdat$year_median <- median(as.numeric(dat$mode))

  hdrdat <- rbind(hdrdat, groupdat)
}

# As the displacement curves are collapsed for the last few centuries
# (see Figure 2 in the main text), this causes some artefacts for the end of
# the probability distributions when plotted. These are therefore cut off for,
# plotting here, but does not impact the numerical result for the HDRs (above).
corshore <- corshore %>%
  filter(cumsum(replace_na(probability, 0)) < 0.99)

# Call to plot
shrplt <- ggplot(data = hdrdat,
                 aes(x = year_median, y = reorder(site_name, -year_median))) +
  geom_segment(data = hdrdat, aes(x = start, xend = end,
                               yend = site_name), col = NA) +
  ggridges::geom_ridgeline(data = corshore,
                           aes(x = years, y = site_name,
                               height = probability * 200),
                           col = "grey", fill = "grey") +
  ggridges::geom_ridgeline(data = corgroupsn2,
                           aes(x = dates, y = site_name,
                               height = probabilities * 50),
                           col = "blue", fill = NA) +
  ggridges::geom_ridgeline(data = corgroupsn3,
                           aes(x = dates, y = site_name,
                               height = probabilities * 50),
                           col = "darkorange", fill = NA) +
  ggridges::geom_ridgeline(data = corgroupsn4,
                           aes(x = dates, y = site_name,
                               height = probabilities * 50),
                           col = "forestgreen", fill = NA) +
  ggridges::geom_ridgeline(data = corsites,
                           aes(x = dates, y = site_name,
                               height = probabilities * 50),
                           col = "red", fill = NA) +
  geom_linerange(data = hdrdat, aes(xmin = start, xmax = end,
                                      y = site_name), size = 0.5, col = "black",
                 position = position_dodge(width = 0.3, preserve = 'single'),
                 inherit.aes = FALSE) +
  labs(x = "BCE/CE", y = "") +
  scale_x_continuous(breaks = c(seq(-10000, -4000, 2000), -2500, 0, 2000)) +
  theme_bw()

# Save plot
ggsave(file = here("analysis/figures/shoredate.png"), shrplt,
       width = 200, height = 250, units = "mm")

# Next section finds synchroneity between radiocarbon- and shoreline dates
# Grab the probabilities again without filtering for plotting
shorelinedates <- bind_rows(shorelinedates) %>%
  filter(site_name %in% unique(corsites$site_name))

agedif <- age_diff(shorelinedates, corsites)
agedif2 <- age_diff(shorelinedates, corgroupsn2)
agedif3 <- age_diff(shorelinedates, corgroupsn3)
agedif4 <- age_diff(shorelinedates, corgroupsn4)

phasedt <- agedif %>% group_by(site_name) %>%
  group_by(phase, cross) %>%  summarise(n = n() / 100000) %>%
  pivot_wider(names_from = phase, values_from = n)

phasedt2 <- agedif2 %>% group_by(site_name) %>%
  group_by(phase, cross) %>%  summarise(n = n() / 100000) %>%
  pivot_wider(names_from = phase, values_from = n)

phasedt3 <- agedif3 %>% group_by(site_name) %>%
  group_by(phase, cross) %>%  summarise(n = n() / 100000) %>%
  pivot_wider(names_from = phase, values_from = n)

phasedt4 <- agedif4 %>% group_by(site_name) %>%
  group_by(phase, cross) %>%  summarise(n = n() / 100000) %>%
  pivot_wider(names_from = phase, values_from = n) %>%
  mutate(meso = NA, ln = NA)

synchdt <- rbind(phasedt, phasedt2, phasedt3, phasedt4)

phasedata <- synchdt %>% summarise_if(is.numeric, sum, na.rm = TRUE)

synchdata <- synchdt  %>% group_by(cross) %>%
  summarise_all(sum, na.rm = TRUE) %>%
  ungroup() %>%
  mutate(row_total = rowSums(across(where(is.numeric))))

save(synchdata,
     file = here("analysis/data/derived_data/08agediff.RData"))

hdr_func <- function(x) {
  r <- hdrcde::hdr(x, prob = 95)
  r <- data.frame("ymin" = min(r$hdr), "median" = r$mode,
                  "ymax" =  max(r$hdr))
  return(r)
}

# Create base for age difference plot
# Use shoreline date data to the get the ordering correct first
agedfplt <- ggplot(data = hdrdat,
       aes(x = year_median, y = reorder(site_name, -year_median))) +
  geom_segment(data = hdrdat, aes(x = start, xend = end,
                               yend = site_name), col = NA) +
  geom_vline(xintercept = 0, colour = "grey", linetype = "solid") +
  labs(x = "Age difference", y = "") +
  scale_x_continuous(breaks = c(-5000, 0, 5000)) +
  theme_bw()

agedfplt1 <- agedfplt + stat_summary(data = agedif, aes(y = site_name,
                                           x = age_diff, colour = phase,
                                          width = 0),
                        fun.data = hdr_func, geom = 'errorbar', size = 0.7) +
  stat_summary(data = filter(agedif, cross == "false"), aes(y = site_name,
                                  x = age_diff, colour = phase,
                                  width = 0.7),
               fun.data = hdr_func, geom = 'errorbar', size = 0.5) +
  geom_hline(yintercept = Inf, col = "red", size = 2) +
  xlab("Age difference") +
  theme(legend.position = "none")

agedfplt2 <-  agedfplt + stat_summary(data = agedif2, aes(y = site_name,
                                                         x = age_diff,
                                                         colour = phase,
                                                         width = 0),
                                      fun.data = hdr_func, geom = 'errorbar',
                                      size = 0.7) +
  stat_summary(data = filter(agedif2, cross == "false"), aes(y = site_name,
                                   x = age_diff,
                                   colour = phase,
                                   width = 0.7),
               fun.data = hdr_func, geom = 'errorbar',
               size = 0.7) +
  scale_color_manual(name = "", labels = c("Mesolithic",
                                           "Early and Middle Neolithic",
                                           "Late Neolithic"),
                     values = c("black", "#00BA38", "deeppink1"),
                     guide = "none") +
  geom_hline(yintercept = Inf, col = "blue", size = 2) +
  theme(axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "none")

agedfplt3 <- agedfplt + stat_summary(data = agedif3, aes(y = site_name,
                                                         x = age_diff,
                                                         colour = phase,
                                                         width = 0),
                                     fun.data = hdr_func, geom = 'errorbar',
                                     size = 0.7) +
  stat_summary(data = filter(agedif3, cross == "false"), aes(y = site_name,
                                   x = age_diff,
                                   colour = phase,
                                   width = 0.7),
               fun.data = hdr_func, geom = 'errorbar',
               size = 0.7) +
  scale_color_manual(name = "", labels = c("Mesolithic",
                                           "Early and Middle Neolithic",
                                           "Late Neolithic"),
                     values = c("black", "#00BA38", "deeppink1"),
                     guide = "none") +
  geom_hline(yintercept = Inf, col = "darkorange", size = 2) +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "none")

agedfplt4 <- agedfplt + stat_summary(data = agedif4, aes(y = site_name,
                                                        x = age_diff,
                                                        colour = phase,
                                                        width = 0),
                                    fun.data = hdr_func, geom = 'errorbar',
                                    size = 0.7) +
  scale_color_manual(name = "", labels = c("Mesolithic",
                                           "Early and Middle Neolithic",
                                           "Late Neolithic"),
                     values = c("#00BA38", "deeppink1", "black"),
                     guide = "none") +
  geom_hline(yintercept = Inf, col = "forestgreen", size = 2) +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "none")


splt <- agedfplt1 +
  plot_layout(nrow = 1, guides = "collect") &
  scale_color_manual(name = "", labels = c(paste0("Mesolithic (",
                                                  phasedata$meso, ")"),
                                          paste0("Early and Middle Neolithic (",
                                                  phasedata$pre_ln, ")"),
                                          paste0("Late Neolithic (",
                                                  phasedata$ln, ")")),
                     values = c("black", "#00BA38", "deeppink1")) &
  theme(legend.position = 'bottom')

difplt <- splt +  agedfplt2 + agedfplt3 + agedfplt4 + plot_layout(nrow = 1)

ggsave(file = here("analysis/figures/shoredate2.png"), difplt,
       width = 200, height = 200, units = "mm")

# Retrieve non-correlating dates as well, for reporting in the main text
posterior_list <- list()

for(i in 1:length(datfiles)){
  # Load date data
  load(file.path(here("analysis/data/derived_data", datfiles[i])))
  posterior_list[[i]] <- output[!(names(output) %in%
                                    c("", "sitel", "sitecurve"))]$datedat %>%
    # Make sure site name is included on all rows (sums have NA)
    mutate(site_name = unique(na.omit(site_name)))
}

# Unpack list of lists
ardates <- posterior_list %>% bind_rows()
ardates <- ardates %>%  group_by(site_name, group) %>%
  # If there is a sum, exclude the rest
  filter(class == 'sum' | all(class != 'sum')) %>%
  ungroup() %>%
  filter(class != "aprior") %>%
  filter(!(site_name %in% c("Dybdalshei 2", "Lunaveien", "Frebergsvik C")))

save(rdates, ardates,
     file = here("analysis/data/derived_data/08posterior_rcarbon.RData"))
