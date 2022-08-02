library(here)
library(dplyr)
library(purrr)
library(stringr)
library(ggplot2)
library(grid)
library(gridExtra)
library(gtable)
library(cowplot)
library(MASS)

# List all files except those starting with a number.
# That is, all data files resulting from analysis in 05script.R
datfiles <- grep("^[0-9]", list.files(here("analysis/data/derived_data")),
     invert = TRUE, value = TRUE)

# Create empty list to hold results
results_list <- list()

# Loop over, load results and assign them to the list
for(i in 1:length(datfiles)){
  # Load results
  load(file.path(here("analysis/data/derived_data", datfiles[i])))
  results_list[[i]] <- output[!(names(output) %in%
                             c("datedat", "sitel", "sitecurve"))]
}

# Sites and phases to label as agricultural
agr <- c("Nordby 1", "Kvastad A2", "Nauen A")

# results_list is a list of lists of variable length that need to be unpacked.
# (Probably a smoother purrr solution to do all of this?)
results <- list()
for(i in 1:length(results_list)){
  dat <- results_list[[i]]
  tmp <- list()
  for(j in 1:length(dat)){
    pdat <- dat[[j]]["results"][[1]]
    if(is.element(unique(pdat$sitename), agr) & any(pdat$year > -4000)){
        pdat$type <- "agricultural"
      } else{
        pdat$type <- "forager"
      }
    tmp[[j]] <- pdat
  }
  results[[i]] <- tmp
}

# Collapse list of lists into a single data frame
distances <- results %>% bind_rows()
# distances <-  distances$results

# Lunaveien and Dybdalshei 2 are excluded
# (see site evaluations in the supplementary material)

# Negative values for Gunnarsrød 5 and Pjonkerød R1 are to be set to zero
# (see site evaluations in the supplementary material)
distances <- distances %>%
  mutate(vertdist = ifelse(sitename %in% c("Gunnarsrød 5", "Pjonkerød R1")
                           & vertdist < 0,  0, vertdist),
         hordist = ifelse(sitename %in% c("Gunnarsrød 5", "Pjonkerød R1")
                          & hordist < 0,  0, hordist),
         topodist = ifelse(sitename %in% c("Gunnarsrød 5", "Pjonkerød R1")
                           & topodist < 0,  0, topodist))

sum_all <- distances %>%
  summarise_at(c("hordist", "topodist", "vertdist"),
               c(mean, median, sd, IQR)) %>%
  round()

sum_cor <- distances %>% filter(rcarb_cor == "t") %>%
  summarise_at(c("hordist", "topodist", "vertdist"),
               c(mean, median, sd, IQR)) %>%
   round()

# Colour for horisontal distance
colsh = c("forager" = "#00ba38", "agricultural" = "darkgrey")
# Colour for topographic distance
colst = c("forager" = "#fad510", "agricultural" = "darkgrey") # "yellow3"
# Colour for vertical distance
colsv = c("forager" = "steelblue", "agricultural" = "darkgrey") # "#046c9a"

# Horizontal distance, all dates to the Stone Age
h1 <- ggplot(distances, aes(x = year, y = hordist,
                            colour = type)) +
  geom_point(shape = 16, alpha = 0.01) +
  scale_colour_manual(values = colsh) +
  labs(y = "Horizontal distance (m)", x = "BCE",
       title = "All dates") +
  scale_x_continuous(breaks = seq(-8000, -1000, by = 2000)) +
  theme_bw() +
  theme(legend.position = "none")

th1 <- tableGrob(sum_all[c("hordist_fn1", "hordist_fn2",
                         "hordist_fn3", "hordist_fn4")],
                 cols = c("Mean", "Median",
                          "SD", "IQR"),
                 rows = NULL, theme = ttheme_minimal())

# Find table height for plotting
th <- sum(th1$heights)
# Arrange plot and table
hor1 <- grid.arrange(h1, th1, heights = unit.c(unit(1, "null"), th))

# Horizontal distance, dates corresponding to inventory
h2 <- distances %>%  filter(rcarb_cor == "t") %>%
  ggplot(aes(x = year, y = hordist,
                        colour = type)) +
  geom_point(shape = 16, alpha = 0.01) +
  scale_colour_manual(values = colsh) +
  labs(y = "Horizontal distance (m)", x = "BCE",
       title = "Corresponding dates") +
  scale_x_continuous(breaks = seq(-8000, -1000, by = 2000)) +
  theme_bw() +
  theme(legend.position = "none")

th2 <- tableGrob(sum_cor[c("hordist_fn1", "hordist_fn2",
                           "hordist_fn3", "hordist_fn4")],
                 cols = c("Mean", "Median",
                          "SD", "IQR"),
                 rows = NULL, theme = ttheme_minimal())

hor2 <- grid.arrange(h2, th2, heights = unit.c(unit(1, "null"), th))

# Repeat for topographic distance
t1 <- ggplot(distances, aes(x = year, y = topodist,
                            colour = type)) +
  geom_point(shape = 16, alpha = 0.01) +
  scale_colour_manual(values = colst) +
  labs(y = "Topographic distance (m)", x = "BCE", title = "") +
  scale_x_continuous(breaks = seq(-8000, -1000, by = 2000)) +
  theme_bw() +
  theme(legend.position = "none")

tt1 <- tableGrob(sum_all[c("topodist_fn1", "topodist_fn2",
                           "topodist_fn3", "topodist_fn4")],
                 cols = c("Mean", "Median",
                          "SD", "IQR"),
                 rows = NULL, theme = ttheme_minimal())

topo1 <- grid.arrange(t1, tt1, heights = unit.c(unit(1, "null"), th))

# Topograhic distance, corresponding dates
t2 <- distances %>%  filter(rcarb_cor == "t") %>%
  ggplot(aes(x = year, y = topodist,
                        colour = type)) +
  geom_point(shape = 16, alpha = 0.01) +
  scale_colour_manual(values = colst) +
  labs(y = "Topographic distance (m)", x = "BCE", title = "") +
   scale_x_continuous(breaks = seq(-8000, -1000, by = 2000)) +
  theme_bw() +
  theme(legend.position = "none")

tt2 <- tableGrob(sum_cor[c("topodist_fn1", "topodist_fn2",
                           "topodist_fn3", "topodist_fn4")],
                 cols = c("Mean", "Median",
                          "SD", "IQR"),
                 rows = NULL, theme = ttheme_minimal())

topo2 <- grid.arrange(t2, tt2, heights = unit.c(unit(1, "null"), th))

# Repeat for vertical distance
v1 <- ggplot(distances, aes(x = year, y = vertdist,
                            colour = type)) +
  geom_point(shape = 16, alpha = 0.01) +
  scale_colour_manual(values = colsv) +
  labs(y = "Vertical distance (m)", x = "BCE", title = "") +
  scale_x_continuous(breaks = seq(-8000, -1000, by = 2000)) +
  ylim(-30, 90) +
  theme_bw() +
  theme(legend.position = "none")

vt1 <- tableGrob(sum_all[c("vertdist_fn1", "vertdist_fn2",
                           "vertdist_fn3", "vertdist_fn4")],
                 cols = c("Mean", "Median",
                          "SD", "IQR"),
                 rows = NULL, theme = ttheme_minimal())

vert1 <- grid.arrange(v1, vt1, heights = unit.c(unit(1, "null"), th))

# Vertical distance, corresponding dates
v2 <- distances %>%  filter(rcarb_cor == "t") %>%
  ggplot(aes(x = year, y = vertdist,
                        colour = type)) +
  geom_point(shape = 16, alpha = 0.01) +
  scale_colour_manual(values = colsv) +
  scale_x_continuous(breaks = seq(-8000, -1000, by = 2000)) +
  labs(y = "Vertical distance (m)", x = "BCE", title = "") +
  ylim(-30, 90) +
  theme_bw() +
  theme(legend.position = "none")

vt2 <- tableGrob(sum_cor[c("vertdist_fn1", "vertdist_fn2",
                           "vertdist_fn3", "vertdist_fn4")],
                 cols = c("Mean", "Median",
                          "SD", "IQR"),
                 rows = NULL, theme = ttheme_minimal())

vert2 <- grid.arrange(v2, vt2, heights = unit.c(unit(1, "null"), th))

grd1 <- grid.arrange(hor1, topo1, vert1, nrow = 1)
grd2 <- grid.arrange(hor2, topo2, vert2, nrow = 1)

grd <- plot_grid(grd1, grd2, labels = c("A", "B"), ncol = 1, scale = 0.9)

ggsave(file = here("analysis/figures/results.png"), grd,
       bg = "white", width = 250, height = 200, units = "mm")

# Measures for vertical distance to inform recommendations for shoreline dating.
# Using only corresponding dates, including dates to agricultural activity

# Dates older than 2500
distances2 <- distances %>%
  filter(year <= -2500 & rcarb_cor == "t") %>%
  # mutate(hordist = ifelse(str_detect(sitename, "Løvås"), 0, hordist),
  #      topodist = ifelse(str_detect(sitename, "Løvås"), 0, topodist),
  #      vertdist = ifelse(str_detect(sitename, "Løvås"), 0, vertdist)) %>%
  mutate(hordist = ifelse(hordist < 0, 0, hordist),
         topodist = ifelse(topodist < 0, 0, topodist),
         vertdist = ifelse(vertdist < 0, 0, vertdist))

# Dates older than 4000 (which also excludes agricultural sites)
distances3 <- distances %>%
  filter(year <= -4000 & rcarb_cor == "t") %>%
  mutate(hordist = ifelse(hordist < 0, 0, hordist),
         topodist = ifelse(topodist < 0, 0, topodist),
         vertdist = ifelse(vertdist < 0, 0, vertdist))

# Dates older than 7000, to accommodate for more uncertaint relationship
# between sites and shorelines for the earliest sites.
distances4 <- distances %>%
  filter(year <= -7000 & rcarb_cor == "t") %>%
  mutate(hordist = ifelse(hordist < 0, 0, hordist),
         topodist = ifelse(topodist < 0, 0, topodist),
         vertdist = ifelse(vertdist < 0, 0, vertdist))

# Fit exponential functions
expfit <- fitdistr(distances2$vertdist, "exponential")

expdat <- data.frame(
  x = distances2$vertdist,
  px = dexp(as.numeric(distances2$vertdist), rate = expfit$estimate))

# Fit exponential functions
expfit2 <- fitdistr(distances3$vertdist, "exponential")

expdat2 <- data.frame(
  x = distances3$vertdist,
  px = dexp(as.numeric(distances3$vertdist), rate = expfit2$estimate))

expfit3 <- fitdistr(distances4$vertdist, "exponential")

# # Visualise in base R
# hist(distances2$vertdist, freq =  FALSE)
# curve(dexp(x, rate = expfit$estimate, log = FALSE), from = 0, add = TRUE)

sum_cor2 <- distances2 %>%
  summarise_at(c("hordist", "topodist", "vertdist"),
               c(mean, median, sd, IQR)) %>%
  round()

sum_cor3 <- distances3 %>%
  summarise_at(c("hordist", "topodist", "vertdist"),
               c(mean, median, sd, IQR)) %>%
  round()

# Vertical distance,
# Vertical distance, dates older than 4000 BCE (i.e. Mesolithic)

histhor <- distances2 %>%
  ggplot() + geom_histogram(aes(hordist, y = ..density..), binwidth = 20,
                            fill = "#00ba38",
                            colour = "black",
                            alpha = 0.5) +
  theme_bw() +
  labs(title= "Horizontal distance (m) \nBinwidth = 20m",
       x = "", y = "Density")

# Table for Horizontal statistics after exclusion of younger sites
tabhor3 <- tableGrob(sum_cor2[c("hordist_fn1", "hordist_fn2",
                                "hordist_fn3", "hordist_fn4")],
                     cols = c("Mean", "Median",
                              "SD", "IQR"),
                     rows = NULL, theme = ttheme_minimal())

th <- sum(tabhor3$heights)

# Assemble Horizontal histogram and table
hor3 <- grid.arrange(histhor, tabhor3, heights = unit.c(unit(1, "null"), th))

histtopo <- distances2 %>%
  ggplot() + geom_histogram(aes(topodist, y = ..density..), binwidth = 20,
                            fill = "#fad510",
                            colour = "black",
                            alpha = 0.5) +
  theme_bw() +
  labs(title= "Topographic distance (m) \nBinwidth = 20m",
       x = "", y = "Density")

tabtopo3 <- tableGrob(sum_cor2[c("topodist_fn1", "topodist_fn2",
                                 "topodist_fn3", "topodist_fn4")],
                      cols = c("Mean", "Median",
                               "SD", "IQR"),
                      rows = NULL, theme = ttheme_minimal())

topo3 <- grid.arrange(histtopo, tabtopo3, heights = unit.c(unit(1, "null"), th))

label_text <- paste("lambda ==" , as.character(round(expfit$estimate, 3)))

histvert <- distances2 %>%
  ggplot() + geom_histogram(aes(vertdist, y = ..density..), binwidth = 1,
                            fill = "#046c9a",
                            colour = "black",
                            alpha = 0.5) +
  geom_line(data = expdat, aes(x = x, y = px), col = "red") +
  geom_text(aes(x = Inf, y = Inf, hjust = 1.2,
            vjust = 5,
            label = label_text), col = "red", parse = TRUE) +
  theme_bw() +
  labs(title= "Vertical distance (m) \nBinwidth = 1m", x = "", y = "Density")

tabvert3 <- tableGrob(sum_cor2[c("vertdist_fn1", "vertdist_fn2",
                                 "vertdist_fn3", "vertdist_fn4")],
                      cols = c("Mean", "Median",
                               "SD", "IQR"),
                      rows = NULL, theme = ttheme_minimal())

vert3 <- grid.arrange(histvert, tabvert3, heights = unit.c(unit(1, "null"), th))

grd3 <- grid.arrange(hor3, topo3, vert3, nrow = 1)

ggsave(file = here("analysis/figures/temp.png"), grd3,
       width = 250, height = 100, units = "mm")

# Repeat for Mesolithic

histhor2 <- distances3 %>%
  ggplot() + geom_histogram(aes(hordist, y = ..density..), binwidth = 20,
                            fill = "#00ba38",
                            colour = "black",
                            alpha = 0.5) +
  theme_bw() +
  labs(title= "Horizontal distance (m) \nBinwidth = 20m",
       x = "", y = "Density")

tabhor4 <- tableGrob(sum_cor3[c("hordist_fn1", "hordist_fn2",
                                "hordist_fn3", "hordist_fn4")],
                     cols = c("Mean", "Median",
                              "SD", "IQR"),
                     rows = NULL, theme = ttheme_minimal())

th <- sum(tabhor4$heights)

hor4 <- grid.arrange(histhor2, tabhor4, heights = unit.c(unit(1, "null"), th))

histtopo2 <- distances3 %>%
  ggplot() + geom_histogram(aes(topodist, y = ..density..), binwidth = 20,
                            fill = "#fad510",
                            colour = "black",
                            alpha = 0.5) +
  theme_bw() +
  labs(title= "Topographic distance (m) \nBinwidth = 20m",
       x = "", y = "Density")

tabtopo4 <- tableGrob(sum_cor3[c("topodist_fn1", "topodist_fn2",
                                 "topodist_fn3", "topodist_fn4")],
                      cols = c("Mean", "Median",
                               "SD", "IQR"),
                      rows = NULL, theme = ttheme_minimal())

topo4 <- grid.arrange(histtopo2, tabtopo4,
                      heights = unit.c(unit(1, "null"), th))

label_text2 <- paste("lambda ==" , as.character(round(expfit2$estimate, 3)))

histvert2 <- distances3 %>%
  ggplot() + geom_histogram(aes(vertdist, y = ..density..), binwidth = 1,
                            fill = "#046c9a",
                            colour = "black",
                            alpha = 0.5) +
  geom_line(data = expdat2, aes(x = x, y = px), col = "red") +
  geom_text(aes(x = Inf, y = Inf, hjust = 1.2,
                vjust = 5,
                label = label_text2), col = "red", parse = TRUE) +
  theme_bw() +
  labs(title= "Vertical distance (m) \nBinwidth = 1m", x = "", y = "Density")

tabvert4 <- tableGrob(sum_cor3[c("vertdist_fn1", "vertdist_fn2",
                                 "vertdist_fn3", "vertdist_fn4")],
                      cols = c("Mean", "Median",
                               "SD", "IQR"),
                      rows = NULL, theme = ttheme_minimal())

vert4 <- grid.arrange(histvert2, tabvert4,
                      heights = unit.c(unit(1, "null"), th))

grd4 <- grid.arrange(hor4, topo4, vert4, nrow = 1)

grd <- plot_grid(grd3, grd4, labels = c("A", "B"), ncol = 1, scale = 0.9)

# Save plot
ggsave(file = here("analysis/figures/results2.png"), grd,
       bg = "white", width = 250, height = 200, units = "mm")

# Save results and exponential functions
save(distances, expfit, expfit2, expfit3,
     file = here("analysis/data/derived_data/06data.RData"))
