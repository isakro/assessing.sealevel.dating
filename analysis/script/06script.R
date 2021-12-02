library(here)
library(dplyr)
library(purrr)
library(ggplot2)
library(grid)
library(gridExtra)
library(gtable)

# List all files except those starting with a number.
# That is, all data files associated with sites
datfiles <- grep("^[0-9]", list.files(here("analysis/data/derived_data")),
     invert = TRUE, value = TRUE)

# Create empty list to hold results
results_list <- list()

# loop over, load results and assign them to the list
for(i in 1:length(datfiles)){
  # Load results
  load(file.path(here("analysis/data/derived_data", datfiles[i])))
  results_list[[i]] <- output[!(names(output) %in%
                             c("datedat", "sitel", "sitecurve"))]
}

# results_list is a list of lists of variable length...
# Anyway, these need to be unpacked. Might be a smoother purrr solution to
# all of this.
results <- list()
for(i in 1:length(results_list)){
  dat <- results_list[[i]]
  for(j in 1:length(dat)){
    results[[i]] <- dat[[j]]["results"]
  }
}

# Collapse list of lists into a single data frame
distances <- results %>% bind_rows()
distances <-  distances$results

sum_all <- distances %>%
  summarise_at(c("hordist", "topodist", "vertdist"),
               c(mean, median, sd, IQR)) %>%
  round()

sum_cor <- distances %>% filter(rcarb_cor == "t") %>%
  summarise_at(c("hordist", "topodist", "vertdist"),
               c(mean, median, sd, IQR)) %>%
   round()

# Horisontal distance, all dates to the Stone Age
h1 <- ggplot(distances) + geom_point(aes(x = year, y = hordist),
                                     shape = 21,
                                     alpha = 0.01,
                                     fill = scales::viridis_pal()(3)[1],
                                     colour = scales::viridis_pal()(3)[1]) +
  labs(y = "Horisontal distance (m)", x = "cal BCE",
       title = "All dates") +
  theme_bw()

th1 <- tableGrob(sum_all[c("hordist_fn1", "hordist_fn2",
                         "hordist_fn3", "hordist_fn4")],
                 cols = c("Mean", "Median",
                          "SD", "IQR"),
                 rows = NULL, theme = ttheme_minimal())

# # Line under column names
# th1 <- gtable_add_grob(th1,
#         grobs =  segmentsGrob(
#           x0 = unit(0,"npc"),
#           y0 = unit(0,"npc"),
#           x1 = unit(1,"npc"),
#           y1 = unit(0,"npc"),
#           gp = gpar(lwd = 1)),
#          t = 1, b =1, l = 1, r = ncol(th1))

# Find table height for plotting
th <- sum(th1$heights)
# Arrange plot and table
hor1 <- grid.arrange(h1, th1, heights = unit.c(unit(1, "null"), th))

h2 <- distances %>%  filter(rcarb_cor == "t") %>%
  ggplot() +
  geom_point(aes(x = year, y = hordist),
                        shape = 21,
                        alpha = 0.01,
                        fill = scales::viridis_pal()(3)[1],
                        colour = scales::viridis_pal()(3)[1]) +
  labs(y = "Horisontal distance (m)", x = "cal BCE",
       title = "Corresponding dates") +
  theme_bw()

th2 <- tableGrob(sum_cor[c("hordist_fn1", "hordist_fn2",
                           "hordist_fn3", "hordist_fn4")],
                 cols = c("Mean", "Median",
                          "SD", "IQR"),
                 rows = NULL, theme = ttheme_minimal())

# th2 <- gtable_add_grob(th2,
#                        grobs =  segmentsGrob(
#                          x0 = unit(0,"npc"),
#                          y0 = unit(0,"npc"),
#                          x1 = unit(1,"npc"),
#                          y1 = unit(0,"npc"),
#                          gp = gpar(lwd = 1)),
#                        t = 1, b =1, l = 1, r = ncol(th2))

hor2 <- grid.arrange(h2, th2, heights = unit.c(unit(1, "null"), th))

# Repeat for topographic distance
t1 <- ggplot(distances) +
  geom_point(aes(x = year, y = topodist),
             shape = 21,
             alpha = 0.01,
             fill = scales::viridis_pal()(3)[2],
             colour = scales::viridis_pal()(3)[2]) +
  labs(y = "Topographic distance (m)", x = "cal BCE") +
  theme_bw()

tt1 <- tableGrob(sum_all[c("topodist_fn1", "topodist_fn2",
                           "topodist_fn3", "topodist_fn4")],
                 cols = c("Mean", "Median",
                          "SD", "IQR"),
                 rows = NULL, theme = ttheme_minimal())

# tt1 <- gtable_add_grob(tt1,
#                        grobs =  segmentsGrob(
#                          x0 = unit(0,"npc"),
#                          y0 = unit(0,"npc"),
#                          x1 = unit(1,"npc"),
#                          y1 = unit(0,"npc"),
#                          gp = gpar(lwd = 1)),
#                        t = 1, b =1, l = 1, r = ncol(tt1))

topo1 <- grid.arrange(t1, tt1, heights = unit.c(unit(1, "null"), th))

# Topograhic distance, corresponding dates
t2 <- distances %>%  filter(rcarb_cor == "t") %>%
  ggplot() +
  geom_point(aes(x = year, y = topodist),
             shape = 21,
             alpha = 0.01,
             fill = scales::viridis_pal()(3)[2],
             colour = scales::viridis_pal()(3)[2]) +
  labs(y = "Topographic distance (m)", x = "cal BCE") +
  theme_bw()

tt2 <- tableGrob(sum_cor[c("topodist_fn1", "topodist_fn2",
                           "topodist_fn3", "topodist_fn4")],
                 cols = c("Mean", "Median",
                          "SD", "IQR"),
                 rows = NULL, theme = ttheme_minimal())

# tt2 <- gtable_add_grob(tt2,
#                        grobs =  segmentsGrob(
#                          x0 = unit(0,"npc"),
#                          y0 = unit(0,"npc"),
#                          x1 = unit(1,"npc"),
#                          y1 = unit(0,"npc"),
#                          gp = gpar(lwd = 1)),
#                        t = 1, b =1, l = 1, r = ncol(tt2))

topo2 <- grid.arrange(t2, tt2, heights = unit.c(unit(1, "null"), th))

# Repeat for vertical distance
v1 <- ggplot(distances) +
  geom_point(aes(x = year, y = vertdist), shape = 21,
             alpha = 0.01,
             fill = scales::viridis_pal()(3)[3],
             colour = scales::viridis_pal()(3)[3]) +
  labs(y = "Vertical distance (m)", x = "cal BCE") +
  ylim(-30, 90) +
  theme_bw()

vt1 <- tableGrob(sum_all[c("vertdist_fn1", "vertdist_fn2",
                           "vertdist_fn3", "vertdist_fn4")],
                 cols = c("Mean", "Median",
                          "SD", "IQR"),
                 rows = NULL, theme = ttheme_minimal())

# vt1 <- gtable_add_grob(vt1,
#                        grobs =  segmentsGrob(
#                          x0 = unit(0,"npc"),
#                          y0 = unit(0,"npc"),
#                          x1 = unit(1,"npc"),
#                          y1 = unit(0,"npc"),
#                          gp = gpar(lwd = 1)),
#                        t = 1, b =1, l = 1, r = ncol(vt1))

vert1 <- grid.arrange(v1, vt1, heights = unit.c(unit(1, "null"), th))

# Vertical distance, corresponding dates
v2 <- distances %>%  filter(rcarb_cor == "t") %>%
  ggplot() +
  geom_point(aes(x = year, y = vertdist),
             shape = 21,
             alpha = 0.01,
             fill = scales::viridis_pal()(3)[3],
             colour = scales::viridis_pal()(3)[3]) +
  labs(y = "Vertical distance (m)", x = "cal BCE") +
  ylim(-30, 90) +
  theme_bw()

vt2 <- tableGrob(sum_cor[c("vertdist_fn1", "vertdist_fn2",
                           "vertdist_fn3", "vertdist_fn4")],
                 cols = c("Mean", "Median",
                          "SD", "IQR"),
                 rows = NULL, theme = ttheme_minimal())

# vt2 <- gtable_add_grob(vt2,
#                        grobs =  segmentsGrob(
#                          x0 = unit(0,"npc"),
#                          y0 = unit(0,"npc"),
#                          x1 = unit(1,"npc"),
#                          y1 = unit(0,"npc"),
#                          gp = gpar(lwd = 1)),
#                        t = 1, b =1, l = 1, r = ncol(vt2))

vert2 <- grid.arrange(v2, vt2, heights = unit.c(unit(1, "null"), th))

grd1 <- grid.arrange(hor1, topo1, vert1, nrow = 1)
# gb <- rectGrob(height = 0.999, width = 0.999, gp = gpar(fill = NA))
# grd1 <- gTree(children = gList(grd, gb))

grd2 <- grid.arrange(hor2, topo2, vert2, nrow = 1)
# gb <- rectGrob(height = 0.999, width = 0.999, gp = gpar(fill = NA))
# grd2 <- gTree(children = gList(grd, gb))

grd <- plot_grid(grd1, grd2, ncol = 1, scale = 0.9)

ggsave(file = here("analysis/figures/results.png"), grd,
       bg = "white")

# ggsave(file = here("analysis/figures/results_cor.png"), grd2,
#        width = 275, height = 120, units = "mm")
#
# grd <- arrangeGrob(grd1, grd2, nrow = 2)
# ggsave(file = here("analysis/figures/results.png"), grd)

# Exclude data younger than 2500 BCE and only use corresponding 14C-dates.

# Find statistics again
sum_cor2 <- distances %>% filter(year <= -2500 & rcarb_cor == "t") %>%
  summarise_at(c("hordist", "topodist", "vertdist"),
               c(mean, median, sd, IQR)) %>%
  round()

# Histogram for horisontal distance
histhor <- distances %>%  filter(year <= -2500 & rcarb_cor == "t") %>%
  ggplot() + geom_histogram(aes(hordist), binwidth = 20,
                            fill = scales::viridis_pal()(3)[1],
                            colour = "black",
                            alpha = 0.5) +
  theme_bw() +
  labs(title= "Horisontal distance (m) \nBinwidth = 20m", x = "", y = "Frequency")

# Table for horisontal statistics after exclusion of younger sites
tabhor3 <- tableGrob(sum_cor2[c("hordist_fn1", "hordist_fn2",
                           "hordist_fn3", "hordist_fn4")],
                 cols = c("Mean", "Median",
                          "SD", "IQR"),
                 rows = NULL, theme = ttheme_minimal())

# Assemble horisontal histogram and table
hor3 <- grid.arrange(histhor, tabhor3, heights = unit.c(unit(1, "null"), th))

histtopo <- distances %>%  filter(year <= -2500 & rcarb_cor == "t") %>%
  ggplot() + geom_histogram(aes(topodist), binwidth = 20,
                            fill = scales::viridis_pal()(3)[2],
                            colour = "black",
                            alpha = 0.5) +
  theme_bw() +
  labs(title= "Topographic distance (m) \nBinwidth = 20m", x = "", y = "Frequency")

tabtopo3 <- tableGrob(sum_cor2[c("topodist_fn1", "topodist_fn2",
                            "topodist_fn3", "topodist_fn4")],
                  cols = c("Mean", "Median",
                           "SD", "IQR"),
                  rows = NULL, theme = ttheme_minimal())

topo3 <- grid.arrange(histtopo, tabtopo3, heights = unit.c(unit(1, "null"), th))

histvert <- distances %>%  filter(year <= -2500 & rcarb_cor == "t") %>%
  ggplot() + geom_histogram(aes(vertdist), binwidth = 1,
                            fill = scales::viridis_pal()(3)[3],
                            colour = "black",
                            alpha = 0.5) +
  theme_bw() +
  labs(title= "Vertical distance (m) \nBinwidth = 1m", x = "", y = "Frequency")

tabvert3 <- tableGrob(sum_cor2[c("vertdist_fn1", "vertdist_fn2",
                            "vertdist_fn3", "vertdist_fn4")],
                  cols = c("Mean", "Median",
                           "SD", "IQR"),
                  rows = NULL, theme = ttheme_minimal())

vert3 <- grid.arrange(histvert, tabvert3, heights = unit.c(unit(1, "null"), th))

grd3 <- grid.arrange(hor3, topo3, vert3, nrow = 1)

ggsave(file = here("analysis/figures/results2.png"), grd3,
       width = 250, height = 100, units = "mm")

save(distances,
     file = here::here("analysis/data/derived_data/06data.RData"))

