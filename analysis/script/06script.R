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
results <- list()

# loop over, load results and assign them to the list
for (i in 1:length(datfiles)){
  # Load results
  load(file.path(here("analysis/data/derived_data", datfiles[i])))
  results[[i]] <- output
}

# Collapse list of lists into a single data frame
distances <- results %>%  map(1) %>% map(1) %>% bind_rows()

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
                                     alpha = 0.01) +
  labs(y = "Horisontal distance (m)", x = "Cal BCE",
       title = "All dates") +
  theme_classic()

th1 <- tableGrob(sum_all[c("hordist_fn1", "hordist_fn2",
                         "hordist_fn3", "hordist_fn4")],
                 cols = c("Mean", "Median",
                          "SD", "IQR"),
                 rows = NULL, theme = ttheme_minimal())

# Line under column names
th1 <- gtable_add_grob(th1,
        grobs =  segmentsGrob(
          x0 = unit(0,"npc"),
          y0 = unit(0,"npc"),
          x1 = unit(1,"npc"),
          y1 = unit(0,"npc"),
          gp = gpar(lwd = 1)),
         t = 1, b =1, l = 1, r = ncol(th1))

# Find table height for plotting
th <- sum(th1$heights)
# Arrange plot and table
hor1 <- grid.arrange(h1, th1, heights = unit.c(unit(1, "null"), th))

h2 <- distances %>%  filter(rcarb_cor == "t") %>%
  ggplot() + geom_point(aes(x = year, y = hordist), alpha = 0.01) +
  labs(y = "Horisontal distance (m)", x = "Cal BCE",
       title = "Corresponding dates") +
  theme_classic()

th2 <- tableGrob(sum_cor[c("hordist_fn1", "hordist_fn2",
                           "hordist_fn3", "hordist_fn4")],
                 cols = c("Mean", "Median",
                          "SD", "IQR"),
                 rows = NULL, theme = ttheme_minimal())

th2 <- gtable_add_grob(th2,
                       grobs =  segmentsGrob(
                         x0 = unit(0,"npc"),
                         y0 = unit(0,"npc"),
                         x1 = unit(1,"npc"),
                         y1 = unit(0,"npc"),
                         gp = gpar(lwd = 1)),
                       t = 1, b =1, l = 1, r = ncol(th2))

hor2 <- grid.arrange(h2, th2, heights = unit.c(unit(1, "null"), th))

# Repeat for topographic distance
t1 <- ggplot(distances) + geom_point(aes(x = year, y = topodist),
                                     alpha = 0.01) +
  labs(y = "Topographic distance (m)", x = "Cal BCE") +
  theme_classic()

tt1 <- tableGrob(sum_all[c("topodist_fn1", "topodist_fn2",
                           "topodist_fn3", "topodist_fn4")],
                 cols = c("Mean", "Median",
                          "SD", "IQR"),
                 rows = NULL, theme = ttheme_minimal())

tt1 <- gtable_add_grob(tt1,
                       grobs =  segmentsGrob(
                         x0 = unit(0,"npc"),
                         y0 = unit(0,"npc"),
                         x1 = unit(1,"npc"),
                         y1 = unit(0,"npc"),
                         gp = gpar(lwd = 1)),
                       t = 1, b =1, l = 1, r = ncol(tt1))

topo1 <- grid.arrange(t1, tt1, heights = unit.c(unit(1, "null"), th))

# Topograhic distance, corresponding dates
t2 <- distances %>%  filter(rcarb_cor == "t") %>%
  ggplot() + geom_point(aes(x = year, y = topodist),
                                     alpha = 0.01) +
  labs(y = "Topographic distance (m)", x = "Cal BCE") +
  theme_classic()

tt2 <- tableGrob(sum_cor[c("topodist_fn1", "topodist_fn2",
                           "topodist_fn3", "topodist_fn4")],
                 cols = c("Mean", "Median",
                          "SD", "IQR"),
                 rows = NULL, theme = ttheme_minimal())

tt2 <- gtable_add_grob(tt2,
                       grobs =  segmentsGrob(
                         x0 = unit(0,"npc"),
                         y0 = unit(0,"npc"),
                         x1 = unit(1,"npc"),
                         y1 = unit(0,"npc"),
                         gp = gpar(lwd = 1)),
                       t = 1, b =1, l = 1, r = ncol(tt2))

topo2 <- grid.arrange(t2, tt2, heights = unit.c(unit(1, "null"), th))

# Repeat for vertical distance
v1 <- ggplot(distances) +
  geom_point(aes(x = year, y = vertdist), alpha = 0.01) +
  labs(y = "Vertical distance (m)", x = "Cal BCE") +
  ylim(-30, 90) +
  theme_classic()

vt1 <- tableGrob(sum_all[c("vertdist_fn1", "vertdist_fn2",
                           "vertdist_fn3", "vertdist_fn4")],
                 cols = c("Mean", "Median",
                          "SD", "IQR"),
                 rows = NULL, theme = ttheme_minimal())

vt1 <- gtable_add_grob(vt1,
                       grobs =  segmentsGrob(
                         x0 = unit(0,"npc"),
                         y0 = unit(0,"npc"),
                         x1 = unit(1,"npc"),
                         y1 = unit(0,"npc"),
                         gp = gpar(lwd = 1)),
                       t = 1, b =1, l = 1, r = ncol(vt1))

vert1 <- grid.arrange(v1, vt1, heights = unit.c(unit(1, "null"), th))

# Vertical distance, corresponding dates
v2 <- distances %>%  filter(rcarb_cor == "t") %>%
  ggplot() + geom_point(aes(x = year, y = vertdist), alpha = 0.01) +
  labs(y = "Vertical distance (m)", x = "Cal BCE") +
  ylim(-30, 90) +
  theme_classic()

vt2 <- tableGrob(sum_cor[c("vertdist_fn1", "vertdist_fn2",
                           "vertdist_fn3", "vertdist_fn4")],
                 cols = c("Mean", "Median",
                          "SD", "IQR"),
                 rows = NULL, theme = ttheme_minimal())

vt2 <- gtable_add_grob(vt2,
                       grobs =  segmentsGrob(
                         x0 = unit(0,"npc"),
                         y0 = unit(0,"npc"),
                         x1 = unit(1,"npc"),
                         y1 = unit(0,"npc"),
                         gp = gpar(lwd = 1)),
                       t = 1, b =1, l = 1, r = ncol(vt2))

vert2 <- grid.arrange(v2, vt2, heights = unit.c(unit(1, "null"), th))

grd <- grid.arrange(hor1, topo1, vert1, nrow = 1)
gb <- rectGrob(height = 0.999, width = 0.999, gp = gpar(fill = NA))
grd1 <- gTree(children = gList(grd, gb))

grd <- grid.arrange(hor2, topo2, vert2, nrow = 1)
gb <- rectGrob(height = 0.999, width = 0.999, gp = gpar(fill = NA))
grd2 <- gTree(children = gList(grd, gb))


grd <- arrangeGrob(grd1, grd2, nrow = 2)
ggsave(file = here("analysis/figures/results.png"), grd)


