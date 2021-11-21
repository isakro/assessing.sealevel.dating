library(here)
library(dplyr)
library(purrr)
library(ggplot2)
library(patchwork)

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

# Horisontal distance, all dates to the Stone Age
h1 <- ggplot(distances) + geom_point(aes(x = year, y = hordist),
                                     alpha = 0.01) +
  labs(y = "Horisontal distance", x = "",
       title = "All dates") +
  theme_classic()

h2 <- distances %>%  filter(rcarb_cor == "t") %>%
  ggplot() + geom_point(aes(x = year, y = hordist), alpha = 0.01) +
  labs(y = "Horisontal distance", x = "",
       title = "Corresponding dates") +
  theme_classic()

t1 <- ggplot(distances) + geom_point(aes(x = year, y = topodist),
                                     alpha = 0.01) +
  labs(y = "Topographic distance", x = "") +
  theme_classic()

t2 <- distances %>%  filter(rcarb_cor == "t") %>%
  ggplot() + geom_point(aes(x = year, y = topodist),
                                     alpha = 0.01) +
  labs(y = "Topographic distance", x = "") +
  theme_classic()

v1 <- ggplot(distances) +
  geom_point(aes(x = year, y = vertdist), alpha = 0.01) +
  labs(y = "Vertical distance", x = "Cal BCE") +
  ylim(-30, 90) +
  theme_classic()

v2 <- distances %>%  filter(rcarb_cor == "t") %>%
  ggplot() + geom_point(aes(x = year, y = vertdist), alpha = 0.01) +
  labs(y = "Vertical distance", x = "Cal BCE") +
  ylim(-30, 90) +
  theme_classic()

(h1 + h2) /
(t1 + t2) /
(v1 + v2)

