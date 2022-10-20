library(tidyverse)
library(MASS)
library(poweRlaw)

# Load simulation results
load(here::here("analysis/data/derived_data/06data.RData"))

# Vertical distances, including negative values
vertdistneg <- distances %>%
  filter(year <= -2500 & rcarb_cor == "t")

# Vertical distances, forcing negative values to zero
vertdist <- distances %>%
  filter(year <= -2500 & rcarb_cor == "t") %>%
  mutate(vertdist = ifelse(vertdist < 0, 0, vertdist))

# Vertical distances, excluding negative values
vertdistna <- distances %>%
  filter(year <= -2500 & rcarb_cor == "t") %>%
  mutate(vertdist = ifelse(vertdist < 0, NA, vertdist)) %>%
  filter(!is.na(vertdist))

# Fit normal distribution to all values
normfit <- MASS::fitdistr(vertdistneg$vertdist, "normal")
normdat <- data.frame(
  x = vertdistneg$vertdist,
  px = dnorm(vertdistneg$vertdist,
             mean =  normfit$estimate[1],
             sd =  normfit$estimate[1]))

# Fit exponential functions to only positive values
expfit <- fitdistr(vertdist$vertdist, "exponential")
expdat <- data.frame(
  x = vertdist$vertdist,
  px = dexp(as.numeric(vertdist$vertdist), rate = expfit$estimate))

# Fit exponential functions to only positive values
expfitna <- fitdistr((vertdistna$vertdist+0.1), "exponential")
expdatna <- data.frame(
  x = vertdistna$vertdist,
  px = dexp(as.numeric(vertdistna$vertdist), rate = expfitna$estimate))

# Fit a continuous power law
pl <- conpl$new((vertdist$vertdist + 0.01))
pl$setXmin(min((vertdist$vertdist + 0.01)))
pl$setPars(estimate_pars(pl))
pl_xy <- lines(pl)
pl_xy$px <- knots(pl_xy$y)

dst <- density((vertdist$vertdist + 0.01))
plot(dst)
plot(pl_xy, type = "line", add = TRUE)

vertdistneg %>%
ggplot() + geom_histogram(aes(vertdist, y = ..density..), binwidth = 1,
                          fill = NA,
                          colour = "darkgrey",
                          alpha = 0.5) +
  geom_line(data = expdat, aes(x = x, y = px), col = "red") +
  geom_line(data = expdatna, aes(x = x, y = px), col = "gold1") +
  geom_line(data = normdat, aes(x = x, y = px), col = "black") +
  # geom_line(data = long_pl_xy, aes(x = x, y = pxnorm),
  #           col = "darkgreen") +
  # geom_text(aes(x = Inf, y = Inf, hjust = 1.2,
  #               vjust = 5,
  #               label = label_text), col = "red", parse = TRUE) +
  theme_bw() +
  labs(title= "Vertical distance (m) \nBinwidth = 1m", x = "", y = "Density")
