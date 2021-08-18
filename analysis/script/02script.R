library(tidyverse)
library(tmap)

# Load radiocarbon and site data from the first script
# load(here::here("analysis/data/derived_data/data.RData"))

# Read in displacement curves
arendal <- read.csv("bjornebu.csv")
tvedestrand <- read.csv("hanto.csv")
larvik <- read.csv("gunnarsrød.csv")
skoppum
