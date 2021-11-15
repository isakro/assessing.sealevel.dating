library(here)
library(dplyr)
library(purrr)

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
