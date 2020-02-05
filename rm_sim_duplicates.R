## Identify and remove duplicates from the simulations run on the cluster.

library("tidyverse")  # ggplot2, dplyr, readr, purrr, etc

# Load the file with simulations to remove duplicates from
sims <- read_csv("sims_etc/power_simulation_results_append_all-relevant_200205.csv")
head(sims)

# We consider duplicates if same rseed, Nsubj and sim_type appears more than 14 times
# What it looks like mostly:
sims %>% group_by(rseed, Nsubj, sim_type) %>%
  summarise(N=n())
# How many duplicates?
sims %>% group_by(rseed, Nsubj, sim_type) %>%
  summarise(N=n()) %>%
  group_by(N) %>%
  summarise(count = n())


# Function to remove duplicates
rm_dupl <- function (df) {
  # Select every 14th row and only relevant columns
  d14 <- df[1 + 14 * 0 : (nrow(sims)/14 - 1), ] %>%
    select(estimate, model, BF_BIC, rseed, sim_type, Nsubj)

  # row ids of unique rows
  unique_rows <- as.numeric(rownames(unique(d14)))

  # Function that takes an integer and makes a vector of that number and the
  # n following naturals
  seq_n <- function(z, n = 3) {
    z : (z + n)
  }
  # return the original dataframe without duplicates
  df[unlist(lapply(unique_rows, seq_n, n = 13)), ]
}

sims_unique <- rm_dupl(sims)

# Check
# Any duplicates left?
sims_unique %>% group_by(rseed, Nsubj, sim_type) %>%
  summarise(N=n()) %>%
  group_by(N) %>%
  summarise(count = n())

# Save the data frame without duplicates
write_csv(sims_unique, "sims_etc/power_simulation_results_append_all-relevant_200205_noduplic.csv")
