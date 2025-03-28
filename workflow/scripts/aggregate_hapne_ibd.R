library(tidyverse)
args <- commandArgs(trailingOnly = TRUE)
infiles <- args[1:(length(args) - 1)]
names(infiles) <- infiles 
output_file <- args[length(args)]

# Process files
infiles |>
  map(read_csv, show_col_types = FALSE) |>
  bind_rows(.id = "File") |>
  rename(
    generation = TIME,
    ne = `Q0.5`,
    lower_conf95 = `Q0.025`,
    upper_conf95 = `Q0.975`,
  ) |>
  # Extract seed and sample_size from the filename, if possible
  mutate(
    seed = as.character(str_extract(File, "(?<=s)\\d+")),
    sample_size = as.numeric(str_extract(File, "(?<=n)\\d+")),
    ne = ne / 2,
    lower_conf95 = lower_conf95 / 2,
    upper_conf95 = upper_conf95 / 2,
  ) |>
  mutate(method = "HapNe-IBD") |>
  select(
    generation, ne, method, sample_size,
    seed, lower_conf95, upper_conf95
  ) |>
  write_csv(output_file)