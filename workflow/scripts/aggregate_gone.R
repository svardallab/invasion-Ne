library(tidyverse)
args <- commandArgs(trailingOnly = TRUE)
infiles <- args[1:(length(args) - 1)]
names(infiles) <- infiles 
output_file <- args[length(args)]

# Process files
infiles |>
  map(read_tsv, show_col_types = FALSE) |>
  bind_rows(.id = "File") |>
  rename(
    generation = Generation,
    ne = Ne_diploids
  ) |>
  # Extract seed and sample_size from the filename, if possible
  mutate(
    seed = as.character(str_extract(File, "(?<=s)\\d+")),
    sample_size = as.numeric(str_extract(File, "(?<=n)\\d+"))
  ) |>
  mutate(method = "GONE2") |>
  select(generation, ne, method, sample_size, seed) |>
  write_csv(output_file)