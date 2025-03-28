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
    generation = GEN,
    ne = NE,
    lower_conf95 = `LWR-95%CI`,
    upper_conf95 = `UPR-95%CI`,
  ) |>
  # Extract seed and sample_size from the filename, if possible
  mutate(
    seed = as.character(str_extract(File, "(?<=s)\\d+")),
    sample_size = as.numeric(str_extract(File, "(?<=n)\\d+"))
  ) |>
  mutate(method = "IBDNE") |>
  select(
    generation, ne, method, sample_size,
    seed, lower_conf95, upper_conf95
  ) |>
  write_csv(output_file)