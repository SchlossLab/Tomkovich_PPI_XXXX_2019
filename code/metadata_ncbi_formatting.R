#Load packages
library(tidyverse)
library(lubridate)
library(janitor)

metadata <- read.table('data/process/ppi_metadata.txt', header = T, sep = '\t', stringsAsFactors = F) %>% 
  #Unite 2 tube label columns into a new column "description". This conveys group name, mouse ID, and collection date.
  unite(sample_title, Tube.label, Tube.label__1) %>% 
  select(shared_names, sample_title)

# Import ppi.files into data frame
ncbi_metadata <- read.table('data/raw/formatting_metadata_for_ncbi.txt', header = T, sep = '\t', stringsAsFactors = F)  

# Add metadata columns to mimarks file to finish filling out the information for all samples
ncbi_metadata_merge <- left_join(metadata, ncbi_metadata, by = c("shared_names" = "shared.names"))
write_tsv(ncbi_metadata_merge, 'data/raw/ncbi_metadata_merge.tsv', na = " ")

# Paste into appropriate columns of ppi_SRA_metadata file to upload to NCBI


