#Load packages
library(tidyverse)
library(lubridate)
library(janitor)

# Import metadata into data frame
metadata <- read.table('data/process/ppi_metadata.txt', header = T, sep = '\t', stringsAsFactors = F) %>% 
  #Unite 2 tube label columns into a new column "description". This conveys group name, mouse ID, and collection date.
  unite(sample_title, Tube.label, Tube.label__1) %>% 
  mutate(description = Group) %>% 
  mutate(collection_date = excel_numeric_to_date(Date.Collected)) %>% #Convert excel numeric values to date values with excel_numeric_to_date function from the janitor package
  mutate(experiment_day = day, #Reformat names of additional columns since they'll be included in mimarks file so that each sample will pass NCBI's uniqueness check.
         sex = M.F,
         exp_cage_number = Exp..cage,
         mouse_id = Mouse.ID) %>% 
  select(sample_title, description, collection_date, sex, experiment_day, mouse_id, exp_cage_number, shared_names) #Select columns that are missing from mimarks tsv file


# Import MIMARKS ppi.tsv file into data frame
mimarks <- read_tsv('data/raw/ppi.tsv')
# Add metadata columns to mimarks file to finish filling out the information for all samples
mimarks <- left_join(mimarks, metadata, by = c("#This is a tab-delimited file. Additional Documentation can be found at http://www.mothur.org/wiki/MIMarks_Data_Packages." = "shared_names"))
write_tsv(mimarks, 'data/raw/ppi.tsv', na = " ")

# To prepare for NCBI submission. Open up file in excel and move additional column names to their appropriate location on row 10
          
                     