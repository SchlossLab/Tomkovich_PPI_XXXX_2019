#####################
# 
# Experiment: Treat mice with clindamycin and/or PPI (omeprazole) and challenge with Clostridium difficile
# Merge sequence file labels to corresponding metadata
# 
# 
# input:
#	data/mothur/ppi.opti_mcc.0.03.subsample.shared
#	data/raw/PPI_metadata.xlsx
#	data/raw/SarahPPiNov2018_plateMap.xlsx
#
#
# output:
#	data/process/ppi_metadata.txt
#		
#
#####################

library(tidyverse)
library(readxl)

shared_sample_names <- read.table('data/mothur/ppi.opti_mcc.0.03.subsample.shared', 
		sep = '\t', header = T, stringsAsFactors = F) %>% 
	select(Group) # select column with sample names of shared file which we will need to match our metadata to

# load in plate map label from metadata file. Rename Plate map label column and Collection D column.
# Get rid of the D before day number in the Collection D column.
#Convert Plate map labels to way they are labeled in shared file.
metadata <- read_xlsx('data/raw/PPI_metadata.xlsx') %>% 
  separate(Plate_number, c("plate_", "Plate_number")) %>% 
  rename("day" = "Collection D") %>% 
  mutate(day = as.numeric(gsub('D', '', day))) %>% 
  rename("Plate_label" = "Plate map label") %>% 
  # convert names of seq to match the format of the shared samples
  mutate(shared_names = gsub('^_{1,3}','', Plate_label),
         shared_names = gsub('_', '+', shared_names)) %>% 
  mutate(Group = case_when(Group == 'O+' ~ 'PPI', #make group names more readable
                           Group == 'C+' ~ 'Clindamycin',
                           Group == 'CO+' ~ 'Clindamycin + PPI')) 
  # select only the samples that are present in the shared data
matching <- inner_join(metadata, shared_sample_names, by = c('shared_names' = 'Group')) 

# What's missing from metadata shared_names compared to shared_sample_names
missing <- anti_join(shared_sample_names, metadata, by = c('Group' = 'shared_names'))
# COpos+D6+unk and Opos+M5+D9 on shared_sample_names did not match to fixed_shared_names
metadata <- metadata %>% mutate(shared_names = replace(shared_names, shared_names == "COpos+M10+D6+Unk", "COpos+M10+D6+unk")) %>%  
  mutate(shared_names = replace(shared_names, shared_names == "Opos+M5+D9+", "Opos+M5+D9"))

# Now check to see if all 131 shared_sample_names join to fixed_shared_names
final_metadata <- inner_join(metadata, shared_sample_names, by = c('shared_names' = 'Group'))

#Tidy C. difficile CFU data
ppi_cdiff_cfu <- metadata %>%
  select(Group, "Mouse ID", contains('difficile')) %>%
  gather(day, cfu, -Group, -"Mouse ID") %>%
  mutate(day = as.numeric(str_match(day, '\\d{1,2}'))) %>% 
  unique() %>% 
  mutate(Group = replace(Group, Group == "NA", NA)) %>% drop_na() %>%  #Get rid of rows that were placeholders for mock & control sequencing
  select("Mouse ID", day, cfu)
  
#Join CFU data to metadata file
output_df <- full_join(metadata,
                       ppi_cdiff_cfu, by = c('Mouse ID', 'day'))

#Change day to match experimental setup where C. difficile challenge day equals day 0. Implement by subtracting 7 from each day.
output_df <- output_df %>% 
  mutate(day = day - 7)

write.table(output_df, 'data/process/ppi_metadata.txt', 
  row.names = F, quote = F, sep = '\t')

