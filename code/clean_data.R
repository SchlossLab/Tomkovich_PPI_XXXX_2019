#####################
# 
# clean up labeling issues with data - Treat mice with abx/PPI and challenge with C difficile
# 
# input:
#	data/mothur/stability.opti_mcc.0.03.subsample.shared
#	PPI_Exp_Sample_List_2018.xlsx
#	data/raw/Extraction Template_Plate #1#2#3.xlsx
#	data/raw/Plate_4_Sarah_Omep2Box1.xlsx
#	data/raw/SarahPPiNov2018_plateMap.xlsx
#
#
# output:
#	clean_data.txt
#		
#
#####################

library(tidyverse)
library(readxl)

# get a list of all the samples names from sequencing
shared_sample_names <- read.table('data/mothur/stability.opti_mcc.0.03.subsample.shared', 
		sep = '\t', header = T, stringsAsFactors = F) %>% 
	select(Group) # select column with sample names of shared file which we will need to match our metadata to

# load the metafile which has all of the original information 
# with a column of what each sample was listed on the plate as
sample_id_df <- read_xlsx('data/raw/PPI_Exp_Sample_List_2018.xlsx', sheet = 'Sheet1') %>% 
	rename(ExpNumber = 'Exp. #', CollectionDay = 'Collection D', CageNumber = 'Cage #', 
		M_F = 'M/F', OmpDailyDose = 'Omp. Daily Dose', MouseID = 'Mouse ID', TubeDayLabel = 'Tube label__1') %>% 
	filter(ExpNumber == 2) %>% 
	mutate(day = as.numeric(gsub('D', '', CollectionDay)),
 		treatment = gsub('\\+', '', Group),
 		sample_id = paste(treatment, MouseID, day, sep = '_'),
 		OmpDailyDose = as.numeric(OmpDailyDose),
 		mouse = as.numeric(MouseID)) %>% 
	select(sample_id, day, treatment, mouse, plate_label = `Plate map label`, 
		sex = M_F, cage = CageNumber, ear_mark = `Ear Mark`)

# data was recorded in wide format, so tidy cdiff infection data so that one column per type of data
exp_2_data <- read_xlsx('data/raw/08_20_18_Omep._#2_exp.xlsx', 
		range = 'A1:AE15', sheet = 'cfu_final') %>%
	rename(treatment = Group, mouse = 'Mouse ID') %>%
	mutate(treatment = gsub('\\+', '', treatment))
exp_2_wt <- exp_2_data %>%
	select(treatment, mouse, one_of(paste0('D', -1:16))) %>%
	gather(day, weight, -treatment, -mouse)	%>%
	mutate(day = as.numeric(gsub('D', '', day)))
exp_2_cdiff <- exp_2_data %>%
	select(treatment, mouse, ear_mark = `Ear Mark`, contains('difficile')) %>%
	gather(day, cfu, -treatment, -mouse, -ear_mark) %>%
	mutate(day = as.numeric(str_match(day, '\\d{1,2}')))

# join metadata to the exp data
output_df <- full_join(sample_id_df,
	exp_2_cdiff, by = c('treatment', 'mouse', 'day', 'ear_mark'))

# convert shared group names to plate labels
shared_names <- shared_sample_names %>%  
  select(Group) %>%  
  mutate(day = str_extract(Group, '(?<=D)\\d{1,2}') %>% as.numeric, 
    mouse = str_extract(Group, '(?<=M)\\d{1,2}') %>% as.numeric, 
    treatment = str_extract(Group, '[CO]+(?=POS)'), 
    unknown = ifelse(grepl('UNK', Group), T, F)) 
 
plate_names <- output_df %>%  
  select(plate_label) %>%  
  mutate(day = str_extract(plate_label, '(?<=D)\\d{1,2}') %>% as.numeric, 
    mouse = str_extract(plate_label, '(?<=M)\\d{1,2}') %>% as.numeric, 
    treatment = str_extract(plate_label, '[CO]+(?=pos)'), 
    unknown = ifelse(grepl('Unk|unk', plate_label), T, F)) 

output_df <- full_join(plate_names, shared_names, 
		by = c('treatment', 'mouse', 'day', 'unknown')) %>% 
	select(plate_label, Group) %>% 
	full_join(output_df, by = 'plate_label')

write.table(exp_2_wt, 'data/process/exp_2_weight.txt', 
	row.names = F, quote = F, sep = '\t')
write.table(output_df, 'data/process/metadata.txt', 
	row.names = F, quote = F, sep = '\t')