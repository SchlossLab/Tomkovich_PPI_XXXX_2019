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
#	data/ppi_merge_data.txt
#		
#
#####################

library(tidyverse)
library(readxl)

shared_sample_names <- read.table('data/mothur/ppi.opti_mcc.0.03.subsample.shared', 
		sep = '\t', header = T, stringsAsFactors = F) %>% 
	select(Group) # select column with sample names of shared file which we will need to match our metadata to

# load in all the names and positions of the sequence sample names
seq_plate_map <- 'data/raw/SarahPPiNov2018_plateMap.xlsx'
seq_sample_names <- bind_rows(
		mutate(read_xlsx(seq_plate_map, sheet = 'Final Plate Map', range = 'C16:O24'), plate = 2),
		mutate(read_xlsx(seq_plate_map, sheet = 'Final Plate Map', range = 'C27:O35'), plate = 3),
		mutate(read_xlsx(seq_plate_map, sheet = 'Final Plate Map', range = 'C39:O47'), plate = 4)) %>%
	rename(row = X__1) %>%
	gather(column, seq_sample, -row, -plate) %>% 
  unite(Plate_location, row, column, sep = "") %>% 
  rename(Plate_number = plate) %>% 
  mutate(Plate_number = as.character(Plate_number))

# load in plate map label from metadata file
metadata <- read_xlsx('data/raw/PPI_metadata.xlsx')
meta_sample_names <- metadata %>% separate(Plate_number, c("plate_", "Plate_number")) %>% rename("Plate_label" = "Plate map label")
  
# check to see if all meta and seq samples match
fixed_shared_names <- full_join(seq_sample_names, meta_sample_names, by = c('Plate_number', 'Plate_location')) %>%
  select(Plate_label) %>%
  # convert names of seq to match the format of the shared samples
  mutate(shared_names = gsub('^_{1,3}','', Plate_label),
         shared_names = gsub('_', '+', shared_names)) 
  # select only the samples that are present in the shared data
matching <- inner_join(fixed_shared_names, shared_sample_names, by = c('shared_names' = 'Group')) %>%
  mutate(new_shared_names = paste0(ext_treatment, 'POS+M', ext_mouse, '+D', ext_day)) %>%
  select(shared_names, new_shared_names)

# What's missing from fixed_shared_names compared to shared_sample_names
missing <- anti_join(shared_sample_names, fixed_shared_names, by = c('Group' = 'shared_names'))
# COpos+D6+unk and Opos+M5+D9 on shared_sample_names did not match to fixed_shared_names
fixed_shared_names <- fixed_shared_names %>% mutate(shared_names = replace(shared_names, shared_names == "COpos+M10+D6+Unk", "COpos+M10+D6+unk")) %>%  
  mutate(shared_names = replace(shared_names, shared_names == "Opos+M5+D9+", "Opos+M5+D9"))

# Now check to see if all 131 shared_sample_names join to fixed_shared_names
final_matching <- inner_join(fixed_shared_names, shared_sample_names, by = c('shared_names' = 'Group'))

	

# join shared sample names with corrected sample names
shared_sample_names <- shared_sample_names %>%
	full_join(fixed_shared_names, by = c('Group' = 'shared_names')) %>%
	mutate(shared_samples = ifelse(is.na(new_shared_names), Group, new_shared_names))

sample_id_df <- shared_sample_names %>% 
	mutate(day = as.numeric(str_match(shared_samples, '(?<=D)\\d{1,2}')),
		mouse = as.numeric(str_match(shared_samples, '(?<=M)\\d{1,2}')),
		treatment = as.character(str_match(shared_samples, '[:alpha:]{1,2}(?=POS)')),
		sample_id = paste(treatment, mouse, day, sep = '_')) %>%
	filter(treatment %in% c('C', 'CO', 'O'))

#######
## below working to identify any other mis-labeled samples
#######
#
## load the metafile originally recording the data to check out sequence labels against
true_data <- read_xlsx('data/raw/PPI_Exp_Sample_List_2018.xlsx', sheet = 'Sheet1') %>% 
	rename(ExpNumber = 'Exp. #', CollectionDay = 'Collection D', CageNumber = 'Cage #', 
		M_F = 'M/F', OmpDailyDose = 'Omp. Daily Dose', MouseID = 'Mouse ID', TubeDayLabel = 'Tube label__1') %>% 
	# Collection Day is shifted for Exp 2 from Day 4 to 7, so use tube label for day
	mutate(day = as.numeric(gsub('D', '', CollectionDay)),
 		treatment = gsub('\\+', '', Group),
 		sample_id = paste(treatment, MouseID, day, sep = '_'))
#	# if going by collection day, two entries for C_M14_D7, looks to be a shifted copy/paste entry error beginning day 4 of exp 2
#	# in column Collection Day, D6 begins a line early (C14, whereas all others outside of D4-7 starts on O1)
#
#mutate(true_data, mice_id = paste(treatment, MouseID, sep = '_')) %>%
#		pull(mice_id) %>% table %>% data.frame %>% filter(Freq < 10)
#	# all mice have 10 samples except O_1, O_2, O_3 (missing day 9, 4, 10)
#mutate(true_data, mice_id = paste(treatment, MouseID, sep = '_')) %>% filter(mice_id %in% c('O_1', 'O_2', 'O_3')) %>% select(mice_id, day)
## what days were the samples from
## samples were collected on day 0 2 4 6 7 8 9 10 11 12 14 16 17 (day 6 clinda IP and day 7 cdi)
#table(sample_id_df$day)
## 0  2  4  6  7  8  9 10 12 14 16 
##31 28 29 28 29 28 28 26 28 13 14 
#table(true_data$day)
## 0  2  4  6  7  8  9 10 12 14  
##32 32 31 32 32 32 31 31 32 32  
## No sample collected for:
## O_M1_D9, O_M2_D4, O_M3_D10, 
## last day depends on exp - Exp 1 (C_1-6, O_7-12, CO_13-18)_D14, Exp 2 (O_1-5, CO_6-10, C_11-14)_D16
experiment_number <- true_data %>%
	mutate(experiment = OmpDailyDose/20) %>%
	select(treatment, mouse = MouseID, experiment) %>%
	unique
#
## what samples are duplicated
#table(sample_id_df$sample_id) %>% data.frame %>% filter(Freq > 1)

#	# 			Freq
#	#	C_14_7    2 <- two positions on extraction 4 plate marked same
########	looking at the picture of the tubes 
#				plate 4 row 3 column 5(E) C+ M13 D7
#				plate 4 row 6 column 9(I) C+ M14 D7

#	#	CO_8_0    2 <- two positions on extraction 4 plate marked same
########	looking at the picture of the tubes 
#				plate 4 row 5 column 5(E) CO+ M8 D0
#				plate 4 row 6 column 9(I) CO+ M9 D0

#
## which mice were in each treatment
#mutate(sample_id_df, mice_id = paste(treatment, mouse, sep = '_')) %>%
#		pull(mice_id) %>% table %>% data.frame %>% filter(Freq > 2)
#	# clinda (C) - 1 2 3 4 5 6 11 12 13 14 (10 mice total)
#	# clinda PPI (CO) - 6 7 8 9 10 13 14 15 16 17 18 (11 mice total)
#	# PPI (O) - 1 2 3 4 5 7 8 9 10 11 12 (11 mice total)
#mutate(sample_id_df, mice_id = paste(treatment, mouse, sep = '_')) %>%
#		pull(mice_id) %>% table %>% data.frame %>% filter(Freq < 5)
##	CO_5 1 <- either C or O OR 6-10,13-18
##		on plate 4 (D9) could be O_5 but already same day exists
##		but other position the text is red
##	O_13 1 <- most likely mislabeled O
##		on plate 3 (C7) in row with day 10 samples, missing O_3 and C13, 
##		but no O_3 collected so this sample should be C_13
mislabeled_samples <- mutate(sample_id_df, mice_id = paste(treatment, mouse, sep = '_')) %>%
		filter(mice_id %in% c('CO_5', 'O_13')) %>%
		pull(Group)

sample_id_df <- sample_id_df %>%
	filter(!Group %in% mislabeled_samples) %>%
	select(-new_shared_names) %>%
	left_join(experiment_number, by = c('treatment', 'mouse'))

### tidy cdiff infection data
exp_2_data <- read_xlsx('data/raw/08_20_18_Omep._#2_exp.xlsx', 
		range = 'A1:AE15', sheet = 'cfu_final') %>%
	rename(treatment = Group, mouse = 'Mouse ID') %>%
	mutate(treatment = gsub('\\+', '', treatment))
exp_2_wt <- exp_2_data %>%
	select(treatment, mouse, one_of(paste0('D', -1:16))) %>%
	gather(day, weight, -treatment, -mouse)	%>%
	mutate(day = as.numeric(gsub('D', '', day)))
exp_2_cdiff <- exp_2_data %>%
	select(treatment, mouse, contains('difficile')) %>%
	gather(day, cfu, -treatment, -mouse) %>%
	mutate(day = as.numeric(str_match(day, '\\d{1,2}')),
		experiment = 2)

output_df <- full_join(sample_id_df,
	exp_2_cdiff, by = c('treatment', 'mouse', 'day', 'experiment'))

write.table(exp_2_wt, 'data/process/exp_2_weight.txt', 
	row.names = F, quote = F, sep = '\t')
write.table(output_df, 'data/process/metadata.txt', 
	row.names = F, quote = F, sep = '\t')