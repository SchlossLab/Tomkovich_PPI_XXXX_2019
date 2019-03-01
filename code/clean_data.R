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

shared_sample_names <- read.table('data/mothur/stability.opti_mcc.0.03.subsample.shared', 
		sep = '\t', header = T, stringsAsFactors = F) %>% 
	select(Group) # select column with sample names of shared file which we will need to match our metadata to

# load in all the names and positions of the sequence sample names
seq_plate_map <- 'data/raw/SarahPPiNov2018_plateMap.xlsx'
seq_sample_names <- bind_rows(
		mutate(read_xlsx(seq_plate_map, sheet = 'Final Plate Map', range = 'C5:O13'), plate = 1),
		mutate(read_xlsx(seq_plate_map, sheet = 'Final Plate Map', range = 'C16:O24'), plate = 2),
		mutate(read_xlsx(seq_plate_map, sheet = 'Final Plate Map', range = 'C27:O35'), plate = 3),
		mutate(read_xlsx(seq_plate_map, sheet = 'Final Plate Map', range = 'C39:O47'), plate = 4)) %>%
	rename(row = X__1) %>%
	gather(column, seq_sample, -row, -plate)
# load in all the names and positions of the extraction plates
extraction123 <- 'data/raw/Extraction Template_Plate #1#2#3.xlsx'
extraction4 <- 'data/raw/Plate_4_Sarah_Omep2Box1.xlsx'
ext_sample_names <- bind_rows(
		mutate(read_xlsx(extraction123, sheet = 'Plate 1', range = 'A2:L9', col_names = F), row = LETTERS[1:8], plate = 1),
		mutate(read_xlsx(extraction123, sheet = 'Plate 2', range = 'A3:L10', col_names = F), row = LETTERS[1:8], plate = 2),
		mutate(read_xlsx(extraction123, sheet = 'Plate 3', range = 'B4:M11', col_names = F), row = LETTERS[1:8], plate = 3),
		mutate(read_xlsx(extraction4, sheet = 'Sheet1', range = 'A1:L8', col_names = F), row = LETTERS[1:8], plate = 4)) %>%
	gather(column, ext_sample, -row, -plate) %>%
	mutate(column = gsub('X__', '', column), 
		ext_day = as.numeric(str_match(ext_sample, '(?<=D)\\d{1,2}')),
		ext_mouse = as.numeric(str_match(ext_sample, '(?<=M)\\d{1,2}')),
		ext_treatment = str_match(ext_sample, '[:alpha:]{1,2}(?=\\+)'),
		ext_sample_id = paste(ext_treatment, ext_mouse, ext_day, sep = '_'))
# check to see if all ext and seq samples match
fixed_shared_names <- full_join(seq_sample_names, ext_sample_names, by = c('row', 'column', 'plate')) %>%
	# find which samples don't match between their seqencing and their extraction (or were labeled with unknown)
	mutate(match = seq_sample == ext_sample_id,
		match = ifelse(grepl('U|u', seq_sample), FALSE, match)) %>%
	filter(match == F) %>% 
	select(seq_sample, ext_sample, ext_treatment, ext_mouse, ext_day) %>%
	# convert names of seq to match the format of the shared samples
	mutate(shared_names = gsub('^_{1,3}','', seq_sample),
		shared_names = gsub('_', '+', shared_names),
		shared_names = gsub('pos', 'POS', shared_names),
		shared_names = gsub('Unk|unk', 'UNK', shared_names)) %>% 
	# select only the samples that are present in the shared data
	inner_join(shared_sample_names, by = c('shared_names' = 'Group')) %>%
	mutate(new_shared_names = paste0(ext_treatment, 'POS+M', ext_mouse, '+D', ext_day)) %>%
	select(shared_names, new_shared_names)
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
#	#	CO_8_0    2 <- two positions on extraction 4 plate marked same
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