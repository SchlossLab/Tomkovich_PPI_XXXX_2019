#####################
# 
# clean up labeling issues with data - Treat mice with abx/PPI and challenge with C difficile
# 
# input:
#	stability.opti_mcc.0.03.subsample.shared
#	PPI_Exp_Sample_List_2018.xlsx
#
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
	gather(column, seq_sample, -row, -plate) %>%
	mutate(seq_day = as.numeric(str_match(seq_sample, '(?<=D)\\d{1,2}')),
		seq_mouse = as.numeric(str_match(seq_sample, '(?<=M)\\d{1,2}')),
		seq_treatment = str_match(seq_sample, '[:alpha:]{1,2}(?=pos)'),
		seq_sample_id = paste(seq_treatment, seq_mouse, seq_day, sep = '_'))
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
	mutate(match = seq_sample_id == ext_sample_id,
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

shared_sample_names <- shared_sample_names %>%
	full_join(fixed_shared_names, by = c('Group' = 'shared_names')) %>%
	mutate(shared_samples = ifelse(is.na(new_shared_names), Group, new_shared_names))

sample_id_df <- shared_sample_names %>% 
	mutate(sample_id = gsub("\\+", "\\_", group), # change separator, don't necessarily need to
		sample_id = gsub('^M', 'NAPOS_M', sample_id)) %>% # add NA to samples missing label treatment
	filter(grepl('POS', sample_id)) %>% # remove sequencing or unlabeled samples
	mutate(label_error = grepl('UNK', sample_id), # create column to mark samples that have unknown parts of their label
		sample_id = gsub('_UNK', '', sample_id)) %>% # remove UNK since moved unknonwn label to new column
 	split(., grepl('^D', .$sample_id)) # create label to differentiate the samples with different ordering of info
# separated the piping because we have to split the data and process each group differently based on the order of the group id
sample_id_df <- rbind(sample_id_df$`TRUE` %>% 
		separate(sample_id, c('day', 'treatment', 'mouse')) %>% 
		select(group, mouse, day, treatment, label_error), # make sure columns in both sections are in the same order
	sample_id_df$`FALSE` %>% 
		separate(sample_id, c('treatment', 'mouse', 'day')) %>% 
		select(group, mouse, day, treatment, label_error) # make sure columns in both sections are in the same order 
		) %>% 
	mutate(mouse = as.numeric(gsub('M', '', mouse)), # remove prefix and convert to numeric for day and mouse number
		day = as.numeric(gsub('D', '', day)),
		treatment = gsub('POS', '', treatment),
		treatment = gsub('NA', NA,treatment)) %>%  # remove POS label since all samples are were challenged with C diff (pos)
	### remove all samples with known labeling issues (remove once addressed)
	# identify rows which didnt have a day, mouse or treatment in the label
	mutate(label_error = apply(., 1, function(x)any(is.na(x))) == T | label_error == T) # %>% 
	# filter(label_error == F) 

# load the metafile originally recording the data to check out sequence labels against
true_data <- readxl::read_xlsx('data/raw/PPI_Exp_Sample_List_2018.xlsx', sheet = 'Sheet1') %>% 
	rename(ExpNumber = 'Exp. #', CollectionDay = 'Collection D', CageNumber = 'Cage #', 
		M_F = 'M/F', OmpDailyDose = 'Omp. Daily Dose', MouseID = 'Mouse ID', TubeDayLabel = 'Tube label__1') %>% 
	# Collection Day is shifted for Exp 2 from Day 4 to 7, so use tube label for day
	mutate(day = as.numeric(gsub('D', '', TubeDayLabel)),
 		treatment = gsub('\\+', '', Group),
 		sample_id = paste(treatment, MouseID, day, sep = '_'))
	# if going by collection day, two entries for C_M14_D7, looks to be a shifted copy/paste entry error beginning day 4 of exp 2
	# in column Collection Day, D6 begins a line early (C14, whereas all others outside of D4-7 starts on O1)

######
# below I started to work to identify the mis-labeled samples
######
mutate(true_data, mice_id = paste(treatment, MouseID, sep = '_')) %>%
		pull(mice_id) %>% table %>% data.frame %>% filter(Freq < 10)
	# all mice have 10 samples except O_1, O_2, O_3 (missing day 9, 4, 10)
mutate(true_data, mice_id = paste(treatment, MouseID, sep = '_')) %>% filter(mice_id %in% c('O_1', 'O_2', 'O_3')) %>% select(mice_id, day)
# what days were the samples from
# samples were collected on day 0 2 4 6 7 8 9 10 11 12 14 16 17 (day 6 clinda IP and day 7 cdi)
table(sample_id_df$day)
#	 0  2  3  4  6  7  8  9 10 12 14 16 
#	31 25  1 28 28 30 24 29 26 28 13 14 
table(true_data$day)
#	 0  2  4  6  7  8  9 10 12 14 16 
#	32 32 30 32 33 32 31 31 32 18 14 
# No sample collected for:
# O_M1_D9, O_M2_D4, O_M3_D10, 
# last day depends on exp - Exp 1 (C_1-6, O_7-12, CO_13-18)_D14, Exp 2 (O_1-5, CO_6-10, C_11-14)_D16

# # errors with following samples
filter(sample_id_df, label_error == T)
	# COPOS+M10+D6+UNK < D8
	# COPOS+M7+D0+UNK
	# COPOS+M7+D7+UNK < D8
	# COPOS+M9+D
	# CPOS+M11+D
	# CPOS+M12+D < D4
	# CPOS+M13+D
	# CPOS+M5+D
	# M11+D4
	# M14+D7
	# M6+D6
	# OPOS+M3+D
# # according to the corrected plate file (Sarah_Omep2Box1.xlsx), 
# # COPOS+M10+D6+UNK, COPOS+M7+D0+UNK, COPOS+M7+D7+UNK, M11+D4, M14+D7, M6+D6 
# # should be COPOS_M10_D8, COPOS_M7_D0, COPOS_M7_D8, CPOS_M11_D4, CPOS_M14_D7, COPOS_M6_D6 
# # might still be incorrect, need to confirm

# what samples are duplicated
mutate(sample_id_df, ids = paste(treatment, mouse, day, sep = '_')) %>%
	pull(ids) %>% table %>% data.frame %>% filter(Freq > 1)
	# 			Freq
	# CO_10_6	2 < with samples marked UNK
	# CO_7_0	2 < with samples marked UNK
	# CO_7_7	2 < with samples marked UNK
	# O_10_2	2
	# O_11_6	2
mutate(sample_id_df, ids = paste(treatment, mouse, day, sep = '_')) %>%
	filter(ids %in% c('O_10_2', 'O_11_6'))
filter(sample_id_df, treatment == 'O', mouse %in% c(10, 11)) %>%
	arrange(mouse, day)
	#	D2+OPOS+M10    10   2         O       FALSE
 	#	OPOS+M10+D2    10   2         O       FALSE <- mislabel (mismatch from individual plate to plate master)
 	OPOS+M10+D2 <- COPOS+M10+D2
# which mice were in each treatment
mutate(sample_id_df, mice_id = paste(treatment, mouse, sep = '_')) %>%
		pull(mice_id) %>% table %>% data.frame %>% filter(Freq > 2)
	# clinda (C) - 1 2 3 4 5 6 11 12 13 14 (10 mice total)
	# clinda PPI (CO) - 6 7 8 9 10 13 14 15 16 17 18 (11 mice total)
	# PPI (O) - 1 2 3 4 5 7 8 9 10 11 12 (11 mice total)

mutate(sample_id_df, mice_id = paste(treatment, mouse, sep = '_')) %>%
		pull(mice_id) %>% table %>% data.frame %>% filter(Freq < 5)
#	CO_5 1 <- either C or O OR 6-10,13-18
#	NA_11 1 <- either O or C
#	NA_14 1 <- either C or CO
#	NA_6 1 <- either C or CO
#	O_13 1 <- either CO or C OR 1-5,7-12
#	O_14 1 <- either CO or C OR 1-5,7-12
nonexistent_mouse <- mutate(sample_id_df, mice_id = paste(treatment, mouse, sep = '_')) %>%
		filter(mice_id %in% c('CO_5', 'O_13', 'O_14')) %>% 
		pull(group)
sample_id_df <- mutate(sample_id_df, label_error = ifelse(group %in% nonexistent_mouse, TRUE, label_error))


mutate(true_data, mice_id = paste(treatment, MouseID, sep = '_')) %>% 
filter(mice_id == 'C_12') %>%
data.frame

sample_id_df[which(sample_id_df$day == 3), 'label_error'] <- TRUE # use which to ignore rows with NA
# CPOS+M12+D3
filter(sample_id_df, mouse == 12, treatment == 'C')
# create a label 
# sample_id_df[sample_id_df$group == 'M11+D4', 'treatment'] <- 'C'
# sample_id_df[sample_id_df$group == 'M14+D7', 'treatment'] <- 'C'
# sample_id_df[sample_id_df$group == 'M6+D6', 'treatment'] <- 'CO'
# sample_id_df[sample_id_df$group == 'COPOS+M10+D6+UNK', 'day'] <- 8
# sample_id_df[sample_id_df$group == 'COPOS+M7+D0+UNK', ''] <- 
# sample_id_df[sample_id_df$group == 'COPOS+M7+D7+UNK', ''] <- 
# 
# sample_id_df %>% 
# 	mutate(new_id = paste(treatment, mouse, day, sep = '_')) %>% 
# 	pull(new_id) %>% 
# 	table %>% 
# 	data.frame %>% 
# 	filter(Freq > 1)
# 
# filter(sample_id_df, mouse == 7) 
# sample_id_df %>% 
# 	arrange(treatment, mouse, day)
# 
# table(sample_id_df$day)
# table(sample_id_df$treatment)
# 
#

#table(nmds$mouse_id)
# read in 3d nmds data
nmds <- read.table('data/mothur/stability.opti_mcc.thetayc.0.03.lt.ave.nmds.axes', 
		sep = '\t', header = T, stringsAsFactors = F) %>% 
	right_join(sample_id_df) %>% # join with cleaned up metadata for day, mouse, treatment
	mutate(sample_id = paste(treatment, mouse, day, sep = '_'), 
		mouse_id = paste(treatment, mouse, sep = '_')) %>% # create a mouse id so samples can be grouped by individual
	left_join(select(true_data, ExpNumber, CageNumber, M_F, OmpDailyDose, sample_id)) 

# plot for before abx
nmds %>% 
	filter(ExpNumber == 1) %>% 
	# convert treatment letter labels to more descriptive labels
	mutate(treatment = case_when(treatment == 'O' ~ 'PPI', 
		treatment == 'C' ~ 'Clindamycin',
		treatment == 'CO' ~ 'Clindamycin + PPI'),
	# modify day number for plotting by size
		size = case_when(day <= 6 ~ 1,
		day == 7 ~ 2,
		day > 7 ~ day - 4)) %>% 
	# look at days before 6 to only include pre-clinda treatment
	filter(day < 6) %>% 
	plot_ly(x = ~axis1, y = ~axis2, z = ~axis3, type = 'scatter3d', mode = 'markers',
		color = ~treatment, size = ~size, 
		marker = list(symbol = 'circle', sizemode = 'diameter'), sizes = c(10,20)) %>% # sizes determines the range of sizes c(smallest, largest)
		layout(title = 'Community by Treatment for Experiment 1') # dont forget to change title based on filtering of experiment

# plot for after abx
nmds %>% 
	filter(ExpNumber == 2) %>% 
	# convert treatment letter labels to more descriptive labels
	mutate(treatment = case_when(treatment == 'O' ~ 'PPI',
		treatment == 'C' ~ 'Clindamycin',
		treatment == 'CO' ~ 'Clindamycin + PPI'),
		# modify day number for plotting by size
		size = case_when(day <= 6 ~ 1,
		day == 7 ~ 2,
		day > 7 ~ day - 4)) %>% 
	# look at days after to compare abx treatment + cdiff challenge
	filter(day > 7) %>% 
	plot_ly(x = ~axis1, y = ~axis2, z = ~axis3, type = 'scatter3d', mode = 'markers',
		color = ~treatment, size = ~size, 
		marker = list(symbol = 'circle', sizemode = 'diameter'), sizes = c(10,20)) %>% # sizes determines the range of sizes c(smallest, largest)
		layout(title = 'Community by Treatment for Experiment 2') # dont forget to change title based on filtering of experiment


# plot for ppi
nmds %>% 
	filter(ExpNumber == 2) %>% 
	mutate(treatment = case_when(treatment == 'O' ~ 'PPI',
		treatment == 'C' ~ 'Clindamycin',
		treatment == 'CO' ~ 'Clindamycin + PPI'),
		size = case_when(day <= 6 ~ 1,
		day == 7 ~ 2,
		day > 7 ~ day - 4)) %>% 
	filter(treatment == 'PPI') %>% 
	plot_ly(x = ~axis1, y = ~axis2, z = ~axis3, type = 'scatter3d', mode = 'markers',
		color = ~treatment, size = ~size, 
		marker = list(symbol = 'circle', sizemode = 'diameter'), sizes = c(10,20)) %>% 
		layout(title = 'Community by Treatment for Experiment 2')

#####
# if we wanted to create a line for time instead of sizing by time
#####
#mouse_list <- split(nmds, nmds$mouse_id)
#nmds %>% 
#	plot_ly(x = ~axis1, y = ~axis2, z = ~axis3) %>% 
#		add_trace(mouse_list[[1]], x = mouse_list[[1]]$axis1, y = mouse_list[[1]]$axis2, z = mouse_list[[1]]$axis3)

