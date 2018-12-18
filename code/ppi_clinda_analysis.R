#####################
# 
# analysis of experiment - Treat mice with abx/PPI and challenge with C difficile
# 
# input:
#	stability.opti_mcc.0.03.subsample.shared
#	stability.opti_mcc.thetayc.0.03.lt.ave.nmds.axes
#	PPI_Exp_Sample_List_2018.txt
#
# output:
#	plotly html plot
#		running code will open plot in browser
#		plot can be modified and saved from browser
#
#####################

library(tidyverse)
library(plotly)

shared <- read.table('stability.opti_mcc.0.03.subsample.shared', sep = '\t', header = T, stringsAsFactors = F) %>% 
	rename(group = Group) # set column with sample name to match


sample_id_df <- data.frame(group = shared$group, stringsAsFactors = F) %>% 
	mutate(sample_id = gsub("\\+", "\\_", group), # change separator, don't necessarily need to
		sample_id = gsub('^M', 'NAPOS_M', sample_id)) %>% # add NA to samples missing label treatment
	filter(grepl('POS', sample_id)) %>% # remove sequencing or unlabeled samples
	mutate(label_error = grepl('NA|UNK', sample_id), # create column to mark samples that have unknown parts of their label
		sample_id = gsub('POS', '', sample_id),  # remove POS label since all samples are were challenged with C diff (pos)
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
		day = as.numeric(gsub('D', '', day))) %>% 
	### remove all samples with known labeling issues (remove once addressed)
	filter(label_error == F) 

# load the metafile originally recording the data to check out sequence labels against
true_data <- read.table('PPI_Exp_Sample_List_2018.txt', sep = '\t', header = T) %>% 
	mutate(day = as.numeric(gsub('D', '', CollectionDay)),
 		treatment = gsub('\\+', '', Group),
 		sample_id = paste(treatment, MouseID, day, sep = '_'))

######
# below I started to work to identify the mis-labeled samples
######
# # errors with COPOS+M10+D6+UNK, COPOS+M7+D0+UNK, COPOS+M7+D7+UNK, M11+D4, M14+D7, M6+D6 
# # according to the corrected plate file (Sarah_Omep2Box1.xlsx), might still be incorrect, need to confirm
# # those should be COPOS_M10_D8, COPOS_M7_D0, COPOS_M7_D7, CPOS_M11_D4, CPOS_M14_D7, COPOS_M6_D6 
# # double of following samples
# #  C_14_7    2
# # CO_10_6    2
# #  CO_7_0    2
# #  CO_7_7    2
# #  O_10_2    2
# #  O_11_6    2
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
# head(true_data)
# table(true_data$day)
# table(true_data$treatment)
# true_data %>% 
# 	mutate(new_id = paste(treatment, MouseID, day, sep = '_')) %>% 
# 	pull(new_id) %>% 
# 	table %>% 
# 	data.frame %>% 
# 	filter(Freq > 1)
# # input error with C_14_7 (looks to be shifted sample from D4)
# table(nmds$mouse_id)
# table(mutate(true_data, mouse_id = paste(treatment, MouseID, sep = '_'))$mouse_id)

# read in 3d nmds data
nmds <- read.table('stability.opti_mcc.thetayc.0.03.lt.ave.nmds.axes', 
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

