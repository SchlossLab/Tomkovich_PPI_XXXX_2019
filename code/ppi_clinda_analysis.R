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

shared <- read.table('data/mothur/stability.opti_mcc.0.03.subsample.shared', sep = '\t', header = T, stringsAsFactors = F) %>% 
	rename(group = Group) # set column with sample name to match

metadata <- read.table('data/process/metadata.txt', header = T, sep = '\t', stringsAsFactors = F)

# read in 3d nmds data
nmds <- read.table('data/mothur/stability.opti_mcc.thetayc.0.03.lt.ave.nmds.axes', 
		sep = '\t', header = T, stringsAsFactors = F) %>% 
	right_join(metadata, by = c('group' = 'Group')) %>% # join with cleaned up metadata for day, mouse, treatment
	mutate(mouse_id = paste(treatment, mouse, sep = '_')) # create a mouse id so samples can be grouped by individual 

# plot for before abx
nmds %>% 
	filter(experiment == 1) %>% 
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

