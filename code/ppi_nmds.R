#####################
# 
# analysis of experiment - Treat mice with abx/PPI and challenge with C difficile
# 
# input:
#	data/mothur/ppi.opti_mcc.0.03.subsample.shared
#	data/mothur/ppi.opti_mcc.thetayc.0.03.lt.ave.nmds.axes
#	data/process/ppi_metadata.txt
#
# output:
#	plotly html plot
#		running code will open plot in browser
#		plot can be modified and saved from browser
#
#####################

library(tidyverse)
library(plotly)

shared <- read.table('data/mothur/ppi.opti_mcc.0.03.subsample.shared', sep = '\t', header = T, stringsAsFactors = F) %>% 
	rename(shared_names = Group) # set column with sample name to matching

metadata <- read.table('data/process/ppi_metadata.txt', header = T, sep = '\t', stringsAsFactors = F)

# read in 3d nmds data
nmds <- read.table('data/mothur/ppi.opti_mcc.thetayc.0.03.lt.ave.nmds.axes', 
		sep = '\t', header = T, stringsAsFactors = F) %>% rename(shared_names = group) %>% 
	right_join(metadata, by = "shared_names") # join with cleaned up metadata for day, mouse, treatment
	

# plot for before abx
nmds %>% 
	# convert Group letter labels to more descriptive labels
	mutate(Group = case_when(Group == 'O+' ~ 'PPI', 
		Group == 'C+' ~ 'Clindamycin',
		Group == 'CO+' ~ 'Clindamycin + PPI'),
	# modify day number for plotting by size
		size = case_when(day <= 6 ~ 1,
		day == 7 ~ 2,
		day > 7 ~ day - 4)) %>% 
	# look at days before 7 to only include pre-clinda treatment. Mice were treated with clinda on D6, so stool
  # sample from D6 is still pre-clinda.
	filter(day < 7) %>% 
	plot_ly(x = ~axis1, y = ~axis2, z = ~axis3, type = 'scatter3d', mode = 'markers',
		color = ~Group, size = ~size, 
		marker = list(symbol = 'circle', sizemode = 'diameter'), sizes = c(10,20)) %>% # sizes determines the range of sizes c(smallest, largest)
		layout(title = 'Community by Treatment')

# plot for before abx
nmds %>% 
	# convert Group letter labels to more descriptive labels
	mutate(Group = case_when(Group == 'O+' ~ 'PPI',
		Group == 'C+' ~ 'Clindamycin',
		Group == 'CO+' ~ 'Clindamycin + PPI'),
		# modify day number for plotting by size
		size = case_when(day <= 6 ~ 1,
		day == 7 ~ 2,
		day > 7 ~ day - 4)) %>% 
	# look at days after to compare abx Group + cdiff challenge
	filter(day < 7) %>% 
	plot_ly(x = ~axis1, y = ~axis2, z = ~axis3, type = 'scatter3d', mode = 'markers',
		color = ~Group, size = ~size, 
		marker = list(symbol = 'circle', sizemode = 'diameter'), sizes = c(10,20)) %>% # sizes determines the range of sizes c(smallest, largest)
		layout(title = 'Community by Treatment before antibiotics') # dont forget to change title based on filtering of experiment


# plot day after abx
nmds %>% 
  # convert Group letter labels to more descriptive labels
  mutate(Group = case_when(Group == 'O+' ~ 'PPI',
                           Group == 'C+' ~ 'Clindamycin',
                           Group == 'CO+' ~ 'Clindamycin + PPI'),
         # modify day number for plotting by size
         size = case_when(day <= 6 ~ 1,
                          day == 7 ~ 2,
                          day > 7 ~ day - 4)) %>% 
	filter(day == 7) %>% 
	plot_ly(x = ~axis1, y = ~axis2, z = ~axis3, type = 'scatter3d', mode = 'markers',
		color = ~Group, size = ~size, 
		marker = list(symbol = 'circle', sizemode = 'diameter'), sizes = c(10,20)) %>% 
		layout(title = 'Community by Treatment after antibiotics')

# plot after abx & C. diff
nmds %>% 
  # convert Group letter labels to more descriptive labels
  mutate(Group = case_when(Group == 'O+' ~ 'PPI',
                           Group == 'C+' ~ 'Clindamycin',
                           Group == 'CO+' ~ 'Clindamycin + PPI'),
         # modify day number for plotting by size
         size = case_when(day <= 6 ~ 1,
                          day == 7 ~ 2,
                          day > 7 ~ day - 4)) %>% 
  filter(day > 7) %>% 
  plot_ly(x = ~axis1, y = ~axis2, z = ~axis3, type = 'scatter3d', mode = 'markers',
          color = ~Group, size = ~size, 
          marker = list(symbol = 'circle', sizemode = 'diameter'), sizes = c(10,20)) %>% 
  layout(title = 'Community by Treatment after antibiotics and C. difficile infection')

#####
# if we wanted to create a line for time instead of sizing by time
#####
#mouse_list <- split(nmds, nmds$mouse_id)
#nmds %>% 
#	plot_ly(x = ~axis1, y = ~axis2, z = ~axis3) %>% 
#		add_trace(mouse_list[[1]], x = mouse_list[[1]]$axis1, y = mouse_list[[1]]$axis2, z = mouse_list[[1]]$axis3)

