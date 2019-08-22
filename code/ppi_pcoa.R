#####################
# 
# analysis of experiment - Treat mice with abx and/or PPI and challenge with C difficile
# 
# input:
#	data/mothur/ppi.opti_mcc.braycurtis.0.03.lt.ave.pcoa.axes
# data/mothur/ppi.opti_mcc.braycurtis.0.03.lt.ave.dist
#	data/process/ppi_metadata.txt
#
# output:
#	pcoa plots
#		
#
#####################

library(tidyverse)
library(cowplot)

metadata <- read.table('data/process/ppi_metadata.txt', header = T, sep = '\t', stringsAsFactors = F) %>% 
  mutate(c.diff_colonized = if_else(D9.C..difficile.CFU.g > 0, "colonized", "resistant")) %>% # Create a column to differentiate mice that were colonized with C. difficile from mice that were resistant. Based off of D9 (2 day post challenge) CFU counts.
  mutate(Mouse.ID = as.factor(Mouse.ID)) %>% # Make sure mouse.ID is treated as a factor to use a discrete color scale.
  mutate(Group=factor(Group, levels=c("Clindamycin", "Clind. + Omep.", "Omeprazole"))) #Make sure Group is treated as a factor

# read in pcoa data
pcoa <- read_tsv('data/mothur/ppi.opti_mcc.braycurtis.0.03.lt.ave.pcoa.axes') %>%
  select(group, axis1, axis2) %>% rename(shared_names = group) %>% 
	right_join(metadata, by = "shared_names") # join with cleaned up metadata for day, mouse, treatment
  # Select O+, C+, and CO+ groups. Convert Group letter labels to more descriptive labels
pcoa <- pcoa %>% filter(Group != "NA")

# Define color palette:
color_scheme <- c("#d95f02", "#1b9e77",  "#7570b3")
color_ppi <-  c("#7570b3") # Use for graphs looking at just the PPI group over time
color_cppi <- c("#1b9e77") # Use for graphs looking at just the Clindamycin + PPI group over time
color_c <- c("#d95f02") # Use for graphs looking at just the Clindamycin group over time
color_mouse <- c("#7fc97f", "#beaed4", "#fdc086", "#386cb0", "#f0027f") # To differentiate between 5 individual mice
shape_mouse <- c(16, 17, 15, 3, 7)

#plot all samples
pcoa %>%  
  ggplot(aes(x=axis1, y=axis2, color=Group, alpha = day)) +
  scale_colour_manual(name=NULL, 
                      values=color_scheme, 
                      breaks=c("Clindamycin", "Clind. + Omep.", "Omeprazole"),
                      labels=c("Clindamycin", "Clind. + Omep.", "Omeprazole"))+
  geom_point() +
#  geom_path() + #Add's lines to plots but looks messy  
  theme_classic()+
  labs(title="All groups, all timepoints") +
  theme(plot.title = element_text(hjust = 0.5))


#plot just PPI samples
pcoa %>% filter(Group == "Omeprazole") %>%  
  ggplot(aes(x=axis1, y=axis2, color=Group, alpha = day)) +
  scale_colour_manual(values=color_ppi) +
  geom_point()+
#  geom_path() + #Add's lines to plots but looks messy
  theme_classic() +
  labs(title="PPI-treated mice") +
  theme(plot.title = element_text(hjust = 0.5))

#plot just PPI samples over the first 7 days
pcoa %>% filter(Group == "Omeprazole") %>%  filter(day < 8) %>% 
  ggplot(aes(x=axis1, y=axis2, color=Group, alpha = day)) +
  scale_colour_manual(values=color_ppi) +
  geom_point()+
#  geom_path() + #Add's lines to plots but looks messy
  theme_classic() +
  labs(title="PPI-treated mice over first 7 days") +
  theme(plot.title = element_text(hjust = 0.5))

#plot just Clindamycin samples
pcoa %>% filter(Group == "Clindamycin") %>%  
  ggplot(aes(x=axis1, y=axis2, color=Group, alpha = day)) +
  scale_colour_manual(values=color_c) +
  geom_point() +
#  geom_path() + #Add's lines to plots but looks messy 
  theme_classic()+
  labs(title="Clindamycin-treated mice") +
  theme(plot.title = element_text(hjust = 0.5))

#plot just Clindamycin + PPI samples
pcoa %>% filter(Group == "Clind. + Omep.") %>%  
  ggplot(aes(x=axis1, y=axis2, color=Group, alpha = day)) +
  scale_colour_manual(values=color_cppi) +
  geom_point() +
#  geom_path() + #Add's lines to plots but looks messy
  theme_classic() +
  labs(title="Clindamycin/PPI treated mice") +
  theme(plot.title = element_text(hjust = 0.5))

# plot for before abx
before_abx <- pcoa %>% filter(day < -1) %>% 
  ggplot(aes(x=axis1, y=axis2, color=Group, alpha = day)) +
  scale_colour_manual(values=color_scheme) +
  geom_point() +
#  geom_path() + #Add's lines to plots but looks messy
  theme_classic() +
  labs(title="Before antibiotic treatment") +
  theme(plot.title = element_text(hjust = 0.5))+
  ggsave("results/figures/before_abx.png")
  
#plot before abx with D0 removed
pcoa %>% filter(day < -1) %>% filter(day > -7) %>% 
  ggplot(aes(x=axis1, y=axis2, color=Group, alpha = day)) +
  scale_colour_manual(values=color_scheme) +
  geom_point() +
#  geom_path() + #Add's lines to plots but looks messy
  theme_classic()+
  labs(title="Before antibiotic treatment without initial baseline day") +
  theme(plot.title = element_text(hjust = 0.5))

# plot day after abx
before_plus_day_after_abx <- pcoa %>%	filter(day < 1) %>% 
  ggplot(aes(x=axis1, y=axis2, color=Group, alpha = day)) +
  scale_colour_manual(values=color_scheme) +
  scale_alpha_continuous(range = c(.2, 1))+
  geom_point() +
#  geom_path() + #Add's lines to plots but looks messy
  theme_classic()+
  labs(title="Before antibiotic treatment + 1 day after") +
  theme(plot.title = element_text(hjust = 0.5))+
  ggsave("results/figures/before_plus_day_after_abx.png")

# Figure 1B----
# plot 1st 7 days: Includes abx treatment day. 
#Get a list of the shared_names from metadata file that are just the omeprazole treated mice before antibiotic treatment
ppi_before_challenge_subset <- metadata %>% 
  filter(day < 1, Group == "Omeprazole") %>% 
  pull(shared_names)
#groups specification for running dist.shared & pcoa commands on just this subset of samples in mothur". See README for instructions on analysis using mothur commands.
# Opos+M1+D0-Opos+M2+D0-Opos+M3+D0-Opos+M4+D0-Opos+M5+D0-Opos+M1+D2-Opos+M2+D2-Opos+M3+D2-Opos+M4+D2-Cpos+M5+D-Opos+M1+D4-Opos+M3+D4-Opos+M4+D4-Opos+M5+D4-Opos+M11+D6-Opos+M2+D6-Opos+M3+D-Opos+M4+D6-Opos+M5+D6-Opos+M1+D7-Opos+M2+D7-Opos+M3+D7-Opos+M4+D7-Opos+M5+D7
#Subsetted distance matrices and pcoa files outputted by mothur are located in data/mothur/subset

# Read in PCoA file for subset of PPI samples taken before C. difficile challenge & merge to metadata
pcoa_subset <- read_tsv('data/mothur/subset/ppi.opti_mcc.braycurtis.0.03.lt.ave.pcoa.axes') %>%
  select(group, axis1, axis2) %>% rename(shared_names = group) %>% 
  right_join(metadata, by = "shared_names") # join with cleaned up metadata for day, mouse, treatment
# Select O+, C+, and CO+ groups. Convert Group letter labels to more descriptive labels
pcoa_subset <- pcoa_subset %>% filter(Group != "NA")
# Plot of samples from omeprazole-treated mice before challenge
pcoa_before_challenge <- pcoa_subset %>% 
  filter(day < 1, Group == "Omeprazole", Mouse.ID == as.factor(Mouse.ID)) %>%	
  ggplot(aes(x=axis1, y=axis2, color=Mouse.ID, shape=Mouse.ID)) +
  labs(color = "Omeprazole Mice")+
  scale_color_manual(name="Omeprazole Mice", values=color_mouse) + 
  scale_shape_manual(name= "Omeprazole Mice", values=shape_mouse)+
  geom_point() +
  geom_path() +
  theme_classic()+
  labs(x = "PCoA 1",
       y = "PCoA 2") +
  xlim(-0.4, 0.4)+
  ylim(-0.4, 0.4)+
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.position = c(0.15, 0.8)) #Move legend position
save_plot("results/figures/before_C._diff_challenge.png", pcoa_before_challenge, base_aspect_ratio = 2) #Use save_plot over ggsave because it works better with cowplot

# Figure 2B----
shape_legend <-c(expression(paste(italic("C. difficile"), "\ status"))) #Expression variable for the title so that bacteria name will be in italics
# plot after abx & C. diff. Colonized mice are represented by x shapes. Resistant mice are represented as circles.
pcoa_after_challenge <- pcoa %>% 
  filter(day > 0) %>% 
  ggplot(aes(x=axis1, y=axis2, color=Group, alpha = day, shape = c.diff_colonized)) +
  scale_colour_manual(name=NULL, 
                       values=color_scheme, 
                       breaks=c("Clindamycin", "Clind. + Omep.", "Omeprazole"),
                       labels=c("Clindamycin", "Clind. + Omep.", "Omeprazole")) +
  scale_shape_manual(name=shape_legend,
                     values=c(4, 19),
                     breaks=c("colonized", "resistant"),
                     labels=c("colonized", "resistant")) +
  scale_alpha_continuous(range = c(.3, 1),
                         breaks=c(1, 3, 6, 9),
                         labels=c(1, 3, 6, 9))+
  labs(alpha = "Day")+
  geom_point() +
  theme_classic() +
  labs(x = "PCoA 1",
       y = "PCoA 2") +
  theme(plot.title = element_text(hjust = 0.5))+
save_plot("results/figures/after_abx_C.diff.png", pcoa_after_challenge, base_aspect_ratio = 2)

