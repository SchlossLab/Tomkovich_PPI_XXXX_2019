#####################
# 
# analysis of experiment - Treat mice with abx and/or PPI and challenge with C difficile
# 
# input:
#	data/mothur/ppi.opti_mcc.braycurtis.0.03.lt.ave.pcoa.axes
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
  mutate(abx_status = if_else(day > -8 & day < 0, "pre", "post")) %>% # Create a column to differentiate between timepoints that are from before or after exposure to the antibiotic clindamycin
  mutate(c.diff_colonized = if_else(D9.C..difficile.CFU.g > 0, "colonized", "resistant")) # Create a column to differentiate mice that were colonized with C. difficile from mice that were resistant. Based off of D9 (2 day post challenge) CFU counts.

# read in pcoa data
pcoa <- read_tsv('data/mothur/ppi.opti_mcc.braycurtis.0.03.lt.ave.pcoa.axes') %>%
  select(group, axis1, axis2) %>% rename(shared_names = group) %>% 
	right_join(metadata, by = "shared_names") # join with cleaned up metadata for day, mouse, treatment
  # Select O+, C+, and CO+ groups. Convert Group letter labels to more descriptive labels
pcoa <- pcoa %>% filter(Group != "NA")

# Define color palette:
color_scheme <- c("#d95f02", "#1b9e77", "#7570b3")
color_ppi <-  c("#7570b3") # Use for graphs looking at just the PPI group over time
color_cppi <- c("#1b9e77") # Use for graphs looking at just the Clindamycin + PPI group over time
color_c <- c("#d95f02") # Use for graphs looking at just the Clindamycin group over time

#plot all samples
pcoa %>%  
  ggplot(aes(x=axis1, y=axis2, color=Group, alpha = day)) +
  scale_colour_manual(values=color_scheme) +
  geom_point() +
#  geom_path() + #Add's lines to plots but looks messy  
  theme_classic()+
  labs(title="All groups, all timepoints") +
  theme(plot.title = element_text(hjust = 0.5))


#plot just PPI samples
pcoa %>% filter(Group == "PPI") %>%  
  ggplot(aes(x=axis1, y=axis2, color=Group, alpha = day)) +
  scale_colour_manual(values=color_ppi) +
  geom_point()+
#  geom_path() + #Add's lines to plots but looks messy
  theme_classic() +
  labs(title="PPI-treated mice") +
  theme(plot.title = element_text(hjust = 0.5))

#plot just PPI samples over the first 7 days
pcoa %>% filter(Group == "PPI") %>%  filter(day < 8) %>% 
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
pcoa %>% filter(Group == "Clindamycin + PPI") %>%  
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

# plot 1st 7 days: Includes abx treatment day. Pre-clindamycin treatment represented by circles.
# Post-clindamycin treatment represented by open diamonds.
pcoa %>% filter(day < 1) %>%	
  ggplot(aes(x=axis1, y=axis2, color=Group, alpha = day, shape = abx_status)) +
  scale_colour_manual(values=color_scheme) +
  scale_shape_manual(values=c(5, 19)) +
  scale_alpha_continuous(range = c(.3, 1))+
  geom_point() +
#  geom_path() + #Add's lines to plots but looks messy
  theme_classic()+
  labs(title="PCoA of timepoints before spore challenge") +
  theme(plot.title = element_text(hjust = 0.5))+
  theme(legend.position="none")
  ggsave("results/figures/before_C._diff_challenge.png")

# plot after abx & C. diff. Colonized mice are represented by x shapes. Resistant mice are represented as circles.
pcoa %>% filter(day > 1) %>% 
  ggplot(aes(x=axis1, y=axis2, color=Group, alpha = day, shape = c.diff_colonized)) +
  scale_colour_manual(values=color_scheme) +
  scale_shape_manual(values=c(4, 19)) +
  scale_alpha_continuous(range = c(.3, 1))+
  geom_point() +
#  geom_path() + #Add's lines to plots but looks messy
  theme_classic() +
  labs(title="PCoA after spore challenge") +
  theme(plot.title = element_text(hjust = 0.5))+
  theme(legend.position="none")
  ggsave("results/figures/after_abx_C.diff.png")

