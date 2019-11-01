#####################
#####################
#Setup----
library(tidyverse)
library(broom)
library(cowplot)

# Define color palette:
color_scheme <- c("#d95f02", "#1b9e77",  "#7570b3")
color_ppi <-  c("#7570b3") # Use for graphs looking at just the PPI group over time
color_cppi <- c("#1b9e77") # Use for graphs looking at just the Clindamycin + PPI group over time
color_c <- c("#d95f02") # Use for graphs looking at just the Clindamycin group over time
color_day <- c('#c7c5e0' ,"#7570b3", "#514e7d") # Tints and Shades of PPI color for graphs showing relative abundances over time

# Import metadata into data frame
metadata <- read.table('data/process/ppi_metadata.txt', header = T, sep = '\t', stringsAsFactors = F) %>% 
  filter(Group != "NA") %>% #Exclude the mock community
  mutate(Group=factor(Group, levels=c("Clindamycin", "Clind. + Omep.", "Omeprazole"))) # make sure Group is treated as a factor

  
# Import taxonomy into data frame and clean up taxonomy names
taxonomy <- read_tsv(file="data/mothur/ppi.taxonomy") %>%
  rename_all(tolower) %>% #remove uppercase from column names
# Split taxonomic information into separate columns for each taxonomic level  
  mutate(taxonomy=str_replace_all(taxonomy, c("\\(\\d*\\)" = "", #drop digits with parentheses around them
                                              ';$' = "", #removes semi-colon at end of line
                                              'Bacteria_unclassified' = 'Unclassified',
                                              "Clostridium_" = "Clostridium ", #Remove underscores after Clostridium
                                              "_unclassified" = " Unclassified"))) %>% 
# Separate taxonomic levels into separate columns according to semi-colon.
  separate(taxonomy, into=c("kingdom", "phylum", "class", "order", "family", "genus"), sep=';')

# Import otu_data for PPI experiment samples
otu_data <- read_tsv("data/mothur/ppi.opti_mcc.0.03.subsample.shared", col_types=cols(Group=col_character())) %>% 
  select(-label, -numOtus) %>% 
  rename(shared_names=Group) %>% 
  gather(-shared_names, key="otu", value="count") %>% 
  mutate(rel_abund=count/3000) #Use 3000, because this is the subsampling parameter chosen.

# Merge otu_data to taxonomy data frame
agg_taxa_data <- inner_join(otu_data, taxonomy) 

# Summarize relative abundance data for genus level
agg_genus_data <- agg_taxa_data %>% 
  group_by(shared_names, genus) %>% 
  summarize(agg_rel_abund=sum(rel_abund)) %>% 
  # Merge relative abundance data to phyla data
  inner_join(., metadata, by = "shared_names") %>% 
  ungroup() 

# Summarize relative abundance data for family level
agg_family_data <- agg_taxa_data %>% 
  group_by(shared_names, family) %>% 
  summarize(agg_rel_abund=sum(rel_abund)) %>% 
  # Merge relative abundance data to phyla data
  inner_join(., metadata, by = "shared_names") %>% 
  ungroup() 


#Hypothesis testing & plotting----

#Kruskal_wallis test for family differences across treatment groups with Benjamini-Hochburg correction Day -7, 0, 2, & 9 of the experiment
family_tests_start_day <- agg_family_data %>% 
  filter(day == -7) %>%
  group_by(family) %>% 
  do(tidy(kruskal.test(agg_rel_abund~factor(Group), data=.))) %>% ungroup() %>% 
  mutate(p.value.adj=p.adjust(p.value, method="BH")) %>% 
  arrange(p.value.adj)

family_tests_day0 <- agg_family_data %>% 
  filter(day == 0) %>%
  group_by(family) %>% 
  do(tidy(kruskal.test(agg_rel_abund~factor(Group), data=.))) %>% ungroup() %>% 
  mutate(p.value.adj=p.adjust(p.value, method="BH")) %>% 
  arrange(p.value.adj)

family_tests_day2 <- agg_family_data %>% 
  filter(day == 2) %>%
  group_by(family) %>% 
  do(tidy(kruskal.test(agg_rel_abund~factor(Group), data=.))) %>% ungroup() %>% 
  mutate(p.value.adj=p.adjust(p.value, method="BH")) %>% 
  arrange(p.value.adj)

family_tests_day9 <- agg_family_data %>% 
  filter(day == 9) %>%
  group_by(family) %>% 
  do(tidy(kruskal.test(agg_rel_abund~factor(Group), data=.))) %>% ungroup() %>% 
  mutate(p.value.adj=p.adjust(p.value, method="BH")) %>% 
  arrange(p.value.adj)

#List the significant family based on treatment groups after Benjamini-Hochburg correction
sig_family_start_day <- family_tests_start_day %>% 
  filter(p.value.adj <= 0.05) %>% 
  pull(family)

sig_family_day0 <- family_tests_day0 %>% 
  filter(p.value.adj <= 0.05) %>% 
  pull(family)

sig_family_day2 <- family_tests_day2 %>% 
  filter(p.value.adj <= 0.05) %>% 
  pull(family)

sig_family_day9 <- family_tests_day9 %>% 
  filter(p.value.adj <= 0.05) %>% 
  pull(family)

# No families are significantly different on day -7, 0, 2, & 9 after FDR correction

# Figure 1C----
#Graph the families associated with human PPI use according to the literature [@Imhann2017]
ppi_family <- c("Enterococcaceae", "Lactobacillaceae", "Micrococcaceae", "Staphylococcaceae", "Streptococcaceae", "Ruminococcaceae")
ppi_family_plot_neg7 <- agg_family_data %>% 
  filter(family %in% ppi_family) %>% 
  filter(day == -7) %>% 
  mutate(family=factor(family, levels=ppi_family)) %>% 
  mutate(agg_rel_abund = agg_rel_abund + 1/6000) %>% 
  ggplot(aes(x= reorder(family, agg_rel_abund), y=agg_rel_abund, color=Group))+
  scale_colour_manual(values=color_scheme) +
  geom_hline(yintercept=1/3000, color="gray")+
  geom_boxplot(outlier.shape = NA, size = 1.2)+
  geom_jitter(shape=19, size=2, alpha=0.6, position=position_jitterdodge(dodge.width=0.7, jitter.width=0.2)) +
  labs(title=NULL, 
       x=NULL,
       y="Relative abundance (%)",
       color = "Day -7")+
  scale_y_log10(breaks=c(1e-4, 1e-3, 1e-2, 1e-1, 1), labels=c(1e-2, 1e-1, 1, 10, 100), limits = c(1e-4, 1))+
  coord_flip()+
  theme_classic()+
  theme(plot.title=element_text(hjust=0.5))+
  theme(axis.text.y = element_text(face = "italic"))+ #Have the families show up as italics
  theme(plot.title=element_text(hjust=0.5))+
  theme(legend.position = c(0.85, 0.2)) + #Move legend position
  theme(text = element_text(size = 16))  # Change font size for entire plot
save_plot("results/figures/families_prev_assoc_w_PPIs_-7.png", ppi_family_plot_neg7, base_aspect_ratio = 2)

# Figure 1D----
ppi_family_plot_day0 <- agg_family_data %>% 
  filter(family %in% ppi_family) %>% 
  filter(day == 0) %>% 
  mutate(family=factor(family, levels=ppi_family)) %>% 
  mutate(agg_rel_abund = agg_rel_abund + 1/6000) %>% 
  ggplot(aes(x= reorder(family, agg_rel_abund), y=agg_rel_abund, color=Group))+
  scale_colour_manual(values=color_scheme) +
  geom_hline(yintercept=1/3000, color="gray")+
  geom_boxplot(outlier.shape = NA, size = 1.2)+
  geom_jitter(shape=19, size=2, alpha=0.6, position=position_jitterdodge(dodge.width=0.7, jitter.width=0.2)) +
  labs(title=NULL, 
       x=NULL,
       y="Relative abundance (%)",
       color = "Day 0")+
  scale_y_log10(breaks=c(1e-4, 1e-3, 1e-2, 1e-1, 1), labels=c(1e-2, 1e-1, 1, 10, 100), limits = c(1e-4, 1))+
  coord_flip()+
  theme_classic()+
  theme(plot.title=element_text(hjust=0.5))+
  theme(plot.title=element_text(hjust=0.5))+
  theme(axis.text.y = element_text(face = "italic"))+ #Have the families show up as italics
  theme(legend.position = c(0.85, 0.2)) + #Move legend position
  theme(text = element_text(size = 16))  # Change font size for entire plot
save_plot("results/figures/families_prev_assoc_w_PPIs_0.png", ppi_family_plot_day0, base_aspect_ratio = 2)

# FDR corrected P values for these families analyzed at day -7 & day 0
ppi_fdr_family_start_day <- family_tests_start_day %>% 
  filter(family %in% ppi_family) %>% 
  select(family, p.value.adj)
# all P > 0.54 and Enterococcaceae is NaN (Not a number)
ppi_fdr_family_day0 <- family_tests_day0 %>% 
  filter(family %in% ppi_family) %>% 
  select(family, p.value.adj)
# all P > 0.2 and Staphylococcaceae, Streptococcaceae, and Enterococcaceae are NaN (Not a number)

#Kruskal_wallis test for genus differences across treatment groups with Benjamini-Hochburg correction 
genus_tests_day0 <- agg_genus_data %>% 
  filter(day == 0) %>% 
  group_by(genus) %>% 
  do(tidy(kruskal.test(agg_rel_abund~factor(Group), data=.))) %>% ungroup() %>% 
  mutate(p.value.adj=p.adjust(p.value, method="BH")) %>% 
  arrange(p.value.adj)

genus_tests_day2 <- agg_genus_data %>% 
  group_by(genus) %>% 
  filter(day == 2) %>% 
  do(tidy(kruskal.test(agg_rel_abund~factor(Group), data=.))) %>% ungroup() %>% 
  mutate(p.value.adj=p.adjust(p.value, method="BH")) %>% 
  arrange(p.value.adj)

genus_tests_day9 <- agg_genus_data %>% 
  group_by(genus) %>% 
  filter(day == 9) %>% 
  do(tidy(kruskal.test(agg_rel_abund~factor(Group), data=.))) %>% ungroup() %>% 
  mutate(p.value.adj=p.adjust(p.value, method="BH")) %>% 
  arrange(p.value.adj)

#List the significant genera based on treatment groups after Benjamini-Hochburg correction
sig_genus_day0 <- genus_tests_day0 %>% 
  filter(p.value.adj <= 0.05) %>% 
  pull(genus)
sig_genus_day2 <- genus_tests_day2 %>% 
  filter(p.value.adj <= 0.05) %>% 
  pull(genus)
sig_genus_day9 <- genus_tests_day9 %>% 
  filter(p.value.adj <= 0.05) %>% 
  pull(genus)
#No genera are significantly different accross treatment groups on day 0, 2, and 8

top_6_genera_across_groups <- top_n(genus_tests_day2, -6, p.value.adj) %>%  #negative to pull rows with the lowest values
  pull(genus)

# Figure 2C----
#Graph the top 6 genera based on treatment groups, no genera are significant after Benjamini-Hochburg correction
group_genera <- agg_genus_data %>% 
  filter(genus %in% top_6_genera_across_groups) %>% 
  filter(day == 2) %>% 
  mutate(genus=factor(genus, top_6_genera_across_groups)) %>% 
  mutate(agg_rel_abund = agg_rel_abund + 1/6000) %>% 
  ggplot(aes(x= reorder(genus, agg_rel_abund), y=agg_rel_abund, color=Group))+
  scale_colour_manual(name="Day 2",
                      values=color_scheme) +
  geom_hline(yintercept=1/3000, color="gray")+
  geom_boxplot(outlier.shape = NA, size = 1.2)+
  geom_jitter(shape=19, size=2, alpha=0.6, position=position_jitterdodge(dodge.width=0.7, jitter.width=0.2)) +
  labs(title=NULL, 
       x=NULL,
       y="Relative abundance (%)")+
  scale_y_log10(breaks=c(1e-4, 1e-3, 1e-2, 1e-1, 1), labels=c(1e-2, 1e-1, 1, 10, 100))+
  coord_flip()+
  theme_classic()+
  theme(axis.text.y = element_text(face = "italic"))+ #Have the genera show up as italics
  theme(plot.title=element_text(hjust=0.5))+
  theme(legend.position = c(0.85, 0.2)) + #Get rid of legend title & move legend position
  theme(text = element_text(size = 16))  # Change font size for entire plot
save_plot("results/figures/genera_assoc_w_treatment.png", group_genera, base_aspect_ratio = 2)
   
#Kruskal_wallis test for family differences across time in the PPI group with Benjamini-Hochburg correction---- 
ppi_family_tests <- agg_family_data %>% 
  group_by(family) %>% 
  filter(Group == "Omeprazole", day %in% c("-7", "0", "9")) %>%
  do(tidy(kruskal.test(agg_rel_abund~factor(day), data=.))) %>% ungroup() %>% 
  mutate(p.value.adj=p.adjust(p.value, method="BH")) %>% 
  arrange(p.value.adj)

#List the significant family across time in the PPI group after Benjamini-Hochburg correction
ppi_sig_family <- ppi_family_tests %>% 
  filter(p.value.adj <= 0.05) %>% 
  pull(family)
#No significant families after Benjamini-Hochburg correction. 

#Kruskal_wallis test for genus differences across time in the PPI group with Benjamini-Hochburg correction---- 
ppi_genus_tests <- agg_genus_data %>% 
  group_by(genus) %>% 
  filter(Group == "Omeprazole", day %in% c("-7", "0", "9")) %>%
  do(tidy(kruskal.test(agg_rel_abund~factor(day), data=.))) %>% ungroup() %>% 
  mutate(p.value.adj=p.adjust(p.value, method="BH")) %>% 
  arrange(p.value.adj)

#List the significant genera across time in the PPI group after Benjamini-Hochburg correction
ppi_sig_genus <- ppi_genus_tests %>% 
  filter(p.value.adj <= 0.05) %>% 
  pull(genus)
#No significant genera after Benjamini-Hochburg correction. 


# Figure S1A----
#Lactobacillaceae family relative abundance over time with mean relative abundance summary line and individual mouse relative abundance points
#Data frame for mean Lactobacillaceae family relative abundance for each group on each day
lacto_family_mean <- agg_family_data %>% 
  filter(family == "Lactobacillaceae") %>% 
  group_by(Group, day) %>% 
  summarize(mean=(mean(agg_rel_abund + 1/6000))) %>% 
  ungroup
# Data frame for individual mouse relative abundance points
lacto_family_mice <-  agg_family_data %>% 
  filter(family == "Lactobacillaceae") %>% 
  mutate(agg_rel_abund = agg_rel_abund + 1/6000) %>%
  select(Group, day, agg_rel_abund, family)

#Plot of Lactobacillaceae family relative abundance over time with mean relative abundance summary line and individual mouse relative abundance points
lacto_family_time <- ggplot(NULL)+
  geom_point(lacto_family_mice, mapping = aes(x=day, y=agg_rel_abund, color=Group, alpha = .2), show.legend = FALSE, size = 2.5)+
  geom_line(lacto_family_mean, mapping = aes(x=day, y=mean, color=Group), size = 1)+
  scale_colour_manual(values=color_scheme) +
  geom_hline(yintercept=1/3000, color="gray")+
  labs(title="Lactobacillaceae",
       x="Day",
       y="Relative abundance (%)") +
  scale_y_log10(breaks=c(1e-4, 1e-3, 1e-2, 1e-1, 1), labels=c(1e-2, 1e-1, 1, 10, 100))+
  theme_classic()+
  theme(plot.title=element_text(hjust=0.5, face = "italic"))+
  theme(legend.title=element_blank(), legend.position = c(0.12, 0.22))+
  theme(text = element_text(size = 16))  # Change font size for entire plot
save_plot("results/figures/lactobacillaceae_time.png", lacto_family_time, base_aspect_ratio = 2)

# Figure S1B---- 
#Ruminococcaceae family relative abundance over time with mean relative abundance summary line and individual mouse relative abundance points
#Data frame for mean Ruminococcaceae family relative abundance for each group on each day
rumino_family_mean <- agg_family_data %>% 
  filter(family == "Ruminococcaceae") %>% 
  group_by(Group, day) %>% 
  summarize(mean=(mean(agg_rel_abund + 1/6000))) %>% 
  ungroup
# Data frame for individual mouse relative abundance points
rumino_family_mice <-  agg_family_data %>% 
  filter(family == "Ruminococcaceae") %>% 
  mutate(agg_rel_abund = agg_rel_abund + 1/6000) %>%
  select(Group, day, agg_rel_abund, family)
#Plot of Ruminococcaceae family relative abundance over time with mean relative abundance summary line and individual mouse relative abundance points
rumino_family_time <- ggplot(NULL)+
  geom_point(rumino_family_mice, mapping = aes(x=day, y=agg_rel_abund, color=Group, alpha = .2), show.legend = FALSE, size = 2.5)+
  geom_line(rumino_family_mean, mapping = aes(x=day, y=mean, color=Group), size = 1)+
  scale_colour_manual(values=color_scheme) +
  geom_hline(yintercept=1/3000, color="gray")+
  labs(title="Ruminococcaceae",
       x="Day",
       y="Relative abundance (%)") +
  scale_y_log10(breaks=c(1e-4, 1e-3, 1e-2, 1e-1, 1), labels=c(1e-2, 1e-1, 1, 10, 100))+
  theme_classic()+
  theme(plot.title=element_text(hjust=0.5, face = "italic"))+
  theme(legend.title=element_blank(), legend.position = c(0.9, 0.2))+
  theme(text = element_text(size = 16))  # Change font size for entire plot
save_plot("results/figures/ruminococcaceae_time.png", rumino_family_time, base_aspect_ratio = 2)



