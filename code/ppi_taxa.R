#####################
# 
# Treat mice with abx and/or PPI and challenge with C difficile and examine how taxa change
# 
# input:
#	data/mothur/ppi.taxonomy
#	data/process/ppi_metadata.txt
# data/mothur/ppi.opti_mcc.groups.rarefaction
# data/mothur/ppi.opti_mcc.0.03.subsample.shared
#
# output:
#	To Be Decided
#		Plots looking at specific taxa levels between treatment groups
#		
#
#####################
#Setup----
library(tidyverse)
library(broom)
library(cowplot)

# Define color palette:
color_scheme <- c("#d95f02", "#1b9e77", "#7570b3")
color_ppi <-  c("#7570b3") # Use for graphs looking at just the PPI group over time
color_cppi <- c("#1b9e77") # Use for graphs looking at just the Clindamycin + PPI group over time
color_c <- c("#d95f02") # Use for graphs looking at just the Clindamycin group over time
color_day <- c('#c7c5e0' ,"#7570b3", "#514e7d") # Tints and Shades of PPI color for graphs showing relative abundances over time

# Import metadata into data frame
metadata <- read.table('data/process/ppi_metadata.txt', header = T, sep = '\t', stringsAsFactors = F) %>% 
  filter(Group != "NA") #Exclude the mock community

  
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


#Look at rarefaction curves, import rarefaction data frame
rarefy <- read_tsv(file="data/mothur/ppi.opti_mcc.groups.rarefaction") %>% 
  select(-contains("lci-"), -contains("hci-")) %>% 
  gather(-numsampled, key=sample, value=sobs) %>% #not including numsampled, gather samples intow two columns
  mutate(sample=str_replace_all(sample, patter="0.03-", replacement="")) %>% #drop prefix 0.03
  drop_na()

#Join metadata to rarefy data frame
metadata_rarefy <- inner_join(metadata, rarefy, by = c("shared_names" = "sample")) %>% 
  filter(Group != "NA") %>% 
  mutate(Group = as.factor(Group)) # make sure Group is treated as a factor


#Plotting----
#Plot rarefaction curves, the more parallel the curves to the x axis, the more confident you can be in the results
ggplot(metadata_rarefy, aes(x=numsampled, y=sobs, group=shared_names, color=Group))+
  scale_colour_manual(values=color_scheme) +
  geom_vline(xintercept=3000, color="gray", size=2) +
  geom_line()+
  coord_cartesian(xlim=c(0, 10000), ylim=c(0,200))+
  labs(subtitle="Vertical line indicates the number of sequences that samples were rarefied to",
       x="Number of sequences sampled per mouse", 
       y="Number of OTUs per mouse") +
  theme_classic()

# Import otu_data for PPI experiment samples
otu_data <- read_tsv("data/mothur/ppi.opti_mcc.0.03.subsample.shared", col_types=cols(Group=col_character())) %>% 
  select(-label, -numOtus) %>% 
  rename(shared_names=Group) %>% 
  gather(-shared_names, key="otu", value="count") %>% 
  mutate(rel_abund=count/3000) #Use 3000, because this is the subsampling parameter chosen.

# Merge otu_data to taxonomy data frame
agg_taxa_data <- inner_join(otu_data, taxonomy) 

# Summarize relative abundance data for phylum level
agg_phylum_data <- agg_taxa_data %>% 
  group_by(shared_names, phylum) %>% 
  summarize(agg_rel_abund=sum(rel_abund)) %>% 
# Merge relative abundance data to phyla data
  inner_join(., metadata, by = "shared_names") %>% 
  ungroup() 

# Use top_n function to determine top 5 phyla in our samples
top_phyla <- agg_phylum_data %>% 
  group_by(phylum) %>% 
  summarize(median=median(agg_rel_abund)) %>% 
  arrange((desc(median))) %>% 
  top_n(5, median) %>% 
  pull(phylum)

# Plot of top phylum data with log-scaled y-axis: 
agg_phylum_data %>% 
  filter(phylum %in% top_phyla) %>% 
  mutate(agg_rel_abund = agg_rel_abund + 1/6000) %>% 
  ggplot(aes(x= reorder(phylum, agg_rel_abund), y=agg_rel_abund, color=Group))+
  scale_colour_manual(values=color_scheme) +
  geom_hline(yintercept=1/3000, color="gray")+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(shape=19, size=1, alpha=0.4, position=position_jitterdodge(dodge.width=0.7, jitter.width=0.2)) +
  labs(title=NULL, 
       x=NULL,
       y="Relative abundance (%)")+
  scale_y_log10(breaks=c(1e-4, 1e-3, 1e-2, 1e-1, 1), labels=c(1e-2, 1e-1, 1, 10, 100))+
  coord_flip()+
  theme_classic()

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

# Use top_n function to determine top 10 genera in our samples
top_genus <- agg_genus_data %>% 
  group_by(genus) %>% 
  summarize(median=median(agg_rel_abund)) %>% 
  arrange((desc(median))) %>% 
  top_n(10, median) %>% 
  pull(genus)

# Plot of top genus data with log-scaled y-axis: 
agg_genus_data %>% 
  filter(genus %in% top_genus) %>% 
  mutate(agg_rel_abund = agg_rel_abund + 1/6000) %>% 
  ggplot(aes(x= reorder(genus, agg_rel_abund), y=agg_rel_abund, color=Group)) +
  scale_colour_manual(values=color_scheme) +
  geom_hline(yintercept=1/3000, color="gray")+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(shape=19, size=1, alpha=0.4, position=position_jitterdodge(dodge.width=0.7, jitter.width=0.1)) +
  labs(title=NULL, 
       x=NULL,
       y="Relative abundance (%)")+
  scale_y_log10(breaks=c(1e-4, 1e-3, 1e-2, 1e-1, 1), labels=c(1e-2, 1e-1, 1, 10, 100))+
  coord_flip()+
  theme_classic()

# Plot of top genus over time
agg_genus_data %>% 
  filter(genus %in% top_genus) %>% 
  mutate(agg_rel_abund = agg_rel_abund + 1/6000) %>%
  select(Group, day, agg_rel_abund, genus) %>% 
  ggplot(aes(x=day, y=agg_rel_abund, group=Group, color=Group))+
  scale_colour_manual(values=color_scheme) +
  facet_wrap("genus")+
  geom_point()+
  geom_line()+
  geom_hline(yintercept=1/3000, color="gray")+
  theme_classic()

# Plot of genus over time for PPI group
agg_genus_data %>% 
  filter(genus %in% top_genus) %>% 
  mutate(agg_rel_abund = agg_rel_abund + 1/6000) %>%
  select(Group, day, agg_rel_abund, genus) %>% 
  filter(Group == "PPI") %>% 
  ggplot(aes(x=day, y=agg_rel_abund, group=Group, color=Group))+
  scale_colour_manual(values=color_ppi) +
  facet_wrap("genus")+
  geom_point()+
  geom_line()+
  geom_hline(yintercept=1/3000, color="gray")+
  theme_classic()

#Hypothesis testing----
#Kruskal_wallis test for phylum differences across treatment groups with Benjamini-Hochburg correction 
phylum_tests <- agg_phylum_data %>% 
  group_by(phylum) %>% 
  do(tidy(kruskal.test(agg_rel_abund~factor(Group), data=.))) %>% ungroup() %>% 
  mutate(p.value.adj=p.adjust(p.value, method="BH")) %>% 
  arrange(p.value.adj)

#List the significant phyla based on treatment groups after Benjamini-Hochburg correction
sig_phyla <- phylum_tests %>% 
  filter(p.value.adj <= 0.05) %>% 
  pull(phylum)

#Pairwise wilcox test to determine which comparisons are significantly different
wilcox_phylum_tests <- agg_phylum_data %>% 
  group_by(phylum) %>% 
  filter(phylum %in% sig_phyla) %>% 
  select(phylum, Group, agg_rel_abund) 
wilcox_phylum_tests <- wilcox_phylum_tests %>% 
  do(tidy(pairwise.wilcox.test(g=wilcox_phylum_tests[["Group"]], x=wilcox_phylum_tests[["agg_rel_abund"]], p.adjust.method = "BH")))
# Warning messages: cannot compute exact p value with ties

#Graph the significant phyla based on treatment groups after Benjamini-Hochburg correction
group_phyla <- agg_phylum_data %>% 
  filter(phylum %in% sig_phyla) %>% 
  mutate(phylum=factor(phylum, levels=sig_phyla)) %>% 
  mutate(agg_rel_abund = agg_rel_abund + 1/6000) %>% 
  ggplot(aes(x= reorder(phylum, agg_rel_abund), y=agg_rel_abund, color=Group))+
  scale_colour_manual(values=color_scheme) +
  geom_hline(yintercept=1/3000, color="gray")+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(shape=19, size=1, alpha=0.4, position=position_jitterdodge(dodge.width=0.7, jitter.width=0.2)) +
  labs(title="Phyla significantly associated with treatment group", 
       x=NULL,
       y="Relative abundance (%)")+
  scale_y_log10(breaks=c(1e-4, 1e-3, 1e-2, 1e-1, 1), labels=c(1e-2, 1e-1, 1, 10, 100))+
  theme_classic()

#Kruskal_wallis test for family differences across treatment groups with Benjamini-Hochburg correction 
family_tests <- agg_family_data %>% 
  group_by(family) %>% 
  do(tidy(kruskal.test(agg_rel_abund~factor(Group), data=.))) %>% ungroup() %>% 
  mutate(p.value.adj=p.adjust(p.value, method="BH")) %>% 
  arrange(p.value.adj)

#List the significant family based on treatment groups after Benjamini-Hochburg correction
sig_family <- family_tests %>% 
  filter(p.value.adj <= 0.05) %>% 
  pull(family)

#Pairwise wilcox test to determine which comparisons are significantly different
wilcox_family_tests <- agg_family_data %>% 
  group_by(family) %>% 
  filter(family %in% sig_family) %>% 
  select(family, Group, agg_rel_abund) 
wilcox_family_tests <- wilcox_family_tests %>% 
  do(tidy(pairwise.wilcox.test(g=wilcox_family_tests[["Group"]], x=wilcox_family_tests[["agg_rel_abund"]], p.adjust.method = "BH")))

#Figure 2D----
#Graph the significant family based on treatment groups after Benjamini-Hochburg correction
group_family <- agg_family_data %>% 
  filter(family %in% sig_family) %>% 
  mutate(family=factor(family, levels=sig_family)) %>% 
  mutate(agg_rel_abund = agg_rel_abund + 1/6000) %>% 
  ggplot(aes(x= reorder(family, agg_rel_abund), y=agg_rel_abund, color=Group, alpha = day))+
  scale_colour_manual(values=color_scheme) +
  geom_hline(yintercept=1/3000, color="gray")+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(shape=19, size=1, alpha=0.4, position=position_jitterdodge(dodge.width=0.7, jitter.width=0.2)) +
  labs(title="Families significantly associated with treatment group", 
       x=NULL,
       y="Relative abundance (%)")+
  scale_y_log10(breaks=c(1e-4, 1e-3, 1e-2, 1e-1, 1), labels=c(1e-2, 1e-1, 1, 10, 100))+
  coord_flip()+
  theme_classic()+
  theme(plot.title=element_text(hjust=.5))+
  theme(legend.title=element_blank(), legend.position = c(0.85, 0.15)) #Get rid of legend title & move legend position
save_plot("results/figures/family_assoc_w_treatment.png", group_family, base_aspect_ratio = 1.5) #Use save_plot instead of ggsave because it's more compatible with cowplot

# Figure 1C----
#Graph the families associated with human PPI use according to the literature [@Imhann2017]
ppi_family <- c("Enterococcaceae", "Lactobacillaceae", "Micrococcaceae", "Staphylococcaceae", "Streptococcaceae", "Ruminococcaceae")
ppi_family_plot <- agg_family_data %>% 
  filter(family %in% ppi_family) %>% 
  mutate(family=factor(family, levels=ppi_family)) %>% 
  mutate(agg_rel_abund = agg_rel_abund + 1/6000) %>% 
  ggplot(aes(x= reorder(family, agg_rel_abund), y=agg_rel_abund, color=Group, alpha = day))+
  scale_colour_manual(values=color_scheme) +
  geom_hline(yintercept=1/3000, color="gray")+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(shape=19, size=1, alpha=0.4, position=position_jitterdodge(dodge.width=0.7, jitter.width=0.2)) +
  labs(title="Families previously associated with human PPI use", 
       x=NULL,
       y="Relative abundance (%)")+
  scale_y_log10(breaks=c(1e-4, 1e-3, 1e-2, 1e-1, 1), labels=c(1e-2, 1e-1, 1, 10, 100))+
  coord_flip()+
  theme_classic()+
  theme(plot.title=element_text(hjust=0.5))+
  theme(plot.title=element_text(hjust=0.5))+
  theme(legend.title=element_blank(), legend.position = c(0.85, 0.15)) #Ger rid of legend title & move legend position
save_plot("results/figures/families_prev_assoc_w_PPIs.png", ppi_family_plot, base_aspect_ratio = 1.5)

#Kruskal_wallis test for genus differences across treatment groups with Benjamini-Hochburg correction 
genus_tests <- agg_genus_data %>% 
  group_by(genus) %>% 
  do(tidy(kruskal.test(agg_rel_abund~factor(Group), data=.))) %>% ungroup() %>% 
  mutate(p.value.adj=p.adjust(p.value, method="BH")) %>% 
  arrange(p.value.adj)

#List the significant genera based on treatment groups after Benjamini-Hochburg correction
sig_genus <- genus_tests %>% 
  filter(p.value.adj <= 0.05) %>% 
  pull(genus)

#Pairwise wilcox test to determine which comparisons are significantly different
wilcox_genus_tests <- agg_genus_data %>% 
  group_by(genus) %>% 
  filter(genus %in% sig_genus) %>% 
  select(genus, Group, agg_rel_abund) 
wilcox_genus_tests <- wilcox_genus_tests %>%
  do(tidy(pairwise.wilcox.test(g=wilcox_genus_tests[["Group"]], x=wilcox_genus_tests[["agg_rel_abund"]], p.adjust.method = "BH")))

# Figure 2C----
#Graph the significant genera based on treatment groups after Benjamini-Hochburg correction
group_genera <- agg_genus_data %>% 
  filter(genus %in% sig_genus) %>% 
  mutate(genus=factor(genus, levels=sig_genus)) %>% 
  mutate(agg_rel_abund = agg_rel_abund + 1/6000) %>% 
  ggplot(aes(x= reorder(genus, agg_rel_abund), y=agg_rel_abund, color=Group))+
  scale_colour_manual(values=color_scheme) +
  geom_hline(yintercept=1/3000, color="gray")+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(shape=19, size=1, alpha=0.4, position=position_jitterdodge(dodge.width=0.7, jitter.width=0.2)) +
  labs(title="Genera significantly associated with treatment group", 
       x=NULL,
       y="Relative abundance (%)")+
  scale_y_log10(breaks=c(1e-4, 1e-3, 1e-2, 1e-1, 1), labels=c(1e-2, 1e-1, 1, 10, 100))+
  coord_flip()+
  theme_classic()+
  theme(axis.text.y = element_text(face = "italic"))+ #Have the genera show up as italics
  theme(plot.title=element_text(hjust=0.5))+
  theme(legend.title=element_blank(), legend.position = c(0.85, 0.2)) #Ger rid of legend title & move legend position
save_plot("results/figures/genera_assoc_w_treatment.png", group_genera, base_aspect_ratio = 1.6)
   
#Kruskal_wallis test with Benjamini-Hochburg correction for phylum differences across time in the PPI treatment group 
ppi_phylum_tests <- agg_phylum_data %>% 
  group_by(phylum) %>% 
  filter(Group == "PPI", day %in% c("-7", "0", "9")) %>%
  do(tidy(kruskal.test(agg_rel_abund~factor(day), data=.))) %>% ungroup() %>% 
  mutate(p.value.adj=p.adjust(p.value, method="BH")) %>% 
  arrange(p.value.adj)

#List the significant phyla over time in the PPI treatment group after Benjamini-Hochburg correction
ppi_sig_phyla <- ppi_phylum_tests %>% 
  filter(p.value.adj <= 0.05) %>% 
  pull(phylum)

#Graph the significant phyla over time in the PPI treatment group after Benjamini-Hochburg correction
ppi_phyla_time <- agg_phylum_data %>% 
  filter(phylum %in% ppi_sig_phyla) %>% 
  filter(Group == "PPI", day %in% c("-7", "0", "9")) %>%
  mutate(phylum=factor(phylum, levels=ppi_sig_phyla)) %>% 
  mutate(agg_rel_abund = agg_rel_abund + 1/6000) %>% 
  ggplot(aes(x= reorder(phylum, agg_rel_abund), y=agg_rel_abund, color=Group))+
  scale_colour_manual(values=color_ppi) +
  geom_hline(yintercept=1/3000, color="gray")+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(shape=19, size=1, alpha=0.4, position=position_jitterdodge(dodge.width=0.7, jitter.width=0.2)) +
  labs(title="Phyla significantly different after 16 days of PPI treatment", 
       x=NULL,
       y="Relative abundance (%)")+
  scale_y_log10(breaks=c(1e-4, 1e-3, 1e-2, 1e-1, 1), labels=c(1e-2, 1e-1, 1, 10, 100))+
  theme_classic()
#No significant phyla so no graph.

#Kruskal_wallis test for family differences across time in the PPI group with Benjamini-Hochburg correction 
ppi_family_tests <- agg_family_data %>% 
  group_by(family) %>% 
  filter(Group == "PPI", day %in% c("-7", "0", "9")) %>%
  do(tidy(kruskal.test(agg_rel_abund~factor(day), data=.))) %>% ungroup() %>% 
  mutate(p.value.adj=p.adjust(p.value, method="BH")) %>% 
  arrange(p.value.adj)

#List the significant family across time in the PPI group after Benjamini-Hochburg correction
ppi_sig_family <- ppi_family_tests %>% 
  filter(p.value.adj <= 0.05) %>% 
  pull(family)

#No significant families after Benjamini-Hochburg correction. Select the families with the lowest corrected p values
top_7_ppi <- top_n(ppi_family_tests, -7, p.value.adj) %>%  #negative to pull rows with the lowest values
  pull(family)

# Figure 1D----
#Graph the 7 families with the lowest Benjamini-Hochburg corrected P-values across time in the PPI group
ppi_family_time <- agg_family_data %>% 
  filter(family %in% top_7_ppi) %>% 
  filter(Group == "PPI", day %in% c("-7", "0", "9")) %>%
  mutate(family=factor(family, levels=top_7_ppi)) %>%
  mutate(agg_rel_abund = agg_rel_abund + 1/6000) %>% 
  ggplot(aes(x= reorder(family, agg_rel_abund), y=agg_rel_abund, color = factor(day, ordered = TRUE, levels =c("-7", "0", "9"))))+ #Turn day into a factor to use set of discrete colors
  scale_colour_manual(name="Day",
                      values=color_day,
                      breaks=c("-7", "0", "9"),
                      labels=c("-7", "0", "9"))+
  geom_hline(yintercept=1/3000, color="gray")+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(shape=19, size=1, position=position_jitterdodge(dodge.width=0.7, jitter.width=0.2), show.legend = NA) +
  labs(title="7 families with lowest adjusted P values (0.45 < P < 0.65) after comparing \n across 3 timepoints", #Break up title into multiple lines with \n
       x=NULL,
       y="Relative abundance (%)")+
  scale_y_log10(breaks=c(1e-4, 1e-3, 1e-2, 1e-1, 1), labels=c(1e-2, 1e-1, 1, 10, 100))+
  coord_flip()+
  theme_classic()+
  theme(legend.position = c(0.9, 0.2))+ #Move legend position
  theme(plot.title=element_text(hjust=.5))
save_plot("results/figures/ppi_family_time.png", ppi_family_time, base_aspect_ratio = 2)

#Kruskal_wallis test for genus differences across time in the PPI group with Benjamini-Hochburg correction 
ppi_genus_tests <- agg_genus_data %>% 
  group_by(genus) %>% 
  filter(Group == "PPI", day %in% c("-7", "0", "9")) %>%
  do(tidy(kruskal.test(agg_rel_abund~factor(day), data=.))) %>% ungroup() %>% 
  mutate(p.value.adj=p.adjust(p.value, method="BH")) %>% 
  arrange(p.value.adj)

#List the significant genera across time in the PPI group after Benjamini-Hochburg correction
ppi_sig_genus <- ppi_genus_tests %>% 
  filter(p.value.adj <= 0.05) %>% 
  pull(genus)

#No significant genera after Benjamini-Hochburg correction. Select the genera with the lowest corrected p values
top_13_ppi <- top_n(ppi_genus_tests, -13, p.value.adj) %>%  #negative to pull rows with the lowest values
  pull(genus)

#Figure S1C----
#Graph the 13 genera with the lowest Benjamini-Hochburg corrected P-values across time in the PPI group
ppi_genus_time <- agg_genus_data %>% 
  filter(genus %in% top_13_ppi) %>% 
  filter(Group == "PPI", day %in% c("-7", "0", "9")) %>%
  mutate(genus=factor(genus, levels=top_13_ppi)) %>% 
  mutate(agg_rel_abund = agg_rel_abund + 1/6000) %>% 
  ggplot(aes(x= reorder(genus, agg_rel_abund), y=agg_rel_abund, color = factor(day, ordered = TRUE, levels =c("-7", "0", "9"))))+
  scale_colour_manual(name="Day",
                      values=color_day,
                      breaks=c("-7", "0", "9"),
                      labels=c("-7", "0", "9"))+
  geom_hline(yintercept=1/3000, color="gray")+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(shape=19, size=1, position=position_jitterdodge(dodge.width=0.7, jitter.width=0.2), show.legend = NA) +
  labs(title="13 genera with lowest adjusted P values (P = 0.568) after comparing across 3 timepoints", 
       x=NULL,
       y="Relative abundance (%)")+
  scale_y_log10(breaks=c(1e-4, 1e-3, 1e-2, 1e-1, 1), labels=c(1e-2, 1e-1, 1, 10, 100))+
  coord_flip()+
  theme_classic()+
  theme(plot.title=element_text(hjust=.5))+ #Adjust plot title to fit graph
  theme(axis.text.y = element_text(face = "italic"))+ #Have the genera show up as italics
  theme(legend.position = c(0.9, 0.2)) #Move legend position
save_plot("results/figures/ppi_genera_time.png", ppi_genus_time, base_aspect_ratio = 2.4)

#Graph Porphyromonadaceae Unclassified over time
porpy_uncl_time <- agg_genus_data %>% 
  filter(genus == "Porphyromonadaceae Unclassified") %>% 
  mutate(agg_rel_abund = agg_rel_abund + 1/6000) %>%
  select(Group, day, agg_rel_abund, genus) %>% 
  ggplot(aes(x=day, y=agg_rel_abund, group=Group, color=Group))+
  scale_colour_manual(values=color_scheme) +
  scale_y_log10(breaks=c(1e-4, 1e-3, 1e-2, 1e-1, 1), labels=c(1e-2, 1e-1, 1, 10, 100))+
  geom_point()+
  geom_line()+
  geom_hline(yintercept=1/3000, color="gray")+
  labs(title="Porphyromonadaceae Unclassified") +
  theme_classic()

#Graph Lachnospiraceae Unclassified over time
lacno_uncl_time <- agg_genus_data %>% 
  filter(genus == "Lachnospiraceae Unclassified") %>% 
  mutate(agg_rel_abund = agg_rel_abund + 1/6000) %>%
  select(Group, day, agg_rel_abund, genus) %>% 
  ggplot(aes(x=day, y=agg_rel_abund, group=Group, color=Group))+
  scale_colour_manual(values=color_scheme) +
  scale_y_log10(breaks=c(1e-4, 1e-3, 1e-2, 1e-1, 1), labels=c(1e-2, 1e-1, 1, 10, 100))+
  geom_point()+
  geom_line()+
  geom_hline(yintercept=1/3000, color="gray")+
  labs(title="Lachnospiraceae Unclassified") +
  theme_classic()

ruminoc_genus_time <- agg_genus_data %>% 
  filter(genus == "Ruminococcus") %>% 
  mutate(agg_rel_abund = agg_rel_abund + 1/6000) %>%
  select(Group, day, agg_rel_abund, genus) %>% 
##  filter(Group == "PPI") %>%
  ggplot(aes(x=day, y=agg_rel_abund, group=Group, color=Group))+
  scale_colour_manual(values=color_scheme) +
  scale_y_log10(breaks=c(1e-4, 1e-3, 1e-2, 1e-1, 1), labels=c(1e-2, 1e-1, 1, 10, 100))+
  geom_point()+
  geom_line()+
  geom_hline(yintercept=1/3000, color="gray")+
  theme_classic()

enteroc_genus_time <- agg_genus_data %>% 
  filter(genus == "Enterococcus") %>% 
  mutate(agg_rel_abund = agg_rel_abund + 1/6000) %>%
  select(Group, day, agg_rel_abund, genus) %>% 
##  filter(Group == "PPI") %>%
  ggplot(aes(x=day, y=agg_rel_abund, group=Group, color=Group))+
  scale_colour_manual(values=color_scheme) +
  scale_y_log10(breaks=c(1e-4, 1e-3, 1e-2, 1e-1, 1), labels=c(1e-2, 1e-1, 1, 10, 100))+
  geom_point()+
  geom_line()+
  geom_line(position=position_dodge(width=0.2))+ #Show lines from all groups
  geom_hline(yintercept=1/3000, color="gray")+
  theme_classic()

strep_genus_time <- agg_genus_data %>% 
  filter(genus == "Streptococcus") %>% 
  mutate(agg_rel_abund = agg_rel_abund + 1/6000) %>%
  select(Group, day, agg_rel_abund, genus) %>% 
  filter(Group == "PPI") %>%
  ggplot(aes(x=day, y=agg_rel_abund, group=Group, color=Group))+
  scale_colour_manual(values=color_ppi) +
  scale_y_log10(breaks=c(1e-4, 1e-3, 1e-2, 1e-1, 1), labels=c(1e-2, 1e-1, 1, 10, 100))+
  geom_point()+
  geom_line()+
  geom_hline(yintercept=1/3000, color="gray")+
  theme_classic()

#Plots of families of bacteria that are associated with PPI use in humans
enteroc_family_time <- agg_family_data %>% 
  filter(family == "Enterococcaceae") %>% 
  mutate(agg_rel_abund = agg_rel_abund + 1/6000) %>%
  select(Group, day, agg_rel_abund, family) %>% 
  filter(Group == "PPI") %>%
  ggplot(aes(x=day, y=agg_rel_abund, group=Group, color=Group))+
  scale_colour_manual(values=color_ppi) +
  scale_y_log10(breaks=c(1e-4, 1e-3, 1e-2, 1e-1, 1), labels=c(1e-2, 1e-1, 1, 10, 100))+
  geom_point()+
  geom_line()+
  geom_hline(yintercept=1/3000, color="gray")+
  theme_classic()

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
# Figure S1A----
#Plot of Lactobacillaceae family relative abundance over time with mean relative abundance summary line and individual mouse relative abundance points
lacto_family_time <- ggplot(NULL)+
  geom_point(lacto_family_mice, mapping = aes(x=day, y=agg_rel_abund, color=Group, alpha = .2), show.legend = FALSE)+
  geom_line(lacto_family_mean, mapping = aes(x=day, y=mean, color=Group))+
  scale_colour_manual(values=color_scheme) +
  geom_hline(yintercept=1/3000, color="gray")+
  labs(title="Lactobacillaceae",
       x="Day",
       y="Relative abundance (%)") +
  scale_y_log10(breaks=c(1e-4, 1e-3, 1e-2, 1e-1, 1), labels=c(1e-2, 1e-1, 1, 10, 100))+
  theme_classic()+
  theme(plot.title=element_text(hjust=0.5))+
  theme(legend.title=element_blank(), legend.position = c(0.12, 0.27))+
  save_plot("results/figures/lactobacillaceae_time.png", lacto_family_time, base_aspect_ratio = 2.4)

 mi_family_time <- agg_family_data %>% 
  filter(family == "Micrococcaceae") %>% 
  mutate(agg_rel_abund = agg_rel_abund + 1/6000) %>%
  select(Group, day, agg_rel_abund, family) %>% 
  filter(Group == "PPI") %>%
  ggplot(aes(x=day, y=agg_rel_abund, group=Group, color=Group))+
  scale_colour_manual(values=color_ppi) +
  scale_y_log10(breaks=c(1e-4, 1e-3, 1e-2, 1e-1, 1), labels=c(1e-2, 1e-1, 1, 10, 100))+
  geom_point()+
  geom_line()+
  geom_hline(yintercept=1/3000, color="gray")+
  theme_classic()

staph_family_time <- agg_family_data %>% 
  filter(family == "Staphylococcaeae") %>% 
  mutate(agg_rel_abund = agg_rel_abund + 1/6000) %>%
  select(Group, day, agg_rel_abund, family) %>% 
  filter(Group == "PPI") %>%
  ggplot(aes(x=day, y=agg_rel_abund, group=Group, color=Group))+
  scale_colour_manual(values=color_ppi) +
  scale_y_log10(breaks=c(1e-4, 1e-3, 1e-2, 1e-1, 1), labels=c(1e-2, 1e-1, 1, 10, 100))+
  geom_point()+
  geom_line()+
  geom_hline(yintercept=1/3000, color="gray")+
  theme_classic()

strep_family_time <- agg_family_data %>% 
  filter(family == "Streptococcaceae") %>% 
  mutate(agg_rel_abund = agg_rel_abund + 1/6000) %>%
  select(Group, day, agg_rel_abund, family) %>% 
  filter(Group == "PPI") %>%
  ggplot(aes(x=day, y=agg_rel_abund, group=Group, color=Group))+
  scale_colour_manual(values=color_ppi) +
  scale_y_log10(breaks=c(1e-4, 1e-3, 1e-2, 1e-1, 1), labels=c(1e-2, 1e-1, 1, 10, 100))+
  geom_point()+
  geom_line()+
  geom_hline(yintercept=1/3000, color="gray")+
  theme_classic()

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
# Figure S1B----
#Plot of Ruminococcaceae family relative abundance over time with mean relative abundance summary line and individual mouse relative abundance points
rumino_family_time <- ggplot(NULL)+
  geom_point(rumino_family_mice, mapping = aes(x=day, y=agg_rel_abund, color=Group, alpha = .2), show.legend = FALSE)+
  geom_line(rumino_family_mean, mapping = aes(x=day, y=mean, color=Group))+
  scale_colour_manual(values=color_scheme) +
  geom_hline(yintercept=1/3000, color="gray")+
  labs(title="Ruminococcaceae",
       x="Day",
       y="Relative abundance (%)") +
  scale_y_log10(breaks=c(1e-4, 1e-3, 1e-2, 1e-1, 1), labels=c(1e-2, 1e-1, 1, 10, 100))+
  theme_classic()+
  theme(plot.title=element_text(hjust=0.5))+
  theme(legend.title=element_blank(), legend.position = c(0.9, 0.2))+
  save_plot("results/figures/ruminococcaceae_time.png", rumino_family_time, base_aspect_ratio = 2.4)

l_family_time <- agg_family_data %>% filter(family == "Lachnospiraceae") %>% 
  mutate(agg_rel_abund = agg_rel_abund + 1/6000) %>%
  select(Group, day, agg_rel_abund, family) %>% 
  filter(Group == "PPI") %>%
  ggplot(aes(x=day, y=agg_rel_abund, group=Group, color=Group))+
  scale_colour_manual(values=color_ppi) +
  scale_y_log10(breaks=c(1e-4, 1e-3, 1e-2, 1e-1, 1), labels=c(1e-2, 1e-1, 1, 10, 100))+
  geom_point()+
  geom_line()+
  geom_hline(yintercept=1/3000, color="gray")+
  theme_classic()

