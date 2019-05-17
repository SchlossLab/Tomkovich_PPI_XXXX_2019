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

# Import metadata into data frame and clean
metadata <- read.table('data/process/ppi_metadata.txt', header = T, sep = '\t', stringsAsFactors = F) %>% 
  filter(Group != "NA") %>% #Exclude the mock community
  mutate(Group = case_when(Group == 'O+' ~ 'PPI', #make group names readable by humans
                           Group == 'C+' ~ 'Clindamycin',
                           Group == 'CO+' ~ 'Clindamycin + PPI'))

# Import taxonomy into data frame and clean up taxonomy names
taxonomy <- read_tsv(file="data/mothur/ppi.taxonomy") %>%
  rename_all(tolower) %>% #remove uppercase from column names
# Split taxonomic information into separate columns for each taxonomic level  
  mutate(taxonomy=str_replace_all(taxonomy, c("\\(\\d*\\)" = "", #drop digits with parentheses around them
                                              ';$' = "", #removes semi-colon at end of line
                                              'Bacteria_unclassified' = 'Unclassified',
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
  filter(Group != "NA") 

#Plotting----
#Plot rarefaction curves, the more parallel the curves to the x axis, the more confident you can be in the results
ggplot(metadata_rarefy, aes(x=numsampled, y=sobs, group=shared_names, color=Group))+
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
  mutate(rel_abund=count/3000) 

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

# Plot of phylum data with log-scaled y-axis: 
agg_phylum_data %>% 
  filter(phylum %in% top_phyla) %>% 
  mutate(agg_rel_abund = agg_rel_abund + 1/6000) %>% 
  ggplot(aes(x= reorder(phylum, agg_rel_abund), y=agg_rel_abund, color=Group))+
  geom_hline(yintercept=1/3000, color="gray")+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(shape=19, size=1, alpha=0.7, position=position_jitterdodge(dodge.width=0.7, jitter.width=0.2)) +
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

# Use top_n function to determine top 10 genera in our samples
top_genus <- agg_genus_data %>% 
  group_by(genus) %>% 
  summarize(median=median(agg_rel_abund)) %>% 
  arrange((desc(median))) %>% 
  top_n(10, median) %>% 
  pull(genus)

# Plot of genus data with log-scaled y-axis: 
agg_genus_data %>% 
  filter(genus %in% top_genus) %>% 
  mutate(agg_rel_abund = agg_rel_abund + 1/6000) %>% 
  ggplot(aes(x= reorder(genus, agg_rel_abund), y=agg_rel_abund, color=Group)) +
  geom_hline(yintercept=1/3000, color="gray")+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(shape=19, size=1, alpha=0.7, position=position_jitterdodge(dodge.width=0.7, jitter.width=0.1)) +
  labs(title=NULL, 
       x=NULL,
       y="Relative abundance (%)")+
  scale_y_log10(breaks=c(1e-4, 1e-3, 1e-2, 1e-1, 1), labels=c(1e-2, 1e-1, 1, 10, 100))+
  coord_flip()+
  theme_classic()

# Plot of genus over time
agg_genus_data %>% 
  filter(genus %in% top_genus) %>% 
  mutate(agg_rel_abund = agg_rel_abund + 1/6000) %>%
  select(Group, day, agg_rel_abund, genus) %>% 
  ggplot(aes(x=day, y=agg_rel_abund, group=Group, color=Group))+
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

#Graph the significant phyla based on treatment groups after Benjamini-Hochburg correction
agg_phylum_data %>% 
  filter(phylum %in% sig_phyla) %>% 
  mutate(phylum=factor(phylum, levels=sig_phyla)) %>% 
  mutate(agg_rel_abund = agg_rel_abund + 1/6000) %>% 
  ggplot(aes(x= reorder(phylum, agg_rel_abund), y=agg_rel_abund, color=Group))+
  geom_hline(yintercept=1/3000, color="gray")+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(shape=19, size=1, alpha=0.7, position=position_jitterdodge(dodge.width=0.7, jitter.width=0.2)) +
  labs(title="Phyla significantly associated with treatment group", 
       x=NULL,
       y="Relative abundance (%)")+
  scale_y_log10(breaks=c(1e-4, 1e-3, 1e-2, 1e-1, 1), labels=c(1e-2, 1e-1, 1, 10, 100))+
  theme_classic()

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

#Graph the significant genera based on treatment groups after Benjamini-Hochburg correction
agg_genus_data %>% 
  filter(genus %in% sig_genus) %>% 
  mutate(genus=factor(genus, levels=sig_genus)) %>% 
  mutate(agg_rel_abund = agg_rel_abund + 1/6000) %>% 
  ggplot(aes(x= reorder(genus, agg_rel_abund), y=agg_rel_abund, color=Group))+
  geom_hline(yintercept=1/3000, color="gray")+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(shape=19, size=1, alpha=0.7, position=position_jitterdodge(dodge.width=0.7, jitter.width=0.2)) +
  labs(title="Genera significantly associated with treatment group", 
       x=NULL,
       y="Relative abundance (%)")+
  scale_y_log10(breaks=c(1e-4, 1e-3, 1e-2, 1e-1, 1), labels=c(1e-2, 1e-1, 1, 10, 100))+
  coord_flip()+
  theme_classic()
   
#Kruskal_wallis test with Benjamini-Hochburg correction for phylum differences across time in the PPI treatment group 
ppi_phylum_tests <- agg_phylum_data %>% 
  group_by(phylum) %>% 
  filter(Group == "PPI", day %in% c("0", "7", "16")) %>%
  do(tidy(kruskal.test(agg_rel_abund~factor(day), data=.))) %>% ungroup() %>% 
  mutate(p.value.adj=p.adjust(p.value, method="BH")) %>% 
  arrange(p.value.adj)

#List the significant phyla over time in the PPI treatment group after Benjamini-Hochburg correction
ppi_sig_phyla <- ppi_phylum_tests %>% 
  filter(p.value.adj <= 0.05) %>% 
  pull(phylum)

#Graph the significant phyla over time in the PPI treatment group after Benjamini-Hochburg correction
agg_phylum_data %>% 
  filter(phylum %in% ppi_sig_phyla) %>% 
  filter(Group == "PPI", day %in% c("0", "7", "16")) %>%
  mutate(phylum=factor(phylum, levels=ppi_sig_phyla)) %>% 
  mutate(agg_rel_abund = agg_rel_abund + 1/6000) %>% 
  ggplot(aes(x= reorder(phylum, agg_rel_abund), y=agg_rel_abund, color=Group))+
  geom_hline(yintercept=1/3000, color="gray")+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(shape=19, size=1, alpha=0.7, position=position_jitterdodge(dodge.width=0.7, jitter.width=0.2)) +
  labs(title="Phyla significantly different after 16 days of PPI treatment", 
       x=NULL,
       y="Relative abundance (%)")+
  scale_y_log10(breaks=c(1e-4, 1e-3, 1e-2, 1e-1, 1), labels=c(1e-2, 1e-1, 1, 10, 100))+
  theme_classic()

#Kruskal_wallis test for genus differences across time in the PPI group with Benjamini-Hochburg correction 
ppi_genus_tests <- agg_genus_data %>% 
  group_by(genus) %>% 
  filter(Group == "PPI", day %in% c("0", "7", "16")) %>%
  do(tidy(kruskal.test(agg_rel_abund~factor(day), data=.))) %>% ungroup() %>% 
  mutate(p.value.adj=p.adjust(p.value, method="BH")) %>% 
  arrange(p.value.adj)

#List the significant genera across time in the PPI group after Benjamini-Hochburg correction
ppi_sig_genus <- ppi_genus_tests %>% 
  filter(p.value.adj <= 0.05) %>% 
  pull(genus)

#Graph genera across time in the PPI group
agg_genus_data %>% 
  filter(genus %in% ppi_genus_tests) %>% 
  filter(Group == "PPI", day %in% c("0", "7", "16")) %>%
  mutate(genus=factor(genus, levels=ppi_genus_tests)) %>% 
  mutate(agg_rel_abund = agg_rel_abund + 1/6000) %>% 
  ggplot(aes(x= reorder(genus, agg_rel_abund), y=agg_rel_abund, color=day))+
  geom_hline(yintercept=1/3000, color="gray")+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(shape=19, size=1, alpha=0.7, position=position_jitterdodge(dodge.width=0.7, jitter.width=0.2)) +
  labs(title="Genera levels over time in PPI treatment", 
       x=NULL,
       y="Relative abundance (%)")+
  scale_y_log10(breaks=c(1e-4, 1e-3, 1e-2, 1e-1, 1), labels=c(1e-2, 1e-1, 1, 10, 100))+
  theme_classic()

agg_genus_data %>% 
  filter(genus == "Ruminococcus") %>% 
  mutate(agg_rel_abund = agg_rel_abund + 1/6000) %>%
  select(Group, day, agg_rel_abund, genus) %>% 
  filter(Group == "PPI") %>%
  ggplot(aes(x=day, y=agg_rel_abund, group=Group, color=Group))+
  geom_point()+
  geom_line()+
  geom_hline(yintercept=1/3000, color="gray")+
  theme_classic()
