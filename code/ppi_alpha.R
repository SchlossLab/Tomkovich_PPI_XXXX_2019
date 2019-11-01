#Setup----
library(tidyverse)
library(broom)
library(cowplot)

# Import alpha diversity measures for each sample into data frame.
alpha <- read_tsv(file="data/mothur/ppi.opti_mcc.groups.ave-std.summary",
                  col_types=cols(group = col_character())) %>%
  filter(method=='ave') %>%
  select(group, sobs, shannon, invsimpson, coverage)

# Import metadata into data frame
metadata <- read.table('data/process/ppi_metadata.txt', header = T, sep = '\t', stringsAsFactors = F) %>% 
  filter(Group != "NA") %>% #Exclude the mock community
  mutate(Group=factor(Group, levels=c("Clindamycin", "Clind. + Omep.", "Omeprazole"))) 


# Combine metadata with alpha data frame
meta_alpha <- inner_join(metadata, alpha, by=c('shared_names'='group')) %>% 
  select(Group, sobs, shannon, day) 

# Define color palette:
color_scheme <- c("#d95f02", "#1b9e77",  "#7570b3")
color_ppi <-  c("#7570b3") # Use for graphs looking at just the PPI group over time

#Test if Shannon diversity and richness is significantly different over D-7, 0, and 9 of the experiment within the omeprazole group----
ppi_day_alpha <- meta_alpha %>% 
  filter( Group == "Omeprazole") %>% # select just PPI group 
  filter(day == -7 | day == 0 | day ==9) %>% # select -7, 0, and -9 timepoints (beginning, middle and end of the experiment)
  mutate(day = as.factor(day))# turn day into factor so we can do post hoc comparisons

# Test if PPI Shannon data are normally distributed using Shapiro-Wilk normality test
shapiro.test(ppi_day_alpha$shannon) # Test for Shannon diversity
shapiro.test(ppi_day_alpha$sobs)
# P > 0.05, means the data are normally distributed and we should do parametric test (aov/ANOVA)

#Hypothesis testing across 3 timepoints in the PPI group by analysis of variance (ANOVA)
#shannon
ppi_day_shannon_aov <- aov(shannon~day, data = ppi_day_alpha) 
summary(ppi_day_shannon_aov)
#Since P = 0.0467, do Tukey's honest significance difference test to determine which comparisons between group are significant.
TukeyHSD(ppi_day_shannon_aov)
#No comparison's significant after adjusting p values. Closest was 9 to -7.

#richness
ppi_day_sobs_aov <- aov(sobs~day, data = ppi_day_alpha) 
summary(ppi_day_sobs_aov)
#Since P = 0.265, do not do Tukey's honest significance difference test to determine which comparisons between group are significant.

#Test if Shannon diversity and richness is significantly different over D-7, 0, and 9 of the experiment within the clindamycin group----
clind_day_alpha <- meta_alpha %>% 
  filter( Group == "Clindamycin") %>% # select just clindamycin group 
  filter(day == -7 | day == 0 | day ==9) %>% # select -7, 0, and -9 timepoints (beginning, middle and end of the experiment)
  mutate(day = as.factor(day))# turn day into factor so we can do post hoc comparisons

# Test if Clindamycin Shannon data are normally distributed using Shapiro-Wilk normality test
shapiro.test(clind_day_alpha$shannon) # Test for Shannon diversity
shapiro.test(clind_day_alpha$sobs)
# P < 0.05, means the data are not normally distributed and we should do nonparametric test (Kruskal-Wallis)

#Hypothesis testing across 3 timepoints in the Clindamycin group
#shannon
clind_shannon_kruskal <- kruskal.test(shannon~day, data=clind_day_alpha) # P-value = 0.02962
pairwise.wilcox.test(g=clind_day_alpha[["day"]], x=clind_day_alpha[["shannon"]], p.adjust.method="BH")
# No pairwise comparisons between timepoints were significant

#richness
clind_sobs_kruskal <- kruskal.test(sobs~day, data=clind_day_alpha) # P-value = 0.02491
pairwise.wilcox.test(g=clind_day_alpha[["day"]], x=clind_day_alpha[["sobs"]], p.adjust.method="BH")
# -7 versus 0 and 9 versus 0 were significant (P = 0.043)

#Test if Shannon diversity and richness is significantly different over D-7, 0, and 9 of the experiment within the clindamycin + PPI group----
clindPPI_day_alpha <- meta_alpha %>% 
  filter( Group == "Clind. + Omep.") %>% # select clindamycin + PPI
  filter(day == -7 | day == 0 | day ==9) %>% # select -7, 0, and -9 timepoints (beginning, middle and end of the experiment)
  mutate(day = as.factor(day))# turn day into factor so we can do post hoc comparisons

# Test if PPI Shannon data are normally distributed using Shapiro-Wilk normality test
shapiro.test(clindPPI_day_alpha$shannon) # Test for Shannon diversity
# P < 0.05, means the data are not normally distributed and we should do nonparametric test (Kruskal-Wallis)
shapiro.test(clindPPI_day_alpha$sobs)
# P > 0.05 means data are normally distrubed and we should do parametric test (aov, analysis of variance)

#Hypothesis testing across 3 timepoints in the Clindamycin + PPI group
#shannon
clindPPI_shannon_kruskal <- kruskal.test(shannon~day, data=clindPPI_day_alpha)
# P = 0.1655, not significant so don't do post hoc comparisons

#richness
clindPPI_day_sobs_aov <- aov(sobs~day, data = clindPPI_day_alpha) 
summary(clindPPI_day_sobs_aov)
# P = 0.0274, do post hoc comparisons
TukeyHSD(clindPPI_day_sobs_aov)
# 9 versus 0 comparison was significant P = 0.022

# Data frame so we can show what Shannon diversity & richness look like over 3 timepoints for all 3 groups
groups_day_alpha <- meta_alpha %>% 
  filter(day == -7 | day == 0 | day ==9) %>% # select -7, 0, and -9 timepoints (beginning, middle and end of the experiment)
  mutate(day = as.factor(day))# turn day into factor

# Figure S2A----
#Boxplots of Shannon diversity for all samples at thetimepoints being compared within PPI group:
shannon <- groups_day_alpha %>% 
  ggplot(aes(x= day, y=shannon, colour=Group)) +
  scale_colour_manual(name=NULL, 
                      values=color_scheme, 
                      breaks=c("Clindamycin", "Clind. + Omep.", "Omeprazole"),
                      labels=c("Clindamycin", "Clind. + Omep.", "Omeprazole")) +
  geom_boxplot(outlier.shape = NA, size = 1.2)+
  geom_jitter(shape=19, size=2, alpha=0.6, position=position_jitterdodge(dodge.width=0.7, jitter.width=0.1), show.legend = FALSE) +
  labs(title=NULL, 
       x="Day",
       y="Shannon Diversity Index")+
  ylim(0, 3.5)+
  theme_classic()+
  theme(legend.position = c(0.8, 0.2))+
  theme(text = element_text(size = 16))  # Change font size for entire plot
save_plot("results/figures/groups_shannon.png", shannon) #Use save_plot instead of ggsave because it works better with cowplot

# Figure S2B----
#Boxplots of sobs (richness) for all samples at the 3 timepoints being compared within PPI group:
richness <- groups_day_alpha %>% 
  ggplot(aes(x= day, y=sobs, colour=Group)) +
  scale_colour_manual(name=NULL, 
                      values=color_scheme, 
                      breaks=c("Clindamycin", "Clind. + Omep.", "Omeprazole"),
                      labels=c("Clindamycin", "Clind. + Omep.", "Omeprazole")) +
  geom_boxplot(outlier.shape = NA, size = 1.2)+
  geom_jitter(shape=19, size=2, alpha=0.6, position=position_jitterdodge(dodge.width=0.7, jitter.width=0.1), show.legend = FALSE) +
  labs(title=NULL, 
       x="Day",
       y="Number of Observed OTUs")+
  ylim(0, 100)+
  theme_classic()+
  theme(legend.position = c(0.8, 0.2))+
  theme(text = element_text(size = 16))  # Change font size for entire plot
save_plot("results/figures/groups_richness.png", richness) #Use save_plot instead of ggsave because it works better with cowplot



