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
  filter(Group != "NA") #Exclude the mock community

# Combine metadata with alpha data frame
meta_alpha <- inner_join(metadata, alpha, by=c('shared_names'='group')) %>% 
  select(Group, sobs, shannon, day) %>% 
  mutate(Group = as.factor(Group)) # make sure Group is treated as a factor

# Define color palette:
color_scheme <- c("#d95f02", "#1b9e77", "#7570b3")
color_ppi <-  c("#7570b3") # Use for graphs looking at just the PPI group over time

# Boxplots of shannon index (measure of community diversity) across treatment groups
a <- meta_alpha %>% 
  ggplot(aes(x= reorder(Group, shannon), y=shannon, colour=Group)) +
  scale_colour_manual(values=color_scheme) +
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(shape=19, size=1, alpha=0.7, position=position_jitterdodge(dodge.width=0.7, jitter.width=0.1)) +
  labs(title=NULL, 
       x=NULL,
       y="Shannon")+
  ylim(0, 3.5)+
  theme_classic()+
  theme(legend.title=element_blank(), legend.position="bottom")

# Test if Shannon and sobs (richness) data are normally distributed using Shapiro-Wilk normality test
shapiro.test(meta_alpha$shannon) # Test for Shannon diversity
# P < 0.05, means the data are not normally distributed and we should do non parametric test (Kruskal-Wallis)

#Hypothesis testing between groups, all timepoints. Since data is not normally distributed use the Kruskal-Wallis test
group_shannon_kruskal <- kruskal.test(shannon~Group, data=meta_alpha)
# P = 6.29e-06, so do post hoc test to determin which comparisons are significant
pairwise.wilcox.test(g=meta_alpha[["Group"]], x=meta_alpha[["shannon"]], p.adjust.method="BH")
## PPI vs clindamycin and PPI vs clindamycin+PPI are significantly different (P = 2.3e-05 and P = 2.3e-05)

# Boxplots of sobs (measure of community richness) across treatment groups
d <- meta_alpha %>% 
  ggplot(aes(x= reorder(Group, sobs), y=sobs, color=Group)) +
  scale_colour_manual(values=color_scheme) +
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(shape=19, size=1, alpha=0.7, position=position_jitterdodge(dodge.width=0.7, jitter.width=0.1)) +
  labs(title=NULL, 
       x=NULL,
       y="richness")+
  ylim(0, 100)+
  theme_classic()+
  theme(legend.position="none")

# Test if sobs (richness) is normally distributed
shapiro.test(meta_alpha$sobs) 
# P > 0.05 means the data are normally distrubuted and we can do a parametric anova test.

#Hypothesis testing between groups, all timepoints analysis of variance (ANOVA)
group_sobs_aov <- aov(sobs~Group, data = meta_alpha)
summary(group_sobs_aov)
#Since P = 5.72e-05, do Tukey's honest significance difference test to determine which comparisons between group are significant.
TukeyHSD(group_sobs_aov)
# PPI vs clindamycin and PPI vs clindamycin+PPI are significantly different (P = 0.003 and P = 0.0000647)

# Linegraph of shannnon index across time within the PPI group
b <- meta_alpha %>% 
  filter(Group == "PPI") %>% 
  select(Group, day, shannon) %>% 
  ggplot(aes(x=day, y=shannon, group=Group, color=Group))+
  scale_colour_manual(values=color_ppi) +
  geom_point(show.legend = FALSE)+
  geom_line(show.legend = FALSE)+
  ylim(0, 3.5)+
  labs(y="Shannon") +
  theme_classic()

# Linegraph of sobs (observed community richness) over time within the PPI group
e <- meta_alpha %>% 
  filter(Group == "PPI") %>% 
  select(Group, day, sobs) %>% 
  ggplot(aes(x=day, y=sobs, group=Group, color=Group))+
  scale_colour_manual(values=color_ppi) +
  geom_point(show.legend = FALSE)+
  geom_line(show.legend = FALSE)+
  ylim(0, 100)+
  labs(y="richness") +
  theme_classic()

#Test if Shannon diversity and richness is significantly different over D-7, 0, and 9 of the experiment
ppi_day_alpha <- meta_alpha %>% 
  filter( Group == "PPI") %>% # select just PPI group 
  filter(day == -7 | day == 0 | day ==9) %>% # select -7, 0, and -9 timepoints (beginning, middle and end of the experiment)
  mutate(day = as.factor(day))# turn day into factor so we can do post hoc comparisons

# Test if PPI Shannon data are normally distributed using Shapiro-Wilk normality test
shapiro.test(ppi_time_alpha$shannon) # Test for Shannon diversity
shapiro.test(ppi_time_alpha$sobs)
# P > 0.05, means the data are normally distributed and we should do parametric test (aov/ANOVA)

#Hypothesis testing between groups, all timepoints analysis of variance (ANOVA)
#shannon
ppi_day_shannon_aov <- aov(shannon~day, data = ppi_day_alpha) 
summary(ppi_day_shannon_aov)
#Since P = 0.0436, do Tukey's honest significance difference test to determine which comparisons between group are significant.
TukeyHSD(ppi_day_shannon_aov)
#No comparison's significant after adjusting p values. Closest was 9 to -7.

#richness
ppi_day_sobs_aov <- aov(sobs~day, data = ppi_day_alpha) 
summary(ppi_day_sobs_aov)
#Since P = 0.0603, do not do Tukeys honest significance difference test to determine which comparisons between group are significant.

#Boxplots of Shannon diversity for timepoints being compared within PPI group:
c <- ppi_day_alpha %>% 
  ggplot(aes(x= day, y=shannon, colour=Group)) +
  scale_colour_manual(values=color_ppi) +
  geom_boxplot(outlier.shape = NA, show.legend = FALSE)+
  geom_jitter(shape=19, size=1, alpha=0.7, position=position_jitterdodge(dodge.width=0.7, jitter.width=0.1), show.legend = FALSE) +
  labs(title=NULL, 
       x=NULL,
       y="Shannon")+
  ylim(0, 3.5)+
  theme_classic()

#Boxplots of sobs (richness) for timepoints being compared within PPI group:
f <- ppi_day_alpha %>% 
  ggplot(aes(x= day, y=sobs, colour=Group)) +
  scale_colour_manual(values=color_ppi) +
  geom_boxplot(outlier.shape = NA, show.legend = FALSE)+
  geom_jitter(shape=19, size=1, alpha=0.7, position=position_jitterdodge(dodge.width=0.7, jitter.width=0.1), show.legend = FALSE) +
  labs(title=NULL, 
       x=NULL,
       y="richness")+
  ylim(0, 100)+
  theme_classic()

shannon_sobs <- plot_grid(a, b, c, d, e, f, labels = "AUTO") +
  ggsave("results/figures/shannon_sobs.pdf", width=11, height=6)
