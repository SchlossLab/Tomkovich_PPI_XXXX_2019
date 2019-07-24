#Load packages
library(tidyverse)
library(vegan)
library(cowplot)

# Define color palette:
color_scheme <- c("#d95f02", "#1b9e77", "#7570b3")

# Import metadata into data frame
metadata <- read.table('data/process/ppi_metadata.txt', header = T, sep = '\t', stringsAsFactors = F) %>% 
  filter(Group != "NA") #Exclude the mock community'

#Function to have y-axis in scientific notation
fancy_scientific <- function(l) {
  # turn in to character string in scientific notation
  l <- format(l, scientific = TRUE)
  # quote the part before the exponent to keep all the digits
  l <- gsub("^(.*)e", "'\\1'e", l)
  # turn the 'e+' into plotmath format
  l <- gsub("e", "%*%10^", l)
  # return this as an expression
  parse(text=l)
}

#Create data frame of CFU data by selecting group column & columns ending with CFU/g
ppi_cfu <- select(metadata, Group, ends_with("CFU.g")) %>% 
  gather(Day, CFU, -Group) %>% 
  mutate(Group=factor(Group, levels=c("Clindamycin", "Clindamycin + PPI", "PPI")),
         Day = as.numeric(str_replace(Day, "^[^.](\\d\\d?)[\\w\\.\\s\\/]+", "\\1")),
         Day = Day - 7, #Modify day notation so that C. difficile challenge day is represented as day 0.
         CFU=as.numeric(CFU) + 1) %>% 
  ungroup

ppi_cfu_summary <- ppi_cfu %>% 
  group_by(Group, Day) %>% filter(!is.na(CFU)) %>% 
  summarise(median=median(CFU, na.rm = TRUE),
            mean=mean(CFU, na.rm = TRUE),
            lci=quantile(CFU, na.rm = TRUE, probs = 0.25),
            uci=quantile(CFU, na.rm = TRUE, probs = 0.75)) %>% 
  ungroup

# Figure 2A----
#Plot of mean CFU line with dots representing CFU of individual mice
title <-c(expression(paste(italic("C. difficile"), " colonization over time"))) #Expression variable for the title so that bacteria name will be in italics
cfu_time <- ggplot(NULL) + 
  geom_point(ppi_cfu, mapping = aes(x= Day, y = CFU, color=Group, fill=Group), alpha = .04, size = 2, shape = 20, show.legend = FALSE, position = position_dodge(width = 0.6))+
  geom_line(ppi_cfu_summary, mapping = aes(x=Day, y=mean, color=Group))+
  scale_colour_manual(values=color_scheme) +
  geom_linerange(show.legend=FALSE)+
  geom_hline(yintercept = 100, linetype=2) +
  labs(x='Days Post-Infection', y='CFU/g Feces')+
  geom_text(x = 28, y = 102, color="black", label="LOD")+
  scale_y_log10(labels=fancy_scientific, breaks = c(10, 100, 10^3, 10^4, 10^5, 10^6, 10^7, 10^8, 10^9))+
  theme_classic()+
  theme(legend.position = c(.9, .8)) +
  labs(title=title) +
  theme(plot.title=element_text(hjust=0.5))
save_plot("results/figures/ppi_cfu.png", cfu_time, base_aspect_ratio = 2)
