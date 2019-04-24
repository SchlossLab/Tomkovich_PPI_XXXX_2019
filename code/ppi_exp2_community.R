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
library(cowplot)

shared_df <- read.table('data/mothur/ppi.opti_mcc.0.03.subsample.shared', sep = '\t', header = T, 
	stringsAsFactors = F, row.names = 'Group') # set column with sample name to match
metadata <- read.table('data/process/ppi_metadata.txt', header = T, sep = '\t', stringsAsFactors = F) %>% 
	# convert Group letter labels to more descriptive labels
	mutate(Group = case_when(Group == 'O+' ~ 'PPI', 
		Group == 'C+' ~ 'Clindamycin',
		Group == 'CO+' ~ 'Clindamycin + PPI'),
		day = day - 7) # C difficile challenge on day 7

mouse_label_df <- metadata %>% 
  select(Mouse.ID, Group) %>% 
  unique %>% 
  group_by(Group) %>% 
  arrange(Mouse.ID) %>% 
  mutate(mouse_position = 1:length(Mouse.ID)) %>% 
  ungroup %>% 
  select(Mouse.ID, mouse_position)

meta_df <- left_join(metadata, mouse_label_df, by = 'Mouse.ID') %>% rename(treatment = Group)

## need to get taxonomy file
#taxonomy_function <- 'code/sum_otu_by_taxa.R'
#source(taxonomy_function)
#taxonomy_df <- read.table(taxonomy_file, sep = '\t', header = T, stringsAsFactors = F)

shared_df <- select(shared_df, -label, -numOtus, )

sample_summed_counts <- apply(shared_df, 1, sum)
rel_abund <- data.frame(100 * shared_df/sample_summed_counts) %>% 
  mutate(shared_names = rownames(shared_df))

pre_post_end_df <- meta_df %>% 
  right_join(select(rel_abund, shared_names), by = "shared_names") %>% 
  group_by(Mouse.ID) %>% 
  mutate(endpoint = max(day), initial = min(day),
    sample_timepoint = case_when(day == initial ~ 'initial',
      day == 0 ~ 'day_0',
      day == endpoint ~ 'endpoint',
      T ~ NA_character_)) %>%
  ungroup %>% 
  filter(!is.na(sample_timepoint)) %>% 
  select(-endpoint, -initial)


colonization_df <- meta_df %>% 
  # to get level of colonization post C. difficile challenge, look of CFU on day 1 
  # but there are some cases day 1 might be deceptive so extending to day 2 
  # to see those that had delayed colonization
  filter(day > 0) %>% 
  group_by(Mouse.ID) %>% 
  mutate(cfu = ifelse(is.na(cfu), 0, cfu)) %>% 
  summarize(colonization = max(cfu))

top_otus <- function(input_dataframe, top_n){
	output_dataframe <- input_dataframe %>% 
	    gather(OTU, abundance, contains('Otu0'))
	    
   top_otus <- output_dataframe %>% 
      group_by(OTU) %>% 
      summarise(median_abundance = mean(abundance, na.rm = T)) %>% 
      top_n(top_n, median_abundance) %>% 
      select(OTU) %>% 
      mutate(top_otus = OTU)

    output_dataframe <- output_dataframe %>%
      full_join(top_otus, by = 'OTU') %>%  
      mutate(OTU = ifelse(is.na(top_otus), 'Other', OTU)) %>% 
      group_by(shared_names, OTU) %>% 
      summarise(abundance = sum(abundance)) %>% 
      mutate(dataset = paste0('top_', top_n, '_by_OTU')) %>% 
      ungroup 
}


pre_post_end_otu_df <- pre_post_end_df %>% 
  split(.$treatment) %>% 
  map_dfr(function(df) left_join(df, 
        top_otus(input_dataframe = filter(rel_abund, shared_names %in% df$shared_names), top_n = 11),
      by = "shared_names"))

all_days_otu_df <- meta_df %>% 
  split(.$treatment) %>% 
  map_dfr(function(df) left_join(df, 
        top_otus(input_dataframe = filter(rel_abund, shared_names %in% df$shared_names), top_n = 11),
      by = "shared_names")) 
	
barwidth <- 1
df <- pre_post_end_otu_df
pre_postabx_end_plot <- df %>% 
  mutate(sample_timepoint = factor(sample_timepoint, c('initial', 'day_0', 'endpoint')),
    taxa = factor(OTU, c('Other', unique(df$OTU)[unique(df$OTU)!='Other']))) %>% 
  ggplot(aes(x = Mouse.ID, y = abundance, fill = taxa)) + 
    geom_bar(stat="identity", position='stack', width = 1, color = "black", size = 0.1) + 
    facet_grid(treatment ~ sample_timepoint) + 
    theme_bw() + labs(x = NULL) + 
    theme(legend.position = 'top', legend.title=element_blank(),
      axis.text.x = element_blank(), axis.ticks.x = element_blank(),
      legend.text=element_text(size=6)) + 
    scale_fill_brewer(palette = 'Paired')

cdiff_plot <- pre_post_end_df %>% 
  filter(sample_timepoint == 'endpoint') %>% 
  left_join(colonization_df, by = 'Mouse.ID') %>% 
  rename(CFU = colonization) %>% 
  ggplot(aes(x = Mouse.ID, y = CFU)) + 
    geom_bar(stat="identity", position='stack', width = barwidth) + 
    facet_grid(treatment ~ .) + 
    scale_y_log10() + 
    theme_bw() + labs(x = NULL, y = 'CFU (log10)', title = '\n\n\n\n\n\nC. difficile Colonization') + 
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

plot_header <- ggdraw() +
 draw_label('Community overview at timepoint before antibiotic treatment, looking at the effect of Clindamycin and PPI', 
 	x = 0.05, hjust = 0)

ggsave(paste0('exploratory/notebook/pre_post_end_plots_top_otus.jpg'),
  plot_grid(plot_header, plot_grid(pre_postabx_end_plot, cdiff_plot, rel_widths = c(2, 1)), 
    ncol = 1, rel_heights = c(0.1, 2)),
  width = 10, height = 10)

daily_plot <- all_days_otu_df %>% 
  mutate(taxa = factor(OTU, c('Other', unique(df$OTU)[unique(df$OTU)!='Other']))) %>% 
  ggplot(aes(x = Mouse.ID, y = abundance, fill = taxa)) + 
    geom_bar(stat="identity", position='stack', width = 1, color = "black", size = 0.1) + 
    facet_grid(treatment ~ day) + 
    theme_bw() + labs(x = NULL) + 
    theme(legend.position = 'top', legend.title=element_blank(),
      axis.text.x = element_blank(), axis.ticks.x = element_blank(),
      legend.text=element_text(size=6)) + 
    scale_fill_brewer(palette = 'Paired')

ggsave(paste0('exploratory/notebook/daily_community_plots_top_otus.jpg'),
  plot_grid(plot_header, plot_grid(daily_plot, cdiff_plot, rel_widths = c(10, 1)), 
    ncol = 1, rel_heights = c(0.1, 2)),
  width = 25, height = 10)


