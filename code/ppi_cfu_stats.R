library(tidyverse)
library(broom)

#Load in metadata file and tidy CFU/g data from each timepoint
ppi_cfu <- select(metadata, Group, ends_with("CFU.g")) %>% 
  gather(Day, CFU, -Group) %>% 
  mutate(Group=factor(Group, levels=c("PPI", "Clindamycin", "Clindamycin + PPI")),
         Day = as.numeric(str_replace(Day, "^[^.](\\d\\d?)[\\w\\.\\s\\/]+", "\\1")),
         Day = Day - 7, #Modify day notation so that C. difficile challenge day is represented as day 0.
         CFU=as.numeric(CFU) + 1) %>% 
  ungroup

#Analysis
get_mean_cfu <- function(x){
      x %>% 
        group_by(Group) %>% 
        summarize(mean=mean(CFU)) %>% 
        spread(key=Group, value=mean)
}

select(ppi_cfu, Group, Day, CFU) %>% 
  group_by(Day) %>% 
  nest() %>% 
  mutate(model=map(data, ~kruskal.test(x=.x$CFU, g=as.factor(.x$Group)) %>% tidy())) %>% 
  mutate(mean= map(data, get_mean_cfu)) %>% 
  unnest(model, mean) %>% 
  mutate(model=map(data,
                   ~wilcox.test(x=.x$CFU, g=as.factor(.x$Group), p.adjust.method="BH") %>% 
                     tidy() %>% 
                     mutate(compare=paste(group1, group2, sep="-")) %>% 
##                     select(-group1, -group2) %>% 
##                     spread(key=compare, value=p.value)
                   )
         ) %>% 
##  unnest(model) %>% 
##  select(-data, -parameter, -statistic) %>% 
  write_tsv("data/process/ppi_cfu_stats.tsv")

cfu_tests <- ppi_cfu %>% 
  group_by(Day) %>% 
  do(tidy(kruskal.test(CFU~Group, data=.))) %>% 
  ungroup() %>% 
  mutate(p.value.adj=p.adjust(p.value, method="BH")) %>% 
  arrange(p.value.adj)

for (i in ppi_cfu$Day){
print(pairwise.wilcox.test(ppi_cfu$CFU, ppi_cfu$Group, p.adjust.method = "BH"))
}
