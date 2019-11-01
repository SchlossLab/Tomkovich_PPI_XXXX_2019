library(tidyverse)
library(magick)
library(cowplot)

a <- ggdraw() + draw_image("results/figures/groups_shannon.png")
b <- ggdraw() + draw_image("results/figures/groups_richness.png")

plot_grid(a, b, labels = "AUTO", label_size = 12)+
  ggsave("results/figures/figure_s2.pdf", width=6.875, height=3)+
  ggsave("submission/figure_s2.pdf", width=6.875, height=3)
