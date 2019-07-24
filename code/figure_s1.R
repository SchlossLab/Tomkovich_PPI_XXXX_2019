library(tidyverse)
library(magick)
library(cowplot)

a <- ggdraw() + draw_image("results/figures/lactobacillaceae_time.png")
b <- ggdraw() + draw_image("results/figures/ruminococcaceae_time.png")
c <- ggdraw() + draw_image("results/figures/ppi_genera_time.png")

plot_grid(a, b, c, labels = "AUTO", label_size = 12)+
  ggsave("results/figures/figure_s1.pdf", width=6.5, height=3)+
  ggsave("submission/figure_s1.pdf", width=6.5, height=3.5)



