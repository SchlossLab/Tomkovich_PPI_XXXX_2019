library(tidyverse)
library(magick)
library(cowplot)

a <- ggdraw() + draw_image("results/figures/ppi_cfu.png")
b <- ggdraw() + draw_image("results/figures/after_abx_C.diff.png")
c <- ggdraw() + draw_image("results/figures/genera_assoc_w_treatment.png")

plot_grid(a, b, c, labels = "AUTO", label_size = 12)+
  ggsave("results/figures/figure_2.pdf", width=6.875, height=4.2)+
  ggsave("submission/figure_2.pdf", width=6.875, height=4.2)

