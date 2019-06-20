library(tidyverse)
library(magick)
library(cowplot)

a <- ggdraw() + draw_image("results/figures/ppi_cfu.png")
b <- ggdraw() + draw_image("results/figures/after_abx_C.diff.png")
c <- ggdraw() + draw_image("results/figures/genera_assoc_w_treatment.png")

plot_grid(a, b, c, labels = "AUTO")+
  ggsave("results/figures/figure_2.pdf", width=6.5, height=6)
