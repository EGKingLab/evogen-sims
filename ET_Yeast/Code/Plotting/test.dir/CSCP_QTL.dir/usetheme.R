usetheme <- library(cowplot)
theme_set(theme_cowplot(font_size = 12))
mytheme <- list(
  theme(axis.title = element_text(face = "bold"),
        axis.text = element_text(face = "italic"),
        plot.title = element_text(hjust = 0.5),
        legend.position = "none"))
