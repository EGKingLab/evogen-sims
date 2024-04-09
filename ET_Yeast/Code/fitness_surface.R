library(tidyverse)
library(cowplot)
theme_set(theme_cowplot())
phenotypes <- seq(-8,8, length.out = 200)

new_opts <- c(-4,-1,0,1,4)

phenotypes <- crossing(phenotypes, new_opts)

scale <- 1/10

phenotypes$fitness <- exp(-(scale) * ((phenotypes$phenotypes - phenotypes$new_opts))^2)

phenotypes |>
  ggplot(aes(phenotypes, fitness, color = as.factor(new_opts))) +
  geom_point() +
  geom_vline(xintercept = seq(-8,8, by =1), color = "grey80")

