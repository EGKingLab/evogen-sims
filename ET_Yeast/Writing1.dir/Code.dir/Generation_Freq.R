
############ Libraries and theme function ########

library(tidyverse)
library(stringr)
library(plotly)
library(patchwork)
library(cowplot)
library(purrr)
library(doParallel)
theme_set(theme_cowplot())
mytheme <- function(){
  theme(axis.title = element_text(face = "bold"),
        axis.text = element_text(face = "italic"),
        plot.title = element_text(hjust = 0.5),
        legend.position = "none")
}

####################
# File: reusable_plotting_function.R

gen_plotting_function <- function(dirpath, pattern, ncol_plots = 2, label_size = 14, label_x = 0, label_y = 1) {
  library(tidyverse)
  library(cowplot)
  
  # Data Processing
  files <- list.files(dirpath, pattern, full.names = TRUE)
  
  combined_data <- map_dfr(files, ~{
    file <- .x
    herit <- as.numeric(str_extract(file, "(?<=H)0\\.\\d+"))
    loci <- as.numeric(str_extract(file, "(?<=_n)\\d+"))
    sd <- as.numeric(str_extract(file, "(?<=SD)\\d+"))
    gen <- as.numeric(str_extract(file, "(?<=Gen)\\d+"))
    
    read.csv(file, header = TRUE) %>% 
      select(Generation, Position, Frequency) %>% 
      mutate(Position = factor(Position),
             herit = rep(herit, n()),
             loci = rep(loci, n()),
             sd = rep(sd, n()),
             gen = rep(gen, n())
      )
  }) %>% 
    mutate(h2_sd = paste("h2 = ", herit, " ", "sd = ", sd, " ", "generations = ", gen, sep = ""))
  
  # Plotting
  plots <- combined_data %>%
    group_by(loci) %>%
    nest() %>%
    mutate(plot = map(data, ~{
      ggplot(.x, aes(Generation, Frequency, group = Position, color = Position)) +
        geom_line(linewidth = 0.1) +
        facet_wrap(~h2_sd, ncol = 4) +
        ylim(min = 0, max = 1) +
        theme_bw() +
        theme(legend.position = "none")
    })) %>%
    pull(plot)
  
  # Combine all plots using plot_grid
  combined_plot <- cowplot::plot_grid(plotlist = plots, ncol = ncol_plots, align = 'v', 
                                      labels = "AUTO", label_size = label_size, 
                                      label_x = label_x, label_y = label_y)
  
  return(combined_plot)
}




# 
# dirpath <- "~/YeastProj.dir/evogen-sims/ET_Yeast/output.dir/Selection_Models/WF.dir/LinFS.dir/"
# pattern <- "^genome1_n(1|10|100|300)_H(0\\.1|0\\.8)SD(1|4)Gen(10|30)\\.csv$"
# files <- list.files(dirpath, pattern, full.names = TRUE)
# 
# combined_data <- map_dfr(files, ~{
#   file <- .x
#   herit <- as.numeric(str_extract(file, "(?<=H)0\\.\\d+"))
#   loci <- as.numeric(str_extract(file, "(?<=_n)\\d+"))
#   sd <- as.numeric(str_extract(file, "(?<=SD)\\d+"))
#   gen <- as.numeric(str_extract(file, "(?<=Gen)\\d+"))
#   
#   read.csv(file, header = TRUE) %>% 
#     select(Generation, Position, Frequency) %>% 
#     mutate(Position = factor(Position),
#            herit = rep(herit, n()),
#            loci = rep(loci, n()),
#            sd = rep(sd, n()),
#            gen = rep(gen, n())
#     )
# }) %>% 
#   mutate(h2_sd = paste("h2 = ", herit, " ", "sd = ", sd, " ", "generations = ", gen, sep = ""))
# 
# plots <- combined_data %>%
#   group_by(loci) %>%
#   nest() %>%
#   mutate(plot = map(data, ~{
#     ggplot(.x, aes(Generation, Frequency, group = Position, color = Position)) +
#       geom_line(linewidth = 0.1) +
#       facet_wrap(~h2_sd, ncol = 4) +
#       ylim(min = 0, max = 1) +
#       theme_bw() +
#       theme(legend.position = "none")
#   })) %>%
#   pull(plot)
# 
# cowplot::plot_grid(plotlist = plots, ncol = 2, align = 'v', 
#                    labels = "AUTO", label_size = 14, 
#                    label_x = 0, label_y = 1)
