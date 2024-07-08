
############ Libraries and theme function ########

library(dplyr)
library(tidyr)
library(stringr)
library(plotly)
library(patchwork)
library(cowplot)
library(purrr)
library(doParallel)
theme_set(theme_cowplot())
mytheme <- function(){
    theme_set(theme_cowplot())+
    theme(axis.text = element_text(size = 12),
          axis.line = element_line(size = 2),
          axis.title = element_text(face = "bold"),
          legend.position = "none")
}

####################
# File: reusable_plotting_function.R

gen_plotting_function <- function(dirpath, pattern, ncol_plots = 2, label_size = 14, label_x = 0, label_y = 1) {
  
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
      select(Generation, Position, Frequency, Effect) %>% 
      mutate(Position = factor(Position),
             herit = rep(herit, n()),
             loci = rep(loci, n()),
             sd = rep(sd, n()),
             gen = rep(gen, n()),
             Pos_Effect = paste("Position = ", Position, " Effect = ", Effect),
             h2_sd = paste("h2 = ", herit, " sd = ", sd, " generations = ", gen, sep = ""))
  })
  
  #   # Select 30 random positions
  set.seed(12345)
  
  selected_positions <- combined_data %>%
    group_by(loci) %>%
    sample_n(size = min(30, n()), replace = F) %>%
    pull(Pos_Effect)
  

    # Filter the data for the selected positions
    combined_data <- combined_data %>%
      filter(Pos_Effect %in% selected_positions)
  
  # Plotting
  plots <- combined_data %>%
    group_by(loci) %>%
    nest() %>%
    mutate(plot = map(data, ~{
      ggplot(.x, aes(Generation, Frequency, group = Pos_Effect, color = Pos_Effect)) +
        geom_line(size = 0.1) +
        facet_wrap(~h2_sd, ncol = 4) +
        ylim(min = 0, max = 1) +
        mytheme()
    })) %>%
    pull(plot)
  
  # Combine all plots using plot_grid
  combined_plot <- cowplot::plot_grid(plotlist = plots, ncol = ncol_plots, align = 'v', 
                                      labels = "AUTO", label_size = label_size, 
                                      label_x = label_x, label_y = label_y)
  
  return(combined_plot)
}




#################################################################
###Function to select only 30 loci for better visualization #####
#################################################################

gen2_plotting_function <- function(dirpath, pattern, ncol_plots = 2, label_size = 14, label_x = 0.5, label_y = 1) {
  
  # Data Processing
  files <- list.files(dirpath, pattern, full.names = TRUE)
  
  combined_data <- map_dfr(files, ~{
    file <- .x
    herit <- as.numeric(str_extract(file, "(?<=H)0\\.\\d+"))
    loci <- as.numeric(str_extract(file, "(?<=_n)\\d+"))
    sd <- as.numeric(str_extract(file, "(?<=SD)\\d+"))
    gen <- as.numeric(str_extract(file, "(?<=Gen)\\d+"))
    
    read.csv(file, header = TRUE) %>% 
      select(Generation, Position, Frequency, Effect) %>% 
      mutate(Position = factor(Position),
             herit = rep(herit, n()),
             loci = rep(loci, n()),
             sd = rep(sd, n()),
             gen = rep(gen, n()),
             Pos_Effect = paste("Position = ", Position, " Effect = ", Effect),
             h2_sd = paste("h2 = ", herit, " sd = ", sd, " generations = ", gen, sep = ""))
  })
  
  # Ensure unique Pos_Effects are sampled correctly
  set.seed(12345)
  
  selected_positions <- combined_data %>%
    group_by(loci) %>%
    sample_n(size = min(30, n()), replace = F) %>%
    pull(Pos_Effect)
  
  
  # Filter the data for the selected positions
  combined_data <- combined_data %>%
    filter(Pos_Effect %in% selected_positions)
  
  # Plotting with ggplotly for interactivity
  plots <- combined_data_filtered %>%
    group_by(loci) %>%
    nest() %>%
    mutate(plot = map(data, ~{
      p <- ggplot(.x, aes(Generation, Frequency, group = Pos_Effect, color = Pos_Effect)) +
        geom_line(size = 0.1) +
        facet_wrap(~h2_sd, ncol = 4) +
        ylim(0, 1) +
        theme_minimal() +
        theme(legend.position = "none") +
        labs(x = "Generation", y = "Frequency")
      ggplotly(p)
    })) %>%
    pull(plot)
  
  # Since ggplotly plots are interactive, they can't be combined using cowplot.
  # Consider returning a list of plots or displaying them individually.
  
  return(plots)
}


# gen2_plotting_function <- 
#   function(dirpath, pattern, ncol_plots = 2, 
#            label_size = 14, label_x = 0, 
#            label_y = 1) {
#   
#   # Data Processing
#   files <- list.files(dirpath, pattern, full.names = TRUE)
#   
#   combined_data <- map_dfr(files, ~{
#     file <- .x
#     herit <- as.numeric(str_extract(file, "(?<=H)0\\.\\d+"))
#     loci <- as.numeric(str_extract(file, "(?<=_n)\\d+"))
#     sd <- as.numeric(str_extract(file, "(?<=SD)\\d+"))
#     gen <- as.numeric(str_extract(file, "(?<=Gen)\\d+"))
#     
#     read.csv(file, header = TRUE) %>% 
#       select(Generation, Position, Frequency, Effect) %>% 
#       mutate(Position = factor(Position),
#              herit = rep(herit, n()),
#              loci = rep(loci, n()),
#              sd = rep(sd, n()),
#              gen = rep(gen, n()),
#              Pos_Effect = paste("Position = ", Position, " Effect = ", Effect)
#       )
#   }) %>% 
#     mutate(h2_sd = paste("h2 = ", herit, " ", "sd = ", sd, " ", 
#                          "generations = ", gen, sep = ""))
#   
#   # Select 30 random positions
#   set.seed(12345)
#   selected_positions <- 
#     sample(combined_data$Pos_Effect, 
#            size = min(30, nrow(combined_data)), 
#            replace = T)
#   
#   # Filter the data for the selected positions
#   combined_data <- combined_data %>% 
#     filter(Pos_Effect %in% selected_positions)
#   
#   # Plotting
#   plots <- combined_data %>%
#     group_by(loci) %>%
#     nest() %>%
#     mutate(plot = map(data, ~{
#       ggplot(.x, aes(Generation, Frequency, 
#                      group = Pos_Effect, color = Pos_Effect)) +
#         geom_line(linewidth = 0.1) +
#         facet_wrap(~h2_sd, ncol = 4) +
#         ylim(min = 0, max = 1) +
#         theme_bw() +
#         theme(legend.position = "none")
#     })) %>%
#     pull((plot))
#   
#   # Combine all plots using plot_grid
#   combined_plot <- 
#     cowplot::plot_grid(plotlist = (plots), ncol = ncol_plots, align = 'v', 
#                                       labels = "AUTO", label_size = label_size, 
#                                       label_x = label_x, label_y = label_y)
#   
#   return(combined_plot)
# }
################################################################################