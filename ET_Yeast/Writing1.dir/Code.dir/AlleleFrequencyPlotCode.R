# # Load necessary libraries
# library(tidyverse)
# library(stringr)
# library(plotly)
# library(patchwork)
# library(cowplot)
# library(purrr)
# library(doParallel)
# 
# # Set theme
# theme_set(theme_cowplot())
# 
# # Define custom theme function
# mytheme <- function(){
#   theme(axis.title = element_text(face = "bold"),
#         axis.text = element_text(face = "italic"),
#         plot.title = element_text(hjust = 0.5),
#         legend.position = "none")
# }
# 
# process_files <- function(dirpath, pattern, plot_type) {
#   files <- list.files(dirpath, pattern, full.names = TRUE)
#   
#   cl <- makeCluster(detectCores() - 1) # leave one core free
#   registerDoParallel(cl)
#   
#   dataframes <- foreach(file = files, .combine = 'rbind', .packages = c("tidyverse", "stringr")) %dopar% {
#     herit <- as.numeric(str_extract(file, "(?<=H)0\\.\\d+"))
#     loci <- as.numeric(str_extract(file, "(?<=_n)\\d+"))
#     sd <- as.numeric(str_extract(file, "(?<=SD)\\d+"))
#     gen <- as.numeric(str_extract(file, "(?<=Gen)\\d+"))
#     
#     h2_sd <- paste("h2 = ", herit," ", "sd = ", sd, sep = "")
#     loci_gen <- paste("loci = ", loci," ", "gen = ", gen, sep = "")
#     
#     data <- read.csv(file, header = TRUE) %>%
#       select(Generation, Position, Frequency, Effect) %>%
#       group_by(Position) %>%
#       mutate(Position = factor(Position),
#              herit = herit,
#              loci = loci,
#              sd = sd,
#              initFreq = Frequency[Generation == 1],
#              h2_sd = h2_sd,
#              loci_gen = loci_gen,
#              postion_effect_init = paste("position = ", Position," ","Effect = ", round(Effect, 2)," ", "Initial Freq = ", round(initFreq, 2) , sep = ""))
#     
#     if(!is.na(gen)) {
#       data$gen <- rep(gen, nrow(data))
#     }
#     
#     data
#   }
#   
#   stopCluster(cl)
#   
#   combined_data <- bind_rows(dataframes)
#   
#   plots <- list()
#   if (plot_type == "loci_gen") {
#     loci_gens <- unique(combined_data$loci_gen)
#     for(loci_geni in loci_gens){
#       locus_data <- combined_data %>%
#         filter(loci_gen == loci_geni)
#       p <- locus_data %>%
#         ggplot(aes(Generation, Frequency, group = postion_effect_init, 
#                    color = postion_effect_init, 
#                    linetype = postion_effect_init))+
#         geom_line(linewidth = 0.5)+
#         facet_wrap(~h2_sd, ncol = 4)+
#         theme(legend.position = "none")+
#         ylim(min = 0, max = 1)+
#         theme_bw()
#       plots[[loci_geni]] <- ggplotly(p)
#     }
#   } else if (plot_type == "loci") {
#     loci <- unique(combined_data$loci)
#     for(locus in loci){
#       locus_data <- combined_data %>% 
#         filter(loci == locus)
#       p <- locus_data %>%
#         ggplot(aes(Generation, Frequency, group = postion_effect_init, 
#                    color = postion_effect_init, 
#                    linetype = postion_effect_init))+
#         geom_line(linewidth = 0.5)+
#         facet_wrap(~h2_sd, ncol = 4)+
#         ylim(min = 0, max = 1)+
#         theme(legend.position = "none")+
#         theme_bw()
#       plots[[locus]] <- ggplotly(p)
#     }
#   }
#   return(list(combined_data = combined_data, plots = plots))
# }

############## This is to help with packages in case theyre not there #####
install_and_load <- function(packages) {
  for (package in packages) {
    if (!require(package, character.only = TRUE)) {
      install.packages(package, dependencies = TRUE)
      library(package, character.only = TRUE)
    }
  }
}

# List of required packages
packages <- c("tidyverse", "stringr", "plotly", "patchwork", "cowplot", "purrr", "doParallel")

# Install and load the packages
#install_and_load(packages)

############ Libraries and theme function ########

#library(tidyverse)
library(dplyr)
library(tidyr)
library(dplyr)
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

#################### Data Processing and plotting ###############

# The first part deals with constant and sinusoidal II
process_files <- function(dirpath, pattern, plot_type) {
  files <- list.files(dirpath, pattern, full.names = TRUE)

  dataframes <- list()
  for(file in files){
    replicate_id <- as.numeric(str_extract(file, "(?<=genome)\\d+"))
    herit <- as.numeric(str_extract(file, "(?<=H)0\\.\\d+"))
    loci <- as.numeric(str_extract(file, "(?<=_n)\\d+"))
    sd <- as.numeric(str_extract(file, "(?<=SD)\\d+"))
    gen <- as.numeric(str_extract(file, "(?<=Gen)\\d+"))

    data <- read.csv(file, header = TRUE) %>%
      select(Generation, Position, Frequency, Effect) %>%
      group_by(Position) %>%
      mutate(Position = factor(Position),
             replicate = replicate_id,
             herit = herit,
             loci = loci,
             sd = sd,
             initFreq = Frequency[Generation == 1],
             h2_sd = paste("h2 = ", herit," ", "sd = ", sd, sep = ""),
             postion_effect_init = paste("position = ", Position," ",
                                         "Effect = ", round(Effect, 2)," ", 
                                         "Initial Freq = ", 
                                         round(initFreq, 2)," ",
                                         "repl = ", " ",replicate, sep = ""))


    # argument for linear and sinusoidal I selections

    if(!is.na(gen)) {
      data$gen <- rep(gen, nrow(data))
    }

    dataframes[[file]] <- data
  }
  combined_data <- bind_rows(dataframes) %>%
    mutate(h2_sd = paste("h2 = ", herit," ", "sd = ", sd, sep = ""),
           loci_gen = paste("loci = ", loci," ", "gen = ", gen, sep = ""))

  # plots for either linear or sinusoidal selections

  plots <- list()
  if (plot_type == "loci_gen") {
    loci_gens <- unique(combined_data$loci_gen)
    for(loci_geni in loci_gens){
      locus_data <- combined_data %>%
        filter(loci_gen == loci_geni)#,
               #Generation %% (gen+1) == 0,
               #Effect > 1.5)
      # Select 10 random positions
      #selected_positions <- sample(locus_data$Position,
    #size = min(10, nrow(locus_data)), replace = F)
      # Filter the data for the selected positions
      #locus_data <- locus_data %>%
        #filter(Position %in% selected_positions)
      
      # Select 10 random positions
      #set.seed(12345)
      #selected_positions <- sample(locus_data$Position, 
       #                            size = min(5, nrow(locus_data)), replace = F)
      
      # Filter the data for the selected positions
      #locus_data <- locus_data %>% 
       # filter(Position %in% selected_positions)
      
      
      p <- locus_data %>%
        ggplot(aes(Generation, Frequency, group = postion_effect_init,
                   color = postion_effect_init))+#ostion_effect_init
        geom_line(linewidth = 0.2)+
        facet_wrap(~h2_sd, ncol = 4)+
        ylim(min = 0, max = 1)+
        theme_bw()+
        theme(legend.position = "none")
      
      plots[[loci_geni]] <- ggplotly(p) #ggplotly
    }

    # plots for constant or sinusoidal II selections

  } else if (plot_type == "loci") {
    loci <- unique(combined_data$loci)
    for(locus in loci){
      locus_data <- combined_data %>%
        filter(loci == locus)#,
               #Generation %% (gen+1) == 0,
               #Effect > 1.5)

      # Select 10 random positions
      #selected_positions <- sample(locus_data$Position,
      #size = min(10, nrow(locus_data)), replace = F)
      # Filter the data for the selected positions
      #locus_data <- locus_data %>%
        #filter(Position %in% selected_positions)
      
      p <- locus_data %>%
        ggplot(aes(Generation, Frequency, group = postion_effect_init, #postion_effect_init 
                   color = postion_effect_init))+
        geom_line(linewidth = 0.2)+
        facet_wrap(~h2_sd, ncol = 4)+
        ylim(min = 0, max = 1)+
        theme_bw()#+
        #theme(legend.position = "none")
      plots[[locus]] <- ggplotly(p) #ggplotly
    }
  }

  return(list(combined_data = combined_data, plots = plots))
}
