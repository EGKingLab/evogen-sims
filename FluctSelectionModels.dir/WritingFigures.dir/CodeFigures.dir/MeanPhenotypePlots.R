############ Libraries and theme function ########
#"^MeanPhenotypes\\d+_n1_H0\\.(1)(SD(4))?(Gen(30))?\\.csv$"
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
    theme_cowplot() +
    theme(legend.position = "none",
          axis.text = element_text(size = 15, face = "bold"),
          axis.line = element_line(size = 2),
          axis.title = element_text(size = 15, face = "bold"),
          strip.text = element_text(size = 15, face = "bold"),
          plot.title = element_text(hjust = 0.5)) #,legend.position = "none"
}

#################### Data Processing and plotting ###############

# The first part deals with constant and sinusoidal II
process_files <- function(dirpath, pattern, plot_type) {
  files <- list.files(dirpath, pattern, full.names = TRUE)

  dataframes <- list()
  for(file in files){
    replicate_id <- as.numeric(str_extract(file, "(?<=MeanPhenotypes)\\d+"))
    herit <- as.numeric(str_extract(file, "(?<=H)0\\.\\d+"))
    loci <- as.numeric(str_extract(file, "(?<=_n)\\d+"))
    sd <- as.numeric(str_extract(file, "(?<=SD)\\d+"))
    gen <- as.numeric(str_extract(file, "(?<=Gen)\\d+"))

    data <- read.csv(file, header = TRUE) %>%
      mutate(herit = herit,loci = loci, sd = sd, replicate = replicate_id,
             Optimum = ifelse(Optimum == max(Optimum), "H", "L"))


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
        filter(loci_gen == loci_geni)
      p <- locus_data %>%
        group_by(Generation, h2_sd) %>%  # Include h2_sd in the grouping
        summarise(mean_phenotype = mean(Phenotype, na.rm = TRUE),
                  se_phenotype = sd(Phenotype, na.rm = TRUE) / sqrt(n()), .groups = 'drop') %>%  # Calculate mean and SE
        ggplot(aes(Generation, mean_phenotype, ymin = mean_phenotype - se_phenotype, ymax = mean_phenotype + se_phenotype)) +
        geom_line(aes()) + #fill = Optimum
        geom_ribbon(alpha = 0.5) +  # Add a ribbon for SE
        facet_wrap(~h2_sd, ncol = 2, scales = "free_y")+
        ylab("Phenotype")+
        mytheme()
      #theme(legend.position = "none")
      plots[[loci_geni]] <- (p) #ggplotly
    }

    # plots for constant or sinusoidal II selections

  } else if (plot_type == "loci") {
    loci <- unique(combined_data$loci)
    for(locus in loci){
      locus_data <- combined_data %>%
        filter(loci == locus)
      p <- locus_data %>%
        group_by(Generation, h2_sd) %>%  # Include h2_sd in the grouping
        summarise(mean_phenotype = mean(Phenotype, na.rm = TRUE),
                  se_phenotype = sd(Phenotype, na.rm = TRUE) / sqrt(n()), .groups = 'drop') %>%  # Calculate mean and SE
        ggplot(aes(Generation, mean_phenotype, ymin = mean_phenotype - se_phenotype, ymax = mean_phenotype + se_phenotype)) +
        geom_line(aes()) + #fill = Optimum
        geom_ribbon(alpha = 0.5) +  # Add a ribbon for SE
        facet_wrap(~h2_sd, ncol = 2, scales = "free_y")+
        ylab("Phenotype")+
        mytheme()
      plots[[locus]] <- (p) #ggplotly
    }
  }

  return(list(combined_data = combined_data, plots = plots))
}

##````````````````````````````````````````

# library(ggplot2)
# library(dplyr)
# library(stringr)
# library(cowplot)
# 
# # Define a custom theme function
# mytheme <- function() {
#   theme_minimal() + 
#     theme(
#       legend.position = "none",
#       axis.text = element_text(size = 17, face = "bold"),
#       axis.line = element_line(size = 2),
#       axis.title = element_text(size = 17, face = "bold"),
#       strip.text = element_text(size = 17, face = "bold")
#     )
# }
# 
# # Function to process files and plot data
# process_files <- function(dirpath, pattern, plot_type) {
#   files <- list.files(dirpath, pattern, full.names = TRUE)
#   
#   dataframes <- list()
#   for (file in files) {
#     replicate_id <- as.numeric(str_extract(file, "(?<=MeanPhenotypes)\\d+"))
#     herit <- as.numeric(str_extract(file, "(?<=H)0\\.\\d+"))
#     loci <- as.numeric(str_extract(file, "(?<=_n)\\d+"))
#     sd <- as.numeric(str_extract(file, "(?<=SD)\\d+"))
#     gen <- as.numeric(str_extract(file, "(?<=Gen)\\d+"))
#     
#     data <- read.csv(file, header = TRUE) %>% 
#       mutate(herit = herit, loci = loci, sd = sd, replicate = replicate_id,
#              Optimum = ifelse(Phenotype == max(Phenotype, na.rm = TRUE), "H", "L"))
#     
#     if (!is.na(gen)) {
#       data$gen <- rep(gen, nrow(data))
#     }
#     
#     dataframes[[file]] <- data
#   }
#   
#   combined_data <- bind_rows(dataframes) %>% 
#     mutate(h2_sd = paste("h2 = ", herit, " ", "sd = ", sd, sep = ""),
#            loci_gen = paste("loci = ", loci, " ", "gen = ", gen, sep = ""))
#   
#   plots <- list()
#   if (plot_type == "loci_gen") {
#     loci_gens <- unique(combined_data$loci_gen)
#     for (loci_geni in loci_gens) {
#       locus_data <- combined_data %>% 
#         filter(loci_gen == loci_geni)
#       p <- locus_data %>% 
#         group_by(Generation, h2_sd) %>% 
#         summarise(mean_phenotype = mean(Phenotype, na.rm = TRUE),
#                   se_phenotype = sd(Phenotype, na.rm = TRUE) / sqrt(n()), .groups = 'drop') %>% 
#         ggplot(aes(Generation, mean_phenotype, ymin = mean_phenotype - se_phenotype, ymax = mean_phenotype + se_phenotype)) +
#         geom_line() +
#         geom_ribbon(alpha = 0.5) +
#         facet_wrap(~h2_sd, ncol = 2, scales = "free_y") +
#         ylab("Phenotype") +
#         mytheme()
#       plots[[loci_geni]] <- p
#     }
#   } else if (plot_type == "loci") {
#     loci <- unique(combined_data$loci)
#     for (locus in loci) {
#       locus_data <- combined_data %>% 
#         filter(loci == locus)
#       p <- locus_data %>% 
#         group_by(Generation, h2_sd) %>% 
#         summarise(mean_phenotype = mean(Phenotype, na.rm = TRUE),
#                   se_phenotype = sd(Phenotype, na.rm = TRUE) / sqrt(n()), .groups = 'drop') %>% 
#         ggplot(aes(Generation, mean_phenotype, ymin = mean_phenotype - se_phenotype, ymax = mean_phenotype + se_phenotype)) +
#         geom_line() +
#         geom_ribbon(alpha = 0.5) +
#         facet_wrap(~h2_sd, ncol = 2, scales = "free_y") +
#         ylab("Phenotype") +
#         mytheme()
#       plots[[locus]] <- p
#     }
#   }
#   
#   return(list(combined_data = combined_data, plots = plots))
# }
