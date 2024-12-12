###################################################
############## Box plots for phenotypes ##########
#################################################

############ Libraries and theme function ########

library(dplyr)
library(stringr)
library(plotly)
library(patchwork)
library(cowplot)
library(purrr)
library(doParallel)

mytheme <- function(){theme_classic() +
    theme(legend.position = "none",
          axis.text.x = element_text(size = 20, face = "bold"),
          axis.text.y = element_text(size = 20, face = "bold"),
          axis.line = element_line(linewidth = 5),
          axis.title.x = element_text(size = 40, face = "bold", margin = margin(t = 25)), # 
          axis.title.y = element_text(size = 40, face = "bold", margin = margin(r = 20)), # 
          strip.text = element_text(size = 40, face = "bold"),
          panel.spacing = unit(3, "lines"),
          panel.grid = element_blank())
}
#################### Data Processing and plotting ###############

# The first part deals with constant and sinusoidal II

process_files <- function(dirpath, pattern, plot_type) {
  files <- list.files(dirpath, pattern, full.names = TRUE)
  
  dataframes <- list()
  for(file in files){
    replicate_id <- as.numeric(str_extract(file, "(?<=MeanPhenotypes)\\d+"))
    H <- as.numeric(str_extract(file, "(?<=H)0\\.\\d+"))
    loci <- as.numeric(str_extract(file, "(?<=_n)\\d+"))
    SD <- as.numeric(str_extract(file, "(?<=SD)\\d+"))
    Gen <- as.numeric(str_extract(file, "(?<=Gen)\\d+"))
    
    data <- read.csv(file, header = TRUE)
    
    data$H <- rep(H, nrow(data))
    data$loci <- rep(loci, nrow(data))
    data$SD <- rep(SD, nrow(data))
    data$replicate <- rep(factor(replicate_id), nrow(data))
    
    # argument for instantaneous and sinusoidal I selections
    
    if(!is.na(Gen)) {
      data$Gen <- rep(Gen, nrow(data))
    }
    
    dataframes[[file]] <- data
  }
  combined_data <- bind_rows(dataframes) %>% 
    mutate(h2_sd = paste("H = ", H," ", "SD = ", SD, sep = ""),
           loci_gen = paste("loci = ", loci," ", "Gen = ", Gen, sep = ""))
  
  # plots for either linear or sinusoidal selections
  
  plots <- list()
  if (plot_type == "loci_gen") {
    loci_gens <- unique(combined_data$loci_gen)
    for(loci_geni in loci_gens){
      locus_data <- combined_data %>% 
        filter(loci_gen == loci_geni,
               Generation == 1 | Generation %% (Gen * 3) == 0)
      p <- locus_data %>% 
        ggplot(aes(factor(Generation), Phenotype, group = factor(Generation)))+ 
        geom_boxplot()+
        facet_wrap(~h2_sd, ncol = 4, scales = "free_y")+
        labs(x = "Generation")+ 
        scale_x_discrete(
          breaks = c("0", "300", "600", "900", "1200", "1500", "1800", "2000"),
          labels = c("0", "300", "600", "900", "1200", "1500", "1800", "2000")
        )+
        mytheme()
      plots[[loci_geni]] <- (p) #ggplotly
    }
    
    # plots for constant or sinusoidal II selections
    
  } else if (plot_type == "loci") {
    loci <- unique(combined_data$loci)
    for(locus in loci){
      locus_data <- combined_data %>% 
        filter(loci == locus)
      p <- locus_data %>% 
        ggplot(aes(factor(Generation), Phenotype, group = factor(replicate)))+
        geom_line(linewidth = 0.1, alpha = 1)+
        facet_wrap(~h2_sd, ncol = 4, scales = "free_y")+
        labs(x = "Generation") + 
        scale_x_discrete(
          breaks = c("0", "300", "600", "900", "1200", "1500", "1800", "2000"),
          labels = c("0", "300", "600", "900", "1200", "1500", "1800", "2000")
        )+
        mytheme()
      plots[[locus]] <- (p) #ggplotly
    }
  }
  
  return(list(combined_data = combined_data, plots = plots))
}

##################################################################
#################### Line plots for phenotypes ###################
##################################################################

# # The first part deals with constant and sinusoidal II
# process_files <- function(dirpath, pattern, plot_type) {
#   files <- list.files(dirpath, pattern, full.names = TRUE)
#   
#   dataframes <- list()
#   for(file in files){
#     replicate_id <- as.numeric(str_extract(file, "(?<=MeanPhenotypes)\\d+"))
#     herit <- as.numeric(str_extract(file, "(?<=H)0\\.\\d+"))
#     loci <- as.numeric(str_extract(file, "(?<=_n)\\d+"))
#     sd <- as.numeric(str_extract(file, "(?<=SD)\\d+"))
#     gen <- as.numeric(str_extract(file, "(?<=Gen)\\d+"))
#     
#     data <- read.csv(file, header = TRUE)
#     
#     data$herit <- rep(herit, nrow(data))
#     data$loci <- rep(loci, nrow(data))
#     data$sd <- rep(sd, nrow(data))
#     data$replicate <- rep(factor(replicate_id), nrow(data))
#     
#     # argument for linear and sinusoidal I selections
#     
#     if(!is.na(gen)) {
#       data$gen <- rep(gen, nrow(data))
#     }
#     
#     dataframes[[file]] <- data
#   }
#   combined_data <- bind_rows(dataframes) %>% 
#     mutate(h2_sd = paste("h2 = ", herit," ", "sd = ", sd, sep = ""),
#            loci_gen = paste("loci = ", loci," ", "gen = ", gen, sep = ""))
#  
#    # plots for either linear or sinusoidal selections
#   
#   plots <- list()
#   if (plot_type == "loci_gen") {
#     loci_gens <- unique(combined_data$loci_gen)
#     for(loci_geni in loci_gens){
#       locus_data <- combined_data %>% 
#         filter(loci_gen == loci_geni)
#       p <- locus_data %>% 
#         ggplot(aes(Generation, Phenotype, group = replicate, color = replicate))+
#         geom_line(linewidth = 0.1, alpha = 1)+
#         #geom_line(aes(y = Optimum), color = "darkred", linewidth = 0.05) +
#         #geom_point(data = filter(locus_data, Generation == 1), aes(y = Optimum), size = 0.1, color = "magenta") +
#         facet_wrap(~h2_sd, ncol = 4, scales = "free_y")+
#         theme_bw()
#         #theme(legend.position = "none")
#       plots[[loci_geni]] <- (p) #ggplotly
#     }
#     
#     # plots for constant or sinusoidal II selections
#     
#   } else if (plot_type == "loci") {
#     loci <- unique(combined_data$loci)
#     for(locus in loci){
#       locus_data <- combined_data %>% 
#         filter(loci == locus)
#       p <- locus_data %>% 
#         ggplot(aes(Generation, Phenotype, group = replicate, color = replicate))+
#         geom_line(linewidth = 0.1, alpha = 1)+
#         #geom_line(aes(y = Optimum), color = "darkred",linewidth = 0.05) +
#         #geom_point(data = filter(locus_data, Generation == 1), aes(y = Optimum), size = 0.1, color = "magenta")+
#         facet_wrap(~h2_sd, ncol = 4, scales = "free_y")+
#         theme_bw()#+
#         #theme(legend.position = "none")
#       plots[[locus]] <- (p) #ggplotly
#     }
#   }
#   
#   return(list(combined_data = combined_data, plots = plots))
# }
# 
# #####################################
# ############ Calling the code #######
# #####################################
# 
