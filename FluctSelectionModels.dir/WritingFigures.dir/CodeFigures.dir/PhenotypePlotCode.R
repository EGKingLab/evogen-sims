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

# mytheme <- function(){theme_classic() +
#     theme(legend.position = "none",
#           axis.text.x = element_text(size = 15, face = "bold"),
#           axis.text.y = element_text(size = 15, face = "bold"),
#           axis.line = element_line(linewidth = 3),
#           axis.title.x = element_text(size = 30, face = "bold", margin = margin(t = 25)), # 
#           axis.title.y = element_text(size = 30, face = "bold", margin = margin(r = 20)), # 
#           strip.text = element_text(size = 30, face = "bold"),
#           panel.spacing = unit(2, "lines"),
#           panel.grid = element_blank())
# }

mythemes <- theme_bw() +
  theme(
    text            = element_text(family = "sans"), 
    legend.position = "none",
    axis.text.x     = element_text(face = "bold", size = 50, angle =25, margin = margin(t = 20), hjust = 1),
    axis.text.y     = element_text(face = "bold", size = 50, angle = 15, margin = margin(r = 10)),
    axis.line       = element_line(size = 3),
    plot.title      = element_text(hjust = 0.01, face = "bold", size = 50, margin = margin(b = 1, unit = "lines")),
    plot.margin     = unit(c(5, 1, 1, 1), "lines"),
    axis.title.x    = element_text(size = 50, face = "bold", margin = margin(t = 30)),
    axis.title.y    = element_text(size = 50, face = "bold", margin = margin(r = 35)),
    strip.text      = element_text(size = 50, face = "bold"),
    panel.spacing   = unit(5, "lines"),
    panel.grid      = element_blank()
  )

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
    
    data <- read.csv(file)
    
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
    
    # Determine the interval based on Gen
    interval <- Gen * 3 
      
    #   if (Gen == 10) {
    #   Gen * 3   
    # # } else if (Gen == 20) {
    # #   as.integer(Gen * 4.5)  
    # } else {
    #   Gen * 3       
    # }
    
    for (loci_geni in loci_gens) {
      locus_data <- combined_data %>% 
        filter(loci_gen == loci_geni)# Generation == 1 | Generation %% (Gen*3) == 0)
      
      p <- locus_data %>% 
        ggplot(aes(Generation, Phenotype, 
                   group = factor(replicate), 
                   color = factor(replicate))) + 
        geom_line() +
        facet_wrap(~h2_sd, ncol = 3, dir = "v", scales = "free_y") +
        labs(x = "Generation") + 
        #scale_x_discrete(
        #  breaks = c("0", "300", "600", "900", "1200", "1500", "1800", "2000"),
        #  labels = c("0", "300", "600", "900", "1200", "1500", "1800", "2000")) +
        mythemes
      
      plots[[loci_geni]] <- p #ggplotly
    }
    
  } else if (plot_type == "loci") {
    loci <- unique(combined_data$loci)
    for(locus in loci) {
      locus_data <- combined_data %>% 
        filter(loci == locus)
      
      p <- locus_data %>% 
        ggplot(aes(Generation, Phenotype, 
                   group = factor(replicate),
                   color = factor(replicate))) +
        geom_line(size = 0.5, alpha = 1) +
        facet_wrap(~h2_sd, ncol = 3, dir = "v", scales = "free_y") +
        labs(x = "Generation") + 
        # scale_x_discrete(
        #   breaks = c("0", "300", "600", "900", "1200", "1500", "1800", "2000"),
        #   labels = c("0", "300", "600", "900", "1200", "1500", "1800", "2000")) +
        mythemes
      
      plots[[locus]] <- p #ggplotly
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
#     data <- read.csv(file)
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
