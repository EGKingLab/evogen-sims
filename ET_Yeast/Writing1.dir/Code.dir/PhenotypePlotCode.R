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
  theme_set(theme_cowplot())+
    theme(axis.title = element_text(face = "bold"),
          axis.text = element_text(face = "italic"),
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
    
    data <- read.csv(file, header = TRUE)
    
    data$herit <- rep(herit, nrow(data))
    data$loci <- rep(loci, nrow(data))
    data$sd <- rep(sd, nrow(data))
    data$replicate <- rep(factor(replicate_id), nrow(data))
    
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
        ggplot(aes(Generation, Phenotype, group = replicate, color = replicate))+
        geom_line(linewidth = 0.1)+
        #geom_line(aes(y = Optimum), color = "darkred", linewidth = 0.05) +
        #geom_point(data = filter(locus_data, Generation == 1), aes(y = Optimum), size = 0.1, color = "magenta") +
        facet_wrap(~h2_sd, ncol = 4, scales = "free_y")+
        theme_bw()
        #theme(legend.position = "none")
      plots[[loci_geni]] <- ggplotly(p)
    }
    
    # plots for constant or sinusoidal II selections
    
  } else if (plot_type == "loci") {
    loci <- unique(combined_data$loci)
    for(locus in loci){
      locus_data <- combined_data %>% 
        filter(loci == locus)
      p <- locus_data %>% 
        ggplot(aes(Generation, Phenotype, group = replicate, color = replicate))+
        geom_line(linewidth = 0.1)+
        #geom_line(aes(y = Optimum), color = "darkred",linewidth = 0.05) +
        #geom_point(data = filter(locus_data, Generation == 1), aes(y = Optimum), size = 0.1, color = "magenta")+
        facet_wrap(~h2_sd, ncol = 4, scales = "free_y")+
        theme_bw()#+
        #theme(legend.position = "none")
      plots[[locus]] <- ggplotly(p)
    }
  }
  
  return(list(combined_data = combined_data, plots = plots))
}

#####################################
############ Calling the code #######
#####################################

# result <- process_files(dirpath, pattern, "loci_gen")
# combined_data <- result$combined_data
# plots <- result$plots
