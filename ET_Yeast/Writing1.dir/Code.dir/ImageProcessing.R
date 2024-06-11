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

# Custom labeller function
facet_labeller <- function(variable, value) {
  labels <- LETTERS[1:length(value)]
  return(labels)
}

#################### Data Processing and plotting ###############

# The first part deals with constant and sinusoidal II
process_files <- function(dirpath, pattern, plot_type, save_path) {
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
        filter(loci_gen == loci_geni)
      
      p <- locus_data %>%
        ggplot(aes(Generation, Frequency, group = postion_effect_init,
                   color = postion_effect_init)) +
        geom_line(linewidth = 0.2) +
        facet_wrap(~h2_sd, ncol = 4, labeller = facet_labeller) +
        ylim(min = 0, max = 1) +
        theme_bw() +
        theme(legend.position = "none")
      
      plots[[loci_geni]] <- p
      
      # Save plot as PNG
      ggsave(file = file.path(save_path, paste0("plot_", loci_geni, ".png")), 
             plot = p, width = 10, height = 7)
    }
    
  } else if (plot_type == "loci") {
    loci <- unique(combined_data$loci)
    for(locus in loci){
      locus_data <- combined_data %>%
        filter(loci == locus)
      
      p <- locus_data %>%
        ggplot(aes(Generation, Frequency, group = postion_effect_init, 
                   color = postion_effect_init)) +
        geom_line(linewidth = 0.2) +
        facet_wrap(~h2_sd, ncol = 4, labeller = facet_labeller) +
        ylim(min = 0, max = 1) +
        theme_bw() +
        theme(legend.position = "none")
      
      plots[[locus]] <- p
      
      # Save plot as PNG
      ggsave(file = file.path(save_path, paste0("plot_", locus, ".png")), 
             plot = p, width = 10, height = 7)
    }
  }
  
  return(list(combined_data = combined_data, plots = plots))
}

# Example of how to call the function with a user-provided file path
dirpath <- "~/YeastProj.dir/evogen-sims/ET_Yeast/output.dir/Selection_Models/WF.dir/LinFS.dir/"
pattern <- "^genome15\\D"
save_path <- "/home/etb68/YeastProj.dir/evogen-sims/ET_Yeast/Writing1.dir/Output.dir/AllPlots/AlleleFreq.dir"
result <- process_files(dirpath, pattern, "loci_gen", save_path)
combined_data <- result$combined_data
plots <- result$plots
plots
