library(tidyverse)
library(purrr)
library(doParallel)
theme_set(theme_cowplot())
################################################################################
############################# Allele Frequency Function ########################
################################################################################

allele_feq_fx <- function(dir_path, pattern, output_dir) {
  n_cores <- detectCores() - 2 
  
  # Initialize parallel processing
  cl <- makeCluster(n_cores)
  registerDoParallel(cl)
  
  file_list <- list.files(dir_path, pattern = pattern, full.names = TRUE)
  
  # Create an empty list to store the data frames
  data_frames <- list()
  
  for (file_name in file_list) {
    repl <- as.integer(str_extract(file_name,  "(?<=genome)[0-9]+"))
    herit <- as.numeric(str_extract(file_name, "(?<=H0\\.)\\d+"))
    loci <- as.numeric(str_extract(file_name, "(?<=_n)\\d+"))
    SD <- as.numeric(str_extract(file_name, "(?<=SD)\\d+"))
    GenDist <- as.numeric(str_extract(file_name, "(?<=Gen)\\d+"))
    
    # Read the file (adjust the reading method as needed)
    data <- read.csv(file_name, header = TRUE)
    
    # Add the "replicate", "herititability", and "loci" columns
    data$repl <- rep(repl, nrow(data))
    data$herit <- rep(herit, nrow(data))
    data$loci <- rep(loci, nrow(data))
    data$SD <- rep(SD, nrow(data))
    data$GenDist <- rep(GenDist, nrow(data))
    
    # Add the data frame to the list
    data_frames[[file_name]] <- data
  }
  
  # Combine all the data frames into one
  combined_data <- bind_rows(data_frames)
  
  allele_data <- combined_data %>% 
    mutate(loci_herit = paste(loci, herit, sep = "_"),
           SD_GenDist = paste(SD, GenDist, sep = "_"),
           hetero = 2*Frequency*(1-Frequency)) 
  
  mytheme <- function(){
    theme_set(theme_cowplot())+
      theme(axis.title = element_text(face = "bold"),
            axis.text = element_text(face = "italic"),
            plot.title = element_text(hjust = 0.5),
            legend.position = "none",
            plot.background = element_rect(fill = "white"))
  }
  
  loci.heritabilities = unique(allele_data$loci_herit)
  SD.GenDistList = unique(allele_data$SD_GenDist)
  plots_list <- list()
  het_plotlist <- list()
  
  for (loci.herit in loci.heritabilities) {
    allele_data2 <- allele_data %>%
      filter(loci_herit == loci.herit) %>%
      group_by(Generation, Position, loci_herit, SD_GenDist) %>%
      summarize(mean_Freq = mean(Frequency), mean_Hetero = mean(hetero), .groups = "keep") 
    plot1 <- allele_data2%>% 
      ggplot(aes(x = Generation, y = mean_Freq, color = as.factor(Position))) +
      geom_line(linewidth = 0.2) +
      facet_wrap(~SD_GenDist, ncol = 3) +
      ggtitle(paste("Frequency plot for Locus Heritability", loci.herit)) +
      labs(y = "Mean Frequency")+
      mytheme()
    
    plot_name <- paste0("Frequency_plot_for_Locus_Herit_", loci.herit, sep = "")
    plots_list[[plot_name]] <- plot1
    plots_list <- plot1
    
    image_file <- file.path(output_dir, paste0(plot_name, ".png"))
    ggsave(image_file, plot1, width = 8, height = 6)
  }
  
  stopCluster(cl)
  
  print(plots_list)
}

################################################################################
####################### Phenotype Function #####################################
################################################################################

pheno_fx <- function(dir_path, pattern, output_dir) {
  n_cores <- detectCores() - 2  
  
  # Initialize parallel processing
  cl <- makeCluster(n_cores)
  registerDoParallel(cl)
  
  file_list <- list.files(dir_path, pattern = pattern, full.names = TRUE)
  
  # Create an empty list to store the data frames
  data_frames <- list()
  
  for (file_name in file_list) {
    repl <- as.integer(str_extract(file_name,  "(?<=MeanPhenotypes)[0-9]+"))
    herit <- as.numeric(str_extract(file_name, "(?<=H0\\.)\\d+"))
    loci <- as.numeric(str_extract(file_name, "(?<=_n)\\d+"))
    SD <- as.numeric(str_extract(file_name, "(?<=SD)\\d+"))
    GenDist <- as.numeric(str_extract(file_name, "(?<=Gen)\\d+"))
    
    # Read the file (adjust the reading method as needed)
    data <- read.csv(file_name, header = TRUE)
    
    # Add the "replicate", "herititability", and "loci" columns
    data$repl <- rep(repl, nrow(data))
    data$herit <- rep(herit, nrow(data))
    data$loci <- rep(loci, nrow(data))
    data$SD <- rep(SD, nrow(data))
    data$GenDist <- rep(GenDist, nrow(data))
    
    # Add the data frame to the list
    data_frames[[file_name]] <- data
  }
  
  # Combine all the data frames into one
  combined_data <- bind_rows(data_frames)
  
  allele_data <- combined_data %>% 
    mutate(loci_herit = paste(loci, herit, sep = "_"),
           SD_GenDist = paste(SD, GenDist, sep = "_")) 
  
  mytheme <- function(){
    theme_set(theme_cowplot())+
      theme(axis.title = element_text(face = "bold"),
            axis.text = element_text(face = "italic"),
            plot.title = element_text(hjust = 0.5),
            legend.position = "none",
            plot.background = element_rect(fill = "white"))
  }
  
  loci.heritabilities = unique(allele_data$loci_herit)
  SD.GenDistList = unique(allele_data$SD_GenDist)
  plots_list <- list()
  
  for (loci.herit in loci.heritabilities) {
    plot <- allele_data %>%
      filter(loci_herit == loci.herit) %>%
      group_by(Generation, loci_herit, SD_GenDist) %>%
      summarize(mean_Pheno = mean(Phenotype), .groups = "keep") %>% 
      ggplot(aes(x = Generation, y = mean_Pheno)) +
      geom_line(linewidth = 0.2) +
      facet_wrap(~SD_GenDist, ncol = 3) +
      ggtitle(paste("Phenotype plot for Locus Heritability", loci.herit)) +
      labs(y = "Mean Phenotype")+
      mytheme()
    
    plot_name <- paste0("Phenotype_plot_for_Locus_Herit_", loci.herit)
    plots_list <- plot
    # 
     image_file <- file.path(output_dir, paste0(plot_name, ".png"))
     ggsave(image_file, plot, width = 8, height = 6)
  }
  stopCluster(cl)
  
  print(plots_list)
}