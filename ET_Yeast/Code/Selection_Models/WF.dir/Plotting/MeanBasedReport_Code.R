################################## packages ###################################
# rm(list = ls())
#install.packages("tidyverse")
library(cowplot)
library(data.table)
library(tidyverse)
library(purrr)
library(doParallel)
theme_set(theme_cowplot())
###############################################################################
##### Allele Frequency Function For Linear and Sin I FS ######################
###############################################################################

allele_feq_fx <- function(dir_path, pattern, output_dir, output_dir2) {
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
    
    plot2 <- allele_data2 %>% 
      ggplot(aes(x = Generation, y = mean_Hetero, color = as.factor(Position))) +
      geom_line(linewidth = 0.2) +
      facet_wrap(~SD_GenDist, ncol = 3) +
      ggtitle(paste("Heterozygosity plot for Locus Heritability", loci.herit, sep=" ", "respectively")) +
      labs(y = "Mean Heterozygosity")+
      mytheme()
    # plot_name <- paste0("Frequency_plot_for_Locus_Herit_", loci.herit)
    # plots_list[[plot_name]] <- plot
    
    plot_name1 <- paste0("Frequency_plot_for_Locus_Herit_", loci.herit, sep = "")
    plot_name2 <- paste0("Heterozygosity_plot_for_Locus_Herit_", loci.herit, sep = "")
    plots_list[[plot_name1]] <- plot1
    het_plotlist[[plot_name2]] <- plot2
    
    image_file <- file.path(output_dir, paste0(plot_name1, ".png"))
    ggsave(image_file, plot1, width = 8, height = 6)
    
    image_file2 <- file.path(output_dir2, paste0(plot_name2, ".png"))
    ggsave(image_file2, plot2, width = 8, height = 6)
  }
  
  stopCluster(cl)
  
  return(list(plots_list, het_plotlist))
}

###############################################################################
################## Phenotype Function For Linear and Sin I FS ################
###############################################################################

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
    plots_list[[plot_name]] <- plot
    
    image_file <- file.path(output_dir, paste0(plot_name, ".png"))
    ggsave(image_file, plot, width = 8, height = 6)
  }
  
  stopCluster(cl)
  
  return(list(plots_list))
}

################################################################################
############################# Spectral Analysis ################################
################################################################################

# Function to process files based on specified patterns
process_files <- function(dir_path, pattern, spectral_span, output_dir) {
  file_list <- list.files(dir_path, pattern = pattern, full.names = TRUE)
  
  process_file <- function(file_name) {
    # Extract information from the file name
    repl <- as.integer(str_extract(file_name, "(?<=genome)\\d+"))
    herit <- as.numeric(str_extract(file_name, "(?<=H0\\.)\\d+"))
    loci <- as.numeric(str_extract(file_name, "(?<=_n)\\d+"))
    SD <- as.numeric(str_extract(file_name, "(?<=SD)\\d+"))
    GenDist <- as.numeric(str_extract(file_name, "(?<=Gen)\\d+"))
    
    data <- read.csv(file_name, header = TRUE)
    
    data$repl <- repl
    data$herit <- herit
    data$loci <- loci
    data$SD <- SD
    data$GenDist <- GenDist
    
    return(data)
  }
  
  combined_data <- purrr::map_dfr(file_list, process_file)
  
  allele_data <- combined_data %>% 
    mutate(loci_herit = paste(loci, herit, sep = "_"),
           SD_GenDist = paste(SD, GenDist, sep = "_")) 
  
  all_results <- list()
  
  loci.heritabilities = unique(allele_data$loci_herit)
  SD.GenDistList = unique(allele_data$SD_GenDist)
  
  # Function to process each combination of loci.herit and SD.Gen
  process_combination <- function(loci.herit, SD.Gen) {
    allele_data2 <- allele_data %>%
      filter(loci_herit == loci.herit, SD_GenDist == SD.Gen) %>%
      group_by(Generation, Position, loci_herit, SD_GenDist) %>%
      summarize(mean_Freq = mean(Frequency), .groups = "keep")
    
    # Pivot the data_wide to a wide format
    data_wide <- pivot_wider(
      data = allele_data2,
      names_from = Position,
      values_from = mean_Freq
    )
    
    data_wide[is.na(data_wide)] <- 0
    
    # Select columns from 4th to the end (assuming columns 1 to 3 are Generation, loci_herit, SD_GenDist)
    ff1 <- data_wide[, 4:ncol(data_wide)]
    
    ttest <- spectrum(ff1, spans = spectral_span, plot = FALSE)
    
    if (!is.matrix(ttest$spec)) {
      out.spect <- ttest$spec
    } else {
      out.spect <- rowMeans(ttest$spec)
    }
    
    dd <- data.frame(
      "Frequency" = rep(ttest$freq),
      "loci_herit" = loci.herit,
      "spec" = out.spect,
      "SD_GenDist" = SD.Gen
    )
    
    return(dd)
  }
  
  all_results2 <- list()
  
  for (loci.herit in loci.heritabilities) {
    for (SD.Gen in SD.GenDistList) {
      result <- process_combination(loci.herit, SD.Gen)
      all_results2[[paste("loci_herit_", loci.herit, "_SD_Gen_", SD.Gen, sep = "")]] <- result
    }
  }
  
  combined_data2 <- bind_rows(all_results2)
  
  # Call the process_plots function to create and save the plots
  process_plots(combined_data2, output_dir)
}

# Function to create and save plots
process_plots <- function(combined_data2, output_dir) {
  PlotsList <- list()
  loci.heritabilities = unique(combined_data2$loci_herit)
  SD.GenDistList = unique(combined_data2$SD_GenDist)
  
  for (loci.herit in loci.heritabilities) {
    freqplot <- combined_data2 %>%
      filter(loci_herit == loci.herit) %>%
      mutate(Frequency = 1/Frequency) %>%
      filter(Frequency < 80) %>%
      ggplot(aes(Frequency, spec)) +
      geom_line(linewidth = 1) +
      facet_wrap( ~ SD_GenDist, ncol = 3) +
      theme(legend.position = "none") +
      xlab("Periodicity") +
      ylab("Spectral Density") +
      ggtitle(paste("Spectrum plot for Locus Heritability", loci.herit))
    
    PlotName <- paste0("Spectral_plot_for_Locus_Herit_", loci.herit, ".png", sep = "")
    PlotsList[[PlotName]] <- freqplot
    
    image_file <- file.path(output_dir, PlotName)
    ggsave(image_file, freqplot, width = 8, height = 6)
  }
  
  return(list(PlotsList))
}

###################### CS and Sin II plots #################

################## Allele frequency plots ###############

CSin_freq_fx <- function(dir_path, pattern, output_dir, output_dir2) {
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
    
    # Read the file (adjust the reading method as needed)
    data <- read.csv(file_name, header = TRUE)
    
    # Add the "replicate", "herititability", and "loci" columns
    
    data$repl <- rep(repl, nrow(data))
    data$herit <- rep(herit, nrow(data))
    data$loci <- rep(loci, nrow(data))
    data$SD <- rep(SD, nrow(data))
    
    # Add the data frame to the list
    data_frames[[file_name]] <- data
  }
  
  # Combine all the data frames into one
  combined_data <- bind_rows(data_frames)
  
  allele_data <- combined_data %>% 
    mutate(loci_herit = paste(loci, herit, sep = "_"),
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
  plots_list <- list()
  het_plotlist <- list()
  
  for (loci.herit in loci.heritabilities) {
    allele_data2 <- allele_data %>%
      filter(loci_herit == loci.herit) %>%
      group_by(Generation, Position, loci_herit, SD) %>%
      summarize(mean_Freq = mean(Frequency), mean_Hetero = mean(hetero), .groups = "keep")
    
    plot1 <- allele_data2 %>% ggplot(aes(x = Generation, y = mean_Freq, color = as.factor(Position))) +
      geom_line(linewidth = 0.2) +
      facet_wrap(~SD, ncol = 2) +
      ggtitle(paste("Frequency plot for Locus Heritability", loci.herit)) +
      labs(y = "Mean Frequency")+
      mytheme()
    
    plot2 <- allele_data2 %>%
      ggplot(aes(x = Generation, y = mean_Hetero, color = as.factor(Position))) +
      geom_line(linewidth = 0.2) +
      facet_wrap(~SD, ncol = 2) +
      ggtitle(paste("Heterozygosity plot for Locus Heritability", loci.herit, sep=" ", "respectively")) +
      labs(y = "Mean Heterozygosity")+
      mytheme()
    
    plot_name1 <- paste0("Frequency_plot_for_Locus_Herit_", loci.herit, sep = "")
    plot_name2 <- paste0("Heterozygosity_plot_for_Locus_Herit_", loci.herit, sep = "")
    plots_list[[plot_name1]] <- plot1
    het_plotlist[[plot_name2]] <- plot2
    
    image_file <- file.path(output_dir, paste0(plot_name1, ".png"))
    ggsave(image_file, plot1, width = 8, height = 6)
    
    image_file2 <- file.path(output_dir2, paste0(plot_name2, ".png"))
    ggsave(image_file2, plot2, width = 8, height = 6)
  }
  
  stopCluster(cl)
  
  return(list(plots_list, het_plotlist))
}

dir_path <- "~/YeastProj.dir/evogen-sims/ET_Yeast/output.dir/Selection_Models/WF.dir/SinFSGen.dir/"
pattern <- "^genome\\d+_n\\d+_H0[.0-9]+SD\\d+\\.csv$"
output_dir <- "~/YeastProj.dir/evogen-sims/ET_Yeast/output.dir/Selection_Models/WF.dir/Plots.dir/SinFSGen.dir/Frequency/"
output_dir2 <- "~/YeastProj.dir/evogen-sims/ET_Yeast/output.dir/Selection_Models/WF.dir/Plots.dir/SinFSGen.dir/Heterozygosity/"
CSin_freq_fx(dir_path, pattern, output_dir, output_dir2)


################## Phenotypes frequency plots ###############

CSin_pheno_fx <- function(dir_path, pattern, output_dir) {
  n_cores <- detectCores() - 2
  cl <- makeCluster(n_cores)
  registerDoParallel(cl)
  
  file_list <- list.files(dir_path, pattern = pattern, full.names = TRUE)
  
  data_frames <- list()
  
  for (file_name in file_list) {
    repl <- as.integer(str_extract(file_name,  "(?<=MeanPhenotypes)[0-9]+"))
    herit <- as.numeric(str_extract(file_name, "(?<=H0\\.)\\d+"))
    loci <- as.numeric(str_extract(file_name, "(?<=_n)\\d+"))
    SD <- as.numeric(str_extract(file_name, "(?<=SD)\\d+"))
    
    # Read the file (adjust the reading method as needed)
    data <- read.csv(file_name, header = TRUE)
    
    # Add the "replicate", "herititability", and "loci" columns
    data$repl <- rep(repl, nrow(data))
    data$herit <- rep(herit, nrow(data))
    data$loci <- rep(loci, nrow(data))
    data$SD <- rep(SD, nrow(data))
    #data$GenDist <- rep(GenDist, nrow(data))
    
    # Add the data frame to the list
    data_frames[[file_name]] <- data
  }
  
  # Combine all the data frames into one
  combined_data <- bind_rows(data_frames)
  
  allele_data <- combined_data %>% 
    mutate(loci_herit = paste(loci, herit, sep = "_")) 
  
  mytheme <- function(){
    theme_set(theme_cowplot())+
      theme(axis.title = element_text(face = "bold"),
            axis.text = element_text(face = "italic"),
            plot.title = element_text(hjust = 0.5),
            legend.position = "none",
            plot.background = element_rect(fill = "white"))
  }
  
  loci.heritabilities = unique(allele_data$loci_herit)
  plots_list <- list()
  
  for (loci.herit in loci.heritabilities) {
    plot <- allele_data %>%
      filter(loci_herit == loci.herit) %>%
      group_by(Generation, loci_herit, SD) %>%
      summarize(mean_Pheno = mean(Phenotype), .groups = "keep") %>% 
      ggplot(aes(x = Generation, y = mean_Pheno)) +
      geom_line(linewidth = 0.2) +
      facet_wrap(~SD, ncol = 2) +
      ggtitle(paste("Phenotype plot for Locus Heritability", loci.herit)) +
      labs(y = "Mean Phenotypes")+
      mytheme()
    
    plot_name <- paste("Pheno_plot_for_Locus_Herit_", loci.herit)
    plots_list[[plot_name]] <- plot
    
    image_file <- file.path(output_dir, paste0(plot_name, ".png"))
    ggsave(image_file, plot, width = 8, height = 6)
  }
  
  stopCluster(cl)
  
  return(unlist(plots_list))
}

dir_path <- "~/YeastProj.dir/evogen-sims/ET_Yeast/output.dir/Selection_Models/WF.dir/SinFSGen.dir/"
pattern <- "^MeanPhenotypes\\d+_n\\d+_H0[.0-9]+SD\\d+\\.csv$"
output_dir <- "~/YeastProj.dir/evogen-sims/ET_Yeast/output.dir/Selection_Models/WF.dir/Plots.dir/SinFSGen.dir/Phenotypes/"
CSin_pheno_fx(dir_path, pattern, output_dir)

########################## Spectral Analysis #######################

data_process_files <- function(dir_path, pattern, spectral_span, output_dir) {
  file_list <- list.files(dir_path, pattern = pattern, full.names = TRUE)
  
  process_file <- function(file_name) {
    # Extract information from the file name
    repl <- as.integer(str_extract(file_name, "(?<=genome)\\d+"))
    herit <- as.numeric(str_extract(file_name, "(?<=H0\\.)\\d+"))
    loci <- as.numeric(str_extract(file_name, "(?<=_n)\\d+"))
    SD <- as.numeric(str_extract(file_name, "(?<=SD)\\d+"))
    GenDist <- as.numeric(str_extract(file_name, "(?<=Gen)\\d+"))
    
    data <- read.csv(file_name, header = TRUE)
    
    data$repl <- repl
    data$herit <- herit
    data$loci <- loci
    data$SD <- SD
    data$GenDist <- GenDist
    
    return(data)
  }
  
  combined_data <- purrr::map_dfr(file_list, process_file)
  
  allele_data <- combined_data %>% 
    mutate(loci_herit = paste(loci, herit, sep = "_"),
           SD_GenDist = paste(SD, GenDist, sep = "_")) 
  
  all_results <- list()
  
  loci.heritabilities = unique(allele_data$loci_herit)
  SD.GenDistList = unique(allele_data$SD_GenDist)
  
  # Function to process each combination of loci.herit and SD.Gen
  process_combination <- function(loci.herit, SD.Gen) {
    allele_data2 <- allele_data %>%
      filter(loci_herit == loci.herit, SD_GenDist == SD.Gen) %>%
      group_by(Generation, Position, loci_herit, SD_GenDist) %>%
      summarize(mean_Freq = mean(Frequency), .groups = "keep")
    
    # Pivot the data_wide to a wide format
    data_wide <- pivot_wider(
      data = allele_data2,
      names_from = Position,
      values_from = mean_Freq
    )
    
    data_wide[is.na(data_wide)] <- 0
    
    # Select columns from 4th to the end (assuming columns 1 to 3 are Generation, loci_herit, SD_GenDist)
    ff1 <- data_wide[, 4:ncol(data_wide)]
    
    ttest <- spectrum(ff1, spans = spectral_span, plot = FALSE)
    
    if (!is.matrix(ttest$spec)) {
      out.spect <- ttest$spec
    } else {
      out.spect <- rowMeans(ttest$spec)
    }
    
    dd <- data.frame(
      "Frequency" = rep(ttest$freq),
      "loci_herit" = loci.herit,
      "spec" = out.spect,
      "SD_GenDist" = SD.Gen
    )
    
    return(dd)
  }
  
  all_results2 <- list()
  
  for (loci.herit in loci.heritabilities) {
    for (SD.Gen in SD.GenDistList) {
      result <- process_combination(loci.herit, SD.Gen)
      all_results2[[paste("loci_herit_", loci.herit, "_SD_Gen_", SD.Gen, sep = "")]] <- result
    }
  }
  
  combined_data2 <- bind_rows(all_results2)
  
  # Call the process_plots function to create and save the plots
  process_plots(combined_data2, output_dir)
}

# Function to create and save plots
process_plots <- function(combined_data2, output_dir) {
  PlotsList <- list()
  loci.heritabilities = unique(combined_data2$loci_herit)
  SD.GenDistList = unique(combined_data2$SD_GenDist)
  
  for (loci.herit in loci.heritabilities) {
    freqplot <- combined_data2 %>%
      filter(loci_herit == loci.herit) %>%
      mutate(Frequency = 1/Frequency) %>%
      filter(Frequency < 40) %>%
      ggplot(aes(Frequency, spec)) +
      geom_line(linewidth = 1) +
      facet_wrap( ~ SD_GenDist, ncol = 2) +
      theme(legend.position = "none") +
      xlab("Periodicity") +
      ylab("Spectral Density") +
      ggtitle(paste("Spectrum plot for Locus Heritability", loci.herit))
    
    PlotName <- paste("Spectral_plot_for_Locus_Herit_", loci.herit, ".png", sep = "")
    PlotsList[[PlotName]] <- freqplot
    
    image_file <- file.path(output_dir, PlotName)
    ggsave(image_file, freqplot, width = 8, height = 6)
  }
  
  return(unlist(PlotsList))
}

dir_path <- "~/YeastProj.dir/evogen-sims/ET_Yeast/output.dir/Selection_Models/WF.dir/SinFSGen.dir/"
pattern <- "^genome\\d+_n\\d+_H0[.0-9]+SD\\d+\\.csv$" #"^genome\\d+_n\\d+_H0\\.\\d+SD\\d+\\.csv$
output_dir <- "~/YeastProj.dir/evogen-sims/ET_Yeast/output.dir/Selection_Models/WF.dir/Plots.dir/SinFSGen.dir/Spectrum/"

spectral_span = 3
data_process_files(dir_path, "^genome\\d+_n\\d+_H0\\.\\d+SD\\d+\\.csv$", spectral_span, output_dir)