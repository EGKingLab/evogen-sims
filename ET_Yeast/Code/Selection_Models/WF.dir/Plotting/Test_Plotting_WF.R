
# Load the necessary packages
rm(list = ls())
library(cowplot)
library(data.table)
library(tidyverse)
library(doParallel)
n_cores <- detectCores() - 2  # Adjust this based on your system

# Initialize parallel processing
cl <- makeCluster(n_cores)
registerDoParallel(cl)

# Define the directory path and file name pattern
dir_path <- "~/YeastProj.dir/evogen-sims/ET_Yeast/output.dir/Selection_Models/WF.dir/LinFS.dir/"
pattern <- "^genome\\d+_n\\d+_H0[.0-9]+SD\\d+Gen\\d+\\.csv$" #"^genome.*_n100_H0\\.\\dSD3Gen20\\.csv$" # 

# allele_feq_fx <- function(dir_path, pattern){
file_list <- list.files(dir_path, pattern = pattern, full.names = TRUE)

# Create an empty list to store the data frames
data_frames <- list()

for (file_name in file_list) {
  # Extract "replicate", "heritability", and "loci" values from the file name using regular expressions
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
         hetero = 2*Frequency*(1 - Frequency)) 

mytheme <- function(){
  theme_set(theme_cowplot())+
    theme(axis.title = element_text(face = "bold"),
          axis.text = element_text(face = "italic"),
          plot.title = element_text(hjust = 0.5),
          legend.position = "none",
          #panel.background = element_rect(fill = "white"), 
          plot.background = element_rect(fill = "white"))
  }

loci.heritabilities = unique(allele_data$loci_herit)
SD.GenDistList = unique(allele_data$SD_GenDist)
plots_list <- list()
het_plotlist <- list()
output_dir <- "~/YeastProj.dir/evogen-sims/ET_Yeast/output.dir/Selection_Models/WF.dir/Plots.dir/CS.dir/"


loci.heritabilities <- unique(allele_data$loci_herit)
SD.GenDistList <- unique(allele_data$SD_GenDist)

for (loci.herit in loci.heritabilities) {
  plot1 <- allele_data %>%
    filter(loci_herit == loci.herit) %>%
    group_by(Generation, Position, loci_herit, SD_GenDist) %>%
    summarize(mean_Freq = mean(Frequency), .groups = "keep") %>% 
    ggplot(aes(x = Generation, y = mean_Freq, color = as.factor(Position))) +
    geom_line(linewidth = 0.2) +
    facet_wrap(~SD_GenDist, ncol = 3) +
    ggtitle(paste("Frequency plot for Locus Heritability", loci.herit, sep=" ", "respectively")) +
    labs(y = "Mean Frequency")+
    mytheme()
  
  plot2 <- allele_data %>%
    filter(loci_herit == loci.herit) %>%
    group_by(Generation, Position, loci_herit, SD_GenDist) %>%
    summarize(mean_Hetero = mean(hetero), .groups = "keep") %>% 
    ggplot(aes(x = Generation, y = mean_Hetero, color = as.factor(Position))) +
    geom_line(linewidth = 0.2) +
    facet_wrap(~SD_GenDist, ncol = 3) +
    ggtitle(paste("Heterozygosity plot for Locus Heritability", loci.herit, sep=" ", "respectively")) +
    labs(y = "Mean Heterozygosity")+
    mytheme()
  
  plot_name1 <- paste("Frequency_plot_for_Locus_Herit_", loci.herit, sep = "")
  plot_name2 <- paste("Heterozygosity_plot_for_Locus_Herit_", loci.herit, sep = "")
  plots_list[[plot_name1]] <- plot1
  het_plotlist[[plot_name2]] <- plot2
  
  # image_file <- file.path(output_dir, paste0(plot_name, ".png"))
  # ggsave(image_file, plot, width = 8, height = 6)
}

unlist(plots_list1)

unlist(het_plotlist)


stopCluster(cl)

###############################################################################
#####################  Spectral Density Plots  ################################
###############################################################################

# Load the necessary packages
rm(list = ls())
library(cowplot)
library(data.table)
library(tidyverse)
library(doParallel)
n_cores <- detectCores() - 2  # Adjust this based on your system

# Initialize parallel processing
cl <- makeCluster(n_cores)
registerDoParallel(cl)

library(dplyr)
library(tidyr)
library(ggplot2)

dir_path <- "~/YeastProj.dir/evogen-sims/ET_Yeast/output.dir/Selection_Models/WF.dir/LinFS.dir/"
pattern <- "^genome\\d+_n\\d+_H0[.0-9]+SD\\d+Gen\\d+\\.csv$" 
spectral_span = 3

file_list <- list.files(dir_path, pattern = pattern, full.names = TRUE)

process_file <- function(file_name) {
  repl <- as.integer(str_extract(file_name,  "(?<=genome)[0-9]+"))
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

# Create facet plots based on loci_herit and SD_GenDist
PlotsList <- list()
for (loci.herit in loci.heritabilities) {
    freqplot <- combined_data2 %>%
      filter(loci_herit == loci.herit) %>%#
      mutate(Frequency = 1/Frequency) %>%
      filter(Frequency < 80) %>%
      ggplot(aes(Frequency, spec)) +
      geom_line(linewidth = 1) +
      facet_wrap( ~ SD_GenDist, ncol = 3) +
      theme(legend.position = "none") +
      xlab("Periodicity") +
      ylab("Spectral Density")+
      ggtitle(paste("Spectrum plot for Locus Heritability", loci.herit)) +
    
    PlotName <- paste("Spectral_plot_for_Locus_Herit_", loci.herit, "_SD_Gen_", SD.Gen, sep = "")
    PlotsList[[PlotName]] <- freqplot
    image_file <- file.path(output_dir, paste0(PlotName, ".png"))
    ggsave(image_file, plot, width = 8, height = 6)

}

PlotsList

stopCluster(cl)

######################################################################

