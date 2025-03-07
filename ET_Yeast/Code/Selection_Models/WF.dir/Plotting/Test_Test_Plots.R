library(dplyr)
library(ggplot2)
library(cowplot)
library(data.table)
library(doParallel)

# Number of CPU cores on the cluster
n_cores <- detectCores()  # Adjust this based on your cluster setup

# Initialize parallel processing on the cluster
cl <- makeCluster(n_cores)
registerDoParallel(cl)

# Define the directory path and file name pattern on the cluster
# Make sure these paths are accessible on the cluster nodes
dir_path <- "/storage/hpc/group/kinglab/etb68/evogen-sims/ET_Yeast/output.dir/Selection_Models/WF.dir/CS.dir/"  # Adjust this
pattern <- "^genome\\d+_n\\d+_H0[.0-9]+SD\\d+Gen\\d+\\.csv$"

# List the files in the directory that match the pattern
file_list <- list.files(dir_path, pattern = pattern, full.names = TRUE)

# Create an empty list to store the data frames
data_frames <- list()

# Loop over the files
foreach(file_name = file_list) %dopar% {
  # Extract "replicate", "heritability", and "loci" values from the file name using regular expressions
  repl <- as.integer(gsub(".*genome([0-9]+).*", "\\1", file_name))
  herit <- as.numeric(gsub(".*H0\\.([0-9]+).*", "\\1", file_name))
  loci <- as.numeric(gsub(".*_n([0-9]+).*", "\\1", file_name))
  SD <- as.numeric(gsub(".*SD([0-9]+).*", "\\1", file_name))
  GenDist <- as.numeric(gsub(".*Gen([0-9]+).*", "\\1", file_name))

  # Read the file (adjust the reading method as needed)
  data <- read.csv(file_name, header = TRUE)

  # Add the "replicate", "herititability", and "loci" columns
  data$repl <- rep(repl, nrow(data))
  data$herit <- rep(herit, nrow(data))
  data$loci <- rep(loci, nrow(data))
  data$SD <- rep(SD, nrow(data))
  data$GenDist <- rep(GenDist, nrow(data))

  # Return the data frame
  return(data)
} -> data_frames

# Combine all the data frames into one
combined_data <- do.call(rbind, data_frames)

# View the first few rows of the combined data
head(combined_data)

# Add new columns
allele_data <- transform(combined_data,
                      loci_herit = paste(loci, herit, sep = "_"),
                      SD_GenDist = paste(SD, GenDist, sep = "_"))

# Define the theme for the plots
mytheme <- function(){
  theme_set(theme_cowplot()) +
    theme(axis.title = element_text(face = "bold"),
          axis.text = element_text(face = "italic"),
          plot.title = element_text(hjust = 0.5),
          legend.position = "none")
}

# Get the unique values of "loci_herit" and "SD_GenDist"
loci.heritabilities <- unique(allele_data$loci_herit)
SD.GenDistList <- unique(allele_data$SD_GenDist)

# Create an empty list to store the plots
plots_list <- list()

# Define the output directory
output_dir <- "/storage/hpc/group/kinglab/etb68/evogen-sims/ET_Yeast/output.dir/Selection_Models/WF.dir/Plots.dir/Test.dir"

# Loop over the "SD_GenDist" and "loci_herit" values
foreach(SD.Gen = SD.GenDistList, .combine = 'c') %:%
foreach(loci.herit = loci.heritabilities, .combine = 'c') %dopar% {
  # Filter the data
  filtered_data <- subset(allele_data, loci_herit == loci.herit & SD_GenDist == SD.Gen)

  # Group the data
  grouped_data <- dplyr::group_by(filtered_data, Generation, Position, loci_herit, SD_GenDist)

  # Summarize the data
  summarized_data <- dplyr::summarize(grouped_data, mean_Freq = mean(Frequency), .groups = "keep")

  # Create the plot
  plot <- ggplot(summarized_data, aes(x = Generation, y = mean_Freq, color = as.factor(Position))) +
    geom_line(linewidth = 0.2) +
    facet_wrap(~SD_GenDist, ncol = 3) +
    ggtitle(paste("Frequency plot for", loci.herit)) +
    mytheme()

  # Save the plot
  plot_name <- paste("Frequency_plot_for_", loci.herit, "_Locus_heritability_values_respectively")
  image_file <- file.path(output_dir, paste("plot_", plot_name, ".png", sep = ""))
  ggsave(image_file, plot, width = 8, height = 6)

  # Return the plot
  return(plot)
} -> plots_list

# View the list of plots
plots_list

