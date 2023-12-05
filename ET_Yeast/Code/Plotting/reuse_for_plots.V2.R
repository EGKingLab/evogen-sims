library(tidyverse)
library(cowplot)

# Define a function to create a custom theme
mytheme <- function() {
  theme_set(theme_cowplot()) +
    theme(axis.title = element_text(face = "bold"),
          axis.text = element_text(face = "italic"),
          plot.title = element_text(hjust = 0.5),
          legend.position = "none")
}

# Define your plotting function
my_plots <- function(df, mygenomefile, mypath) {

  # For Allele Frequency plots
  Freq_plot <- ggplot(df, aes(x = Generation, y = Frequency, color = as.character(Position))) +
    geom_line(linewidth = 0.2) +
    mytheme()

  # For Effect plots
  Effect_plot <- ggplot(df, aes(x = Effect)) +
    geom_density() +
    mytheme()

  # Display the Frequency plot
  print(Freq_plot)

  # Display the Effect plot
  print(Effect_plot)
}

mypath <- "/home/path.dir"

mygenomefiles <- list.files(path = "/home/path.dir/out.dir", pattern = "^genome.*\\.csv$", full.names = TRUE)

for (mygenomefile in mygenomefiles) {
  df <- read.csv(mygenomefile)
  my_plots(df, mygenomefile, mypath)
}

mon_plots <- function(data, myphenofile, mypath) {
  pheno_plot <- ggplot(data, aes(x = Generation, y = Phenotype)) +
    geom_line() +
    mytheme()

  Fitness_plot <- ggplot(data, aes(x = Generation, y = meanFitness)) +
    geom_line() +
    mytheme()

  PhenoDens <- ggplot(data, aes(x = Phenotype)) +
    geom_histogram(bins = 200) +
    mytheme()

  FitDens <- ggplot(data, aes(x = meanFitness)) +
    geom_density() +
    mytheme()

  # Display the Phenotype plot
  print(pheno_plot)

  # Display the Phenotype Density plot
  print(PhenoDens)

  # Display the Fitness plot
  print(Fitness_plot)

  # Display the Fitness Density plot
  print(FitDens)
}

myphenofiles <- list.files(path = "/home/path.dir/out.dir", pattern = "^MeanPhenotypes.*\\.csv$", full.names = TRUE)

for (myphenofile in myphenofiles) {
  df <- read.csv(myphenofile)
  mon_plots(df, myphenofile, mypath)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

process_files <- function(files) {
  all_results <- data.frame()
  
  for (file_name in files) {
    file <- read.csv(file_name) %>% dplyr::select(!c(Origin, Effect))
    dat_wide <- pivot_wider(file[,c("Generation","Position","Frequency")], names_from=Position, values_from=Frequency)
    dat_wide[is.na(dat_wide)] <- 0 
    ff1 <- dat_wide[,2:101]
    ttest <- spectrum(ff1, spans=2, plot=FALSE)
    out.spect <- rowMeans(ttest$spec)
    dd <- data.frame("Frequency" = rep(ttest$freq), 
                     "File" = basename(file_name),
                     "sel" = rep("FS", each = 1000),
                     "spec" = out.spect)
    all_results <- rbind(all_results, dd)
  }
  
  return(all_results)
}

