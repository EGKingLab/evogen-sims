library(tidyverse)

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
  
  freq_filename <- file.path(mypath, paste0(basename(mygenomefile), "_Freq.png"))
  
  # Save the plot Frequency plots
  png(filename = freq_filename, width = 6, height = 4, units = "in", res = 300)
  print(Freq_plot)
  dev.off()
  
  effect_filename <- file.path(mypath, paste0(basename(mygenomefile), "_Effect.png"))
  
  # Save the plot Effect plots
  png(filename = effect_filename, width = 6, height = 4, units = "in", res = 300)
  print(Effect_plot)
  dev.off()
}

# mypath <- "/home/path.dir"
# 
# mygenomefiles <- list.files(path = "/home/path.dir/out.dir", pattern = "^genome.*\\.csv$", full.names = TRUE)

for (mygenomefile in mygenomefiles) {
  df <- read.csv(mygenomefile)
  my_plots(df, mygenomefile, mypath)
  
  mon_plots <- function(data, myphenofile, mypath){
    
    pheno_plot <- 
      ggplot(data, aes(x = Generation, y = Phenotype))+
      geom_line()+
      mytheme()
    
    Fitness_plot <- ggplot(data, aes(x = Generation, y = meanFitness))+
      geom_line()+
      mytheme()
    
    PhenoDens <- ggplot(data, aes(x = Phenotype))+
      geom_histogram(bins = 200)
    mytheme()
    
    FitDens <- ggplot(data, aes(x = meanFitness))+
      geom_density()
    mytheme()
    
    pheno_filename <- file.path(mypath, paste0(basename(myphenofile), "_Phenotype.png"))
    pheno_DensName <- file.path(mypath, paste0(basename(myphenofile), "_PhenoDens.png"))
    Fitness_filename <- file.path(mypath, paste0(basename(myphenofile), "_Fitness.png"))
    Fitness_DensName <- file.path(mypath, paste0(basename(myphenofile), "_FitDens.png"))
    
    # Save the plot Effect plots
    png(filename = pheno_filename, width = 6, height = 4, units = "in", res = 300)
    print(pheno_plot)
    dev.off()
    
    png(filename = pheno_DensName, width = 6, height = 4, units = "in", res = 300)
    print(PhenoDens)
    dev.off()
    
    png(filename = Fitness_filename, width = 6, height = 4, units = "in", res = 300)
    print(Fitness_plot)
    dev.off()
    
    png(filename = Fitness_DensName, width = 6, height = 4, units = "in", res = 300)
    print(FitDens)
    dev.off()
  }
}

#myphenofiles <- list.files(path = "/home/path.dir/out.dir", pattern = "^MeanPhenotypes.*\\.csv$", full.names = TRUE)

# for (myphenofile in myphenofiles) {
#   df <- read.csv(myphenofile)
#   mon_plots(df, myphenofile, mypath)
# }



