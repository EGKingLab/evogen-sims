
theme_set(theme_cowplot())

# Define a function to create a custom theme
mytheme <- function() {
  theme_set(theme_cowplot()) +
    theme(axis.title = element_text(face = "bold"),
          axis.text = element_text(face = "italic"),
          plot.title = element_text(hjust = 0.5),
          legend.position = "none")
}

###### Function to plot allele frequency vs Generation ######

allele_plots <- function(df, mygenomefile) {
  replicate <- str_extract(mygenomefile, "(?<=genome)[0-9]+")
  df$replicate <- replicate
  
  # For Allele Frequency plots
  Freq_plot <- ggplot(df, aes(x = Generation, y = Frequency, color = as.character(Position))) +
    geom_line(linewidth = 0.5) +
    mytheme() +
    ggtitle(paste0(basename(mygenomefile))) +
    # Use facet_wrap with the replicate column
    facet_wrap(~replicate, ncol = 4)
  
  print(Freq_plot)
  
  # # For Effect plots
  # Effect_plot <- ggplot(df, aes(x = Effect)) +
  #   geom_density() +
  #   mytheme() +
  #   ggtitle(paste0((mygenomefile)))  # Added line
  
  # # Display the Effect plot
  # print(Effect_plot)
}

########## Function for Phenotype plots ###########

phenotype_plots <- function(data, myphenofile) {
  
  replicate <- str_extract(myphenofile, "(?<=MeanPhenotypes)[0-9]+")
  data$replicate <- replicate
  
  pheno_plot <- ggplot(data, aes(x = Generation, y = Phenotype)) +
    geom_line() +
    mytheme() +
    #ggtitle(basename(myphenofile)) +
    facet_wrap(~replicate, ncol = 4)
  
  return(list(replicate = replicate, plot = pheno_plot))



  
  #Fitness_plot <- ggplot(data, aes(x = Generation, y = meanFitness)) +
   # geom_line() +
  #  mytheme() +
  #  ggtitle(myphenofile)  # Added line
  
  #PhenoDens <- ggplot(data, aes(x = Phenotype)) +
  #  geom_histogram(bins = 200) +
  #  mytheme() +
   # ggtitle(paste0("Phenotype Density plot for ", basename(myphenofile)))  # Added line
  
  #FitDens <- ggplot(data, aes(x = meanFitness)) +
   # geom_density() +
  #  mytheme() +
  #  ggtitle(paste0("Fitness Density plot for ", basename(myphenofile)))  # Added line
  
  # Display the Phenotype Density plot
  #print(PhenoDens)
  
  # Display the Fitness plot
  #print(Fitness_plot)
  
  # Display the Fitness Density plot
  #print(FitDens)
}

########### Fourier Transform Function ########

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
    
    fft1 <- all_results %>%
      mutate(Frequency = 1/Frequency) %>% 
      filter(Frequency < 40) %>% 
      ggplot(aes(Frequency, spec, color = File)) + 
      geom_line(linewidth = 1.1) +
      theme(legend.position = "none") +
      xlab("Periodicity") +
      ylab("Spectral Density") +
      ggtitle(paste0("Spectral Density for ", basename(file_name)))  # Added line
    
    print(fft1)
  }
  
  return(all_results)
}
