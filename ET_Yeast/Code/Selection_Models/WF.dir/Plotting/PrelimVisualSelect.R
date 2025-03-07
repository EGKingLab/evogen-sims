## Here is the code used for plotting allele frequency, phenotype, and spectral density
mytheme <- function(){
  theme_set(theme_cowplot())+
    theme(axis.title = element_text(face = "bold"),
          axis.text = element_text(face = "italic"),
          plot.title = element_text(hjust = 0.5),
          legend.position = "none")
}

### Allele Frequency

allele_plots <- function(df, mygenomefile){
  replicate <- str_extract(mygenomefile, "(?<=genome)[0-9]+")
  df$replicate <- replicate

  Freq_plot <- ggplot(df, aes(x = Generation, y = Frequency, color = as.character(Position)))+
    geom_line(linewidth = 0.2)+
    mytheme()+
    facet_grid(~replicate)
  return(list(replicate = replicate, plot = Freq_plot))
}

plot_files <- function(dir_path, pattern, ncol){
  mygenomefiles <- list.files(dir_path, pattern = pattern, full.names = TRUE)
  plots_list <- list()
  for (mygenomefile in mygenomefiles){
    df <- read.csv(mygenomefile)
    result <- allele_plots(df, mygenomefile)
    replicate <- result$replicate
    plot <- result$plot
    
    if (!is.null(plots_list[[replicate]])){
      plots_list[[replicate]][[mygenomefile]] <- plot
    }else{
      plots_list[[replicate]] <- list(mygenomefile = plot)
    }
  }
  for (replicate in names(plots_list)){
    print(plot_grid(plotlist = plots_list[[replicate]], ncol = ncol))
  }
}


### Phenotypes

Pheno_plots <- function(df, myphenofile){
  replicate <- str_extract(myphenofile, "(?<=MeanPhenotypes)[0-9]+")
  df$replicate <- replicate
  Pheno_plot <- ggplot(df, aes(x = Generation, y = Phenotype))+
    geom_line(linewidth = 0.2)+
    mytheme()+
    facet_grid(~replicate)
  
  return(list(replicate = replicate, plot = Pheno_plot))
}

Pheno_files <- function(dir_path, pattern, ncol){
  myphenofiles <- list.files(dir_path, pattern = pattern, full.names = TRUE)
  plots_list <- list()
  for (myphenofile in myphenofiles){
    df <- read.csv(myphenofile)
    result <- Pheno_plots(df, myphenofile)
    replicate <- result$replicate
    plot <- result$plot
    
    plots_list[[replicate]] <- c(plots_list[[replicate]], list(plot))
  }
  for(replicate in names(plots_list)){
    if (length(plots_list[[replicate]]) >= ncol){
      print(plot_grid(plotlist = plots_list[[replicate]], ncol = ncol))
    }else{
      next
    }
  }
}
  
### Spectral Density

fft_fx <- function(dir_path, pattern, spectral_span){
  files <- list.files(dir_path, pattern = pattern, full.names = TRUE)
  files_df <- data.frame(file_name = files) %>%
    mutate(myloci_full = str_extract(file_name, "genome\\d+_n\\d+")) %>%
    mutate(myloci = str_remove(myloci_full, "genome\\d+_n"))

  files_by_myloci <- split(files_df, files_df$myloci)

  all_results <- data.frame()

  for(i in 1:nrow(files_df)){
    file_name <- files_df$file_name[i]
    myloci <- files_df$myloci[i]

    file <- read.csv(file_name) %>% dplyr::select(-Origin, -Effect) 

    data_wide <- pivot_wider(file[, c("Generation", "Position", "Frequency")],
                             names_from = Position, values_from = Frequency)
    data_wide[is.na(data_wide)] <- 0

    ff1 <- data_wide[,2:(as.numeric(files_df$myloci[i])+1)]
    ttest <- spectrum(ff1, spans = spectral_span, plot = FALSE)

    if(!is.matrix(ttest$spec)){
      out.spect <- ttest$spec
    }else{
      out.spect <- rowMeans(ttest$spec)
    }

    dd <- data.frame("Frequency" = rep(ttest$freq),
                     "file" = basename(file_name),
                     "spec" = out.spect,
                     "myloci" = myloci)

    all_results <- rbind(all_results, dd)
  }
  return(all_results)
}


