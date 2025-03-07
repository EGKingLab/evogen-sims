### SPECTRAL ANALYSIS

fft_fx <- function(dir_path, pattern, spectrum_span){
  files <- list.files(dir_path, pattern = pattern, full.names = TRUE)

  files_df <- data.frame(file_name = files) %>%
    mutate(myloci_full = str_extract(file_name, "genome\\d+_n\\d+")) %>%
    mutate(myloci = str_remove(myloci_full, "genome\\d+_n"))

  files_by_myloci <- split(files_df, files_df$myloci)

  all_results <- data.frame()

  # Iterate files in the data frame
  for (i in 1:nrow(files_df)) {
    file_name <- files_df$file_name[i]
    myloci <- files_df$myloci[i]

    file <- read.csv(file_name) %>% dplyr::select(-Origin, -Effect)

    dat_wide <- pivot_wider(file[,c("Generation","Position","Frequency")], names_from=Position, values_from=Frequency)
    dat_wide[is.na(dat_wide)] <- 0
    ff1 <- dat_wide[,2:(as.numeric(files_df$myloci[i])+1)]

    # Apply spectrum and calculate row means
    ttest <- spectrum(ff1, spans=spectrum_span, plot=FALSE)

    # If ttest$spec is not a matrix, skip this file
    if (!is.matrix(ttest$spec)) {
      next
    }

    out.spect <- rowMeans(ttest$spec)

    # Prepare data frame for plot
    dd <- data.frame("Frequency" = rep(ttest$freq),
                     "File" = basename(file_name),
                     "spec" = out.spect,
                     "myloci" = myloci)

    # Combine with previous results
    all_results <- rbind(all_results, dd)
  }

  # Return the result
  return(all_results)
}

# fft_fx <- function(dir_path, pattern, spectrum_span){
#   
#   print(paste("Directory path:", dir_path))
#   print(paste("File pattern:", pattern))
#   
#   files <- list.files(dir_path, pattern = pattern, full.names = TRUE)
#   
#   print("Found files:")
#   print(files)
#   
#   files_df <- data.frame(file_name = files) %>%
#     mutate(myloci_full = str_extract(file_name, "genome\\d+_n\\d+")) %>%
#     mutate(myloci = str_remove(myloci_full, "genome\\d+_n"))
#   
#   files_by_myloci <- split(files_df, files_df$myloci)
#   
#   all_results <- data.frame()
#   
#   # Iterate files in the data frame
#   for (i in 1:nrow(files_df)) {
#     file_name <- files_df$file_name[i]
#     
#     print(paste("Reading file:", file_name))
#     
#     myloci <- files_df$myloci[i]
#     
#     file <- read.csv(file_name) %>% dplyr::select(-Origin, -Effect)
#     
#     dat_wide <- pivot_wider(file[,c("Generation","Position","Frequency")], names_from=Position, values_from=Frequency)
#     dat_wide[is.na(dat_wide)] <- 0 
#     ff1 <- dat_wide[,2:(as.numeric(files_df$myloci[i])+1)]
#     
#     # Apply spectrum and calculate row means
#     ttest <- spectrum(ff1, spans=spectrum_span, plot=FALSE)
#     
#     # If ttest$spec is not a matrix, skip this file
#     if (!is.matrix(ttest$spec)) {
#       next
#     }
#     
#     out.spect <- rowMeans(ttest$spec)
#     
#     # Prepare data frame for plot
#     dd <- data.frame("Frequency" = rep(ttest$freq), 
#                      "File" = basename(file_name),
#                      "sel" = rep("FS", each = 1000),
#                      "spec" = out.spect,
#                      "myloci" = myloci)  
#     
#     # Combine with previous results
#     all_results <- rbind(all_results, dd)
#   }
#   
#   # Return the result
#   return(all_results)
# }

### ALLELE FREQUENCE

allele_plots <- function(df, mygenomefile) {
  
  replicate <- str_extract(mygenomefile, "(?<=genome)[0-9]+")
  
  df$replicate <- replicate
  
  Freq_plot <- ggplot(df, aes(x = Generation, y = Frequency, color = as.character(Position))) +
    geom_line(linewidth = 0.5) +
    mytheme() +
    facet_grid(~replicate)
  
  return(list(replicate = replicate, plot = Freq_plot))
}

plot_files <- function(dir_path, pattern, ncol) {
  mygenomefiles <- list.files(dir_path, pattern = pattern, full.names = TRUE)
  
  plots_list <- list()
  
  for (mygenomefile in mygenomefiles){
    df <- read.csv(mygenomefile)
    result <- allele_plots(df, mygenomefile)
    replicate <- result$replicate
    plot <- result$plot
    
    if (!is.null(plots_list[[replicate]])) {
      plots_list[[replicate]][[mygenomefile]] <- plot
    } else {
      plots_list[[replicate]] <- list(mygenomefile = plot)
    }
  }
  
  for (replicate in names(plots_list)) {
    print(plot_grid(plotlist = plots_list[[replicate]], ncol = ncol))  
  }
}

### PHENOTYPES
# Pheno_plots <- function(df, myphenofile) {
#   
#   replicate <- str_extract(myphenofile, "(?<=MeanPhenotypes)[0-9]+")
#   
#   df$replicate <- replicate
#   
#   Pheno_plot <- ggplot(df, aes(x = Generation, y = Phenotype)) +
#     geom_line(linewidth = 0.5) +
#     mytheme() +
#     facet_grid(~replicate)
#   
#   return(list(replicate = replicate, plot = Pheno_plot))
# }
# 
# pheno_files <- function(dir_path, pattern, ncol) {
#   myphenofiles <- list.files(dir_path, pattern = pattern, full.names = TRUE)
#   
#   plots_list <- list()
#   
#   for (myphenofile in myphenofiles){
#     df <- read.csv(myphenofile)
#     result <- Pheno_plots(df, myphenofile)
#     replicate <- result$replicate
#     plot <- result$plot
#     
#     if (!is.null(plots_list[[replicate]])) {
#       plots_list[[replicate]][[paste("plot_", replicate)]] <- plot
#     } else {
#       plots_list[[replicate]] <- list(paste("plot_", replicate) = plot)
#     }
#   }
#   
#   for (replicate in names(plots_list)) {
#     print(plot_grid(plotlist = plots_list[[replicate]], ncol = ncol))  
#   }
# }


Pheno_plots <- function(df, myphenofile) {
  
  replicate <- str_extract(myphenofile, "(?<=MeanPhenotypes)[0-9]+")
  
  df$replicate <- replicate
  
  Pheno_plot <- ggplot(df, aes(x = Generation, y = Phenotype)) +
    geom_line(linewidth = 0.5) +
    mytheme() +
    facet_grid(~replicate)
  
  return(list(replicate = replicate, plot = Pheno_plot))
}

pheno_files <- function(dir_path, pattern, ncol) {
  myphenofiles <- list.files(dir_path, pattern = pattern, full.names = TRUE)
  
  plots_list <- list()
  
  for (myphenofile in myphenofiles){
    df <- read.csv(myphenofile)
    result <- Pheno_plots(df, myphenofile)
    replicate <- result$replicate
    plot <- result$plot
    
    # Add each plot to the plot list
    plots_list[[replicate]] <- c(plots_list[[replicate]], list(plot))
  }
  
  # Print each group of plots
  for (replicate in names(plots_list)) {
    # Check if there are enough plots for the specified number of columns
    if (length(plots_list[[replicate]]) >= ncol) {
      print(plot_grid(plotlist = plots_list[[replicate]], ncol = ncol))
    } else {
      print(paste("Not enough plots for the specified number of columns for replicate", replicate))
    }
  }
}
