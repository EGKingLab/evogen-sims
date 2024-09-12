run2_fft_analysis <- function(dirpath, pattern, spectrum_span = 2) {
  library(tidyr)
  library(forcats)
  library(dplyr)
  library(tidyselect)
  library(ggplot2)
  library(cowplot)
  library(lubridate)
  library(future)
  library(future.apply)
  theme_set(theme_cowplot())
  
  # Define a custom theme function
  mytheme <- theme(axis.title = element_text(face = "bold"),
                   axis.text = element_text(face = "italic"),
                   plot.title = element_text(hjust = 0.5),
                   legend.position = "none",
                   axis.text.x = element_text(face = "bold", size = 12),
                   axis.text.y = element_text(face = "bold", size = 12))
  
  # Create the directory to save images if it doesn't exist
  last_dir <- basename(normalizePath(dirpath))
  save_dir <- file.path("myimages", last_dir)
  if (!dir.exists(save_dir)) {
    dir.create(save_dir, recursive = TRUE)
  }
  
  # Pattern to match all specified file formats
  files <- list.files(dirpath, pattern = pattern, full.names = TRUE)
  
  # Extract information from file names
  files_df <- data.frame(file_name = files) %>%
    mutate(
      myloci_full = str_extract(file_name, "genome\\d+_n\\d+"),
      myGen_full = str_extract(file_name, "Gen\\d*"),
      mySD_full = str_extract(file_name, "SD\\d*"),
      myHert_full = str_extract(file_name, "(?<=H)\\d*\\.?\\d*"),
      myloci = str_remove(myloci_full, "genome\\d+_n"),
      myGen = str_remove(myGen_full, "Gen"),
      mySD = str_remove(mySD_full, "SD"),
      myHert = myHert_full,
      loci_Gen = paste(myloci, myGen, sep="_"),
      SD_H2 = paste("SD =", mySD, "H =", myHert)
    ) %>%
    select(-myloci_full, -myGen_full, -mySD_full, -myHert_full)
  
  all_results <- data.frame()
  
  # Process each file
  for (i in seq_along(files_df$file_name)) {
    file_name <- files_df$file_name[i]
    loci_Gen <- files_df$loci_Gen[i]
    SD_H2 <- files_df$SD_H2[i]
    
    file <- read.csv(file_name) %>%
      select(-Origin, -Effect) %>%
      pivot_wider(names_from = Position, values_from = Frequency, 
                  values_fill = list(Frequency = 0))
    
    ff1 <- file[, -1] # Exclude non-FFT columns
    
    # Apply FFT and calculate spectrum
    ttest <- spectrum(ff1, spans = spectrum_span, plot = FALSE)
    
    out.spect <- if (is.matrix(ttest$spec)) {
      rowMeans(ttest$spec)
    } else {
      ttest$spec
    }
    
    # Prepare results for plotting
    dd <- data.frame(
      Frequency = rep(ttest$freq),
      File = basename(file_name),
      spec = out.spect,
      loci_Gen = loci_Gen,
      SD_H2 = SD_H2
    )
    
    all_results <- rbind(all_results, dd)
  }
  
  # Plot results with conditional filtering based on the gener interval
  nloci <- unique(all_results$loci_Gen)
  fft1 <- list()
  for (loci in nloci) {
    threshold <- ifelse(grepl("_30$", loci), 80, 
                        ifelse(grepl("_20$", loci), 50, 
                               ifelse(grepl("_10$", loci), 40, 80)))
    fft1[[loci]] <- all_results %>%
      filter(loci_Gen == loci) %>%
      mutate(Frequency = 1 / Frequency, spec = round((spec), 5)) %>%
      filter(Frequency < threshold) %>%
      ggplot(aes(Frequency, spec, color = File)) + 
      geom_line(linewidth = 1.1) +
      facet_wrap(~SD_H2) + # scales = "free_y"
      xlab("Periodicity") +
      ylab("Spectral Density") +
      mytheme
    
    # Save the plot
    ggsave(filename = file.path(save_dir, paste0("Spectral_", loci, ".png")), plot = fft1[[loci]])
  }
  
  return(fft1)
}
##############################################################################
#################### Previously used code that does not save plots ###########
##############################################################################
# run2_fft_analysis <- function(dirpath, pattern, spectrum_span = 2) {
#   library(tidyr)
#   library(forcats)
#   library(dplyr)
#   library(tidyselect)
#   library(ggplot2)
#   library(cowplot)
#   library(lubridate)
#   library(future)
#   library(future.apply)
#   theme_set(theme_cowplot())
#   
#   # Define a custom theme function
#   mytheme <- theme(axis.title = element_text(face = "bold"),
#                    axis.text = element_text(face = "italic"),
#                    plot.title = element_text(hjust = 0.5),
#                    legend.position = "none")
#   
#   # Pattern to match all specified file formats
#  #pattern <- "^genome[1-3]_n\\d+_H0\\.(1|5|8)(SD(1|4))?(Gen(10|30|20))?\\.csv$"
#   
#   files <- list.files(dirpath, pattern = pattern, full.names = TRUE)
#   
#   # Extract information from file names
#   files_df <- data.frame(file_name = files) %>%
#     mutate(
#       myloci_full = str_extract(file_name, "genome\\d+_n\\d+"),
#       myGen_full = str_extract(file_name, "Gen\\d*"),
#       mySD_full = str_extract(file_name, "SD\\d*"),
#       myHert_full = str_extract(file_name, "(?<=H)\\d*\\.?\\d*"),
#       myloci = str_remove(myloci_full, "genome\\d+_n"),
#       myGen = str_remove(myGen_full, "Gen"),
#       mySD = str_remove(mySD_full, "SD"),
#       myHert = myHert_full,
#       loci_Gen = paste(myloci, myGen, sep="_"),
#       SD_H2 = paste(mySD, myHert, sep="_")
#     ) %>%
#     select(-myloci_full, -myGen_full, -mySD_full, -myHert_full)
#   
#   all_results <- data.frame()
#   
#   # Process each file
#   for (i in seq_along(files_df$file_name)) {
#     file_name <- files_df$file_name[i]
#     loci_Gen <- files_df$loci_Gen[i]
#     SD_H2 <- files_df$SD_H2[i]
#     
#     file <- read.csv(file_name) %>%
#       select(-Origin, -Effect) %>%
#       pivot_wider(names_from = Position, values_from = Frequency, 
#                   values_fill = list(Frequency = 0))
#     
#     # if (ncol(file) < 3) {
#     #   next
#     # }
#     # 
#     ff1 <- file[, -1] # Exclude non-FFT columns
#     
#     # Apply FFT and calculate spectrum
#     ttest <- spectrum(ff1, spans = spectrum_span, plot = FALSE)
# 
#    # This helps us deal with a single locus vs other loci
#     
#     out.spect <- if (is.matrix(ttest$spec)) {
#       rowMeans(ttest$spec)
#     } else {
#       ttest$spec
#     }
#     
#     # Prepare results for plotting
#     dd <- data.frame(
#       Frequency = rep(ttest$freq),
#       File = basename(file_name),
#       spec = out.spect,
#       loci_Gen = loci_Gen,
#       SD_H2 = SD_H2
#     )
#     
#     all_results <- rbind(all_results, dd)
#   }
#   
#   # Plot results with conditional filtering based on the gener interval
#   nloci <- unique(all_results$loci_Gen)
#   fft1 <- list()
#   for (loci in nloci) {
#     threshold <- ifelse(grepl("_30$", loci), 80, 
#                         ifelse(grepl("_20$", loci), 50, 
#                                ifelse(grepl("_10$", loci), 40, 80)))
#     fft1[[loci]] <- all_results %>%
#       filter(loci_Gen == loci) %>%
#       mutate(Frequency = 1 / Frequency) %>%
#       filter(Frequency < threshold) %>%
#       ggplot(aes(Frequency, spec, color = File)) + 
#       geom_line(size = 1.1) +
#       facet_wrap(~SD_H2, scales = "free_y") +
#       xlab("Periodicity") +
#       ylab("Spectral Density") +
#       mytheme
#   }
#   
#   return(fft1)
# }

##########################################################
# library(tidyverse)
# library(cowplot)
# library(lubridate)
# library(future)
# library(future.apply)
# 
# run3_fft_analysis <- function(dirpath, spectrum_span = 2) {
#   theme_set(theme_cowplot())
#   plan(multisession) # Enable parallel processing
#   
#   # Custom theme
#   mytheme <- theme(axis.title = element_text(face = "bold"),
#                    axis.text = element_text(face = "italic"),
#                    plot.title = element_text(hjust = 0.5),
#                    legend.position = "none")
#   
#   # Pattern to match specified file formats
#   pattern <- "^genome[1-3]_n\\d+_H0\\.(1|5|8)(SD(1|4))?(Gen(10|20|30))?\\.csv$"
#   files <- list.files(dirpath, pattern = pattern, full.names = TRUE)
#   
#   # Parallel processing of files
#   process_file <- function(file_path) {
#     # Extract file info
#     file_info <- str_extract_all(file_path, "genome[1-3]_n\\d+|Gen\\d*|SD\\d*|H\\d*\\.?\\d*") %>% 
#       unlist() %>%
#       set_names(c("loci", "Gen", "SD", "Hertz"))
#     
#     # Adjust names and combine
#     file_info <- map(file_info, ~str_remove(.x, "genome\\d+_n|Gen|SD|H")) %>%
#       as.list()
#     file_info$loci_Gen <- paste(file_info$loci, file_info$Gen, sep = "_")
#     file_info$SD_H2 <- paste(file_info$SD, file_info$Hertz, sep = "_")
#     
#     # Read and prepare data
#     data <- read.csv(file_path) %>%
#       select(-Origin, -Effect) %>%
#       pivot_wider(names_from = Position, values_from = Frequency, values_fill = list(Frequency = 0))
#     
#     # Skip if not enough columns
#     if (ncol(data) < 3) return(NULL)
#     
#     # Exclude non-FFT columns and apply FFT
#     ff1 <- data[, -1]
#     ttest <- spectrum(ff1, spans = spectrum_span, plot = FALSE)
#     
#     if (!is.matrix(ttest$spec)) return(NULL)
#     
#     out.spect <- rowMeans(ttest$spec)
#     
#     # Prepare results for plotting
#     data.frame(
#       Frequency = rep(ttest$freq),
#       File = basename(file_path),
#       spec = out.spect,
#       loci_Gen = file_info$loci_Gen,
#       SD_H2 = file_info$SD_H2
#     )
#   }
#   
#   all_results <- future_lapply(files, process_file) %>% bind_rows
#   
#   # Adjust threshold based on loci_Gen after converting Frequency to 1/Frequency
#   fft_plots <- lapply(unique(all_results$loci_Gen), function(loci) {
#     subset_data <- all_results %>%
#       filter(loci_Gen == loci) %>%
#       mutate(Frequency = 1 / Frequency)
#     
#     threshold <- ifelse(grepl("_30$", loci), 70, ifelse(grepl("_20$", loci), 50, 30))
#     
#     subset_data %>%
#       filter(Frequency < threshold) %>%
#       ggplot(aes(x = Frequency, y = spec, color = File)) +
#       geom_line(linewidth = 1.1) +
#       facet_wrap(~SD_H2, scales = "free_y") +
#       xlab("Periodicity (1/Frequency)") +
#       ylab("Spectral Density") +
#       mytheme
#   })
#   
#   return(fft_plots)
# }
