# Main function to run FFT analysis

run_fft_analysis <-
  function(dirpath,
           pattern = "^genome\\d+_n\\d+_H0.(1|8|5)(SD[1-4])?(Gen\\d+)?\\.csv$",
           spectrum_span = 2) {
  #library(tidyverse)
    library(tidyr)
    library(forcats)
    library(dplyr)
    library(tidyselect)
    library(ggplot2)
  library(future)
  library(future.apply)

    mytheme <- theme(axis.title = element_text(face = "bold"),
                     axis.text = element_text(face = "italic"),
                     plot.title = element_text(hjust = 0.5),
                     legend.position = "none")
  # Set up parallel backend
  plan(multisession) # Use multicore, if available

  # Function to process each file
  fft_fx <- function(file_name, spectrum_span) {
    file <- read.csv(file_name) %>%
      select(-Origin, -Effect)

    dat_wide <- pivot_wider(file, names_from = Position,
                            values_from = Frequency,
                            values_fill = list(Frequency = 0))
    myloci_full <- str_extract(file_name, "genome\\d+_n\\d+")
    myGen_full <- str_extract(file_name, "Gen\\d+") %>% replace_na("")
    myloci <- str_remove(myloci_full, "genome\\d+_n")
    myGen <- str_remove(myGen_full, "Gen")

    ff1 <- dat_wide[, 2:(as.numeric(myloci))] #  + 1

    ttest <- spectrum(ff1, spans = spectrum_span, plot = FALSE)

    out.spect <- if (is.matrix(ttest$spec)) {
      rowMeans(ttest$spec)
    } else {
      ttest$spec
    }

    dd <- data.frame(Frequency = rep(ttest$freq),
                     File = basename(file_name),
                     spec = out.spect,
                     myloci = myloci,
                     myGen = myGen)

    return(dd)
  }

  # Process files
  files <- list.files(dirpath, pattern = pattern, full.names = TRUE)
  all_results <- future_lapply(files, fft_fx, spectrum_span) %>% bind_rows()

  # Plot results
  nloci <- unique(all_results$myloci)
  fft1 <- list()

  for (loci in nloci) {
    fft1[[loci]] <- all_results %>%
      filter(myloci == loci) %>%
      mutate(Frequency = 1 / Frequency) %>%
      filter(Frequency < 30) %>%
      ggplot(aes(Frequency, spec, color = File)) +
      geom_line(linewidth = 1.1) +
      facet_wrap(~myGen, scales = "free_y") +
      mytheme
  }

  return(fft1)
  }

######################################################
######################################################

library(tidyr)
library(cowplot)
library(lubridate)
library(future)
library(future.apply)
theme_set(theme_cowplot())

# Define a custom theme function
mytheme <- theme(
  axis.title = element_text(face = "bold"),
  axis.text = element_text(face = "italic"),
  plot.title = element_text(hjust = 0.5),
  legend.position = "none"
)

run2_fft_analysis <- function(dirpath, pattern, spectrum_span = 2) {
  
  files <- list.files(dirpath, pattern = pattern, full.names = TRUE)
  
  # Extract information from file names
  files_df <- data.frame(file_name = files) %>%
    mutate(
      myloci = str_extract(file_name, "(?<=_n)\\d+"),
      myGen = str_extract(file_name, "(?<=Gen)\\d*"),
      mySD = str_extract(file_name, "(?<=SD)\\d*"),
      myHert = str_extract(file_name, "(?<=H)\\d*\\.?\\d*"),
      loci_Gen = paste(myloci, myGen, sep="_"),
      SD_H2 = paste(mySD, myHert, sep="_")
    )
  
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
    
    # Calculate the mean of the spectrum
    out.spect <- if (is.matrix(ttest$spec)) {
      rowMeans(ttest$spec)
    } else {
      ttest$spec
    }
    
    # Prepare results for plotting
    dd <- data.frame(
      Frequency = ttest$freq,
      File = basename(file_name),
      spec = out.spect,
      loci_Gen = loci_Gen,
      SD_H2 = SD_H2
    )
    
    all_results <- rbind(all_results, dd)
  }
  
  # Plot results with conditional filtering
  nloci <- unique(all_results$loci_Gen)
  fft1 <- list()
  for (loci in nloci) {
    threshold <- ifelse(grepl("_30$", loci), 70, ifelse(grepl("_20$", loci), 50, 40))
    fft1[[loci]] <- all_results %>%
      filter(loci_Gen == loci) %>%
      mutate(Frequency = 1 / Frequency) %>%
      filter(Frequency < threshold) %>%
      ggplot(aes(Frequency, spec, color = File)) + 
      geom_line(linewidth = 1.1) +
      facet_wrap(~SD_H2, scales = "free_y") +
      xlab("Periodicity") +
      ylab("Spectral Density") +
      mytheme
  }
  
  return(fft1)
}



# # Define the main function to process FFT analysis and generate plots
# run2_fft_analysis <- function(dirpath, spectrum_span = 2) {
#   library(tidyverse)
#   library(cowplot)
#   library(lubridate)
#   library(future)
#   library(future.apply)
#   theme_set(theme_cowplot())
#   
#   # Set up parallel backend, if desired
#   # plan(multisession)
#   mytheme <- function(){
#     theme_set(theme_cowplot())+
#       theme(axis.title = element_text(face = "bold"),
#             axis.text = element_text(face = "italic"),
#             plot.title = element_text(hjust = 0.5),
#             legend.position = "none")
#   }
#   # Pattern to match all specified file formats
#   pattern <- "^genome[1-3]_n\\d+_H0\\.(1|8|5)(SD(1|4))?(Gen(10|30))?\\.csv$"
#   
#   files <- list.files(dirpath, pattern = pattern, full.names = TRUE)
#   
#   # Extract information from file names and store in a data frame
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
#     if (ncol(file) < 3) { # Check there are enough columns
#       next
#     }
#     
#     ff1 <- file[, -1] # Exclude non-FFT columns
#     
#     # Apply FFT and calculate spectrum
#     ttest <- spectrum(ff1, spans = spectrum_span, plot = FALSE)
#     
#     if (!is.matrix(ttest$spec)) {
#       next # Skip non-matrix results
#     }
#     
#     out.spect <- rowMeans(ttest$spec)
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
#   # Plot results
#   nloci <- unique(all_results$loci_Gen)
#   fft1 <- list()
#   for (loci in nloci) {
#     fft1[[loci]] <- all_results %>%
#       filter(loci_Gen == loci) %>%
#       mutate(Frequency = 1 / Frequency) %>%
#       filter(Frequency < 40) %>%
#       ggplot(aes(Frequency, spec, color = File)) + 
#       geom_line(linewidth = 1.1) +
#       facet_wrap(~SD_H2, scales = "free_y") + 
#       #theme(legend.position = "none") +
#       xlab("Periodicity") +
#       ylab("Spectral Density")+
#       mytheme()
#   }
#   
#   return(fft1)
# }

