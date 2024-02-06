# Main function to run FFT analysis
run_fft_analysis <- 
  function(dirpath, 
           pattern = "^genome\\d+_n\\d+_H0.(1|8|5)(SD[1-4])?(Gen\\d+)?\\.csv$", 
           spectrum_span = 2) {
  library(tidyverse)
  library(future)
  library(future.apply)
  
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
    
    ff1 <- dat_wide[, 2:(as.numeric(myloci) + 1)]
    
    ttest <- spectrum(ff1, spans = spectrum_span, plot = FALSE)
    
    if (!is.matrix(ttest$spec)) {
      return(NULL) # Skip this file if ttest$spec is not a matrix
    }
    
    out.spect <- rowMeans(ttest$spec)
    
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
      filter(Frequency < 40) %>%
      ggplot(aes(Frequency, spec, color = File)) + 
      geom_line(linewidth = 1.1) +
      facet_wrap(~myGen, scales = "free_y") +
      theme(legend.position = "none") +
      xlab("Periodicity") +
      ylab("Spectral Density")
  }
  
  return(fft1)
}
