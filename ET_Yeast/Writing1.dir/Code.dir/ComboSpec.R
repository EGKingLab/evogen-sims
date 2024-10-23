run_fft_analysis <- function(folder_path, pattern, spectrum_span = 2) {
  # Load necessary libraries
  library(parallel)
  library(doParallel)
  library(foreach)
  library(dplyr)
  library(ggplot2)
  library(cowplot)
  library(stringr)
  
  theme_set(theme_cowplot())
  
  # Define a custom theme function with bold facet labels and centered plot title
  mytheme <- theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1, size = 15, face = "bold"),
    axis.text.y = element_text(size = 15, face = "bold"),
    axis.line = element_line(linewidth = 1),
    axis.title = element_text(size = 15, face = "bold"),
    strip.text = element_text(size = 15, face = "bold"),
    plot.title = element_text(hjust = 0.5, size = 15, face = "bold"),  # Centered title
    panel.spacing = unit(1, "lines")
  ) +
    theme(panel.grid = element_blank())
  
  # Create the directory to save images if it doesn't exist
  last_dir <- basename(normalizePath(folder_path))
  save_dir <- file.path("myimages", last_dir)
  if (!dir.exists(save_dir)) {
    dir.create(save_dir, recursive = TRUE)
  }
  
  # Get a list of all replicate files matching the pattern
  file_list <- list.files(
    path = folder_path,
    pattern = pattern,
    full.names = TRUE
  )
  
  # Check if files are found
  if (length(file_list) == 0) {
    stop("No files found with the given pattern.")
  }
  
  # Create a data frame with filenames and extract parameter information
  files_df <- data.frame(file_name = file_list) %>%
    mutate(
      base_name = basename(file_name),
      genome = str_extract(base_name, "genome\\d+") %>% str_remove("genome"),
      n = str_extract(base_name, "_n\\d+") %>% str_remove("_n"),
      H = str_extract(base_name, "_H\\d+\\.\\d*") %>% str_remove("_H"),
      SD = str_extract(base_name, "SD\\d+(\\.\\d+)?"),  # Adjusted regex
      SD = ifelse(is.na(SD), NA, str_remove(SD, "SD")),
      Gen = str_extract(base_name, "Gen\\d+"),
      Gen = ifelse(is.na(Gen), NA, str_remove(Gen, "Gen"))
    ) %>%
    mutate(
      param_combination = paste0(
        "n_", n, "_H", H,
        ifelse(!is.na(SD), paste0("_SD", SD), ""),
        ifelse(!is.na(Gen), paste0("_Gen", Gen), "")
      )
    )
  
  # Set up parallel backend
  num_cores <- detectCores()
  num_cores_to_use <- max(1, num_cores - 16)  # Leave at least 16 cores unused
  
  cl <- makeCluster(num_cores_to_use)
  registerDoParallel(cl)
  
  # Read and combine data from all replicate files in parallel
  combined_data_list <- foreach(i = seq_along(files_df$file_name), .packages = c("dplyr", "stringr")) %dopar% {
    file <- files_df$file_name[i]
    params <- files_df[i, ]
    df <- read.csv(file)
    data <- df %>%
      mutate(
        allele_id = as.factor(Position),
        generation = as.numeric(Generation),
        allele_frequency = Frequency,
        sample_size = 10000,
        file = params$base_name,
        genome = params$genome,
        n = as.numeric(params$n),
        H = as.numeric(params$H),
        SD = as.numeric(params$SD),
        Gen = as.numeric(params$Gen),
        param_combination = params$param_combination
      ) %>%
      dplyr::select(
        allele_id, generation, allele_frequency, sample_size, file,
        genome, n, H, SD, Gen, param_combination
      )
    data
  }
  combined_data <- bind_rows(combined_data_list)
  
  # Calculate the mean allele frequency for each grouping
  mean_allele_data <- combined_data %>%
    group_by(generation, n, genome, H, SD, Gen, param_combination) %>%
    summarize(mean_allele_frequency = mean(allele_frequency), .groups = "drop")
  
  # Create time series objects
  ts_data_list <- mean_allele_data %>%
    group_by(n, genome, H, SD, Gen, param_combination) %>%
    summarize(ts_data = list(ts(mean_allele_frequency, start = min(generation), frequency = 1)), .groups = 'drop')
  
  # Compute spectral density estimates in parallel
  spec_results_list <- foreach(i = 1:nrow(ts_data_list), .packages = c("stats")) %dopar% {
    ts_data <- ts_data_list$ts_data[[i]]
    spectrum(ts_data, spans = spectrum_span, plot = FALSE)
  }
  
  # Prepare data for plotting in parallel
  plot_data_list <- foreach(i = seq_along(spec_results_list), .packages = c("dplyr", "stringr")) %dopar% {
    spec_result <- spec_results_list[[i]]
    params <- ts_data_list[i, ]
    # Compute period and filter based on threshold
    period <- 1 / spec_result$freq
    # Construct n_Gen variable
    n_Gen <- ifelse(!is.na(params$Gen), paste0("n_", params$n, "_Gen_", params$Gen), paste0("n_", params$n))
    # Set threshold based on n_Gen value using grepl
    threshold <- ifelse(grepl("_30$", n_Gen), 80,
                        ifelse(grepl("_20$", n_Gen), 50,
                               ifelse(grepl("_10$", n_Gen), 40, 80)))
    valid_indices <- which(period < threshold)
    data.frame(
      Period = period[valid_indices],
      Spectrum = spec_result$spec[valid_indices],
      n = params$n,
      genome = params$genome,
      H = params$H,
      SD = params$SD,
      Gen = params$Gen,
      param_combination = params$param_combination,
      n_Gen = n_Gen
    )
  }
  
  # Stop cluster
  stopCluster(cl)
  
  # Combine all plot data
  plot_data <- bind_rows(plot_data_list)
  
  # Prepare data for faceting
  plot_data <- plot_data %>%
    mutate(
      H_label = paste0("H = ", H),
      SD_label = ifelse(!is.na(SD), paste0("SD = ", SD), "No SD"),
      Gen_label = ifelse(!is.na(Gen), paste0("Gen = ", Gen), "No Gen"),
      n_label = paste0("n_", n)
    )
  
  # Generate plots for each unique combination of n and Gen
  unique_n_Gen_values <- unique(plot_data$n_Gen)
  
  fft_plots <- list()
  
  for (n_Gen_value in unique_n_Gen_values) {
    plot_df <- plot_data %>%
      filter(n_Gen == n_Gen_value)
    
    # Determine the faceting variables
    if (all(plot_df$SD_label == "No SD")) {
      # If SD is missing, facet only by H_label
      facet_formula <- ~ H_label
    } else {
      # Facet by SD_label and H_label
      facet_formula <- H_label ~ SD_label
    }
    
    p <- ggplot(plot_df, aes(x = Period, y = Spectrum, color = genome, group = genome)) +
      geom_line(linewidth = 1.1) +
      facet_wrap(facet_formula, scales = "free") +  
      labs(
        #title = paste("Spectral Analysis for", n_Gen_value),
        x = "Period (Generations)",
        y = "Spectral Density"
      ) +
      mytheme
    
    fft_plots[[n_Gen_value]] <- p
    
    # Save the plot with n and Gen in filename
    filename <- paste0("Spectral_", n_Gen_value, ".png")
    ggsave(filename = file.path(save_dir, filename), plot = p)
  }
  
  return(fft_plots)
}

# run_fft_analysis <- function(folder_path, pattern, spectrum_span = 2) {
#   theme_set(theme_cowplot())
#   
#   # Define a custom theme function with bold facet labels
#   mytheme <- theme(
#     legend.position = "none",
#     axis.text.x = element_text(angle = 45, hjust = 1, size = 15, face = "bold"),
#     axis.text.y = element_text(size = 15, face = "bold"),
#     axis.line = element_line(linewidth = 1),
#     axis.title = element_text(size = 15, face = "bold"),
#     strip.text = element_text(size = 15, face = "bold"),
#     panel.spacing = unit(1, "lines")
#   ) +
#     theme(panel.grid = element_blank())
#   
#   # Create the directory to save images if it doesn't exist
#   last_dir <- basename(normalizePath(folder_path))
#   save_dir <- file.path("myimages", last_dir)
#   if (!dir.exists(save_dir)) {
#     dir.create(save_dir, recursive = TRUE)
#   }
#   
#   # Get a list of all replicate files matching the pattern
#   file_list <- list.files(
#     path = folder_path,
#     pattern = pattern,
#     full.names = TRUE
#   )
#   
#   # Check if files are found
#   if (length(file_list) == 0) {
#     stop("No files found with the given pattern.")
#   }
#   
#   # Create a data frame with filenames and extract parameter information
#   files_df <- data.frame(file_name = file_list) %>%
#     mutate(
#       base_name = basename(file_name),
#       genome = str_extract(base_name, "genome\\d+") %>% str_remove("genome"),
#       n = str_extract(base_name, "_n\\d+") %>% str_remove("_n"),
#       H = str_extract(base_name, "_H\\d+\\.\\d*") %>% str_remove("_H"),
#       SD = str_extract(base_name, "SD\\d+(\\.\\d+)?"),  # Adjusted regex
#       SD = ifelse(is.na(SD), NA, str_remove(SD, "SD")),
#       Gen = str_extract(base_name, "Gen\\d+"),
#       Gen = ifelse(is.na(Gen), NA, str_remove(Gen, "Gen"))
#     ) %>%
#     mutate(
#       param_combination = paste0(
#         "n_", n, "_H", H,
#         ifelse(!is.na(SD), paste0("_SD", SD), ""),
#         ifelse(!is.na(Gen), paste0("_Gen", Gen), "")
#       )
#     )
#   
#   combined_data <- data.frame()
#   
#   # Read and combine data from all replicate files
#   for (i in seq_along(files_df$file_name)) {
#     file <- files_df$file_name[i]
#     params <- files_df[i, ]
#     df <- read.csv(file)
#     data <- df %>%
#       mutate(
#         allele_id = as.factor(Position),
#         generation = as.numeric(Generation),
#         allele_frequency = Frequency,
#         sample_size = 10000,
#         file = params$base_name,
#         genome = params$genome,
#         n = as.numeric(params$n),
#         H = as.numeric(params$H),
#         SD = as.numeric(params$SD),
#         Gen = as.numeric(params$Gen),
#         param_combination = params$param_combination
#       ) %>%
#       dplyr::select(
#         allele_id, generation, allele_frequency, sample_size, file,
#         genome, n, H, SD, Gen, param_combination
#       )
#     combined_data <- rbind(combined_data, data)
#   }
#   
#   # Calculate the mean allele frequency for each grouping
#   mean_allele_data <- combined_data %>%
#     group_by(generation, n, genome, H, SD, Gen, param_combination) %>%
#     summarize(mean_allele_frequency = mean(allele_frequency), .groups = "drop")
#   
#   # Create time series objects
#   ts_data_list <- mean_allele_data %>%
#     group_by(n, genome, H, SD, Gen, param_combination) %>%
#     do(ts_data = ts(.$mean_allele_frequency, start = min(.$generation), frequency = 1)) %>%
#     ungroup()
#   
#   # Compute spectral density estimates
#   spec_results <- lapply(ts_data_list$ts_data, function(ts_data) {
#     spectrum(ts_data, spans = spectrum_span, plot = FALSE)
#   })
#   
#   # Prepare data for plotting
#   plot_data_list <- lapply(seq_along(spec_results), function(i) {
#     spec_result <- spec_results[[i]]
#     params <- ts_data_list[i, ]
#     # Compute period and filter based on threshold
#     period <- 1 / spec_result$freq
#     # Set threshold based on n value (adjust as needed)
#     n_value <- params$n
#     threshold <- ifelse(n_value == 30, 80,
#                         ifelse(n_value == 20, 50,
#                                ifelse(n_value == 10, 40, 80)))
#     valid_indices <- which(period < threshold)
#     data.frame(
#       Period = period[valid_indices],
#       Spectrum = spec_result$spec[valid_indices],
#       n = n_value,
#       genome = params$genome,
#       H = params$H,
#       SD = params$SD,
#       Gen = params$Gen,
#       param_combination = params$param_combination
#     )
#   })
#   
#   # Combine all plot data
#   plot_data <- bind_rows(plot_data_list)
#   
#   # Prepare data for faceting
#   plot_data <- plot_data %>%
#     mutate(
#       H_label = paste0("H = ", H),
#       SD_label = ifelse(!is.na(SD), paste0("SD = ", SD), "No SD"),
#       Gen_label = ifelse(!is.na(Gen), paste0("Gen = ", Gen), "No Gen"),
#       n_label = paste0("n_", n)
#     )
#   
#   # Generate plots for each unique combination of n and Gen
#   plot_data <- plot_data %>%
#     mutate(
#       n_Gen = ifelse(!is.na(Gen), paste0("n_", n, "_Gen_", Gen), paste0("n_", n))
#     )
#   
#   unique_n_Gen_values <- unique(plot_data$n_Gen)
#   
#   fft_plots <- list()
#   
#   for (n_Gen_value in unique_n_Gen_values) {
#     plot_df <- plot_data %>%
#       filter(n_Gen == n_Gen_value)
#     
#     # Determine the faceting variables
#     if (all(plot_df$SD_label == "No SD")) {
#       # If SD is missing, facet only by H_label
#       facet_formula <- ~ H_label
#     } else {
#       # Facet by SD_label and H_label
#       facet_formula <- SD_label ~ H_label
#     }
#     
#     p <- ggplot(plot_df, aes(x = Period, y = Spectrum, color = genome, group = genome)) +
#       geom_line(linewidth = 1.1) +
#       facet_grid(facet_formula, scales = "free_y") +
#       labs(
#         title = paste("Spectral Analysis for", n_Gen_value),
#         x = "Period (Generations)",
#         y = "Spectral Density"
#       ) +
#       mytheme
#     
#     fft_plots[[n_Gen_value]] <- p
#     
#     # Save the plot with n and Gen in filename
#     filename <- paste0("Spectral_", n_Gen_value, ".png")
#     ggsave(filename = file.path(save_dir, filename), plot = p)
#   }
#   
#   return(fft_plots)
# }

#########################################################
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
  
  # Define a custom theme function with bold facet labels
  mytheme <- theme(legend.position = "none",
                   axis.text.x = element_text(angle = 45, hjust = 1, size = 15, face = "bold"),
                   axis.text.y = element_text(size = 15, face = "bold"),
                   axis.line = element_line(linewidth = 3),
                   axis.title = element_text(size = 15, face = "bold"),
                   strip.text = element_text(size = 15, face = "bold"),
                   panel.spacing = unit(3, "lines")) +
    theme(panel.grid = element_blank())
    
    # theme(axis.title = element_text(face = "bold"),
    #                axis.text = element_text(face = "italic"),
    #                plot.title = element_text(hjust = 0.5),
    #                legend.position = "none",
    #                axis.text.x = element_text(face = "bold", size = 12),
    #                axis.text.y = element_text(face = "bold", size = 12),
    #                strip.text = element_text(face = "bold", size = 14))  # Bold facet labels
  
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
    dplyr::select(-myloci_full, -myGen_full, -mySD_full, -myHert_full)
  
  all_results <- data.frame()
  
  # Process each file
  for (i in seq_along(files_df$file_name)) {
    file_name <- files_df$file_name[i]
    loci_Gen <- files_df$loci_Gen[i]
    SD_H2 <- files_df$SD_H2[i]
    
    file <- read.csv(file_name) %>%
      dplyr::select(-Origin, -Effect) %>%
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
      mutate(Frequency = 1 / Frequency) %>%
      filter(Frequency < threshold) %>%
      ggplot(aes(Frequency, spec, color = File)) + 
      geom_line(linewidth = 1.1) +
      facet_wrap(~SD_H2, scales = "free_y") + 
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
