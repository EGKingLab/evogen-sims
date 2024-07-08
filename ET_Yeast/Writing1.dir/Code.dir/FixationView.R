############ Libraries and theme function ########
rm(list = ls())
library(dplyr)
library(tidyr)
library(stringr)
library(plotly)
library(patchwork)
library(cowplot)
library(purrr)
library(doParallel)
library(data.table)
library(ggplot2)

theme_set(theme_cowplot())
mytheme <- function(){
  theme(axis.title = element_text(face = "bold"),
        axis.text = element_text(face = "italic"),
        plot.title = element_text(hjust = 0.5),
        legend.position = "none")
}

#################### Data Processing and plotting ###############

# Directory path and pattern
path <- "~/YeastProj.dir/evogen-sims/ET_Yeast/output.dir/Selection_Models/WF.dir/SinFSGen.dir/"
pattern <- "^genome\\d+_n(10|100)_H0.(1|8)(SD[1|4])?(Gen10|30)?\\.csv$"

# Function to process files in chunks
process_files_in_chunks <- function(files, chunk_size, num_cores) {
  combined_data <- list()
  
  for (i in seq(1, length(files), by = chunk_size)) {
    chunk <- files[i:min(i + chunk_size - 1, length(files))]
    
    cl <- makeCluster(num_cores)
    registerDoParallel(cl)
    
    data_list <- foreach(file = chunk, .packages = c("dplyr", "stringr", "data.table"), .errorhandling = 'remove') %dopar% {
      repl <- as.numeric(str_extract(file, "(?<=genome)\\d+"))
      herit <- as.numeric(str_extract(file, "(?<=H)0\\.\\d+"))
      loci <- as.numeric(str_extract(file, "(?<=_n)\\d+"))
      sd <- as.numeric(str_extract(file, "(?<=SD)\\d+"))
      gen <- as.numeric(str_extract(file, "(?<=Gen)\\d+"))
      
      data <- fread(file) %>%
        select(-Origin) %>%
        mutate(Position = factor(Position), repl = repl,
               herit = herit, loci = loci, sd = sd, gen = gen,
               h2_sd = paste("h2 = ", herit, " ", "sd = ", sd, sep = ""),
               loci_gen = paste("loci = ", loci, " ", "gen = ", gen, sep = ""),
               Loci_h_gen_sd = paste("Loci = ", loci, " ",
                                     "Herit = ", round(herit, 2), " ", 
                                     "Gen = ", gen, " ",
                                     "sd = ", sd, sep = ""))
      return(data)
    }
    
    stopCluster(cl)
    closeAllConnections()  # Close all connections
    
    combined_data <- c(combined_data, data_list)
  }
  
  combined_data <- bind_rows(combined_data)
  
  return(combined_data)
}

# Main function to list files and process them in chunks
myheatmaps <- function(path, pattern, chunk_size = 50, num_cores = 8) {
  files <- list.files(path, pattern, full.names = TRUE)
  combined_data <- process_files_in_chunks(files, chunk_size, num_cores)
  return(combined_data)
}

process_data <- function(data) {
  data %>%
    group_by(herit, loci, sd, gen, repl, Position) %>%
    mutate(
      Fixation = case_when(
        Frequency == 1 ~ "Fixed",
        Frequency < 1 & Generation == max(Generation) & max(Generation) < 2000 ~ "Lost",
        TRUE ~ NA_character_
      ),
      InitFreq = case_when(Generation == 1 ~ Frequency, 
                           TRUE ~ NA_real_)
    ) %>%
    ungroup() %>%
    dplyr::select(Generation, Position, Fixation, Frequency, Effect, InitFreq, repl, loci, Loci_h_gen_sd) %>%
    group_by(Loci_h_gen_sd, repl, Position) %>%
    mutate(InitFreq = ifelse(is.na(InitFreq), lag(InitFreq), InitFreq)) %>%
    fill(InitFreq) %>%
    drop_na() %>%
    select(-Generation, -Frequency) %>%
    distinct()
}

create_summary_df <- function(data) {
  data %>%
    group_by(loci, Loci_h_gen_sd) %>%
    summarize(
      num_fixed = sum(Fixation == "Fixed", na.rm = TRUE),
      num_lost = sum(Fixation == "Lost", na.rm = TRUE)
    ) %>%
    mutate(label = sprintf("Fixed: %d\nLost : %d", num_fixed, num_lost))
}

# Create a labeling function
custom_labels <- function(labels) {
  labels <- LETTERS[1:length(labels)]
  return(setNames(labels, names(labels)))
}

# Function to create plots for each loci
create_plot <- function(data, summary_data, loci_value) {
  plot_data <- data %>% filter(loci == loci_value)
  summary_data <- summary_data %>% filter(loci == loci_value)
  
  p <- ggplot(plot_data, aes(x = Effect, y = InitFreq, size = Effect, color = Fixation)) +
    geom_point(alpha = 0.5) +
    geom_text(aes(label = paste0(sprintf("%.2f", Effect))), 
              vjust = 1, hjust = 1, size = 2, check_overlap = TRUE, fontface = "bold", color = "black") +
    labs(title = paste("Relationship among Position, Fixation, Effect, and Initial Frequency for Loci", loci_value),
         x = "Effect",
         y = "Initial Frequency",
         color = "Fixation Status",
         size = "Effect") +
    facet_wrap(~ Loci_h_gen_sd, ncol = 2, scales = "free_y", labeller = labeller(Loci_h_gen_sd = custom_labels)) +
    geom_label(data = summary_data, aes(x = Inf, y = Inf, label = label), 
               vjust = 1.1, hjust = 1.1, size = 4, fontface = "bold", family = "mono", inherit.aes = FALSE, label.padding = unit(0.2, "lines"), label.size = 0.5) +
    theme_classic() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      axis.title = element_text(face = "bold", size = 12),
      axis.text = element_text(face = "bold", size = 10),
      legend.title = element_text(face = "bold", size = 12),
      legend.text = element_text(face = "bold", size = 10),
      legend.position = "bottom",
      legend.justification = "center",
      panel.grid = element_line(size = 1.5),
      strip.text = element_text(face = "bold", size = 12)  # Bold facet labels
    ) +
    scale_color_manual(values = c("Fixed" = "blue", "Lost" = "red"))
  
  return(p)
}

# Main function to execute all steps
run_analysis <- function(path, pattern, chunk_size = 50, num_cores = 8) {
  mydata <- myheatmaps(path, pattern, chunk_size, num_cores)
  mydata4 <- process_data(mydata)
  summary_df <- create_summary_df(mydata4)
  
  loci_values <- unique(mydata4$loci)
  plots <- lapply(loci_values, function(loci) create_plot(mydata4, summary_df, loci))
  
  return(plots)
}

