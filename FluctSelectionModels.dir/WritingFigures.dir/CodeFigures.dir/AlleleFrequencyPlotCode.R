################# Required Packages ################################

install_load <- function(packages) {
  for (package in packages) {
    if (!require(package, character.only = TRUE)) {
      install.packages(package, dependencies = TRUE)
      library(package, character.only = TRUE)
    }
  }
}

packages <- c("dplyr", "patchwork", "cowplot", "purrr", "doParallel")
install_load(packages)

#################################################
############ Setting theme #####################
#################################################

theme_set(theme_cowplot())

mythemes <- theme_bw() +
  theme(
    text            = element_text(family = "sans"), 
    legend.position = "none",
    axis.text.x     = element_text(face = "bold", size = 50, angle = 25, margin = margin(t = 20), hjust = 1),
    axis.text.y     = element_text(face = "bold", size = 50, angle = 15, margin = margin(r = 10)),
    axis.line       = element_line(size = 3),
    plot.title      = element_text(hjust = 0.01, face = "bold", size = 50, margin = margin(b = 1, unit = "lines")),
    plot.margin     = unit(c(5, 1, 1, 1), "lines"),
    axis.title.x    = element_text(size = 50, face = "bold", margin = margin(t = 30)),
    axis.title.y    = element_text(size = 50, face = "bold", margin = margin(r = 35)),
    strip.text      = element_text(size = 50, face = "bold"),
    panel.spacing   = unit(5, "lines"),
    panel.grid      = element_blank()
  )

################################################################################
####################### Data Processing and plotting ##########################
################################################################################


# The first part of the plot deals with constant and gradual II selections ### 
####### while the second is for instantaneous and gradual I selections #######

process_files <- function(dirpath, pattern, plot_type) {
 
##### Get the files from the folder that corresponds to selection type and read each file
  
  files <- list.files(dirpath, pattern, full.names = TRUE)
  
  dataframes <- list()
  for(file in files){
    herit <- as.numeric(sub(".*_H(0\\.\\d+).*", "\\1", basename(file)))
    loci <- as.numeric( sub(".*_n(\\d+)_.*", "\\1", basename(file)))
    sd <- as.numeric(sub(".*SD(\\d+).*", "\\1", basename(file)))
    gen <- as.numeric(sub(".*Gen(\\d+).*", "\\1", basename(file)))
    replicate_id <- as.numeric(sub(".*genome(\\d+)_.*", "\\1", basename(file)))
    
    data <- read.csv(file, header = TRUE) %>%
      dplyr::select(Generation, Position, Frequency, Effect) %>%
      group_by(Position) %>%
      mutate(Position = factor(Position),
             replicate = replicate_id,
             initFreq = Frequency[Generation == 1],
             position_effect_init = paste("position = ", Position, " ",
                                         "Effect = ", round(Effect, 2), " ",
                                         "Initial Freq = ", round(initFreq, 2), " ",
                                         "repl = ", replicate, sep = "")) %>% 
      ungroup() %>%
      mutate(herit = herit,
             loci = loci,
             sd = sd,
             gen = gen,
             h2_sd = paste("h\u00B2 = ", herit, " ", "sd = ", sd, sep = ""),
             h2_facet = ifelse(is.na(gen) & is.na(sd),
                               paste("h\u00B2 = ", herit, sep = ""),
                               h2_sd))
    
    dataframes[[file]] <- data
  }
  
  combined_data <- bind_rows(dataframes) %>%
    mutate(loci_gen = paste("loci = ", loci, " ", "gen = ", gen, sep = ""))
  
  plots <- list()
  if (plot_type == "loci_gen") {
    set.seed(7670986)
    loci_gens <- unique(combined_data$loci_gen)
    for(loci_geni in loci_gens){
      locus_data <- combined_data %>%
        filter(loci_gen == loci_geni)
      
      unique_positions <- unique(locus_data$Position)
      if (length(unique_positions) > 30) {
        selected_positions <- sample(unique_positions, 30)
        locus_data <- locus_data %>%
          filter(Position %in% selected_positions)
      }
      
      p <- locus_data %>%
        ggplot(aes(Generation, Frequency, group = position_effect_init,
                   color = position_effect_init)) +
        geom_line(size = 0.5) +
        facet_wrap(~h2_facet, ncol = 3, dir = "v") +
        ylim(0, 1) +
        mythemes
      
      plots[[loci_geni]] <- p
    }
  } else if (plot_type == "loci") {
    loci <- unique(combined_data$loci)
    for(locus in loci){
      locus_data <- combined_data %>%
        filter(loci == locus)
      
      unique_positions <- unique(locus_data$Position)
      if (length(unique_positions) > 30) {
        selected_positions <- sample(unique_positions, 30)
        locus_data <- locus_data %>%
          filter(Position %in% selected_positions)
      }
      
      p <- locus_data %>%
        ggplot(aes(Generation, Frequency, group = position_effect_init,
                   color = position_effect_init)) +
        geom_line(linewidth = 0.5) +
        facet_wrap(~h2_facet, ncol = 3, dir = "v") +
        ylim(0, 1) +
        mythemes
      
      plots[[locus]] <- p
    }
  }
  
  # Return data and allele frequency plots. We can use combined data for additional exploration
  return(list(combined_data = combined_data, plots = plots))
}
