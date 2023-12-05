########### Allele Frequency ##############

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

#source("~/YeastProj.dir/evogen-sims/ET_Yeast/Code/Selection_Models/WF.dir/Plotting/PrelimVisualSelect.R")
dir_path <- "~/YeastProj.dir/evogen-sims/ET_Yeast/output.dir/Selection_Models/WF.dir/NS.dir/"
pattern <- "^genome.*0\\.1SD3Gen20\\.csv$"
ncol <- 3

plot_files(dir_path, pattern, ncol = ncol)

geno1 <- read.csv("~/YeastProj.dir/evogen-sims/ET_Yeast/output.dir/Selection_Models/WF.dir/SinFS.dir/genome1_n10_H0.1SD3Gen10.csv") %>% 
  mutate(hetero = 2*(Frequency - Frequency^2),
         Fixation = case_when(Frequency == 0 ~ 1,
                              Frequency == 1 ~ 1,
                              TRUE ~ 0)) 

geno2 <- read.csv("~/YeastProj.dir/evogen-sims/ET_Yeast/output.dir/Selection_Models/WF.dir/SinFS.dir/genome2_n10_H0.1SD3Gen10.csv") %>% 
  mutate(hetero = 2*(Frequency - Frequency^2),
         Fixation = case_when(Frequency == 0 ~ 1,
                              Frequency == 1 ~ 1,
                              TRUE ~ 0)) 

geno <- full_join(geno1, geno2) %>% dplyr::select(!Origin) %>% 
  group_by(Generation, Position) %>% 
  summarize(Mean_Freq = mean(Frequency), Mean_Effect = mean(Effect), 
            Mean_Het = mean(hetero), Mean_Fix = mean(Fixation))

geno %>% ggplot(aes(Generation, Mean_Freq, color = as.character(Position)))+
  geom_line(linewidth = 0.2)+
  theme(axis.title = element_text(face = "bold"),
        axis.text = element_text(face = "italic"),
        plot.title = element_text(hjust = 0.5),
        legend.position = "none")

#########################

# Import and merge all CSV files into one data frame
geno <- list.files(path = "~/YeastProj.dir/evogen-sims/ET_Yeast/output.dir/Selection_Models/WF.dir/SinFS.dir/", # Specify your path
                   pattern = "^genome.*_n100_H0\\.5SD3Gen20\\.csv$",
                   full.names = TRUE) %>% 
  lapply(read_csv) %>% 
  # Store all files in a list
  bind_rows # Combine data sets into one data frame

# Select and summarize the data
geno <- geno %>% 
  mutate(hetero = 2*(Frequency - Frequency^2),
         Fixation = case_when(Frequency == 0 ~ 1,
                              Frequency == 1 ~ 1,
                              TRUE ~ 0)) %>% 
  dplyr::select(!Origin) %>% 
  group_by(Generation, Position) %>% 
  summarize(Mean_Freq = mean(Frequency), 
            Mean_Effect = mean(Effect), 
            Mean_Het = mean(hetero), 
            Mean_Fix = mean(Fixation))

# Plot the data

geno %>% 
  ggplot(aes(Generation, Mean_Freq, color = as.character(Position))) + 
  geom_line(linewidth = 0.5)+ # Add line plot
  theme(axis.title = element_text(face = "bold"), 
        axis.text = element_text(face = "italic"), 
        plot.title = element_text(hjust = 0.5), 
        legend.position = "none")

################################# ZOOM IN ############################

x = seq(0, 2000, by = 40)
y = seq(20, 2020, 40) # x + 25
myrect <- data.frame(x,y)

 ggplot()+
  geom_line(data = geno, aes(x = Generation, y = Mean_Freq, 
                             color = as.character(Position)))+ 
   theme(axis.title = element_text(face = "bold"),
         axis.text = element_text(face = "italic"), 
         plot.title = element_text(hjust = 0.5), 
         legend.position = "none") + 
  facet_zoom(xlim = c(800,1000),ylim = c(0.2, 0.8), 
             horizontal = T, zoom.size = 1)+
  geom_rect(data = myrect,
            aes(xmin = x, xmax = y,
                ymin = 0, ymax = 1), 
            fill = "grey", alpha = 0.2)


#####################$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$###########################
test_test <- read.csv("~/YeastProj.dir/evogen-sims/ET_Yeast/output.dir/Selection_Models/WF.dir/SinFS.dir/genome1_n100_H0.5SD3Gen30.csv")
 
 ggplot()+
   geom_line(data = test_test, aes(x = Generation, y = Frequency, 
                              color = as.character(Position)))+ 
   theme(axis.title = element_text(face = "bold"),
         axis.text = element_text(face = "italic"), 
         plot.title = element_text(hjust = 0.5), 
         legend.position = "none") + 
   facet_zoom(xlim = c(200,520),ylim = c(0.0, 1.0), 
              horizontal = T, zoom.size = 1)+
   geom_rect(data = myrect,
             aes(xmin = x, xmax = y,
                 ymin = 0, ymax = 1), 
             fill = "grey", alpha = 0.2)
 
 ###############################################################################
 ############## Creating functions for multiple heritability values ############
 ###############################################################################
 
 mytheme <- function(){
   theme_set(theme_cowplot())+
     theme(axis.title = element_text(face = "bold"),
           axis.text = element_text(face = "italic"),
           plot.title = element_text(hjust = 0.5),
           legend.position = "none")
 }
 
 allele_plots <- function(df, mygenomefile){
  herit_val <- str__match(string = mygenomefile, pattern = "(?<=genome.*_n100_H0\\.)[0-9]+", str_match(string, pattern)[, 3])
  df$herity <- herit_val
  Freq_plot <- ggplot(df, aes(x = Generation, y = Mean_Freq, color = as.character(Position)))+
    geom_line(linewidth = 0.2)+
    mytheme()+
    facet_grid(~herity)
  return(list(herity = herity, plot = Freq_plot))
 }
 
 plot_files <- function(dir_path, pattern, ncol){
   mygenomefiles <- list.files(dir_path, pattern = pattern, full.names = TRUE)
   plots_list <- list()
   for (mygenomefile in mygenomefiles){
     df <- read.csv(mygenomefile)
     df <- df %>% 
       mutate(hetero = 2*(Frequency - Frequency^2),
              Fixation = case_when(Frequency == 0 ~ 1,
                                   Frequency == 1 ~ 1,
                                   TRUE ~ 0)) %>% 
       dplyr::select(!Origin) %>% 
       group_by(Generation, Position, herit_val) %>% 
       summarize(Mean_Freq = mean(Frequency), 
                 Mean_Effect = mean(Effect), 
                 Mean_Het = mean(hetero), 
                 Mean_Fix = mean(Fixation))
     result <- allele_plots(df, mygenomefile)
     herity <- result$herity
     plot <- result$plot
     if (!is.null(plots_list[[herity]])){
       plots_list[[herity]][[mygenomefile]] <- plot
     }else{
       plots_list[[herity]] <- list(mygenomefile = plot)
     }
   }
   for (herity in names(plots_list)){
     print(plot_grid(plotlist = plots_list[[herity]], ncol = ncol))
   }
 }
 ###########################################################################
 
 
 mytheme <- function(){
   theme_set(theme_cowplot())
   theme(axis.title = element_text(face = "bold"),
         axis.text = element_text(face = "italic"),
         plot.title = element_text(hjust = 0.5),
         legend.position = "none")
 }
 
 allele_plots <- function(df, mygenomefile){
   herit_val <- as.numeric(str_match(string = mygenomefile, pattern = "(?<=genome.*_n100_H0\\.)[0-9]+")[,3])
   df$herit_val <- herit_val
   Freq_plot <- ggplot(df, aes(x = Generation, y = Mean_Freq, color = as.factor(Position))) +
     geom_line(linewidth = 0.2) +
     mytheme() +
     facet_grid(~herit_val)
   return(list(herit_val = herit_val, plot = Freq_plot))
 }
 
 plot_files <- function(dir_path, pattern, ncol){
   mygenomefiles <- list.files(dir_path, pattern = pattern, full.names = TRUE)
   plots_list <- list()
   for (mygenomefile in mygenomefiles){
     df <- read.csv(mygenomefile)
     df <- df %>% 
       mutate(hetero = 2 * (Frequency - Frequency^2),
              Fixation = case_when(Frequency == 0 ~ 1,
                                   Frequency == 1 ~ 1,
                                   TRUE ~ 0)) %>% 
       dplyr::select(!Origin) %>% 
       group_by(Generation, Position, herit_val) %>% 
       summarize(Mean_Freq = mean(Frequency), 
                 Mean_Effect = mean(Effect), 
                 Mean_Het = mean(hetero), 
                 Mean_Fix = mean(Fixation))
     result <- allele_plots(df, mygenomefile)
     herit_val <- result$herit_val
     plot <- result$plot
     if (!is.null(plots_list[[as.character(herit_val)]])){
       plots_list[[as.character(herit_val)]][[mygenomefile]] <- plot
     }else{
       plots_list[[as.character(herit_val)]] <- list(mygenomefile = plot)
     }
   }
   for (herity in names(plots_list)){
     print(plot_grid(plotlist = plots_list[[herit_val]], ncol = ncol))
   }
 }
 
 dir_path <- "~/YeastProj.dir/evogen-sims/ET_Yeast/output.dir/Selection_Models/WF.dir/NS.dir/"
 pattern <- "^genome.*_n100_H0\\.\\dSD3Gen20\\.csv$"
 ncol <- 3
 
 plot_files(dir_path, pattern, ncol = ncol)
 
 