library(tidyverse)
library(cowplot)
library(patchwork)

genome <- read.csv("~/YeastProj.dir/evogen-sims/ET_Yeast/output.dir/Selection_Models/WF.dir/LinFS.dir/genome1_n1_H0.1SD3Gen30.csv")

# Select relevant columns and replace NAs with 0
all_freq <- genome %>% 
  select(Generation, Position, Frequency) %>% 
  pivot_wider(names_from = Position, values_from = Frequency) %>% 
  mutate(across(everything(), ~replace_na(., 0)))

all_freq_mean <- rowMeans(all_freq[-1])

my_fft <- fft(all_freq_mean)
n <- length(my_fft)
amplitude <- abs(my_fft)/n
freq <- seq(1, n) * (1 / n)
new_data <- data.frame(amplitude, freq)

p1 <- new_data %>% filter(freq > 0.01 & freq < 0.80) %>% 
  mutate(freq = 1/freq) %>% 
  ggplot(aes(x = freq, y = amplitude)) +
  geom_line() +
  labs(
       x = "Frequency",
       y = "Amplitude") +
  theme_minimal()+
  theme(legend.position = "none")

p1
###########################################################################
# 
# # Create a list to store the FFT results
# myfft <- list()
# 
# # Get the names of the frequency columns
# positions <- colnames(all_freq)[-1]
# 
# # Loop through each position column and calculate the FFT
# for(position in positions){
#   frq <- all_freq[[position]] # Use double brackets to access column by name
#   myfft[[position]] <- fft(frq)
# }
# 
# # Initialize an empty data frame for plotting
# plot_data <- data.frame(Position = character(), Frequency = numeric(), Amplitude = numeric(), stringsAsFactors = FALSE)
# 
# # Loop through each position to prepare data
# for(position in positions){
#   fft_result <- myfft[[position]]
#   n <- length(fft_result)
#   amplitude <- abs(fft_result)/n
#   freq <- seq(0, n-1) * (1 / n)
#   
#   temp_data <- data.frame(Position = rep(position, n), Frequency = freq, Amplitude = amplitude)
#   plot_data <- rbind(plot_data, temp_data)
# }
# 
# plot_data %>% filter(Frequency > 0.02 & Frequency < 0.55) %>% 
#   mutate(Frequency = 1/Frequency) %>% 
#   ggplot(aes(x = Frequency, y = Amplitude, color = Position)) +
#   geom_line() +
#   labs(title = "Combined Plot of Frequency and Amplitude for All Positions",
#        x = "Frequency",
#        y = "Amplitude") +
#   theme_minimal()+
#   theme(legend.position = "none")

############################ SPECTRAL ANALYSIS #############

#head(all_freq_mean)
spects <- spectrum(all_freq_mean, spans = 2, plot = F)
spec <- spects$spec
fre <- spects$freq

newfft <- data.frame(Amplitude = spects$spec, Frequency = fre)

p2 <- newfft %>% filter(Frequency > 0.01) %>% 
  mutate(Frequency = 1/Frequency) %>% 
  ggplot(aes(x = Frequency, y = Amplitude)) +
  geom_line() +
  labs(x = "Frequency",
       y = "Amplitude") +
  theme_minimal()+
  theme(legend.position = "none")
p2

p1/p2
