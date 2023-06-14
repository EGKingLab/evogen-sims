install.packages("WaveletComp")
library(WaveletComp)

file <- read.csv("/home/etb68/YeastProj.dir/evogen-sims/ET_Yeast/output.dir/Select_Press/ConstFluctSelect.dir/FS.dir/genome4_100_H0.5SD4Gen30.csv") %>% dplyr::select(!c(Origin, Effect))
head(file)
dat_wide <- pivot_wider(file[,c("Generation","Position","Frequency")], names_from=Position, values_from=Frequency)
dat_wide[is.na(dat_wide)] <- 0 
ff1 <-dat_wide[,2:101]
# Assume the data is in the 'ff1' data frame
wavelet_results <- analyze.wavelet(ff1, loess.span = 0)

# Global wavelet spectrum (time-averaged wavelet power spectrum)
out.spect <- wavelet_results$global_spectrum
# Corresponding frequencies
frequencies <- wavelet_results$Period

dd <- data.frame("Frequency"= frequencies, 
                 "sel"=rep("CS",each=length(frequencies)),
                 "spec"=out.spect)

simp <- dd %>% filter(Frequency <= 100) %>%
  ggplot(aes(Frequency, spec)) +
  geom_line(size=1.1) +
  geom_vline(xintercept = 60.5, color= 'grey30') +
  scale_color_manual(values=c("coral","grey30","steelblue")) +
  xlab("Periodicity (Generations)") +
  ylab("Spectral Density") +
  annotate("text",label="simulated \ndatasets", x=15, y =0.012)

simp
