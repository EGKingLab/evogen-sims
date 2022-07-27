library(tseries)
library(slider)
library(spectral)
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())


set.sims <- data.frame("pop"=c("C","C","F","F","C"), "sel"=c("C","F","C","F","N"), 
                       "samp_id"=rep("2",5))

out.spect <- matrix(NA,250,5)
out.abs.diff1 <- matrix(NA,499,5)

for(jj in 1:nrow(set.sims))
{
  pop_type <- set.sims$pop[jj]
  sel_type <- set.sims$sel[jj]
  samp_id <- set.sims$samp_id[jj]
  output <- read_csv(file = paste0("../output/SLiM/genome_track_",sel_type,"S",pop_type,"P_",samp_id,".csv"))
  dat_wide <- pivot_wider(output[,c("Generation","Position","Frequency")], names_from=Position, values_from=Frequency)
  dat_wide[is.na(dat_wide)] <- 0 
  ff1 <-dat_wide[,2:101]
  ttest <- spectrum(ff1, spans=2, plot=FALSE)
  out.spect[,jj] <- rowMeans(ttest$spec)
  ddifs <- diff(as.matrix(ff1), lag=1)
  out.abs.diff1[,jj] <- rowMeans(abs(ddifs))
  
}



dd <- data.frame("Frequency"= rep((1/ttest$freq), 5), 
                 "pop"=rep(c("C","C","F","F","C"), each=250),
                 "sel"=rep(c("C","F","C","F","N"),each=250),
                 "spec"=c(out.spect[,1],out.spect[,2],out.spect[,3],out.spect[,4], out.spect[,5]))


simp <- dd %>% filter(Frequency <= 75) %>%
  ggplot(aes(Frequency, spec, color=sel, linetype=pop)) +
  geom_line(size=1.1) +
  geom_vline(xintercept=50, color= 'grey30') +
  scale_color_manual(values=c("coral","grey30","steelblue")) +
  xlab("Periodicity (Generations)") +
  ylab("Spectral Density") +
  annotate("text",label="simulated \ndatasets", x=15, y =0.012)

emp.set.all <- readRDS(file="../output/emp_data_spec.Rds")

empp <- emp.set.all %>% filter(freq <= 200, samp %in% c("maf_3_","maf_7_","maf_10_","maf_11_")) %>%
  ggplot(aes(freq, spec, group=sample_chromosome, color=samp)) +
  geom_line(alpha=1/2) + 
  ylim(0,0.0142) +
  xlab("Periodicity (Generations)") +
  ylab("Spectral Density") +
  theme(legend.position = "none") +
  annotate("text",label="empirical yeast \ndataset", x=90, y =0.012)


pall <- plot_grid(simp, empp, ncol=2, labels=c("a.","b."))
ggsave(pall, file="../output/sim_spectral.pdf", width=8, height=4)




dd1 <- data.frame("Generation"= rep(seq(1,499),5), 
                 "pop"=rep(c("C","C","F","F","C"), each=499),
                 "sel"=rep(c("C","F","C","F","N"),each=499),
                 "AbsoluteDifference"=c(out.abs.diff1[,1],
                                         out.abs.diff1[,2],
                                         out.abs.diff1[,3],
                                         out.abs.diff1[,4],
                                         out.abs.diff1[,5]))

pplot1 <- ggplot(dd1,aes(Generation, AbsoluteDifference, color=sel, pch=pop)) +
  geom_point(alpha=1/2) +
  scale_color_manual(values=c("coral","grey30","steelblue")) +
  ylab("Mean Absolute Difference") +
  annotate("text",label="simulated \ndatasets", x=70, y =0.011)
pplot1

emp.diff <- readRDS(file="../output/emp_data_diff.Rds")
emp.diff$pop <- "Large"
emp.small <- readRDS(file="../output/emp_data_small_diffs.Rds")
emp.small$pop <- "Small"

emp.diff <- rbind(emp.diff, emp.small)
emp.means <- emp.diff %>% 
  group_by(samp, Generations, pop) %>%
  summarise("ABSD" = mean(AbsoluteDifference))

empp <- emp.diff %>% filter(samp %in% c("maf_10_","maf_11_","freq_YEE_0112_04_02_","freq_YEE_0112_04_03_")) %>%
  ggplot(aes(Generations, AbsoluteDifference, group=sample_chromosome, shape=pop, color=samp)) +
  geom_point(alpha=1/2,position = position_jitter(width = 0.25), size=1.8) + 
  xlab("Generation Span") +
  ylab("Mean Absolute Difference") +
  scale_x_continuous(breaks=c(1,2), labels=c("G180_G360","G360_G540")) +
  guides(color="none") +
  annotate("text",label="empirical yeast \ndataset", x=1.2, y =0.23)
empp

pall <- plot_grid(pplot1, empp, ncol=2, labels=c("a.","b."))
ggsave(pall, file="../output/pop_abs_diffs.pdf", width=8, height=3.5)

temps <- spectrum(out.abs.diff1[,3], spans=2, plot=FALSE)
plot(1/temps$freq, temps$spec, type='l')