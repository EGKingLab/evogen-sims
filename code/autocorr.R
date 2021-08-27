library(tseries)
library(slider)
library(spectral)
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())


set.sims <- data.frame("pop"=c("C","C","F","F"), "sel"=c("C","F","C","F"), 
                       "samp_id"=c("S6","S4","S301","S300"))

out.spect <- matrix(NA,250,4)

for(jj in 1:nrow(set.sims))
{
  pop_type <- set.sims$pop[jj]
  sel_type <- set.sims$sel[jj]
  samp_id <- set.sims$samp_id[jj]
  output <- readRDS(file = paste0("../output/Sim_Output_S",sel_type,"_P",pop_type,"_",samp_id,".rds"))
  ff1 <- output$freqs[seq(1,500, by =1),]
  ttest <- spectrum(ff1, spans=2, plot=FALSE)
  out.spect[,jj] <- rowMeans(ttest$spec)
  
}

output.N <- readRDS(file="../output/Null_Sim_Output.rds")
ff1 <- output.N[seq(11,510,by=1),]
ttest <- spectrum(ff1, spans=2, plot=FALSE)
tt <- rowMeans(ttest$spec)

dd <- data.frame("Frequency"= rep((1/ttest$freq), 5), 
                 "pop"=rep(c("C","C","F","F","C"), each=250),
                 "sel"=rep(c("C","F","C","F","N"),each=250),
                 "spec"=c(out.spect[,1],out.spect[,2],out.spect[,3],out.spect[,4], tt))


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

#plot(1/ttest$freq,tt,type='l')

#acor <- acf(ff1,lag.max=499,pl=FALSE)
#aa <- matrix(NA,1000,500)

#for(ii in 1:1000)
#{
#aa[ii,] <- acor$acf[,ii,ii]
#}

#mm <- colMeans(aa, na.rm=TRUE)
#plot(mm)
#abline(v=c(25,50,75,100))

set.sims <- data.frame("pop"=c("C","C","F","F"), "sel"=c("C","F","C","F"), 
                       "samp_id"=c("S6","S4","S301","S300"))

out.diff1 <- matrix(NA,499,5)
out.abs.diff1 <- matrix(NA,499,5)
out.diff25 <- matrix(NA,475,5)
out.abs.diff25 <- matrix(NA,475,5)

out.abs.diff_sub <- matrix(NA,19,5)

for(jj in 1:4)
{
pop_type <- set.sims$pop[jj]
sel_type <- set.sims$sel[jj]
samp_id <- set.sims$samp_id[jj]
output <- readRDS(file = paste0("../output/Sim_Output_S",sel_type,"_P",pop_type,"_",samp_id,".rds"))
ff1 <- output$freqs[seq(1,500, by =1),]

  dd <- diff(ff1, lag=1)
  out.abs.diff1[,jj] <- rowMeans(abs(dd))
  out.diff1[,jj] <- rowMeans(dd)
  dd <- diff(ff1, lag=25)
  out.abs.diff25[,jj] <- rowMeans(abs(dd))
  out.diff25[,jj] <- rowMeans(dd)
  
  ff2 <- output$freqs[seq(1,500, by = 25),]
  dd <- diff(ff2, lag=1)
  out.abs.diff_sub[,jj] <- rowMeans(abs(dd))
}


output.N <- readRDS(file="../output/Null_Sim_Output.rds")
ff1 <- output.N[seq(11,510,by=1),]
dd <- diff(ff1, lag=1)
out.abs.diff1[,5] <- rowMeans(abs(dd))
out.diff1[,5] <- rowMeans(dd)
dd <- diff(ff1, lag=25)
out.abs.diff25[,5] <- rowMeans(abs(dd))
out.diff25[,5] <- rowMeans(dd)
ff2 <- output$freqs[seq(11,510, by = 25),]
dd <- diff(ff2, lag=1)
out.abs.diff_sub[,5] <- rowMeans(abs(dd))


plot(out.abs.diff1[,3])
points(out.abs.diff1[,2], col='blue')
points(out.abs.diff1[,1],col='green')
points(out.abs.diff1[,4],col='red')
points(out.abs.diff1[,5],col='yellow')


plot(out.diff1[,3])
points(out.diff1[,2], col='blue')
points(out.diff1[,1],col='green')
points(out.diff1[,4],col='red')

plot(out.abs.diff25[,3])
points(out.abs.diff25[,2], col='blue')
points(out.abs.diff25[,1],col='green')
points(out.abs.diff25[,4],col='red')

plot(out.diff25[,3])
points(out.diff25[,2], col='blue')
points(out.diff25[,1],col='green')
points(out.diff25[,4],col='red')
points(out.diff25[,5],col='yellow')

plot(out.abs.diff_sub[,4],ylim=c(0,0.05),type='l')
lines(out.abs.diff_sub[,1], col='blue')
lines(out.abs.diff_sub[,2],col='green')
lines(out.abs.diff_sub[,3],col='red')
lines(out.abs.diff_sub[,5],col='yellow')

dd1 <- data.frame("Generation"= rep(seq(1,499),5), 
                 "pop"=rep(c("C","C","F","F","C"), each=499),
                 "sel"=rep(c("C","F","C","F","N"),each=499),
                 "AbsoluteDifference"=c(out.abs.diff1[,1],
                                         out.abs.diff1[,2],
                                         out.abs.diff1[,3],
                                         out.abs.diff1[,4],
                                         out.abs.diff1[,5]))

dd25 <- data.frame("Generation"= rep(seq(1,475),5), 
                  "pop"=rep(c("C","C","F","F","C"), each=475),
                  "sel"=rep(c("C","F","C","F","N"),each=475),
                  "AverageDifference"=c(out.diff25[,1],
                                          out.diff25[,2],
                                          out.diff25[,3],
                                          out.diff25[,4],
                                          out.diff25[,5]))

ddsub <- data.frame("Generation"= rep(seq(1,19),5), 
                  "pop"=rep(c("C","C","F","F","C"), each=19),
                  "sel"=rep(c("C","F","C","F","N"),each=19),
                  "AbsoluteDifference"=c(out.abs.diff_sub[,1],
                                         out.abs.diff_sub[,2],
                                         out.abs.diff_sub[,3],
                                         out.abs.diff_sub[,4],
                                         out.abs.diff_sub[,5]))

pplot1 <- ggplot(dd1,aes(Generation, AbsoluteDifference, color=sel, pch=pop)) +
  geom_point(alpha=1/2) +
  scale_color_manual(values=c("coral","grey30","steelblue")) +
  ylab("Mean Absolute Difference") +
  annotate("text",label="simulated \ndatasets", x=455, y =0.0073)
pplot1

pplot1 <- ggplot(dd1,aes(Generation, AbsoluteDifference, color=pop, linetype=sel)) +
  geom_line() +
  scale_color_manual(values=c("coral","grey30","steelblue")) +
  ylab("Mean Absolute Difference")
pplot1

pplotsub <- ggplot(ddsub,aes(Generation*25, AbsoluteDifference, color=sel, linetype=pop)) +
  geom_line() +
  scale_color_manual(values=c("coral","grey30","steelblue")) +
  ylab("Mean Absolute Difference") 
pplotsub

pplot25 <- ggplot(dd25,aes(Generation, AverageDifference, color=sel, pch=pop)) +
  geom_point(alpha=1/2) +
  scale_color_manual(values=c("coral","grey30","steelblue")) +
  ylab("Mean Difference")
pplot25

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
  annotate("text",label="empirical yeast \ndataset", x=1.2, y =0.25)
empp

empm <- emp.means %>% 
  filter(samp %in% c("freq_YEE_0112_04_01_","freq_YEE_0112_04_02_",
                     "freq_YEE_0112_04_03_","freq_YEE_0112_04_10_",
                     "freq_YEE_0112_04_11_","freq_YEE_0112_04_12_",
                     "maf_3_","maf_7_",
                     "maf_8_","maf_9_",
                     "maf_10_","maf_11_",
                     "maf_12_")) %>%
  ggplot(aes(Generations, ABSD, color=pop)) +
  geom_point(alpha=1/2,position = position_jitter(width = 0.25)) + 
  xlab("Generation") +
  ylab("Mean Absolute Difference") +
  scale_x_continuous(breaks=c(1,2), labels=c("timept1","timept2"))
#theme(legend.position = "none") +
#ylim(0,0.15)
empm

pall <- plot_grid(pplot1, empp, ncol=2, labels=c("a.","b."))
ggsave(pall, file="../output/pop_abs_diffs.pdf", width=8, height=3.5)


empp <- emp.diff %>% filter(samp %in% c("freq_YEE_0112_04_01_","freq_YEE_0112_04_02_","freq_YEE_0112_04_03_","freq_YEE_0112_04_10_")) %>%
  ggplot(aes(Generations, AbsoluteDifference, group=sample_chromosome, color=samp)) +
  geom_line(alpha=1/2) + 
  xlab("Generation") +
  ylab("Mean Absolute Difference") +
  theme(legend.position = "none") +
#ylim(0,0.15)
empp

set.sims <- data.frame("pop"=c("C","C","F","F"), "sel"=c("C","F","C","F"), 
                       "samp_id"=c("S6","S4","S301","S300"))


jj<-2

  pop_type <- set.sims$pop[jj]
  sel_type <- set.sims$sel[jj]
  samp_id <- set.sims$samp_id[jj]
  output <- readRDS(file = paste0("../output/Sim_Output_S",sel_type,"_P",pop_type,"_",samp_id,".rds"))
  ff1 <- output$freqs[seq(1,500, by =1),]
  
  
dd.mat.abs <- vector(mode='list', length=(nrow(ff1)-1))
dd.mat <- vector(mode='list', length=(nrow(ff1)-1))
for(jj in 1:(nrow(ff1)-1))
#dd.mat.abs <- vector(mode='list', length=100)
#dd.mat <- vector(mode='list', length=100)
#for(jj in 1:100)
  
  {

  dd <- diff(ff1, lag=jj)
  dd.mat.abs[[jj]] <- rowMeans(abs(dd))
  dd.mat[[jj]] <- rowMeans(dd)
  
  #dd_sd <- apply(dd, 2, sd)
  #cc <- cor(dd[,-which(dd_sd==0)])
  #cat(jj,"\t",mean(abs(cc[upper.tri(cc)])),"\n")
}

plot(dd.mat.abs[[1]])
plot(dd.mat.abs[[10]])
plot(dd.mat.abs[[25]])

plot(dd.mat[[1]])
plot(dd.mat[[10]])
plot(dd.mat[[25]])


#need likelihood of resulting stat

#slide.out <- matrix(NA,490,ncol(ff1))

#for(gg in 1:ncol(ff1))
#{
#  slide.out[,gg] <- unlist(slide(ff1[,gg],var,.before=10,.complete=TRUE))
#}

#ss <- rowMeans(slide.out)



dd_sd <- apply(dd, 2, sd)
cc <- cor(dd[,-which(dd_sd==0)])
mean(abs(cc[upper.tri(cc)]))

ff1_sd <- apply(ff1, 2, sd)
cc <- cor(ff1[,-which(ff1_sd==0)])
mean(abs(cc[upper.tri(cc)]))
#WGCNA, PCA
