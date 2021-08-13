library(tseries)
library(slider)
library(spectral)
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
