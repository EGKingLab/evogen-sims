library(tidyverse)
library(cowplot)
theme_set(theme_cowplot())

#make pheno plot
ngens <- 500
nreps <- 30
Long_phenos <- data.frame("Type" = numeric(length = ngens*nreps),
                          "Generation" = rep(seq(1,ngens),nreps),
                          "Phenotype"= numeric(length = ngens*nreps),
                          "Optimum"= numeric(length = ngens*nreps),
                          "Sim_run" = character(length = ngens*nreps))
pop_type <- "C"
counter <- 1
row_counter <- 1
gen_os <- c(1,5,10,25,50)
niter <- 5
for(nn in 1:niter)
{
  
  sel_type <- "F"
  
  for(gg in gen_os)
  {
    samp_id <- paste0("S",counter)
    output <- readRDS(file = paste0("../output/Sim_Output_S",sel_type,"_P",pop_type,"_",samp_id,".rds"))
    Long_phenos[row_counter:(row_counter+499),"Type"] <- rep(paste0(sel_type,gg),500)
    Long_phenos[row_counter:(row_counter+499),"Sim_run"] <- rep(samp_id,500)
    Long_phenos[row_counter:(row_counter+499),"Phenotype"] <- output$pheno[1:500,1]
    Long_phenos[row_counter:(row_counter+499),"Optimum"] <- output$optim
    counter <- counter + 1
    row_counter <- row_counter + 500
    
  }
  
  sel_type <- "C"
  samp_id <- paste0("S",counter)
  output <- readRDS(file = paste0("../output/Sim_Output_S",sel_type,"_P",pop_type,"_",samp_id,".rds"))
  Long_phenos[row_counter:(row_counter+499),"Type"] <- rep(sel_type,500)
  Long_phenos[row_counter:(row_counter+499),"Sim_run"] <- rep(samp_id,500)
  Long_phenos[row_counter:(row_counter+499),"Phenotype"] <- output$pheno[1:500,1]
  Long_phenos[row_counter:(row_counter+499),"Optimum"] <- output$optim
  counter <- counter + 1
  row_counter <- row_counter + 500
  
}

Long_phenos <- Long_phenos %>%
  mutate(Type = factor(Type),
         Type = fct_relevel(Type, "C","F1","F5","F10","F25","F50"))

opts <- c(Long_phenos$Optimum[1], Long_phenos$Optimum[2])

g1 <- Long_phenos[1:3000,] %>%
  filter(Type %in% c("C","F5","F25")) %>%
  ggplot(aes(Generation,Phenotype, color=Type, group=Sim_run)) +
  geom_line() + 
  geom_hline(yintercept = Long_phenos$Optimum[1], lty = 3) +
  geom_hline(yintercept = Long_phenos$Optimum[2], lty = 3) +
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank())


g1


ngens <- 500
nreps <- 4
Long_phenos <- data.frame("Type" = numeric(length = ngens*nreps),
                          "Generation" = rep(seq(1,ngens),nreps),
                          "Phenotype"= numeric(length = ngens*nreps),
                          "Optimum"= numeric(length = ngens*nreps),
                          "Sim_run" = character(length = ngens*nreps))

samp_id <- "S300"
output <- readRDS(file = paste0("../output/Sim_Output_SF_PF_S300.rds"))
Long_phenos[1:500,"Type"] <- rep("FS_FP",500)
Long_phenos[1:500,"Sim_run"] <- rep(samp_id,500)
Long_phenos[1:500,"Phenotype"] <- output$pheno[1:500,1]
Long_phenos[1:500,"Optimum"] <- output$optim

output <- readRDS(file = paste0("../output/Sim_Output_SC_PF_S301.rds"))
Long_phenos[501:1000,"Type"] <- rep("CS_FP",500)
Long_phenos[501:1000,"Sim_run"] <- rep(samp_id,500)
Long_phenos[501:1000,"Phenotype"] <- output$pheno[1:500,1]
Long_phenos[501:1000,"Optimum"] <- output$optim

output <- readRDS(file = paste0("../output/Sim_Output_SC_PC_S6.rds"))
Long_phenos[1001:1500,"Type"] <- rep("CS_CP",500)
Long_phenos[1001:1500,"Sim_run"] <- rep(samp_id,500)
Long_phenos[1001:1500,"Phenotype"] <- output$pheno[1:500,1]
Long_phenos[1001:1500,"Optimum"] <- output$optim

output <- readRDS(file = paste0("../output/Sim_Output_SF_PC_S4.rds"))
Long_phenos[1501:2000,"Type"] <- rep("FS_CP",500)
Long_phenos[1501:2000,"Sim_run"] <- rep(samp_id,500)
Long_phenos[1501:2000,"Phenotype"] <- output$pheno[1:500,1]
Long_phenos[1501:2000,"Optimum"] <- output$optim

g1p <- Long_phenos %>%
  ggplot(aes(Generation,Phenotype, color=Type)) +
  geom_line() + 
  geom_hline(yintercept = opts[1], lty = 3) +
  geom_hline(yintercept = opts[2], lty = 3) +
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank())


g1p

phen.all <- plot_grid(g1,g1p, labels=c("a.","b."),ncol=2)
ggsave(phen.all, filename="../output/phenotype_sim.pdf", width=6.5, height=2.5)

#### Alllele freq plots

#fluctuating sel constant pop

ngens <- 500
pop_type <- "C"
gen_os <- 25
sel_type <- "F"

samp_id <- "S4"


output <- readRDS(file = paste0("../output/Sim_Output_S",sel_type,"_P",pop_type,"_",samp_id,".rds"))

ff <- as.data.frame(output$freqs[1:500,])
ff$Generation <- seq(1,500)


Long_genes <- ff %>%
pivot_longer(-Generation, names_to="GeneId", values_to = "AlleleFrequency")

Long_genes$EffectSize <- rep(output$effects,times=500)
Long_genes$Init_freq <- rep(output$freqs[501,],times=500)

  


ex.ef <- Long_genes %>%
  filter(EffectSize > 5.5, 
         #Generation %in% seq(1,500, by = 25),
         Init_freq > 0.3 & Init_freq < 0.7)
ex.ef$Effect <- "Large"

sm.ef <- Long_genes %>%
  filter(EffectSize <0.5, Init_freq > 0.3 & Init_freq < 0.7)
sm.ef$Effect <- "Small"


ex.ef <- rbind(subset(ex.ef, GeneId %in% ex.ef$GeneId[c(1,2)]), subset(sm.ef, GeneId %in% sm.ef$GeneId[c(1,2)]))

FSCP <- Long_genes %>%
  filter(Generation %in% seq(1,500, by = 1),
         EffectSize > 0,
         Init_freq > 0.1 & Init_freq < 0.9
         ) %>%
  ggplot(aes(Generation, AlleleFrequency, group=GeneId, color=EffectSize)) +
  geom_line(alpha=1/2) +
  scale_colour_gradientn(colors=c("white","black")) +
  theme(panel.background = element_rect(fill = "white")) + 
  theme(legend.position = "none") +
  ggtitle("FSCP") +
  theme(plot.title = element_text(hjust = 0.5, size=12,  color="steelblue"))


#  scale_colour_gradientn(colors=c("red","white","blue","black"))
  #geom_line(data=ex.ef, aes(Generation, AlleleFrequency, color=Effect)) +
  #scale_color_manual(values=c("blue","red"))

##### Constant Sel Constant pop

ngens <- 500
pop_type <- "C"
gen_os <- NA
sel_type <- "C"

samp_id <- "S6"


output <- readRDS(file = paste0("../output/Sim_Output_S",sel_type,"_P",pop_type,"_",samp_id,".rds"))

ff <- as.data.frame(output$freqs[1:500,])
ff$Generation <- seq(1,500)


Long_genes <- ff %>%
  pivot_longer(-Generation, names_to="GeneId", values_to = "AlleleFrequency")

Long_genes$EffectSize <- rep(output$effects,times=500)
Long_genes$Init_freq <- rep(output$freqs[501,],times=500)

mm <- min(which(output$pheno>=(output$optim[1]-5)))

ex.ef <- Long_genes %>%
  filter(EffectSize > 5.5, 
         #Generation %in% seq(1,500, by = 25),
         Init_freq > 0.3 & Init_freq < 0.7)
ex.ef$Effect <- "Large"

sm.ef <- Long_genes %>%
  filter(EffectSize <0.5, Init_freq > 0.3 & Init_freq < 0.7)
sm.ef$Effect <- "Small"


ex.ef <- rbind(subset(ex.ef, GeneId %in% ex.ef$GeneId[c(1,2)]), subset(sm.ef, GeneId %in% sm.ef$GeneId[c(1,2)]))

CSCP <- Long_genes %>%
  filter(EffectSize > 0,
         Generation %in% seq(1,500, by = 1),
         Init_freq > 0.1 & Init_freq < 0.9
  ) %>%
  ggplot(aes(Generation, AlleleFrequency, group=GeneId, color=EffectSize)) +
  geom_line(alpha=1/3) +
  scale_colour_gradientn(colors=c("white","black"))+
  #geom_line(data=ex.ef, aes(Generation, AlleleFrequency, color=Effect)) +
  #scale_color_manual(values=c("blue","red")) +
  geom_vline(xintercept=mm, color='yellow') + 
  theme(legend.position = "none") +
  theme(axis.title.x=element_blank()) +
  ggtitle("CSCP") +
  theme(plot.title = element_text(hjust = 0.5, size=12,  color="steelblue"))





##### Fluctuating Sel Fluctuating pop

ngens <- 500
pop_type <- "F"
gen_os <- NA
sel_type <- "F"

samp_id <- "S300"


output <- readRDS(file = paste0("../output/Sim_Output_S",sel_type,"_P",pop_type,"_",samp_id,".rds"))

ff <- as.data.frame(output$freqs[1:500,])
ff$Generation <- seq(1,500)


Long_genes <- ff %>%
  pivot_longer(-Generation, names_to="GeneId", values_to = "AlleleFrequency")

Long_genes$EffectSize <- rep(output$effects,times=500)
Long_genes$Init_freq <- rep(output$freqs[501,],times=500)



ex.ef <- Long_genes %>%
  filter(EffectSize > 5.5, 
         #Generation %in% seq(1,500, by = 25),
         Init_freq > 0.3 & Init_freq < 0.7)
ex.ef$Effect <- "Large"

sm.ef <- Long_genes %>%
  filter(EffectSize <0.5, Init_freq > 0.3 & Init_freq < 0.7)
sm.ef$Effect <- "Small"


ex.ef <- rbind(subset(ex.ef, GeneId %in% ex.ef$GeneId[c(1,2)]), subset(sm.ef, GeneId %in% sm.ef$GeneId[c(1,2)]))

FSFP <- Long_genes %>%
  filter(EffectSize > 0,
         Generation %in% seq(1,500, by = 1),
         Init_freq > 0.1 & Init_freq < 0.9
  ) %>%
  ggplot(aes(Generation, AlleleFrequency, group=GeneId, color=EffectSize)) +
  geom_line(alpha=1/2) +
  scale_colour_gradientn(colors=c("white","black")) + 
  theme(legend.position = "none") +
  theme(axis.title.y=element_blank()) +
  ggtitle("FSFP") +
  theme(plot.title = element_text(hjust = 0.5, size=12,  color="steelblue"))


  
  #geom_line(data=ex.ef, aes(Generation, AlleleFrequency, color=Effect)) +
  #scale_color_manual(values=c("blue","red")) 



##### Constant Sel Fluctuating pop

ngens <- 500
pop_type <- "F"
gen_os <- NA
sel_type <- "C"

samp_id <- "S301"


output <- readRDS(file = paste0("../output/Sim_Output_S",sel_type,"_P",pop_type,"_",samp_id,".rds"))

ff <- as.data.frame(output$freqs[1:500,])
ff$Generation <- seq(1,500)


Long_genes <- ff %>%
  pivot_longer(-Generation, names_to="GeneId", values_to = "AlleleFrequency")

Long_genes$EffectSize <- rep(output$effects,times=500)
Long_genes$Init_freq <- rep(output$freqs[501,],times=500)



ex.ef <- Long_genes %>%
  filter(EffectSize > 5.5, 
         #Generation %in% seq(1,500, by = 25),
         Init_freq > 0.3 & Init_freq < 0.7)
ex.ef$Effect <- "Large"

sm.ef <- Long_genes %>%
  filter(EffectSize <0.5, Init_freq > 0.3 & Init_freq < 0.7)
sm.ef$Effect <- "Small"


ex.ef <- rbind(subset(ex.ef, GeneId %in% ex.ef$GeneId[c(1,2)]), subset(sm.ef, GeneId %in% sm.ef$GeneId[c(1,2)]))

CSFP <- Long_genes %>%
  filter(EffectSize > 0,
         Generation %in% seq(1,500, by = 1),
         Init_freq > 0.1 & Init_freq < 0.9
  ) %>%
  ggplot(aes(Generation, AlleleFrequency, group=GeneId, color=EffectSize)) +
  geom_line(alpha=1/3) +
  scale_colour_gradientn(colors=c("white","black"))+
  geom_vline(xintercept=mm, color='yellow') + 
  theme(legend.position = "none") +
  theme(axis.title.x=element_blank(),axis.title.y=element_blank()) +
  ggtitle("CSFP") +
  theme(plot.title = element_text(hjust = 0.5, size=12,  color="steelblue"))


allele.all <- plot_grid(CSCP,CSFP,FSCP,FSFP, ncol=2, labels=c("a.","b.","c.","d."), 
                        rel_heights = c(1,1.1,1,1.1), rel_widths = c(1.1,1,1.1,1),
                        label_x=c(0,-0.06,0,-0.06))

ggsave(allele.all, filename="../output/afreq_sim.pdf",width=6.5, height=5)






#genome wide heterozygosity
#sliding window heterozygosity

slideH <- function(aa, step, ww)
{
  sts <- seq(1, (length(aa)-ww), by = step)
  oo <- vector(mode="numeric", length=length(sts))
  counter <- 1
  for(ii in sts) 
  {
   oo[counter] <- sum(aa[ii:(ii+(ww-1))])/ww
   counter <- counter+1
  }
  return(oo)
}

ngens <- 500
pop_type <- "C"
gen_os <- NA
sel_type <- "C"

samp_id <- "S6"


output <- readRDS(file = paste0("../output/Sim_Output_S",sel_type,"_P",pop_type,"_",samp_id,".rds"))

pp <- 2*output$freqs*(1-output$freqs)

gen.w <- rowSums(pp)/1000

plot(gen.w)

H1 <- slideH(pp[1,],1,20)
H2 <- slideH(pp[500,],1,20)

plot(H1, type='l', ylim=c(0.1,0.5))
lines(H2, col='blue')
abline(v=740)
abline(v=230)




