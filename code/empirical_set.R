library(tidyverse)
library(cowplot)
theme_set(theme_cowplot())

MM <- read.table("../output/SNP_table_Yeast.txt",header=TRUE)
#So the “best” replicates are: 3, 7, 8, 9, 10, 11 and 12

samps <- paste0("maf_",seq(1,12),"_")
arms <- unique(MM$chr)

emp.set.all <- data.frame("samp"=character(length=0),
                          "chr" = character(length=0),
                          "freq"=numeric(length=0), "spec"= numeric(length=0))
for(aa in samps)
{
  for(bb in arms)
  {

M_rep <- MM %>% 
  select("chr","pos",starts_with("ancmaf")|starts_with(aa)) %>%
  filter(chr==bb)

Mspec <- as.matrix(t(M_rep[3:19]))

ss <- spectrum(Mspec[,sample(seq(1,ncol(Mspec)),1000)],  plot=FALSE)

ss1 <- rowMeans(ss$spec)

emp.set <- data.frame("samp"=rep(aa, length=length(ss1)),
                      "chr" = rep(bb,length=length(ss1)),
                      "freq"=(1/ss$freq)*30, "spec"= ss1)
emp.set.all <- rbind(emp.set.all, emp.set)
cat(aa,"\t",bb,"\n")
  }
  
}

emp.set.all$sample_chromosome <- paste0(emp.set.all$samp,emp.set.all$chr)

saveRDS(emp.set.all, file="../output/emp_data_spec.Rds")




emp.set.all <- data.frame("samp"=character(length=0),
                          "chr" = character(length=0),
                          "Generations"=numeric(length=0), "AbsoluteDifference"= numeric(length=0))
for(aa in samps)
{
  for(bb in arms)
  {
    
    M_rep <- MM %>% 
      select("chr","pos",starts_with("ancmaf")|starts_with(aa)) %>%
      filter(chr==bb)
    
    Mspec <- as.matrix(t(M_rep[3:19]))
    
    dd <- diff(Mspec, lag=1)
    
emp.set <- data.frame("samp"=rep(aa, length=nrow(dd)),
                      "chr" = rep(bb,length=nrow(dd)),
                      "Generations"=seq(1:16), "AbsoluteDifference"= rowMeans(abs(dd)))
emp.set.all <- rbind(emp.set.all, emp.set)
cat(aa,"\t",bb,"\n")
  }
  
}

emp.set.all$sample_chromosome <- paste0(emp.set.all$samp,emp.set.all$chr)

saveRDS(emp.set.all, file="../output/emp_data_diff.Rds")


    




emp.set.all %>% filter(freq <= 200) %>%
  ggplot(aes(freq, spec, group=sample_chromosome, color=samp)) +
  geom_line(alpha=1/2) #+ 
  #ylim(0,0.01) +
  #theme(legend.position = "none")
  

emp.set.all %>% filter(freq <= 200, samp %in% c("maf_3_","maf_7_","maf_10_","maf_11_")) %>%
  ggplot(aes(freq, spec, group=sample_chromosome, color=samp)) +
  geom_line(alpha=1/2) + 
  ylim(0,0.01) +
  xlab("Periodicity (Generations)") +
  ylab("Spectral Density") +
  theme(legend.position = "none")


dd.mat <- vector(mode='list', length=100)
for(jj in 1:16)
{
  
  dd <- diff(t(M_rep), lag=jj)
  dd.mat[[jj]] <- rowMeans((dd))
  
}

plot(dd.mat[[1]])

plot(dd.mat[[4]])

plot(dd.mat[[6]])


M_rep <- MM %>% 
  select("chr","pos",starts_with("ancmaf")|starts_with("maf_12_")) %>%
  filter(chr=="C02")

Long_genes <- M_rep %>%
  pivot_longer(!c(chr,pos), names_to="Generation", values_to = "AlleleFrequency")

Long_genes %>%
  ggplot(aes(Generation, AlleleFrequency, group=pos)) +
  geom_line(alpha=1/20) +
  scale_x_discrete(labels=as.character(seq(0,16)))



#match small set

samps <- paste0("maf_",seq(1,12),"_")
arms <- unique(MM$chr)

emp.set.all <- data.frame("samp"=character(length=0),
                          "chr" = character(length=0),
                          "Generations"=numeric(length=0), "AbsoluteDifference"= numeric(length=0))
for(aa in samps)
{
  for(bb in arms)
  {
    
    M_rep <- MM %>% 
      select("chr","pos",starts_with("ancmaf")|starts_with(aa)) %>%
      filter(chr==bb)
    
    foc.cols <- paste0(aa, c("6","12","18"))
    
    
    Mspec <- as.matrix(t(M_rep[,foc.cols]))
    
    dd <- diff(Mspec, lag=1)
    
    emp.set <- data.frame("samp"=rep(aa, length=nrow(dd)),
                          "chr" = rep(bb,length=nrow(dd)),
                          "Generations"=c(1,2), "AbsoluteDifference"= rowMeans(abs(dd)))
    emp.set.all <- rbind(emp.set.all, emp.set)
    cat(aa,"\t",bb,"\n")
  }
  
}

emp.set.all$sample_chromosome <- paste0(emp.set.all$samp,emp.set.all$chr)

saveRDS(emp.set.all, file="../output/emp_data_diff.Rds")





load(file="../output/small_sex.rda")
MS <- small_sex %>% 
  select("CHROM","POS",starts_with("freq_YEE_0112_04"))


samps <- c("freq_YEE_0112_04_01_",
           "freq_YEE_0112_04_02_",
           "freq_YEE_0112_04_03_",
           "freq_YEE_0112_04_10_",
           "freq_YEE_0112_04_11_",
           "freq_YEE_0112_04_12_")

arms <- unique(MS$CHROM)


emp.set.all <- data.frame("samp"=character(length=0),
                          "chr" = character(length=0),
                          "Generations"=numeric(length=0), "AbsoluteDifference"= numeric(length=0))
for(aa in samps)
{
  for(bb in arms)
  {
    
    M_rep <- MS %>% 
      select("CHROM","POS",starts_with(aa)) %>%
      filter(CHROM==bb)
    foc.cols <- paste0(aa, c("06","12","18"))
    Mspec <- as.matrix(t(M_rep[,foc.cols]))
    
    dd <- diff(Mspec, lag=1)
    
    Gens <- unlist(lapply(str_split(colnames(M_rep[,foc.cols]),"_"), function(x) x[6]))
    
    emp.set <- data.frame("samp"=rep(aa, length=nrow(dd)),
                          "chr" = rep(bb,length=nrow(dd)),
                          "Generations"=c(1,2), "AbsoluteDifference"= rowMeans(abs(dd)))
    emp.set.all <- rbind(emp.set.all, emp.set)
    cat(aa,"\t",bb,"\n")
  }
  
}

emp.set.all$sample_chromosome <- paste0(emp.set.all$samp,emp.set.all$chr)

saveRDS(emp.set.all, file="../output/emp_data_small_diffs.Rds")
