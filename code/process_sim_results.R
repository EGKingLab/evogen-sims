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
pop_type <- "C"
gen_os <- 25
sel_type <- "F"

samp_id <- "S4"


output <- readRDS(file = paste0("../output/Sim_Output_S",sel_type,"_P",pop_type,"_",samp_id,".rds"))

ff <- as.data.frame(output$freqs[1:100,1:50])
ff$Generation <- seq(1,100)

Long_genes <- pivot_longer(ff, -Generation, names_to="GeneId", values_to = "AlleleFrequency")



ggplot(Long_genes, aes(Generation, AlleleFrequency, group=GeneId)) +
  geom_line(alpha=1/6)


plot(output$freqs[1:50,1], type="l", ylim=c(0,1))
lines()
