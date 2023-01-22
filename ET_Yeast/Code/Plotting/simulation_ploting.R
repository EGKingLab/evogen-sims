library(tidyverse)
library(cowplot)

empty_geno = read.csv("/Users/ezra/OneDrive - University of Missouri/Academics Mizzou (etb68@umsystem.edu)/Biologycal_Sciences/Research_Projects/Yeast_GitHub/evogen-sims/ET_Yeast/output.dir/Empty_geno/genome_CSCP.csv") #%>% 
  #filter(Frequency < 1 & Frequency > 0)

empty_geno_ph = read.csv("/Users/ezra/OneDrive - University of Missouri/Academics Mizzou (etb68@umsystem.edu)/Biologycal_Sciences/Research_Projects/Yeast_GitHub/evogen-sims/ET_Yeast/output.dir/Empty_geno/MeanPhenotypes_CSCP.csv")
plot(empty_geno_ph$Generation, empty_geno_ph$Phenotype)

empty_geno %>%  ggplot(aes(Generation, Frequency, col = as.character(Effect)))+
  geom_point()+
  theme_cowplot(12)+
  theme(legend.position = "none")

dummy = read.csv("/Users/etb68/OneDrive - University of Missouri/Academics Mizzou (etb68@umsystem.edu)/Biologycal_Sciences/Research_Projects/Yeast_GitHub/evogen-sims/ET_Yeast/output.dir/DummyPop/genome_CSCP_MutRate.csv") #%>% 
  #filter(Frequency < 1 & Frequency > 0)

dummy_ph = read.csv("/Users/etb68/OneDrive - University of Missouri/Academics Mizzou (etb68@umsystem.edu)/Biologycal_Sciences/Research_Projects/Yeast_GitHub/evogen-sims/ET_Yeast/output.dir/DummyPop/MeanPhenotypes_CSCP_MutRate.csv")

plot(dummy_ph$Generation, dummy_ph$Phenotype)


dummy %>%  ggplot(aes(Generation, Frequency, col = as.character(Effect)))+
  geom_point()+
  theme_cowplot(12)+
  theme(legend.position = "none")

dummy_norate = read.csv("/Users/ezra/OneDrive - University of Missouri/Academics Mizzou (etb68@umsystem.edu)/Biologycal_Sciences/Research_Projects/Yeast_GitHub/evogen-sims/ET_Yeast/output.dir/DummyPop/genome_CSCP_NoMut.csv")%>% 
  filter(Frequency < 1 & Frequency > 0)
dummy_norate_ph = read.csv("/Users/ezra/OneDrive - University of Missouri/Academics Mizzou (etb68@umsystem.edu)/Biologycal_Sciences/Research_Projects/Yeast_GitHub/evogen-sims/ET_Yeast/output.dir/DummyPop/MeanPhenotypes_CSCP_MutRate.csv")

plot(dummy_norate_ph$Generation, dummy_norate_ph$Phenotype)

dummy_norate %>%  ggplot(aes(Generation, Frequency, col = as.character(Effect)))+
  geom_point()+
  theme_cowplot(12)+
  theme(legend.position = "none")

test = read.csv("/Users/ezra/OneDrive - University of Missouri/Academics Mizzou (etb68@umsystem.edu)/Biologycal_Sciences/Research_Projects/Yeast_GitHub/evogen-sims/ET_Yeast/code/test.dir/genome_test.csv") %>% 
  filter(Frequency < 1 & Frequency > 0)

test_pheno = read.csv("/Users/ezra/OneDrive - University of Missouri/Academics Mizzou (etb68@umsystem.edu)/Biologycal_Sciences/Research_Projects/Yeast_GitHub/evogen-sims/ET_Yeast/code/test.dir/MeanPhenotypes_test.csv") %>% 
  mutate(fitnessScaling = (1 + (Phenotype - Optimum)/200))

plot(test_pheno$Generation, test_pheno$Phenotype)
hist(test_pheno$Phenotype)
hist(test_pheno$Optimum)

test %>%  ggplot(aes(Generation, Frequency, col = as.character(Effect)))+
  geom_point()+
  theme_cowplot(12)+
  theme(legend.position = "none")

plot(test$Generation, test$Frequency)


#####################
library(tidyverse)
library(cowplot)
trial = read.csv("/Users/etb68/OneDrive - University of Missouri/Academics Mizzou (etb68@umsystem.edu)/Biologycal_Sciences/Research_Projects/Yeast_GitHub/evogen-sims/ET_Yeast/output.dir/DummyPop/WF_models/genome_CSFP_WF.csv") #%>% 
  #filter(Frequency < 1 & Frequency > 0)

trial_pheno = read.csv("/Users/etb68/OneDrive - University of Missouri/Academics Mizzou (etb68@umsystem.edu)/Biologycal_Sciences/Research_Projects/Yeast_GitHub/evogen-sims/ET_Yeast/output.dir/DummyPop/WF_models/MeanPhenotypes_CSFP_WF.csv") 

plot(trial_pheno$Generation, trial_pheno$Phenotype)
hist(trial_pheno$Phenotype, breaks = 300)
#hist(trial_pheno$Optimum)

cor.test(trial_pheno$Phenotype,trial_pheno$Optimum)
var(trial_pheno$Phenotype)

trial %>%  ggplot(aes(Generation, Frequency, col = as.character(Effect)))+
  geom_point()+
  theme_cowplot(12)+
  theme(legend.position = "none")

hist(trial$Effect, breaks = 100000)

plot(trial$Generation, trial$Frequency)

