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

############################################################
############################################################
############################################################

########## CSCP ########
CSCP_Genom = read.csv("/Users/etb68/OneDrive - University of Missouri/Academics Mizzou (etb68@umsystem.edu)/Biologycal_Sciences/Research_Projects/Yeast_GitHub/evogen-sims/ET_Yeast/output.dir/DummyPop/WF_models/genome_CSCP_WF.csv")

CSCP_Pheno = read.csv("/Users/etb68/OneDrive - University of Missouri/Academics Mizzou (etb68@umsystem.edu)/Biologycal_Sciences/Research_Projects/Yeast_GitHub/evogen-sims/ET_Yeast/output.dir/DummyPop/WF_models/MeanPhenotypes_CSCP_WF.csv")

plot(CSCP_Pheno$Generation, CSCP_Pheno$Phenotype)
hist(CSCP_Pheno$Phenotype, breaks = 300)
CSCP_Genom %>%  ggplot(aes(Generation, Frequency, col = as.character(Effect)))+
  geom_point()+
  theme_cowplot(12)+
  theme(legend.position = "none")

CSCP_Pheno %>% ggplot(aes(Generation, Phenotype)) + geom_point() +  theme_cowplot(12)
CSCP_Pheno %>% ggplot(aes(Phenotype)) + geom_histogram(bins = 150) +  theme_cowplot(12)

########## FSCP1 ########

FSCP1_Genom = read.csv("/Users/etb68/OneDrive - University of Missouri/Academics Mizzou (etb68@umsystem.edu)/Biologycal_Sciences/Research_Projects/Yeast_GitHub/evogen-sims/ET_Yeast/output.dir/DummyPop/WF_models/genome_FSCP_WF.v1.csv")

FSCP1_Pheno = read.csv("/Users/etb68/OneDrive - University of Missouri/Academics Mizzou (etb68@umsystem.edu)/Biologycal_Sciences/Research_Projects/Yeast_GitHub/evogen-sims/ET_Yeast/output.dir/DummyPop/WF_models/MeanPhenotypes_FSCP_WF.v1.csv")

plot1_phenoGen = FSCP1_Pheno %>% ggplot(aes(Generation, Phenotype)) + geom_point() +  theme_cowplot(12)
plot1_phenoGen
plot1_pheno = FSCP1_Pheno %>% ggplot(aes(Phenotype)) + geom_histogram(bins = 500) +  theme_cowplot(12)
plot1_pheno

plot1_Genom = FSCP1_Genom %>%  ggplot(aes(Generation, Frequency, col = as.character(Effect)))+
  geom_point()+
  theme_cowplot(12)+
  theme(legend.position = "none")

########## FSCP2 ########

FSCP2_Genom = read.csv("/Users/etb68/OneDrive - University of Missouri/Academics Mizzou (etb68@umsystem.edu)/Biologycal_Sciences/Research_Projects/Yeast_GitHub/evogen-sims/ET_Yeast/output.dir/DummyPop/WF_models/genome_FSCP_WF.v2.1.csv")

FSCP2_Pheno = read.csv("/Users/etb68/OneDrive - University of Missouri/Academics Mizzou (etb68@umsystem.edu)/Biologycal_Sciences/Research_Projects/Yeast_GitHub/evogen-sims/ET_Yeast/output.dir/DummyPop/WF_models/MeanPhenotypes_FSCP_WF.v2.1.csv")

plot2_phenoGen = FSCP2_Pheno %>% ggplot(aes(Generation, Phenotype)) + geom_point() +  theme_cowplot(12)
plot2_pheno = FSCP2_Pheno %>% ggplot(aes(Phenotype)) + geom_histogram(bins = 500) +  theme_cowplot(12)

plot2_phenoGen
plot2_pheno

plot2_Genom = FSCP2_Genom %>%  ggplot(aes(Generation, Frequency, col = as.character(Effect)))+
  geom_point()+
  theme_cowplot(12)+
  theme(legend.position = "none")

######## FSCP3 ######
FSCP3_Genom = read.csv("/Users/etb68/OneDrive - University of Missouri/Academics Mizzou (etb68@umsystem.edu)/Biologycal_Sciences/Research_Projects/Yeast_GitHub/evogen-sims/ET_Yeast/output.dir/DummyPop/WF_models/genome_FSCP_WF.v2.2.csv")

FSCP3_Pheno = read.csv("/Users/etb68/OneDrive - University of Missouri/Academics Mizzou (etb68@umsystem.edu)/Biologycal_Sciences/Research_Projects/Yeast_GitHub/evogen-sims/ET_Yeast/output.dir/DummyPop/WF_models/MeanPhenotypes_FSCP_WF.v2.2.csv")

plot3_phenoGen = FSCP3_Pheno %>% ggplot(aes(Generation, Phenotype)) + geom_point() +  theme_cowplot(12)
plot3_pheno = FSCP3_Pheno %>% ggplot(aes(Phenotype)) + geom_histogram(bins = 500) +  theme_cowplot(12)

plot3_phenoGen
plot3_pheno

plot3_Genom = FSCP1_Genom %>%  ggplot(aes(Generation, Frequency, col = as.character(Effect)))+
  geom_point()+
  theme_cowplot(12)+
  theme(legend.position = "none")

################### FSCP4 ##############################

FSCP4_Genom = read.csv("/Users/etb68/OneDrive - University of Missouri/Academics Mizzou (etb68@umsystem.edu)/Biologycal_Sciences/Research_Projects/Yeast_GitHub/evogen-sims/ET_Yeast/output.dir/DummyPop/WF_models/genome_FSCP_WF.v3.csv")

FSCP4_Pheno = read.csv("/Users/etb68/OneDrive - University of Missouri/Academics Mizzou (etb68@umsystem.edu)/Biologycal_Sciences/Research_Projects/Yeast_GitHub/evogen-sims/ET_Yeast/output.dir/DummyPop/WF_models/MeanPhenotypes_FSCP_WF.v3.csv")

plot4_phenoGen = FSCP4_Pheno %>% ggplot(aes(Generation, Phenotype)) + geom_point() +  theme_cowplot(12)
plot4_pheno = FSCP4_Pheno %>% ggplot(aes(Phenotype)) + geom_histogram(bins = 250) +  theme_cowplot(12)

plot4_phenoGen
plot4_pheno

plot4_Genom = FSCP4_Genom %>%  ggplot(aes(Generation, Frequency, col = as.character(Effect)))+
  geom_point()+
  theme_cowplot(12)+
  theme(legend.position = "none")

################### FSCP5 ##############################

FSCP5_Genom = read.csv("/Users/etb68/OneDrive - University of Missouri/Academics Mizzou (etb68@umsystem.edu)/Biologycal_Sciences/Research_Projects/Yeast_GitHub/evogen-sims/ET_Yeast/output.dir/DummyPop/WF_models/genome_FSCP_WF.v4.csv")

FSCP5_Pheno = read.csv("/Users/etb68/OneDrive - University of Missouri/Academics Mizzou (etb68@umsystem.edu)/Biologycal_Sciences/Research_Projects/Yeast_GitHub/evogen-sims/ET_Yeast/output.dir/DummyPop/WF_models/MeanPhenotypes_FSCP_WF.v4.csv")

plot5_phenoGen = FSCP5_Pheno %>% ggplot(aes(Generation, Phenotype)) + geom_point() +  theme_cowplot(12)
plot5_pheno = FSCP5_Pheno %>% ggplot(aes(Phenotype)) + geom_histogram(bins = 300) +  theme_cowplot(12)

plot5_phenoGen
plot5_pheno

plot5_Genom = FSCP5_Genom %>%  ggplot(aes(Generation, Frequency, col = as.character(Effect)))+
  geom_point()+
  theme_cowplot(12)+
  theme(legend.position = "none")
plot5_Genom

############################################################################################
############################################################################################
############################################################################################

####### CSFP ###

