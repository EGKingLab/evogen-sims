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

### Constant Selection & fluctuating pop between 10000 and 1000

CSFPo_Genom = read.csv("/Users/etb68/OneDrive - University of Missouri/Academics Mizzou (etb68@umsystem.edu)/Biologycal_Sciences/Research_Projects/Yeast_GitHub/evogen-sims/ET_Yeast/output.dir/DummyPop/WF_models/genome_CSFP_WF.csv")

CSFPo_Pheno = read.csv("/Users/etb68/OneDrive - University of Missouri/Academics Mizzou (etb68@umsystem.edu)/Biologycal_Sciences/Research_Projects/Yeast_GitHub/evogen-sims/ET_Yeast/output.dir/DummyPop/WF_models/MeanPhenotypes_CSFP_WF.csv")

CSFPo_PhGen_plot = CSFPo_Pheno %>% ggplot(aes(Generation, Phenotype)) + geom_point() +  theme_cowplot(12)
CSFPo_Ph_plot = CSFPo_Pheno %>% ggplot(aes(Phenotype)) + geom_histogram(bins = 150) +  theme_cowplot(12)

CSFPo_PhGen_plot
CSFPo_Ph_plot

CSFPo_Geno_plot = CSFPo_Genom %>%  ggplot(aes(Generation, Frequency, col = as.character(Effect)))+
  geom_point()+
  theme_cowplot(12)+
  theme(legend.position = "none")

CSFPo_Geno_plot

### Constant Selection & fluctuating pop following uniform distribution


CSFPu_Genom = read.csv("/Users/etb68/OneDrive - University of Missouri/Academics Mizzou (etb68@umsystem.edu)/Biologycal_Sciences/Research_Projects/Yeast_GitHub/evogen-sims/ET_Yeast/output.dir/DummyPop/WF_models/Uniform_Pop/genome_CSFPu_WF.csv")

CSFPu_Pheno = read.csv("/Users/etb68/OneDrive - University of Missouri/Academics Mizzou (etb68@umsystem.edu)/Biologycal_Sciences/Research_Projects/Yeast_GitHub/evogen-sims/ET_Yeast/output.dir/DummyPop/WF_models/Uniform_Pop/MeanPhenotypes_CSFPu_WF.csv")

CSFPu_PhGen_plot = CSFPu_Pheno %>% ggplot(aes(Generation, Phenotype)) + geom_point() +  theme_cowplot(12) + ggtitle("Unif CSFP Pheno vs Generation Unif")
CSFPu_Ph_plot = CSFPu_Pheno %>% ggplot(aes(Phenotype)) + geom_histogram(bins = 150) +  theme_cowplot(12) + ggtitle("Unif CSFP Pheno")

CSFPu_PhGen_plot
CSFPu_Ph_plot

CSFPu_Geno_plot = CSFPu_Genom %>%  ggplot(aes(Generation, Frequency, col = as.character(Effect)))+
  geom_point()+
  theme_cowplot(12)+
  theme(legend.position = "none") + 
  ggtitle("Unif CSFP Allele Freq")

CSFPu_Geno_plot

### Constant Selection & fluctuating pop following Poisson distribution


CSFPp_Genom = read.csv("/Users/etb68/OneDrive - University of Missouri/Academics Mizzou (etb68@umsystem.edu)/Biologycal_Sciences/Research_Projects/Yeast_GitHub/evogen-sims/ET_Yeast/output.dir/DummyPop/WF_models/Pois_Pop/genome_CSFPp_WF.csv")

CSFPp_Pheno = read.csv("/Users/etb68/OneDrive - University of Missouri/Academics Mizzou (etb68@umsystem.edu)/Biologycal_Sciences/Research_Projects/Yeast_GitHub/evogen-sims/ET_Yeast/output.dir/DummyPop/WF_models/Pois_Pop/MeanPhenotypes_CSFPp_WF.csv")

CSFPp_PhGen_plot = CSFPp_Pheno %>% ggplot(aes(Generation, Phenotype)) + geom_point() +  theme_cowplot(12) + ggtitle("Pois CSFP Pheno vs Generation ")
CSFPp_Ph_plot = CSFPp_Pheno %>% ggplot(aes(Phenotype)) + geom_histogram(bins = 150) +  theme_cowplot(12) + ggtitle("Pois CSFP Pheno")

CSFPp_PhGen_plot
CSFPp_Ph_plot

CSFPp_Geno_plot = CSFPp_Genom %>%  ggplot(aes(Generation, Frequency, col = as.character(Effect)))+
  geom_point()+
  theme_cowplot(12)+
  theme(legend.position = "none") + 
  ggtitle("Pois CSFP Allele Freq")

CSFPp_Geno_plot

################################################################
################ FSFP1 ###########################
################################################################

######Initial Pop (10000 and 1000) ###### 

FSFPo1_Genom = read.csv("/Users/etb68/OneDrive - University of Missouri/Academics Mizzou (etb68@umsystem.edu)/Biologycal_Sciences/Research_Projects/Yeast_GitHub/evogen-sims/ET_Yeast/output.dir/DummyPop/WF_models/genome_FSFP_WF.v1.csv")

FSFPo1_pheno = read.csv("/Users/etb68/OneDrive - University of Missouri/Academics Mizzou (etb68@umsystem.edu)/Biologycal_Sciences/Research_Projects/Yeast_GitHub/evogen-sims/ET_Yeast/output.dir/DummyPop/WF_models/MeanPhenotypes_FSFP_WF.v1.csv")

FSFPo1_PhGen_plot = FSFPo1_pheno %>% ggplot(aes(Generation, Phenotype)) + geom_point() +  theme_cowplot(12)+ ggtitle("FSFPo1 Pheno vs Generation")
FSFPo1_Ph_plot = FSFPo1_pheno %>% ggplot(aes(Phenotype)) + geom_histogram(bins = 150) +  theme_cowplot(12) + ggtitle("FSFPo1 Pheno")

FSFPo1_Geno_plot = FSFPo1_Genom %>%  ggplot(aes(Generation, Frequency, col = as.character(Effect)))+
  geom_point()+
  theme_cowplot(12)+
  theme(legend.position = "none")+
  ggtitle("FSFPo1 Genome")

FSFPo1_PhGen_plot
FSFPo1_Ph_plot
FSFPo1_Geno_plot

###### Uniform Pop (Min = 100, Max = 10,000) ###### 


FSFPu1_Genom = read.csv("/Users/etb68/OneDrive - University of Missouri/Academics Mizzou (etb68@umsystem.edu)/Biologycal_Sciences/Research_Projects/Yeast_GitHub/evogen-sims/ET_Yeast/output.dir/DummyPop/WF_models/Uniform_Pop/genome_FSFPu_WF.v1.csv")

FSFPu1_pheno = read.csv("/Users/etb68/OneDrive - University of Missouri/Academics Mizzou (etb68@umsystem.edu)/Biologycal_Sciences/Research_Projects/Yeast_GitHub/evogen-sims/ET_Yeast/output.dir/DummyPop/WF_models/Uniform_Pop/MeanPhenotypes_FSFPu_WF.v1.csv")

FSFPu1_PhGen_plot = FSFPu1_pheno %>% ggplot(aes(Generation, Phenotype)) + geom_point() +  theme_cowplot(12)+ ggtitle("FSFPu1 Pheno vs Generation")
FSFPu1_Ph_plot = FSFPu1_pheno %>% ggplot(aes(Phenotype)) + geom_histogram(bins = 150) +  theme_cowplot(12) + ggtitle("FSFPu1 Pheno")

FSFPu1_Geno_plot = FSFPu1_Genom %>%  ggplot(aes(Generation, Frequency, col = as.character(Effect)))+
  geom_point()+
  theme_cowplot(12)+
  theme(legend.position = "none")+
  ggtitle("FSFPu1 Genome")

FSFPu1_PhGen_plot
FSFPu1_Ph_plot
FSFPu1_Geno_plot


###### Poisson Pop (with mean = 10,000) ###### 


FSFPp1_Genom = read.csv("/Users/etb68/OneDrive - University of Missouri/Academics Mizzou (etb68@umsystem.edu)/Biologycal_Sciences/Research_Projects/Yeast_GitHub/evogen-sims/ET_Yeast/output.dir/DummyPop/WF_models/Pois_Pop/genome_FSFPp_WF.v1.csv")

FSFPp1_pheno = read.csv("/Users/etb68/OneDrive - University of Missouri/Academics Mizzou (etb68@umsystem.edu)/Biologycal_Sciences/Research_Projects/Yeast_GitHub/evogen-sims/ET_Yeast/output.dir/DummyPop/WF_models/Pois_Pop/MeanPhenotypes_FSFPp_WF.v1.csv")

FSFPp1_PhGen_plot = FSFPp1_pheno %>% ggplot(aes(Generation, Phenotype)) + geom_point() +  theme_cowplot(12)+ ggtitle("FSFPp1 Pheno vs Generation")
FSFPp1_Ph_plot = FSFPp1_pheno %>% ggplot(aes(Phenotype)) + geom_histogram(bins = 150) +  theme_cowplot(12) + ggtitle("FSFPp1 Pheno")

FSFPp1_Geno_plot = FSFPp1_Genom %>%  ggplot(aes(Generation, Frequency, col = as.character(Effect)))+
  geom_point()+
  theme_cowplot(12)+
  theme(legend.position = "none")+
  ggtitle("FSFPp1 Genome")

FSFPp1_PhGen_plot
FSFPp1_Ph_plot
FSFPp1_Geno_plot

################################################################
##################### FSFP2 ###################################
################################################################

###### Initial Pop (10000 and 1000) ###### 

FSFPo2_Genom = read.csv("/Users/etb68/OneDrive - University of Missouri/Academics Mizzou (etb68@umsystem.edu)/Biologycal_Sciences/Research_Projects/Yeast_GitHub/evogen-sims/ET_Yeast/output.dir/DummyPop/WF_models/genome_FSFP_WF.v2.1.csv")

FSFPo2_pheno = read.csv("/Users/etb68/OneDrive - University of Missouri/Academics Mizzou (etb68@umsystem.edu)/Biologycal_Sciences/Research_Projects/Yeast_GitHub/evogen-sims/ET_Yeast/output.dir/DummyPop/WF_models/MeanPhenotypes_FSFP_WF.v2.1.csv")

FSFPo2_PhGen_plot = FSFPo2_pheno %>% ggplot(aes(Generation, Phenotype)) + geom_point() +  theme_cowplot(12)+ ggtitle("FSFPo2 Pheno vs Generation")
FSFPo2_Ph_plot = FSFPo2_pheno %>% ggplot(aes(Phenotype)) + geom_histogram(bins = 150) +  theme_cowplot(12) + ggtitle("FSFPo2 Pheno")

FSFPo2_Geno_plot = FSFPo2_Genom %>%  ggplot(aes(Generation, Frequency, col = as.character(Effect)))+
  geom_point()+
  theme_cowplot(12)+
  theme(legend.position = "none")+
  ggtitle("FSFPo2 Genome")

FSFPo2_PhGen_plot
FSFPo2_Ph_plot
FSFPo2_Geno_plot

###### Uniform Pop (with mean = 10,000) ###### 


FSFPu2_Genom = read.csv("/Users/etb68/OneDrive - University of Missouri/Academics Mizzou (etb68@umsystem.edu)/Biologycal_Sciences/Research_Projects/Yeast_GitHub/evogen-sims/ET_Yeast/output.dir/DummyPop/WF_models/Uniform_Pop/genome_FSFPu_WF.v2.1.csv")

FSFPu2_pheno = read.csv("/Users/etb68/OneDrive - University of Missouri/Academics Mizzou (etb68@umsystem.edu)/Biologycal_Sciences/Research_Projects/Yeast_GitHub/evogen-sims/ET_Yeast/output.dir/DummyPop/WF_models/Uniform_Pop/MeanPhenotypes_FSFPu_WF.v2.1.csv")

FSFPu2_PhGen_plot = FSFPu2_pheno %>% ggplot(aes(Generation, Phenotype)) + geom_point() +  theme_cowplot(12)+ ggtitle("FSFPu2 Pheno vs Generation")
FSFPu2_Ph_plot = FSFPu2_pheno %>% ggplot(aes(Phenotype)) + geom_histogram(bins = 150) +  theme_cowplot(12) + ggtitle("FSFPu2 Pheno")

FSFPu2_Geno_plot = FSFPu2_Genom %>%  ggplot(aes(Generation, Frequency, col = as.character(Effect)))+
  geom_point()+
  theme_cowplot(12)+
  theme(legend.position = "none")+
  ggtitle("FSFPu2 Genome")

FSFPu2_PhGen_plot
FSFPu2_Ph_plot
FSFPu2_Geno_plot

###### Poisson Pop (with mean = 10,000) ###### 


FSFPp2_Genom = read.csv("/Users/etb68/OneDrive - University of Missouri/Academics Mizzou (etb68@umsystem.edu)/Biologycal_Sciences/Research_Projects/Yeast_GitHub/evogen-sims/ET_Yeast/output.dir/DummyPop/WF_models/Pois_Pop/genome_FSFPp_WF.v2.1.csv")

FSFPp2_pheno = read.csv("/Users/etb68/OneDrive - University of Missouri/Academics Mizzou (etb68@umsystem.edu)/Biologycal_Sciences/Research_Projects/Yeast_GitHub/evogen-sims/ET_Yeast/output.dir/DummyPop/WF_models/Pois_Pop/MeanPhenotypes_FSFPp_WF.v2.1.csv")

FSFPp2_PhGen_plot = FSFPp1_pheno %>% ggplot(aes(Generation, Phenotype)) + geom_point() +  theme_cowplot(12)+ ggtitle("FSFPp2 Pheno vs Generation")
FSFPp2_Ph_plot = FSFPp2_pheno %>% ggplot(aes(Phenotype)) + geom_histogram(bins = 150) +  theme_cowplot(12) + ggtitle("FSFPp2 Pheno")

FSFPp2_Geno_plot = FSFPp2_Genom %>%  ggplot(aes(Generation, Frequency, col = as.character(Effect)))+
  geom_point()+
  theme_cowplot(12)+
  theme(legend.position = "none")+
  ggtitle("FSFPp2 Genome")

FSFPp2_PhGen_plot
FSFPp2_Ph_plot
FSFPp2_Geno_plot

################################################################
##################### FSFP3 ###################################
################################################################

######Initial Pop (10000 and 1000) ###### 

FSFPo3_Genom = read.csv("/Users/etb68/OneDrive - University of Missouri/Academics Mizzou (etb68@umsystem.edu)/Biologycal_Sciences/Research_Projects/Yeast_GitHub/evogen-sims/ET_Yeast/output.dir/DummyPop/WF_models/genome_FSFP_WF.v2.2.csv")

FSFPo3_pheno = read.csv("/Users/etb68/OneDrive - University of Missouri/Academics Mizzou (etb68@umsystem.edu)/Biologycal_Sciences/Research_Projects/Yeast_GitHub/evogen-sims/ET_Yeast/output.dir/DummyPop/WF_models/MeanPhenotypes_FSFP_WF.v2.2.csv")

FSFPo3_PhGen_plot = FSFPo3_pheno %>% ggplot(aes(Generation, Phenotype)) + geom_point() +  theme_cowplot(12)+ ggtitle("FSFPo3 Pheno vs Generation")
FSFPo3_Ph_plot = FSFPo3_pheno %>% ggplot(aes(Phenotype)) + geom_histogram(bins = 150) +  theme_cowplot(12) + ggtitle("FSFPo3 Pheno")

FSFPo3_Geno_plot = FSFPo3_Genom %>%  ggplot(aes(Generation, Frequency, col = as.character(Effect)))+
  geom_point()+
  theme_cowplot(12)+
  theme(legend.position = "none")+
  ggtitle("FSFPo3 Genome")

FSFPo3_PhGen_plot
FSFPo3_Ph_plot
FSFPo3_Geno_plot

###### Uniform Pop (with mean = 10,000) ###### 


FSFPu3_Genom = read.csv("/Users/etb68/OneDrive - University of Missouri/Academics Mizzou (etb68@umsystem.edu)/Biologycal_Sciences/Research_Projects/Yeast_GitHub/evogen-sims/ET_Yeast/output.dir/DummyPop/WF_models/Uniform_Pop/genome_FSFPu_WF.v2.2.csv")

FSFPu3_pheno = read.csv("/Users/etb68/OneDrive - University of Missouri/Academics Mizzou (etb68@umsystem.edu)/Biologycal_Sciences/Research_Projects/Yeast_GitHub/evogen-sims/ET_Yeast/output.dir/DummyPop/WF_models/Uniform_Pop/MeanPhenotypes_FSFPu_WF.v2.2.csv")

FSFPu3_PhGen_plot = FSFPu3_pheno %>% ggplot(aes(Generation, Phenotype)) + geom_point() +  theme_cowplot(12)+ ggtitle("FSFPu3 Pheno vs Generation")
FSFPu3_Ph_plot = FSFPu3_pheno %>% ggplot(aes(Phenotype)) + geom_histogram(bins = 150) +  theme_cowplot(12) + ggtitle("FSFPu3 Pheno")

FSFPu3_Geno_plot = FSFPu3_Genom %>%  ggplot(aes(Generation, Frequency, col = as.character(Effect)))+
  geom_point()+
  theme_cowplot(12)+
  theme(legend.position = "none")+
  ggtitle("FSFPu3 Genome")

FSFPu3_PhGen_plot
FSFPu3_Ph_plot
FSFPu3_Geno_plot


###### Poisson Pop (with mean = 10,000) ###### 


FSFPp3_Genom = read.csv("/Users/etb68/OneDrive - University of Missouri/Academics Mizzou (etb68@umsystem.edu)/Biologycal_Sciences/Research_Projects/Yeast_GitHub/evogen-sims/ET_Yeast/output.dir/DummyPop/WF_models/Pois_Pop/genome_FSFPp_WF.v2.2.csv")

FSFPp3_pheno = read.csv("/Users/etb68/OneDrive - University of Missouri/Academics Mizzou (etb68@umsystem.edu)/Biologycal_Sciences/Research_Projects/Yeast_GitHub/evogen-sims/ET_Yeast/output.dir/DummyPop/WF_models/Pois_Pop/MeanPhenotypes_FSFPp_WF.v2.2.csv")

FSFPp3_PhGen_plot = FSFPp1_pheno %>% ggplot(aes(Generation, Phenotype)) + geom_point() +  theme_cowplot(12)+ ggtitle("FSFPp3 Pheno vs Generation")
FSFPp3_Ph_plot = FSFPp3_pheno %>% ggplot(aes(Phenotype)) + geom_histogram(bins = 150) +  theme_cowplot(12) + ggtitle("FSFPp3 Pheno")

FSFPp3_Geno_plot = FSFPp3_Genom %>%  ggplot(aes(Generation, Frequency, col = as.character(Effect)))+
  geom_point()+
  theme_cowplot(12)+
  theme(legend.position = "none")+
  ggtitle("FSFPp3 Genome")

FSFPp3_PhGen_plot
FSFPp3_Ph_plot
FSFPp3_Geno_plot

################################################################
##################### FSFP4 ###################################
################################################################

######Initial Pop (10000 and 1000) ###### 

FSFPo4_Genom = read.csv("/Users/etb68/OneDrive - University of Missouri/Academics Mizzou (etb68@umsystem.edu)/Biologycal_Sciences/Research_Projects/Yeast_GitHub/evogen-sims/ET_Yeast/output.dir/DummyPop/WF_models/genome_FSFP_WF.v3.csv")

FSFPo4_pheno = read.csv("/Users/etb68/OneDrive - University of Missouri/Academics Mizzou (etb68@umsystem.edu)/Biologycal_Sciences/Research_Projects/Yeast_GitHub/evogen-sims/ET_Yeast/output.dir/DummyPop/WF_models/MeanPhenotypes_FSFP_WF.v3.csv")

FSFPo4_PhGen_plot = FSFPo4_pheno %>% ggplot(aes(Generation, Phenotype)) + geom_point() +  theme_cowplot(12)+ ggtitle("FSFPo4 Pheno vs Generation")
FSFPo4_Ph_plot = FSFPo4_pheno %>% ggplot(aes(Phenotype)) + geom_histogram(bins = 150) +  theme_cowplot(12) + ggtitle("FSFPo4 Pheno")

FSFPo4_Geno_plot = FSFPo4_Genom %>%  ggplot(aes(Generation, Frequency, col = as.character(Effect)))+
  geom_point()+
  theme_cowplot(12)+
  theme(legend.position = "none")+
  ggtitle("FSFPo4 Genome")

FSFPo4_PhGen_plot
FSFPo4_Ph_plot
FSFPo4_Geno_plot

###### Uniform Pop (with mean = 10,000) ###### 


FSFPu4_Genom = read.csv("/Users/etb68/OneDrive - University of Missouri/Academics Mizzou (etb68@umsystem.edu)/Biologycal_Sciences/Research_Projects/Yeast_GitHub/evogen-sims/ET_Yeast/output.dir/DummyPop/WF_models/Uniform_Pop/genome_FSFPu_WF.v3.csv")

FSFPu4_pheno = read.csv("/Users/etb68/OneDrive - University of Missouri/Academics Mizzou (etb68@umsystem.edu)/Biologycal_Sciences/Research_Projects/Yeast_GitHub/evogen-sims/ET_Yeast/output.dir/DummyPop/WF_models/Uniform_Pop/MeanPhenotypes_FSFPu_WF.V3.csv")

FSFPu4_PhGen_plot = FSFPu4_pheno %>% ggplot(aes(Generation, Phenotype)) + geom_point() +  theme_cowplot(12)+ ggtitle("FSFPu4 Pheno vs Generation")
FSFPu4_Ph_plot = FSFPu4_pheno %>% ggplot(aes(Phenotype)) + geom_histogram(bins = 150) +  theme_cowplot(12) + ggtitle("FSFPu4 Pheno")

FSFPu4_Geno_plot = FSFPu4_Genom %>%  ggplot(aes(Generation, Frequency, col = as.character(Effect)))+
  geom_point()+
  theme_cowplot(12)+
  theme(legend.position = "none")+
  ggtitle("FSFPu4 Genome")

FSFPu4_PhGen_plot
FSFPu4_Ph_plot
FSFPu4_Geno_plot



###### Poisson Pop (with mean = 10,000) ###### 


FSFPp4_Genom = read.csv("/Users/etb68/OneDrive - University of Missouri/Academics Mizzou (etb68@umsystem.edu)/Biologycal_Sciences/Research_Projects/Yeast_GitHub/evogen-sims/ET_Yeast/output.dir/DummyPop/WF_models/Pois_Pop/genome_FSFP_WFp.v3.csv")

FSFPp4_pheno = read.csv("/Users/etb68/OneDrive - University of Missouri/Academics Mizzou (etb68@umsystem.edu)/Biologycal_Sciences/Research_Projects/Yeast_GitHub/evogen-sims/ET_Yeast/output.dir/DummyPop/WF_models/Pois_Pop/MeanPhenotypes_FSFPp_WF.v3.csv")

FSFPp4_PhGen_plot = FSFPp4_pheno %>% ggplot(aes(Generation, Phenotype)) + geom_point() +  theme_cowplot(12)+ ggtitle("FSFPp4 Pheno vs Generation")
FSFPp4_Ph_plot = FSFPp4_pheno %>% ggplot(aes(Phenotype)) + geom_histogram(bins = 150) +  theme_cowplot(12) + ggtitle("FSFPp4 Pheno")

FSFPp4_Geno_plot = FSFPp4_Genom %>%  ggplot(aes(Generation, Frequency, col = as.character(Effect)))+
  geom_point()+
  theme_cowplot(12)+
  theme(legend.position = "none")+
  ggtitle("FSFPp4 Genome")

FSFPp4_PhGen_plot
FSFPp4_Ph_plot
FSFPp4_Geno_plot

################################################################
##################### FSFP5 ###################################
################################################################

######Initial Pop (10000 and 1000) ###### 

FSFPo5_Genom = read.csv("/Users/etb68/OneDrive - University of Missouri/Academics Mizzou (etb68@umsystem.edu)/Biologycal_Sciences/Research_Projects/Yeast_GitHub/evogen-sims/ET_Yeast/output.dir/DummyPop/WF_models/genome_FSFP_WF.v4.csv")

FSFPo5_pheno = read.csv("/Users/etb68/OneDrive - University of Missouri/Academics Mizzou (etb68@umsystem.edu)/Biologycal_Sciences/Research_Projects/Yeast_GitHub/evogen-sims/ET_Yeast/output.dir/DummyPop/WF_models/MeanPhenotypes_FSFP_WF.v4.csv")

FSFPo5_PhGen_plot = FSFPo5_pheno %>% ggplot(aes(Generation, Phenotype)) + geom_point() +  theme_cowplot(12)+ ggtitle("FSFPo5 Pheno vs Generation")
FSFPo5_Ph_plot = FSFPo5_pheno %>% ggplot(aes(Phenotype)) + geom_histogram(bins = 150) +  theme_cowplot(12) + ggtitle("FSFPo5 Pheno")

FSFPo5_Geno_plot = FSFPo5_Genom %>%  ggplot(aes(Generation, Frequency, col = as.character(Effect)))+
  geom_point()+
  theme_cowplot(12)+
  theme(legend.position = "none")+
  ggtitle("FSFPo5 Genome")

FSFPo5_PhGen_plot
FSFPo5_Ph_plot
FSFPo5_Geno_plot

###### Uniform Pop (with mean = 10,000) ###### 


FSFPu5_Genom = read.csv("/Users/etb68/OneDrive - University of Missouri/Academics Mizzou (etb68@umsystem.edu)/Biologycal_Sciences/Research_Projects/Yeast_GitHub/evogen-sims/ET_Yeast/output.dir/DummyPop/WF_models/Uniform_Pop/genome_FSFPu_WF.v4.csv")

FSFPu5_pheno = read.csv("/Users/etb68/OneDrive - University of Missouri/Academics Mizzou (etb68@umsystem.edu)/Biologycal_Sciences/Research_Projects/Yeast_GitHub/evogen-sims/ET_Yeast/output.dir/DummyPop/WF_models/Uniform_Pop/MeanPhenotypes_FSFPu_WF.v4.csv")

FSFPu5_PhGen_plot = FSFPu5_pheno %>% ggplot(aes(Generation, Phenotype)) + geom_point() +  theme_cowplot(12)+ ggtitle("FSFPu5 Pheno vs Generation")
FSFPu5_Ph_plot = FSFPu5_pheno %>% ggplot(aes(Phenotype)) + geom_histogram(bins = 150) +  theme_cowplot(12) + ggtitle("FSFPu5 Pheno")

FSFPu5_Geno_plot = FSFPu5_Genom %>%  ggplot(aes(Generation, Frequency, col = as.character(Effect)))+
  geom_point()+
  theme_cowplot(12)+
  theme(legend.position = "none")+
  ggtitle("FSFPu5 Genome")

FSFPu5_PhGen_plot
FSFPu5_Ph_plot
FSFPu5_Geno_plot


###### Poisson Pop (with mean = 10,000) ###### 


FSFPp5_Genom = read.csv("/Users/etb68/OneDrive - University of Missouri/Academics Mizzou (etb68@umsystem.edu)/Biologycal_Sciences/Research_Projects/Yeast_GitHub/evogen-sims/ET_Yeast/output.dir/DummyPop/WF_models/Pois_Pop/genome_FSFPp_WF.v4.csv")

FSFPp5_pheno = read.csv("/Users/etb68/OneDrive - University of Missouri/Academics Mizzou (etb68@umsystem.edu)/Biologycal_Sciences/Research_Projects/Yeast_GitHub/evogen-sims/ET_Yeast/output.dir/DummyPop/WF_models/Pois_Pop/MeanPhenotypes_FSFPp_WF.v4.csv")

FSFPp5_PhGen_plot = FSFPp5_pheno %>% ggplot(aes(Generation, Phenotype)) + geom_point() +  theme_cowplot(12)+ ggtitle("FSFPp5 Pheno vs Generation")
FSFPp5_Ph_plot = FSFPp5_pheno %>% ggplot(aes(Phenotype)) + geom_histogram(bins = 150) +  theme_cowplot(12) + ggtitle("FSFPp5 Pheno")

FSFPp5_Geno_plot = FSFPp5_Genom %>%  ggplot(aes(Generation, Frequency, col = as.character(Effect)))+
  geom_point()+
  theme_cowplot(12)+
  theme(legend.position = "none")+
  ggtitle("FSFPp5 Genome")

FSFPp5_PhGen_plot
FSFPp5_Ph_plot
FSFPp5_Geno_plot
