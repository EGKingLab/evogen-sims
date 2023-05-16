########################
## This is for test ###
library(tidyverse)
library(cowplot)
library(ggforce)
theme_set(theme_cowplot())
######################

############## GENOTYPE ##############

setwd("/home/etb68/YeastProj.dir/evogen-sims/ET_Yeast/Code/Plotting/")
genome_file1 <- read.csv("../../output.dir/Selection_Models/FSCP2sd.dir/genome2sd1_100_0.5.csv")

genome_file1 %>% ggplot(aes(x = Generation, y = Frequency, group = Position))+
  geom_line()+
  facet_zoom(xlim = c(400,600),ylim = c(0, 1), 
             horizontal = T, zoom.size = 1)+
  theme(legend.position = "none", 
        axis.title = element_text(face = "bold"))

genome_file1 %>% ggplot(aes(x = Generation, y = Frequency, color = factor(Position)))+
  theme(legend.position = "none",
        axis.title = element_text(face = "bold"))+
  geom_line()

genome_file2 <- read.csv("../../output.dir/Selection_Models/CSFP.1.dir/genomeCSFP_A1_100_0.5.csv")

genome_file2 %>% ggplot(aes(x = Generation, y = Frequency, color = factor(Position)))+
  geom_line()+
  facet_zoom(xlim = c(400,600),ylim = c(0, 1), 
             horizontal = T, zoom.size = 1)+
  theme(legend.position = "none", 
        axis.title = element_text(face = "bold"))

genome_file2 %>% ggplot(aes(x = Generation, y = Frequency, color = factor(Position)))+
  theme(legend.position = "none",
        axis.title = element_text(face = "bold"))+
  geom_line()

############## PHENOTYPE ##############

pheno1 <- read.csv("../../output.dir/Selection_Models/FSCP2sd.dir/MeanPhenotypes2sd1_100_0.5.csv")
pheno2 <- read.csv("../../output.dir/Selection_Models/CSFP.1.dir/MeanPhenotypesCSFP_A1_100_0.5.csv")

ggplot()+
  geom_line(data = pheno1, aes(x = Generation, y = Phenotype))+
  theme(axis.title = element_text(face = "bold"))+
  facet_zoom(xlim = c(700,1025),ylim = c(109, 112), 
             horizontal = F, zoom.size = 1)+
  geom_rect(data = myrect,
            aes(xmin = x, xmax = y,
                ymin = 109, ymax = Inf), 
            fill = "lightgrey", alpha = 0.2)

pheno1 %>% ggplot(aes(Phenotype))+
  geom_histogram(bins = 50)+
  theme(axis.title = element_text(face = "bold"))

###########################
p1.1 <- pheno1 %>% ggplot(aes(Phenotype))+
  geom_histogram(bins = 30)+
  theme(axis.title = element_text(face = "bold"))

p1.2 <- pheno1 %>% ggplot(aes(x = Generation, y = Phenotype))+
  geom_line()+
  theme(axis.title = element_text(face = "bold"))

plot_grid(p1.2, p1.1, ncol = 1)
##########
x = seq(0, 2000, by = 50)
y = seq(25, 2025, 50) # x + 25
myrect <- data.frame(x,y)

########## More Visual ##########
ggplot()+
  geom_line(data = pheno2, aes(x = Generation, y = Phenotype))+
  theme(axis.title = element_text(face = "bold"))+
  facet_zoom(xlim = c(425,800),ylim = c(109, 114), 
             horizontal = F, zoom.size = 1)+
  geom_rect(data = myrect,
            aes(xmin = x, xmax = y,
                ymin = 109, ymax = Inf), 
            fill = "lightgrey", alpha = 0.2)

pheno2 %>% ggplot(aes(Phenotype))+
  geom_histogram(bins = 50)+
  theme(axis.title = element_text(face = "bold"))

########################## My 2nd Plots ######

p2.1 <- pheno2 %>% ggplot(aes(Phenotype))+
  geom_histogram(bins = 30)+
  theme(axis.title = element_text(face = "bold"))

p2.2 <- pheno2 %>% ggplot(aes(x = Generation, y = Phenotype))+
  geom_line()+
  theme(axis.title = element_text(face = "bold"))

plot_grid(p2.2, p2.1, ncol = 1)

################################ TRANSFORMATION #############

file <- read.csv("../../output.dir/Selection_Models/FSCP2sd.dir/genome2sd1_100_0.5.csv")

mtw <- dwt(ile$Frequency, boundary="periodic", filter = "la14")
mtw_coef <- unlist(lapply(mtw@W, function(x) rep(x, each = 2^(mtw@level - 1))))
length(mtw_coef) # should be equal to nrow(ile)

trans_gen <- ile %>% 
  mutate(Trans_Freq = mtw_coef) #Mod(fft(Frequency)))

trans_gen %>% ggplot(aes(x = Generation, y = Trans_Freq, group = Position))+
  geom_line()+
  theme(legend.position = "none", 
        axis.title = element_text(face = "bold"))


        