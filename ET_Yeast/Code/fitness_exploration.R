library(tidyverse)
library(cowplot)
theme_set(theme_cowplot())

rm(list = ls())
set.seed(7711)

#how many loci
nloci <- 100

#how many individuals
ninds <- 10000

#heritability
h2 <- 0.5

#get effects
effects <- rexp(nloci)
#effects <- rgamma(nloci, 8, 7)
#effects <- rbeta(nloci, 1, 1)
#get freqs
markers <- runif(nloci,0.01,0.99)

additive <- numeric(length = ninds)

for(ii in 1:ninds){
  #get genotype for 10 loci for one individual
  geno <- rbinom(rep(1,10),rep(1,10), markers)
  effs <- geno*effects
  additive[ii] <- sum(effs)
}

V_A <- var(additive)
V_E <- (V_A/h2)- V_A
env <- rnorm(ninds, 0.0, sqrt(V_E))

phenotypes <- tibble(phenotypes = additive + env)
#phenotypes$phenotypes <- (phenotypes$phenotypes/ sd(phenotypes$phenotypes))  + 100
#phenotypes$phenotypes <- ((phenotypes$phenotypes - mean(phenotypes$phenotypes))/sd(phenotypes$phenotypes)) + 100
#phenotypes$phenotypes <- (phenotypes$phenotypes/ max(phenotypes$phenotypes))  + 100
#phenotypes$phenotypes <- quantile(phenotypes$phenotypes, rank(phenotypes$phenotypes)/ length(phenotypes$phenotypes))
 #phenotypes$phenotypes <- (phenotypes$phenotypes - mean(phenotypes$phenotypes))/ sd(phenotypes$phenotypes) 
 phenotypes$phenotypes <- ((phenotypes$phenotypes - min(phenotypes$phenotypes))/ (max(phenotypes$phenotypes) - min(phenotypes$phenotypes)))

#define new optimum
new_opt <- mean(phenotypes$phenotypes) + sd(phenotypes$phenotypes)*2

phenotypes |>
  ggplot(aes(x=phenotypes)) +
  geom_histogram(bins = 50) +
  geom_vline(xintercept = new_opt, color = "steelblue")


scaling_f <- mean(phenotypes$phenotypes)
  
#phenotypes$fitness <- 1 - ((phenotypes$phenotypes - new_opt)^2/phenotypes$phenotypes) #Quadratic
phenotypes$fitness <- exp(-(1/10) * (phenotypes$phenotypes - new_opt)^2)# we can alsp try 1/
#phenotypes$fitness <- exp(-(1/(phenotypes$phenotypes)) * ((phenotypes$phenotypes - new_opt))^2)

#phenotypes$fitness <- exp(-(3/(scaling_f)) * (phenotypes$phenotypes - new_opt)^2)
#phenotypes$fitness <- (dnorm(new_opt - phenotypes$phenotypes, 10))

phenotypes |>
  ggplot(aes(x=fitness)) +
  geom_histogram(bins = 50)


phenotypes |>
  ggplot(aes(phenotypes, fitness)) +
  geom_point() +
  #ylim(0,1) +
  geom_vline(xintercept = new_opt, color = "steelblue")

set.seed(1234)
hist(rexp(50))

hist(rgamma(50, 2, 2))

hist(0.1/rbeta(50, 1, 1))
