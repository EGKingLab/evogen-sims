library(tidyverse)
library(cowplot)
theme_set(theme_cowplot())

set.seed(7711)


#how many loci
nloci <- 100

#how many individuals
ninds <- 10000

#heritability
h2 <- 0.5

#get effects
effects <- rexp(nloci)

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
phenotypes$phenotypes <- (phenotypes$phenotypes/ sd(phenotypes$phenotypes)) + 100
#phenotypes$phenotypes <- ((phenotypes$phenotypes - mean(phenotypes$phenotypes))/sd(phenotypes$phenotypes)) + 10

#define new optimum
new_opt <- mean(phenotypes$phenotypes) + sd(phenotypes$phenotypes)*2

phenotypes |>
  ggplot(aes(x=phenotypes)) +
  geom_histogram(bins = 50) +
  geom_vline(xintercept = new_opt, color = "steelblue")


scaling_f <- 100
  
phenotypes$fitness <- 1 - (abs(phenotypes$phenotypes - new_opt)/scaling_f)

#phenotypes$fitness <- (dnorm(new_opt - phenotypes$phenotypes, 10))

phenotypes |>
  ggplot(aes(x=fitness)) +
  geom_histogram(bins = 50)


phenotypes |>
  ggplot(aes(phenotypes, fitness)) +
  geom_point() +
  #ylim(0,1) +
  geom_vline(xintercept = new_opt, color = "steelblue")

