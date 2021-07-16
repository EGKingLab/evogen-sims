library(tidyverse)
library(cowplot)
theme_set(theme_cowplot())

MM <- read.table("../output/SNP_table_Yeast.txt",header=TRUE)
#So the “best” replicates are: 3, 7, 8, 9, 10, 11 and 12

M_rep <- MM %>% 
  select(starts_with("ancmaf")|starts_with("maf_12_"))

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
