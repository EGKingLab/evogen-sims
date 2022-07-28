library(tidyverse)
library(cowplot)
theme_set(theme_cowplot())

ppall <- vector(mode="list",length=4)

types <- c("CSCP","CSFP","FSCP","FSFP")

for(ii in 1:4)
{

type <- types[ii]
dd <- read_csv(paste0("../output/SLiM/genome_track_",type,"_2.csv"))
#tester <- subset(dd, Generation == 500)
#cat(nrow(tester),"\n")
#ftest <- tester$Frequency
#ftest <- c(ftest, rep(0,(100-length(ftest))))
#cat(mean(2*ftest*(1-ftest)), "\n")
#}
dd$GeneID <- as.character(dd$Position)

#add mean
mm <- read_csv(paste0("../output/SLiM/MeanPhenos_",type,"_2.csv"))

mopt <- min(which(mm$Pheno >= (mm$Optim[1]-10)))


pp <- dd %>%
  ggplot(aes(Generation, Frequency, group=GeneID, color=Effect)) +
  geom_line(alpha=1/3) +
  scale_colour_gradientn(colors=c("white","black"))+
  theme(legend.position = "none") +
  ggtitle(type) +
  theme(plot.title = element_text(hjust = 0.5, size=12,  color="steelblue"))

if(type %in% c("CSCP", "CSFP"))
  {
  pp <- pp +   geom_vline(xintercept=mopt, color='#FDE725FF', size=1.5) 
}

if(type %in% c("CSFP", "FSFP"))
{
  pp <- pp +    theme(axis.title.y=element_blank()) 
  
}

if(type %in% c("CSCP", "CSFP"))
{
  pp <- pp +    theme(axis.title.x=element_blank()) 
  
}

ppall[[ii]] <- pp
}


allele.all <- plot_grid(plotlist = ppall, ncol=2, labels=c("a.","b.","c.","d."), 
                        rel_heights = c(1, 1.1,1,1.1), rel_widths = c(1.1,1,1.1,1),
                        label_x=c(0,-0.06,0,-0.06))

ggsave(allele.all, filename="../output/afreq_SLIM.pdf",width=6.5, height=5)


C_dat <- read_csv("../output/SLiM/MeanPhenos_CSCP_2.csv")
C_dat$Type <- "C"
F25_dat <- read_csv("../output/SLiM/MeanPhenos_FSCP_2.csv")
F25_dat$Type <- "F25"
F5_dat <- read_csv("../output/SLiM/MeanPhenos_FSCP_F5_2.csv")
F5_dat$Type <- "F5"

Long_phenos <- rbind(C_dat, F5_dat, F25_dat)

Long_phenos$Type <- factor(Long_phenos$Type, levels = c("C","F5","F25"))

opts <- c(F5_dat$Optim[1], F5_dat$Optim[6])

g1 <- ggplot(Long_phenos, aes(Generation,Pheno, color=Type)) +
  geom_line() + 
  geom_hline(yintercept = opts[1], lty = 3) +
  geom_hline(yintercept = opts[2], lty = 3) +
  ylab("Phenotype") +
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank())


g1

CSCP <- read_csv("../output/SLiM/MeanPhenos_CSCP_2.csv")
CSCP$Type <- "CSCP"
CSFP <- read_csv("../output/SLiM/MeanPhenos_CSFP_2.csv")
CSFP$Type <- "CSFP"
FSCP <- read_csv("../output/SLiM/MeanPhenos_FSCP_2.csv")
FSCP$Type <- "FSCP"
FSFP <- read_csv("../output/SLiM/MeanPhenos_FSFP_2.csv")
FSFP$Type <- "FSFP"

Long_phenos <- rbind(CSCP, FSCP, CSFP, FSFP)

g1p <- ggplot(Long_phenos, aes(Generation,Pheno, color=Type)) +
  geom_line() + 
  geom_hline(yintercept = opts[1], lty = 3) +
  geom_hline(yintercept = opts[2], lty = 3) +
  ylab("Phenotype") +
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank())


g1p

phen.all <- plot_grid(g1,g1p, labels=c("a.","b."),ncol=2)
ggsave(phen.all, filename="../output/phenotype_sim_SLiM.pdf", width=6.5, height=2.5)




