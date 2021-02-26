
#check that mimiceer package

# starting population - number of haplotypes
#randomly generate 18 effects?

# just chr lengths would be fine?
# just one chr would be fine?
# snp based instead?

# do some generations of crossing
# ignore genetic map for now? one recomb per gen?

# build model for trait
# genomic prediction/selection

# choose heritability and distribution/number of effects

# generate individuals & phenotypes

# assign fitness

# next generation

set.seed(1631)

tic <- Sys.time()

sel_type <- "F"
gen_os <- 20

#number of gens of selection
Nsel <- 120

Nhaps <- 18
haps <- paste0("H",seq(1,Nhaps))
Nmark <- 1000
Nchr <- 16
chroms <- paste0("C",seq(1,Nchr))

poslist <- data.frame("CHROM"=sample(chroms, Nmark,replace=T),"POS"=character(length=Nmark), stringsAsFactors = FALSE)
poslist <- poslist[order(poslist$CHROM),]
poslist$POS <- paste0("M",seq(1,Nmark))

Ngens.cross <- 5

recomb.rate <- 1

pop.size.i <- 1000

h2 <- 0.5

dist.to.trait <- 4 

i.mark.fr <- runif(Nmark, 0,1)

snp.c <- t(sapply(i.mark.fr, function(x) rbinom(Nhaps,1,x)))

#first mating done
#with haplotypes
#pop.i <- array(data=rep(sample(paste0("H",seq(1,Nhaps)),2*pop.size.i,replace=T),Nmark), dim=c(pop.size.i,2,Nmark))

#with snps
rand.haps <- sample(seq(1:Nchr),2*pop.size.i,replace=T)
pop.i <- array(data=as.numeric(snp.c[,rand.haps]), dim=c(Nmark,2,pop.size.i))
pop.i <- aperm(pop.i, c(3,2,1))
dimnames(pop.i)[[2]] <- c("A","B")
dimnames(pop.i)[[3]] <- paste0(poslist$CHROM,"_",poslist$POS)

#Nmat rounds of mating

pop.set <- pop.i
pop.size <- pop.size.i

#for loop here

Nmat <- 10

for(gg in 1:Nmat)
{

#how many offspring from each?
#make more than need - then select randomly
N.offspring <- rep(10,pop.size)

#make gamete pool

gamete.pool <- array(data=NA, dim=c(sum(N.offspring),Nmark))
dimnames(gamete.pool)[[2]] <- paste0(poslist$CHROM,"_",poslist$POS)

counter <- 1
for (ii in 1:length(N.offspring))
{
  if(N.offspring[ii] != 0)
    {
      for(jj in 1:N.offspring[ii])
        {
          r.ind <- sample(seq(2,(Nmark -1)),1)
          gamete.pool[counter,] <- c(pop.set[ii,'A',1:(r.ind -1)],pop.set[ii,'B',r.ind:Nmark])
          counter <- counter + 1
        }
    }
}

#need to extend to multiple chromosomes - now just one
#track types & sexes?

#randomly pair to make new individuals

#check it is an even number - 
#pop.size <- ifelse((sum(N.offspring) %% 2) == 0, sum(N.offspring), (sum(N.offspring) - 1))/2

#shuffle gamete pool
gamete.pool <- gamete.pool[sample(1:nrow(gamete.pool)),]

pop.set <- array(t(gamete.pool[1:(pop.size*2),]),dim=c(Nmark,2,pop.size))
pop.set <- aperm(pop.set, c(3,2,1))
dimnames(pop.set)[[2]] <- c("A","B")
dimnames(pop.set)[[3]] <- paste0(poslist$CHROM,"_",poslist$POS)

cat("F",gg,"\n")

}




#selection

mark.effs <- rexp(Nmark)
BVs <- colSums(mark.effs*(t(pop.set[,"A",]) + t(pop.set[,"B",])))
gen.var <- var(BVs)
env.var <- (gen.var - h2*gen.var)/h2
env.eff <- rnorm(pop.size,0,sqrt(env.var))
phenos <- BVs + env.eff
init.mean <- mean(phenos)
init.sd <- sd(phenos)

if(sel_type == "C")
{
  new.opt.all <- rep(init.mean + dist.to.trait*init.sd,
                     times = Nsel)
}else{
  if(sel_type == "F")
    {
   
    #make vector of new.opt
    new.opt.all <- rep(c(init.mean - dist.to.trait*init.sd, init.mean + dist.to.trait*init.sd),
                       each = gen_os, 
                       times = ((Nsel/2)/gen_os))
    
  }else{
    cat("invalid sel_type")
  }
}


dist.p <- abs(phenos - new.opt.all[1])
seld <- (dist.p/init.sd)/10
rel.w <- 1-seld
rel.w[rel.w < 0] <- 0
N.offspring <- round(rel.w*10)

output <- data.frame("gen" = seq(1,Nsel),"freq" = numeric(length=Nsel),"pheno"= numeric(length=Nsel))
freqs.pos <- colMeans(pop.set, dims=2)
focal.pos <- sample(which(freqs.pos < 0.6 & freqs.pos > 0.4 & mark.effs >2),1)



for(ss in 1:Nsel)
{
  
  #make gamete pool
  
  gamete.pool <- array(data=NA, dim=c(sum(N.offspring),Nmark))
  dimnames(gamete.pool)[[2]] <- paste0(poslist$CHROM,"_",poslist$POS)
  
  counter <- 1
  for (ii in 1:length(N.offspring))
  {
    if(N.offspring[ii] != 0)
    {
      for(jj in 1:N.offspring[ii])
      {
        r.ind <- sample(seq(2,(Nmark -1)),1)
        gamete.pool[counter,] <- c(pop.set[ii,'A',1:(r.ind -1)],pop.set[ii,'B',r.ind:Nmark])
        counter <- counter + 1
      }
    }
  }
  
  #shuffle gamete pool
  gamete.pool <- gamete.pool[sample(1:nrow(gamete.pool)),]
  
  pop.set <- array(t(gamete.pool[1:(pop.size*2),]),dim=c(Nmark,2,pop.size))
  pop.set <- aperm(pop.set, c(3,2,1))
  dimnames(pop.set)[[2]] <- c("A","B")
  dimnames(pop.set)[[3]] <- paste0(poslist$CHROM,"_",poslist$POS)
  
  BVs <- colSums(mark.effs*(t(pop.set[,"A",]) + t(pop.set[,"B",])))
  env.eff <- rnorm(pop.size,0,sqrt(env.var))
  phenos <- BVs + env.eff
  
  new.opt <- new.opt.all[(ss+1)]
  
  dist.p <- abs(phenos - new.opt)
  seld <- (dist.p/init.sd)/10
  rel.w <- 1-seld
  N.offspring <- round(rel.w*10)
  
  output[ss,"pheno"] <- mean(phenos)
  output[ss,"freq"] <- mean(as.numeric(pop.set[,,focal.pos]))
  
  cat("S",ss,"\t",mean(phenos),"\t",var(BVs),"\n")
}

toc<- Sys.time()

toc-tic

out.dat <- write.csv(output,file = "../OutputData/small_test.csv")

#what to keep?!

library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())

ggplot(output, aes(gen,freq)) +
  geom_point()

