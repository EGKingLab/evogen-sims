
seed.s <- 376294

set.seed(seed.s)

#number of gens of selection
Nsel <- 500

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

pop.size.i <- 10000

h2 <- 0.5

#dist.to.trait <- 4 

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

Nmat <- 10 + Nsel
output.freqs <- matrix(NA,Nmat,Nmark)

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

output.freqs[gg,] <- colMeans(pop.set, dims=2)

cat("F",gg,"\n")

}

saveRDS(output.freqs, file = paste0("../output/Null_Sim_Output.rds"))

