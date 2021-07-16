library(tseries)
library(slider)

ngens <- 500
pop_type <- "F"
sel_type <- "C"

samp_id <- "S301"
#s4 PC SF
#s6 PC SC
#s301 PF SC
#300 PF SF

output <- readRDS(file = paste0("../output/Sim_Output_S",sel_type,"_P",pop_type,"_",samp_id,".rds"))

output.N <- readRDS(file="../output/Null_Sim_Output.rds")

ff1 <- output.N[seq(11,510,by=25),]

#ff1 <- output$freqs[seq(1,500, by =25),]


#acor <- acf(ff1,lag.max=499,pl=FALSE)
#aa <- matrix(NA,1000,500)

#for(ii in 1:1000)
#{
#aa[ii,] <- acor$acf[,ii,ii]
#}

#mm <- colMeans(aa, na.rm=TRUE)
#plot(mm)
#abline(v=c(25,50,75,100))

dd.mat.abs <- vector(mode='list', length=(nrow(ff1)-1))
dd.mat <- vector(mode='list', length=(nrow(ff1)-1))
for(jj in 1:(nrow(ff1)-1))
#dd.mat.abs <- vector(mode='list', length=100)
#dd.mat <- vector(mode='list', length=100)
#for(jj in 1:100)
  
  {

  dd <- diff(ff1, lag=jj)
  dd.mat.abs[[jj]] <- rowMeans(abs(dd))
  dd.mat[[jj]] <- rowMeans(dd)
  
  #dd_sd <- apply(dd, 2, sd)
  #cc <- cor(dd[,-which(dd_sd==0)])
  #cat(jj,"\t",mean(abs(cc[upper.tri(cc)])),"\n")
}

plot(dd.mat.abs[[1]])
plot(dd.mat.abs[[10]])
plot(dd.mat.abs[[25]])

plot(dd.mat[[1]])
plot(dd.mat[[10]])
plot(dd.mat[[25]])


#need likelihood of resulting stat

#slide.out <- matrix(NA,490,ncol(ff1))

#for(gg in 1:ncol(ff1))
#{
#  slide.out[,gg] <- unlist(slide(ff1[,gg],var,.before=10,.complete=TRUE))
#}

#ss <- rowMeans(slide.out)



dd_sd <- apply(dd, 2, sd)
cc <- cor(dd[,-which(dd_sd==0)])
mean(abs(cc[upper.tri(cc)]))

ff1_sd <- apply(ff1, 2, sd)
cc <- cor(ff1[,-which(ff1_sd==0)])
mean(abs(cc[upper.tri(cc)]))
#WGCNA, PCA
