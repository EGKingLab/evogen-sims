
set.seed(60941)

rand_seeds <- round(runif(1000,1,1000000))

counter <- 1

gen_os <- c(1,5,10,20,50)

pop_type <- c("C","F")

niter <- 5

cat("", file="cmds_to_run.txt")

for(nn in 1:niter)
{
  
sel_type <- "F"

for(gg in gen_os)
{
  cmd <- paste0("Rscript Initial_model_QGPG.R ",rand_seeds[counter]," ",sel_type," ",gg," ","S",counter," >temp",counter,".txt 2>error",counter,".txt & \n")
  counter <- counter + 1
  
  cat(cmd, file="cmds_to_run.txt", append=TRUE)
}

sel_type <- "C"
gg <- NA

cmd <- paste0("Rscript Initial_model_QGPG.R ",rand_seeds[counter]," ",sel_type," ",gg," ","S",counter," >temp",counter,".txt 2>error",counter,".txt & \n")
counter <- counter + 1

cat(cmd, file="cmds_to_run.txt", append=TRUE)

}