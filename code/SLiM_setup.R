
set.seed(60941)

rand_seeds <- round(runif(100,1,100000))

niter <- 5

#slim -d seed=7 -d rep=1 QG_SLiM.slim >test.txt 2>error.txt


cat("", file="slim_cmds_to_run.txt")

for(ii in 1:niter)
{
  cat(paste0("slim -d seed=",rand_seeds[ii]," -d rep=", ii," QG_SLiM_CSCP.slim >testCSCP_",ii,".txt 2>errorCSCP_",ii,".txt & \n"), file="slim_cmds_to_run.txt", append=TRUE)
  cat(paste0("slim -d seed=",rand_seeds[ii]," -d rep=", ii," QG_SLiM_CSFP.slim >testCSFP_",ii,".txt 2>errorCSFP_",ii,".txt & \n"), file="slim_cmds_to_run.txt", append=TRUE)
  cat(paste0("slim -d seed=",rand_seeds[ii]," -d rep=", ii," QG_SLiM_FSCP.slim >testFSCP_",ii,".txt 2>errorFSCP_",ii,".txt & \n"), file="slim_cmds_to_run.txt", append=TRUE)
  cat(paste0("slim -d seed=",rand_seeds[ii]," -d rep=", ii," QG_SLiM_FSFP.slim >testFSFP_",ii,".txt 2>errorFSFP_",ii,".txt & \n"), file="slim_cmds_to_run.txt", append=TRUE)
  
}


