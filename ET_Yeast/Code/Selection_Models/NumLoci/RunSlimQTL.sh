#!/bin/bash
##################### Running Different QTL Numbers on Constant Selection and Populaton  #######################################
#SBATCH -p BioCompute 
#SBATCH -A kinglab
#SBATCH -J SLiM_qtl 
#SBATCH -c 4
#SBATCH --mem 120G
#SBATCH -o ../../output.dir/Selection_Models/SLiM_QTL_CSCP-%j.log
#SBATCH -e ../../output.dir/Selection_Models/SLiM_QTL_CSCP--%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=etb68@mail.missouri.edu
##########SCIENCE FOLLOWS HERE ########################

#module load SLiM

echo -e "=== Bigining of SLiM run with different QTLs > $(date) ==="

output="/storage/hpc/group/kinglab/etb68/evogen-sims/ET_Yeast/output.dir/Selection_Models/genome5_300_0.5.csv"
if [[ ! -f "$output" ]]
then
  echo "My file ${output} doesn't exist. Running SLiM QTLs now."

slim -d seed=2345 -d repl=1 -d loci=1 -d region=6104140 -d h=0.5 CSCP.QTL.slim  > ../../output.dir/Selection_Models/CSCP1_1.txt 1>errorCSCP1_1.txt & 
slim -d seed=2345 -d repl=1 -d loci=10 -d region=110934 -d h=0.5 CSCP.QTL.slim  > ../../output.dir/Selection_Models/CSCP1_1O.txt 1>errorCSCP1_10.txt &
slim -d seed=2345 -d repl=1 -d loci=70 -d region=17186 -d h=0.5 CSCP.QTL.slim  > ../../output.dir/Selection_Models/CSCP1_70.txt 1>errorCSCP1_70.txt & 
slim -d seed=2345 -d repl=1 -d loci=100 -d region=12081 -d h=0.5 CSCP.QTL.slim  > ../../output.dir/Selection_Models/CSCP1_100.txt 1>error_CSCP1_100.txt & 
slim -d seed=2345 -d repl=1 -d loci=300 -d region=4053 -d h=0.5 CSCP.QTL.slim  > ../../output.dir/Selection_Models/CSCP1_300.txt 1>error_CSCP1_300.txt & 

else

  echo "All is well, Boss.  The ${output} file is there."
fi
echo "=== Fini finito! End of SLiM QTLs CSCP run >" $(date)

