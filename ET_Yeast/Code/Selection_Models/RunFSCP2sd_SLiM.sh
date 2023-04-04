#!/bin/bash
##################### Running Different QTL Numbers on Constant Selection and Populaton  #######################################
#SBATCH -p BioCompute 
#SBATCH -A kinglab
#SBATCH -J SLiM_qtl 
#SBATCH -c 4
#SBATCH --mem 120G
#SBATCH -o ../../output.dir/Selection_Models/FSCP2sd.dir/SLiM_QTL_FSCP2sd-%j.log
#SBATCH -e ../../output.dir/Selection_Models/FSCP2sd.dir/SLiM_QTL_FSCP2sd--%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=etb68@mail.missouri.edu
##########SCIENCE FOLLOWS HERE ########################

#module load SLiM

echo -e "=== Bigining of SLiM run with different QTLs > $(date) ==="

output="/storage/hpc/group/kinglab/etb68/evogen-sims/ET_Yeast/output.dir/Selection_Models/FSCP2sd.dir/genome2sd5_300_0.5.csv"
if [[ ! -f "$output" ]]
then
  echo "My file ${output} doesn't exist. Running SLiM QTLs now."

slim -d seed=2345 -d repl=1 -d loci=1 -d region=610140 -d h=0.5 FSCP2sd.QTL.slim
slim -d seed=2345 -d repl=1 -d loci=10 -d region=61014 -d h=0.5 FSCP2sd.QTL.slim
slim -d seed=2345 -d repl=1 -d loci=70 -d region=8716 -d h=0.5 FSCP2sd.QTL.slim
slim -d seed=2345 -d repl=1 -d loci=100 -d region=6101 -d h=0.5 FSCP2sd.QTL.slim
slim -d seed=2345 -d repl=1 -d loci=300 -d region=2033 -d h=0.5 FSCP2sd.QTL.slim

slim -d seed=78344 -d repl=2 -d loci=1 -d region=610140 -d h=0.5 FSCP2sd.QTL.slim
slim -d seed=78344 -d repl=2 -d loci=10 -d region=61014 -d h=0.5 FSCP2sd.QTL.slim
slim -d seed=78344 -d repl=2 -d loci=70 -d region=8716 -d h=0.5 FSCP2sd.QTL.slim
slim -d seed=78344 -d repl=2 -d loci=100 -d region=6101 -d h=0.5 FSCP2sd.QTL.slim 
slim -d seed=78344 -d repl=2 -d loci=300 -d region=2033 -d h=0.5 FSCP2sd.QTL.slim

slim -d seed=11349 -d repl=3 -d loci=1 -d region=610140 -d h=0.5 FSCP2sd.QTL.slim 
slim -d seed=11349 -d repl=3 -d loci=10 -d region=61014 -d h=0.5 FSCP2sd.QTL.slim
slim -d seed=11349 -d repl=3 -d loci=70 -d region=8716 -d h=0.5 FSCP2sd.QTL.slim
slim -d seed=11349 -d repl=3 -d loci=100 -d region=6101 -d h=0.5 FSCP2sd.QTL.slim 
slim -d seed=11349 -d repl=3 -d loci=300 -d region=2033 -d h=0.5 FSCP2sd.QTL.slim

slim -d seed=85732 -d repl=4 -d loci=1 -d region=610140 -d h=0.5 FSCP2sd.QTL.slim 
slim -d seed=85732 -d repl=4 -d loci=10 -d region=61014 -d h=0.5 FSCP2sd.QTL.slim
slim -d seed=85732 -d repl=4 -d loci=70 -d region=8716 -d h=0.5 FSCP2sd.QTL.slim
slim -d seed=85732 -d repl=4 -d loci=100 -d region=6101 -d h=0.5 FSCP2sd.QTL.slim 
slim -d seed=85732 -d repl=4 -d loci=300 -d region=2033 -d h=0.5 FSCP2sd.QTL.slim

slim -d seed=65741 -d repl=5 -d loci=1 -d region=610140 -d h=0.5 FSCP2sd.QTL.slim 
slim -d seed=65741 -d repl=5 -d loci=10 -d region=61014 -d h=0.5 FSCP2sd.QTL.slim
slim -d seed=65741 -d repl=5 -d loci=70 -d region=8716 -d h=0.5 FSCP2sd.QTL.slim
slim -d seed=65741 -d repl=5 -d loci=100 -d region=6101 -d h=0.5 FSCP2sd.QTL.slim 
slim -d seed=65741 -d repl=5 -d loci=300 -d region=2033 -d h=0.5 FSCP2sd.QTL.slim



else

  echo "All is well, Boss.  The ${output} file is there."
fi
echo "=== Fini finito! End of SLiM QTLs FSCP2sd run >" $(date)
