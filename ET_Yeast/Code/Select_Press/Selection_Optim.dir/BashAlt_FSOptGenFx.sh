#!/bin/bash
##################### Running Different QTL Numbers on Constant Selection and Population  #######################################
#SBATCH -p BioCompute,Lewis,hpc4 
#SBATCH -A kinglab
#SBATCH -J AltOptGenFx 
#SBATCH -c 4
#SBATCH -t 1-03:00:00
#SBATCH --mem 120G
#SBATCH -o ../../../output.dir/Select_Press/ConstFluctSelect.dir/Selection_OptimGen.dir/AltFSOptGenFx-%j.log
#SBATCH -e ../../../output.dir/Select_Press/ConstFluctSelect.dir/Selection_OptimGen.dir/AltFSOptGenFx-%j.err
##########SCIENCE FOLLOWS HERE ########################

#module load SLiM

echo -e "=== Beginning of SLiM run with different QTLs > $(date) ==="

output="/storage/hpc/group/kinglab/etb68/evogen-sims/ET_Yeast/output.dir/Select_Press/ConstFluctSelect.dir/Selection_Optim.dir/genome5_100_0.5.csv"
if [[ ! -f "$output" ]]
then
  echo "My file ${output} doesn't exist. Running SLiM QTLs now."
seeds=(2345 78344 11349 85732 65741)
replicates=(1 2 3 4 5)
heritabilities=(0.1 0.5 0.8)
stdvs=(1 2 3 4)

for seed in ${seeds[@]}
do
  for repl in ${replicates[@]}
  do
    for h in ${heritabilities[@]}
    do
      for SD in ${stdvs[@]}
      do
        slim -d seed=$seed -d repl=$repl -d loci=100 -d region=6101 -d h=$h -d SD=$SD FS_OptGenFx.slim
      done
    done
  done
done

else
  echo "All is well, Boss.  The ${output} file is there."
fi
echo "=== Fini finito! End of SLiM QTLs Constant Selection run >" $(date)

