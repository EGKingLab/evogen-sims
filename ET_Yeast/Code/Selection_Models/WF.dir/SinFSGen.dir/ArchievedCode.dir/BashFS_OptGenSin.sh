#!/bin/bash
##################### Running Different QTL Numbers on Constant Selection and Population  #######################################
#SBATCH -p BioCompute 
#SBATCH -A deckerlab
#SBATCH -J FS_OptGenSin.slim 
#SBATCH -c 4
#SBATCH -t 2-00:00:00
#SBATCH --mem 120G
#SBATCH --array=0-4
#SBATCH -o ../../../../output.dir/Selection_Models/WF.dir/SinFSGen.dir/FS_OptGenSin-%A%a.log
#SBATCH -e ../../../../output.dir/Selection_Models/WF.dir/SinFSGen.dir/FS_OptGenSin-%A%a.err
##########SCIENCE FOLLOWS HERE ########################

#module load SLiM

echo -e "=== Beginning of SLiM run with different QTLs > $(date) ===">>../../../../output.dir/Selection_Models/WF.dir/SinFSGen.dir/FS_OptGenSin-%A%a.log
output="/storage/hpc/group/kinglab/etb68/evogen-sims/ET_Yeast/output.dir/Selection_Models/WF.dir/SinFSGen.dir/genome5_100_0.5.csv"
if [[ ! -f "$output" ]]
then
  echo "My file ${output} doesn't exist. Running SLiM QTLs now."

#seeds=(2345 78344 11349 85732 65741)
#replicates=(1 2 3 4 5)
heritabilities=(0.1 0.5 0.8)
stdvs=(1 2 3 4)

declare -A seeds_to_replicates
seeds_to_replicates=( [2345]=1 [78344]=2 [11349]=3 [85732]=4 [65741]=5 )

declare -A loci_to_regions
loci_to_regions=( [1]=610140 [10]=61014 [70]=8716 [100]=6101 [300]=2033  )


#creating an associative array where keys are the loci and the values are regions
index=${SLURM_ARRAY_TASK_ID}
seed_keys=(${!seeds_to_replicates[@]})
seed=${seed_keys[$index]}
repl=${seeds_to_replicates[$seed]}

#for seed in "${!seeds_to_replicates[@]}"
#do
for h in ${heritabilities[@]}
do
  for SD in ${stdvs[@]}
  do
    for loci in "${!loci_to_regions[@]}"
    do
        repl=${seeds_to_replicates[$seed]}
        region=${loci_to_regions[$loci]}
        slim -d seed=$seed -d repl=$repl -d loci=$loci -d region=$region -d h=$h -d SD=$SD FS_OptGenSin.slim
    done
  done
done
#done
wait
else
  echo "All is well, Boss.  The ${output} file is there."
fi
echo "=== Fini finito! End of SLiM QTLs Constant Selection run >" $(date)>>../../../../output.dir/Selection_Models/WF.dir/SinFSGen.dir/FS_OptGenSin-%A%a.log
