#!/bin/bash
#module load SLiM

echo -e "=== Bigining of SLiM run with different QTLs > $(date) ===" >> ../../../../output.dir/Selection_Models/WF.dir/NS.dir/NS_par4-${SLURM_JOBID}_${SLURM_ARRAY_TASK_ID}.log

output="/storage/hpc/group/kinglab/etb68/evogen-sims/ET_Yeast/output.dir/Selection_Models/WF.dir/NS.dir/genome5_100_0.5.csv"
if [[ ! -f "$output" ]]
then
  echo "My file ${output} doesn't exist. Running SLiM QTLs now."

heritabilities=(0.1 0.5 0.8)
stdvs=(1 2 3 4)

declare -A seeds_to_replicates
seeds_to_replicates=( [2345]=1 [78344]=2 [11349]=3 [85732]=4 [65741]=5 [49831]=6 [49826]=7 [49914]=8 [49969]=9 [49719]=10 [49849]=11 [50022]=12 [50172]=13 [50346]=14 [49970]=15 [50007]=16 [50103]=17 [49991]=18 [49876]=19 [50084]=20 [49993]=21 [50123]=22 [50079]=23 [49801]=24 [49909]=25 [50064]=26 [49950]=27 [50207]=28 [50164]=29 [50196]=30)

declare -A loci_to_regions
loci_to_regions=( [1]=610140 [10]=61014 [70]=8716 [100]=6101 [300]=2033  )

declare -A generations_to_ranges
generations_to_ranges=( [10]=101 [20]=51 [30]=34 )

for seed in "${!seeds_to_replicates[@]}"; do
  repl=${seeds_to_replicates[$seed]}
  for h in "${heritabilities[@]}"; do
    for loci in "${!loci_to_regions[@]}"; do
       region=${loci_to_regions[$loci]}
           {
             slim -d seed=$seed -d repl=$repl -d loci=$loci -d region=$region -d h=$h NS_Par4_Copy.slim
           } &
    done
  done
done

wait

else
  echo "All is well, Boss.  The ${output} file is there."
fi
echo "=== Fini finito! End of SLiM QTLs Constant Selection run >" $(date) >> ../../../../output.dir/Selection_Models/WF.dir/NS.dir/NS_par4--$SLURM_JOBID_${SLURM_ARRAY_TASK_ID}.log

