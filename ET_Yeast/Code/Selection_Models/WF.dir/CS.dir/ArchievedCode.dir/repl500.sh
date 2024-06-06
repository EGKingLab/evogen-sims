#!/bin/bash

# Set the seed for reproducibility
RANDOM=42

echo -e "=== Beginning of SLiM run with different QTLs > $(date) ===" >> ../../../../output.dir/Selection_Models/WF.dir/CS.dir/CS_par4.slim-$SLURM_JOBID.log

output="/storage/hpc/group/kinglab/etb68/evogen-sims/ET_Yeast/output.dir/Selection_Models/WF.dir/CS.dir/genome5_100_0.5.csv"
if [[ ! -f "$output" ]]; then
  echo "My file ${output} doesn't exist. Running SLiM QTLs now."

  heritabilities=(0.1 0.5 0.8)
  stdvs=(1 2 3 4)

  #creating an associative array where keys are the loci and the values are regions
  declare -A loci_to_regions
  loci_to_regions=( [1]=610140 [10]=61014 [70]=8716 [100]=6101 [300]=2033)

  # Run jobs in parallel
  for (( repl=1; repl<=500; repl++ )); do
    seed=$RANDOM
    for h in "${heritabilities[@]}"; do
      for SD in "${stdvs[@]}"; do
        for loci in "${!loci_to_regions[@]}"; do
          region=${loci_to_regions[$loci]}
          {
            start_time=$(date +%s)
            slim -d seed=$seed -d repl=$repl -d loci=$loci -d region=$region -d h=$h -d SD=$SD CS_Par4_Copy.slim
            end_time=$(date +%s)
            runtime=$((end_time - start_time))
            echo "=== Finished replicate $repl. Runtime: $runtime seconds ===" >> ../../../../output.dir/Selection_Models/WF.dir/CS.dir/CS_par4.slim-$SLURM_JOBID.log
          } &
        done
      done
    done
  done
  wait

else
  echo "All is well, Boss. The ${output} file is there."
fi
total_end_time=$(date +%s)
total_runtime=$((total_end_time - total_start_time))
echo "=== Finished! End of SLiM QTLs Constant Selection run > $(date). Total runtime: $total_runtime seconds ===" >> ../../../../output.dir/Selection_Models/WF.dir/CS.dir/CS_par4.slim-$SLURM_JOBID.log

