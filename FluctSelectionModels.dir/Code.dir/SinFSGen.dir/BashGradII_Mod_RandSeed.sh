#!/bin/bash
max_jobs=60

echo -e "=== Beginning of SLiM run with different QTLs > $(date) ===">>../../Output.dir/SinFSGen.dir/FS_OptGenSin${SLURM_JOBID}_${SLURM_ARRAY_TASK_ID}.log
output="/storage/hpc/group/kinglab/etb68/evogen-sims/ET_Yeast/Output.dir/SinFSGen.dir/genome5_100_0.5.csv"
if [[ ! -f "$output" ]]
then
  echo "My file ${output} doesn't exist. Running SLiM QTLs now."

heritabilities=(0.1 0.5 0.8)
stdvs=(1 2 3 4)

declare -A loci_to_regions
loci_to_regions=( [1]=610140 [10]=61014 [70]=8716 [100]=6101 [300]=2033 )

  # Generate 30 unique seeds using awk for a uniform distribution (range: 4000 to 7000)
  #seed="$fixed_seed" That is, I can replace this line with seed variable in the next line
  seeds=($(awk -v seed=12345 'BEGIN { srand(seed); for (i=1; i<=30; i++) printf "%d ", int(4000 + rand() * (7000 - 4000)); }'))

  # Mapping seeds to replicates
  declare -A seeds_to_replicates
  for i in {0..29}; do
    seeds_to_replicates[${seeds[$i]}]=$((i+1))
  done

for seed in "${!seeds_to_replicates[@]}"
do
  repl=${seeds_to_replicates[$seed]}
  for h in ${heritabilities[@]}
  do
    for SD in ${stdvs[@]}
    do
      for loci in "${!loci_to_regions[@]}"
      do
        repl=${seeds_to_replicates[$seed]}
        region=${loci_to_regions[$loci]}
          
          # Wait for free slot
          while (( $(jobs | wc -l) >= max_jobs )); do sleep 1; done
          
	{
	 slim -d seed=$seed -d repl=$repl -d loci=$loci -d region=$region -d h=$h -d SD=$SD GradII_Mod.slim
        } & 
      done
    done
  done
done

wait

else
  echo "All is well, Boss.  The ${output} file is there."
fi
total_end_time=$(date +%s)
total_runtime=$((total_end_time - total_start_time))
echo "=== Finished! End of SLiM QTLs fluctuating 4 sinusoidal run > $(date). Total runtime: $total_runtime seconds ===" >>../../Output.dir/SinFSGen.dir/FS_OptGenSin${SLURM_JOBID}_${SLURM_ARRAY_TASK_ID}.log
