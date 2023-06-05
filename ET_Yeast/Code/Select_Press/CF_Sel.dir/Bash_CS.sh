#!/bin/bash
##################### Running Different QTL Numbers on Constant Selection and Populaton  #######################################
#SBATCH -p BioCompute,Lewis,hpc4 
#SBATCH -A kinglab
#SBATCH -J SLiM_CS 
#SBATCH -c 4
#SBATCH -t 1-00:00:00
#SBATCH --mem 120G
#SBATCH -o ../../../output.dir/Select_Press/ConstFluctSelect.dir/CS.dir/SLiM_QTL_CS-%j.log
#SBATCH -e ../../../output.dir/Select_Press/ConstFluctSelect.dir/CS.dir/SLiM_QTL_CS--%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=etb68@mail.missouri.edu
##########SCIENCE FOLLOWS HERE ########################

#module load SLiM

echo -e "=== Bigining of SLiM run with different QTLs > $(date) ==="

output="/storage/hpc/group/kinglab/etb68/evogen-sims/ET_Yeast/output.dir/Select_Press/ConstFluctSelect.dir/CS.dir/genome5_100_0.5.csv"
if [[ ! -f "$output" ]]
then
  echo "My file ${output} doesn't exist. Running SLiM QTLs now."

slim -d seed=2345 -d repl=1 -d loci=100 -d region=6101 -d h=0.5 -d SD=1 -d rang=101 -d gen=10 CS.slim 
slim -d seed=78344 -d repl=2 -d loci=100 -d region=6101 -d h=0.5 -d SD=1 -d rang=101 -d gen=10 CS.slim 
slim -d seed=11349 -d repl=3 -d loci=100 -d region=6101 -d h=0.5 -d SD=1 -d rang=101 -d gen=10 CS.slim 
slim -d seed=85732 -d repl=4 -d loci=100 -d region=6101 -d h=0.5 -d SD=1 -d rang=101 -d gen=10 CS.slim 
slim -d seed=65741 -d repl=5 -d loci=100 -d region=6101 -d h=0.5 -d SD=1 -d rang=101 -d gen=10 CS.slim

slim -d seed=2345 -d repl=1 -d loci=100 -d region=6101 -d h=0.5 -d SD=1 -d rang=51 -d gen=20 CS.slim
slim -d seed=78344 -d repl=2 -d loci=100 -d region=6101 -d h=0.5 -d SD=1 -d rang=51 -d gen=20 CS.slim
slim -d seed=11349 -d repl=3 -d loci=100 -d region=6101 -d h=0.5 -d SD=1 -d rang=51 -d gen=20 CS.slim
slim -d seed=85732 -d repl=4 -d loci=100 -d region=6101 -d h=0.5 -d SD=1 -d rang=51 -d gen=20 CS.slim
slim -d seed=65741 -d repl=5 -d loci=100 -d region=6101 -d h=0.5 -d SD=1 -d rang=51 -d gen=20 CS.slim

slim -d seed=2345 -d repl=1 -d loci=100 -d region=6101 -d h=0.5 -d SD=1 -d rang=34 -d gen=30 CS.slim
slim -d seed=78344 -d repl=2 -d loci=100 -d region=6101 -d h=0.5 -d SD=1 -d rang=34 -d gen=30 CS.slim
slim -d seed=11349 -d repl=3 -d loci=100 -d region=6101 -d h=0.5 -d SD=1 -d rang=34 -d gen=30 CS.slim
slim -d seed=85732 -d repl=4 -d loci=100 -d region=6101 -d h=0.5 -d SD=1 -d rang=34 -d gen=30 CS.slim
slim -d seed=65741 -d repl=5 -d loci=100 -d region=6101 -d h=0.5 -d SD=1 -d rang=34 -d gen=30 CS.slim



slim -d seed=2345 -d repl=1 -d loci=100 -d region=6101 -d h=0.5 -d SD=2 -d rang=101 -d gen=10 CS.slim
slim -d seed=78344 -d repl=2 -d loci=100 -d region=6101 -d h=0.5 -d SD=2 -d rang=101 -d gen=10 CS.slim
slim -d seed=11349 -d repl=3 -d loci=100 -d region=6101 -d h=0.5 -d SD=2 -d rang=101 -d gen=10 CS.slim
slim -d seed=85732 -d repl=4 -d loci=100 -d region=6101 -d h=0.5 -d SD=2 -d rang=101 -d gen=10 CS.slim
slim -d seed=65741 -d repl=5 -d loci=100 -d region=6101 -d h=0.5 -d SD=2 -d rang=101 -d gen=10 CS.slim

slim -d seed=2345 -d repl=1 -d loci=100 -d region=6101 -d h=0.5 -d SD=2 -d rang=51 -d gen=20 CS.slim
slim -d seed=78344 -d repl=2 -d loci=100 -d region=6101 -d h=0.5 -d SD=2 -d rang=51 -d gen=20 CS.slim
slim -d seed=11349 -d repl=3 -d loci=100 -d region=6101 -d h=0.5 -d SD=2 -d rang=51 -d gen=20 CS.slim
slim -d seed=85732 -d repl=4 -d loci=100 -d region=6101 -d h=0.5 -d SD=2 -d rang=51 -d gen=20 CS.slim
slim -d seed=65741 -d repl=5 -d loci=100 -d region=6101 -d h=0.5 -d SD=2 -d rang=51 -d gen=20 CS.slim

slim -d seed=2345 -d repl=1 -d loci=100 -d region=6101 -d h=0.5 -d SD=2 -d rang=34 -d gen=30 CS.slim
slim -d seed=78344 -d repl=2 -d loci=100 -d region=6101 -d h=0.5 -d SD=2 -d rang=34 -d gen=30 CS.slim
slim -d seed=11349 -d repl=3 -d loci=100 -d region=6101 -d h=0.5 -d SD=2 -d rang=34 -d gen=30 CS.slim
slim -d seed=85732 -d repl=4 -d loci=100 -d region=6101 -d h=0.5 -d SD=2 -d rang=34 -d gen=30 CS.slim
slim -d seed=65741 -d repl=5 -d loci=100 -d region=6101 -d h=0.5 -d SD=2 -d rang=34 -d gen=30 CS.slim

 

slim -d seed=2345 -d repl=1 -d loci=100 -d region=6101 -d h=0.5 -d SD=3 -d rang=101 -d gen=10 CS.slim
slim -d seed=78344 -d repl=2 -d loci=100 -d region=6101 -d h=0.5 -d SD=3 -d rang=101 -d gen=10 CS.slim
slim -d seed=11349 -d repl=3 -d loci=100 -d region=6101 -d h=0.5 -d SD=3 -d rang=101 -d gen=10 CS.slim
slim -d seed=85732 -d repl=4 -d loci=100 -d region=6101 -d h=0.5 -d SD=3 -d rang=101 -d gen=10 CS.slim
slim -d seed=65741 -d repl=5 -d loci=100 -d region=6101 -d h=0.5 -d SD=3 -d rang=101 -d gen=10 CS.slim

slim -d seed=2345 -d repl=1 -d loci=100 -d region=6101 -d h=0.5 -d SD=3 -d rang=51 -d gen=20 CS.slim
slim -d seed=78344 -d repl=2 -d loci=100 -d region=6101 -d h=0.5 -d SD=3 -d rang=51 -d gen=20 CS.slim
slim -d seed=11349 -d repl=3 -d loci=100 -d region=6101 -d h=0.5 -d SD=3 -d rang=51 -d gen=20 CS.slim
slim -d seed=85732 -d repl=4 -d loci=100 -d region=6101 -d h=0.5 -d SD=3 -d rang=51 -d gen=20 CS.slim
slim -d seed=65741 -d repl=5 -d loci=100 -d region=6101 -d h=0.5 -d SD=3 -d rang=51 -d gen=20 CS.slim

slim -d seed=2345 -d repl=1 -d loci=100 -d region=6101 -d h=0.5 -d SD=3 -d rang=34 -d gen=30 CS.slim
slim -d seed=78344 -d repl=2 -d loci=100 -d region=6101 -d h=0.5 -d SD=3 -d rang=34 -d gen=30 CS.slim
slim -d seed=11349 -d repl=3 -d loci=100 -d region=6101 -d h=0.5 -d SD=3 -d rang=34 -d gen=30 CS.slim
slim -d seed=85732 -d repl=4 -d loci=100 -d region=6101 -d h=0.5 -d SD=3 -d rang=34 -d gen=30 CS.slim
slim -d seed=65741 -d repl=5 -d loci=100 -d region=6101 -d h=0.5 -d SD=3 -d rang=34 -d gen=30 CS.slim



slim -d seed=2345 -d repl=1 -d loci=100 -d region=6101 -d h=0.5 -d SD=4 -d rang=101 -d gen=10 CS.slim      
slim -d seed=78344 -d repl=2 -d loci=100 -d region=6101 -d h=0.5 -d SD=4 -d rang=101 -d gen=10 CS.slim
slim -d seed=11349 -d repl=3 -d loci=100 -d region=6101 -d h=0.5 -d SD=4 -d rang=101 -d gen=10 CS.slim
slim -d seed=85732 -d repl=4 -d loci=100 -d region=6101 -d h=0.5 -d SD=4 -d rang=101 -d gen=10 CS.slim      
slim -d seed=65741 -d repl=5 -d loci=100 -d region=6101 -d h=0.5 -d SD=4 -d rang=101 -d gen=10 CS.slim

slim -d seed=2345 -d repl=1 -d loci=100 -d region=6101 -d h=0.5 -d SD=4 -d rang=51 -d gen=20 CS.slim
slim -d seed=78344 -d repl=2 -d loci=100 -d region=6101 -d h=0.5 -d SD=4 -d rang=51 -d gen=20 CS.slim
slim -d seed=11349 -d repl=3 -d loci=100 -d region=6101 -d h=0.5 -d SD=4 -d rang=51 -d gen=20 CS.slim
slim -d seed=85732 -d repl=4 -d loci=100 -d region=6101 -d h=0.5 -d SD=4 -d rang=51 -d gen=20 CS.slim
slim -d seed=65741 -d repl=5 -d loci=100 -d region=6101 -d h=0.5 -d SD=4 -d rang=51 -d gen=20 CS.slim

slim -d seed=2345 -d repl=1 -d loci=100 -d region=6101 -d h=0.5 -d SD=4 -d rang=34 -d gen=30 CS.slim
slim -d seed=78344 -d repl=2 -d loci=100 -d region=6101 -d h=0.5 -d SD=4 -d rang=34 -d gen=30 CS.slim
slim -d seed=11349 -d repl=3 -d loci=100 -d region=6101 -d h=0.5 -d SD=4 -d rang=34 -d gen=30 CS.slim
slim -d seed=85732 -d repl=4 -d loci=100 -d region=6101 -d h=0.5 -d SD=4 -d rang=34 -d gen=30 CS.slim
slim -d seed=65741 -d repl=5 -d loci=100 -d region=6101 -d h=0.5 -d SD=4 -d rang=34 -d gen=30 CS.slim

else

  echo "All is well, Boss.  The ${output} file is there."
fi
echo "=== Fini finito! End of SLiM QTLs Constant Selection run >" $(date)
