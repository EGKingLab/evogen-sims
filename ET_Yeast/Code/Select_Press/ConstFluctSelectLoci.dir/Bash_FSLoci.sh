#!/bin/bash
##################### Running Different QTL Numbers on Constant Selection and Populaton  #######################################
#SBATCH -p BioCompute,Lewis,hpc4 
#SBATCH -A kinglab
#SBATCH -J SLiM_FS 
#SBATCH -c 4
#SBATCH -t 1-03:00:00
#SBATCH --mem 120G
#SBATCH -o ../../../output.dir/Select_Press/ConstFluctSelectLoci.dir/FS_Loci.dir/SLiM_QTL_FS-%j.log
#SBATCH -e ../../../output.dir/Select_Press/ConstFluctSelectLoci.dir/FS_Loci.dir/SLiM_QTL_FS--%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=etb68@mail.missouri.edu
##########SCIENCE FOLLOWS HERE ########################

#module load SLiM

echo -e "=== Bigining of SLiM run with different QTLs > $(date) ==="

output="/storage/hpc/group/kinglab/etb68/evogen-sims/ET_Yeast/output.dir/Select_Press/ConstFluctSelectLoci.dir/FS_Loci.dir/genome5_300_0.5.csv"
if [[ ! -f "$output" ]]
then
  echo "My file ${output} doesn't exist. Running SLiM QTLs now."

slim -d seed=2515 -d repl=1 -d loci=300 -d region=2033 -d mysample=10 -d h=0.5 -d SD=4 -d rang=51 -d gen=20 FS_Loci.slim 
slim -d seed=78514 -d repl=2 -d loci=300 -d region=2033 -d mysample=10 -d h=0.5 -d SD=4 -d rang=51 -d gen=20 FS_Loci.slim 
slim -d seed=11519 -d repl=3 -d loci=300 -d region=2033 -d mysample=10 -d h=0.5 -d SD=4 -d rang=51 -d gen=20 FS_Loci.slim 
slim -d seed=85732 -d repl=4 -d loci=300 -d region=2033 -d mysample=10 -d h=0.5 -d SD=4 -d rang=51 -d gen=20 FS_Loci.slim 
slim -d seed=65741 -d repl=5 -d loci=300 -d region=2033 -d mysample=10 -d h=0.5 -d SD=4 -d rang=51 -d gen=20 FS_Loci.slim

slim -d seed=2515 -d repl=1 -d loci=300 -d region=2033 -d mysample=10 -d h=0.1 -d SD=4 -d rang=51 -d gen=20 FS_Loci.slim
slim -d seed=78514 -d repl=2 -d loci=300 -d region=2033 -d mysample=10 -d h=0.1 -d SD=4 -d rang=51 -d gen=20 FS_Loci.slim
slim -d seed=11519 -d repl=3 -d loci=300 -d region=2033 -d mysample=10 -d h=0.1 -d SD=4 -d rang=51 -d gen=20 FS_Loci.slim
slim -d seed=85732 -d repl=4 -d loci=300 -d region=2033 -d mysample=10 -d h=0.1 -d SD=4 -d rang=51 -d gen=20 FS_Loci.slim
slim -d seed=65741 -d repl=5 -d loci=300 -d region=2033 -d mysample=10 -d h=0.1 -d SD=4 -d rang=51 -d gen=20 FS_Loci.slim

slim -d seed=2515 -d repl=1 -d loci=300 -d region=2033 -d mysample=10 -d h=0.8 -d SD=4 -d rang=51 -d gen=20 FS_Loci.slim
slim -d seed=78514 -d repl=2 -d loci=300 -d region=2033 -d mysample=10 -d h=0.8 -d SD=4 -d rang=51 -d gen=20 FS_Loci.slim
slim -d seed=11519 -d repl=3 -d loci=300 -d region=2033 -d mysample=10 -d h=0.8 -d SD=4 -d rang=51 -d gen=20 FS_Loci.slim
slim -d seed=85732 -d repl=4 -d loci=300 -d region=2033 -d mysample=10 -d h=0.8 -d SD=4 -d rang=51 -d gen=20 FS_Loci.slim
slim -d seed=65741 -d repl=5 -d loci=300 -d region=2033 -d mysample=10 -d h=0.8 -d SD=4 -d rang=51 -d gen=20 FS_Loci.slim



slim -d seed=2515 -d repl=1 -d loci=300 -d region=2033 -d mysample=70 -d h=0.5 -d SD=4 -d rang=51 -d gen=20 FS_Loci.slim
slim -d seed=78514 -d repl=2 -d loci=300 -d region=2033 -d mysample=70 -d h=0.5 -d SD=4 -d rang=51 -d gen=20 FS_Loci.slim
slim -d seed=11519 -d repl=3 -d loci=300 -d region=2033 -d mysample=70 -d h=0.5 -d SD=4 -d rang=51 -d gen=20 FS_Loci.slim
slim -d seed=85732 -d repl=4 -d loci=300 -d region=2033 -d mysample=70 -d h=0.5 -d SD=4 -d rang=51 -d gen=20 FS_Loci.slim
slim -d seed=65741 -d repl=5 -d loci=300 -d region=2033 -d mysample=70 -d h=0.5 -d SD=4 -d rang=51 -d gen=20 FS_Loci.slim

slim -d seed=2515 -d repl=1 -d loci=300 -d region=2033 -d mysample=70 -d h=0.1 -d SD=4 -d rang=51 -d gen=20 FS_Loci.slim
slim -d seed=78514 -d repl=2 -d loci=300 -d region=2033 -d mysample=70 -d h=0.1 -d SD=4 -d rang=51 -d gen=20 FS_Loci.slim
slim -d seed=11519 -d repl=3 -d loci=300 -d region=2033 -d mysample=70 -d h=0.1 -d SD=4 -d rang=51 -d gen=20 FS_Loci.slim
slim -d seed=85732 -d repl=4 -d loci=300 -d region=2033 -d mysample=70 -d h=0.1 -d SD=4 -d rang=51 -d gen=20 FS_Loci.slim
slim -d seed=65741 -d repl=5 -d loci=300 -d region=2033 -d mysample=70 -d h=0.1 -d SD=4 -d rang=51 -d gen=20 FS_Loci.slim

slim -d seed=2515 -d repl=1 -d loci=300 -d region=2033 -d mysample=70 -d h=0.8 -d SD=4 -d rang=51 -d gen=20 FS_Loci.slim
slim -d seed=78514 -d repl=2 -d loci=300 -d region=2033 -d mysample=70 -d h=0.8 -d SD=4 -d rang=51 -d gen=20 FS_Loci.slim
slim -d seed=11519 -d repl=3 -d loci=300 -d region=2033 -d mysample=70 -d h=0.8 -d SD=4 -d rang=51 -d gen=20 FS_Loci.slim
slim -d seed=85732 -d repl=4 -d loci=300 -d region=2033 -d mysample=70 -d h=0.8 -d SD=4 -d rang=51 -d gen=20 FS_Loci.slim
slim -d seed=65741 -d repl=5 -d loci=300 -d region=2033 -d mysample=70 -d h=0.8 -d SD=4 -d rang=51 -d gen=20 FS_Loci.slim

 

slim -d seed=2515 -d repl=1 -d loci=300 -d region=2033 -d mysample=100 -d h=0.5 -d SD=4 -d rang=51 -d gen=20 FS_Loci.slim
slim -d seed=78514 -d repl=2 -d loci=300 -d region=2033 -d mysample=100 -d h=0.5 -d SD=4 -d rang=51 -d gen=20 FS_Loci.slim
slim -d seed=11519 -d repl=3 -d loci=300 -d region=2033 -d mysample=100 -d h=0.5 -d SD=4 -d rang=51 -d gen=20 FS_Loci.slim
slim -d seed=85732 -d repl=4 -d loci=300 -d region=2033 -d mysample=100 -d h=0.5 -d SD=4 -d rang=51 -d gen=20 FS_Loci.slim
slim -d seed=65741 -d repl=5 -d loci=300 -d region=2033 -d mysample=100 -d h=0.5 -d SD=4 -d rang=51 -d gen=20 FS_Loci.slim

slim -d seed=2515 -d repl=1 -d loci=300 -d region=2033 -d mysample=100 -d h=0.1 -d SD=4 -d rang=51 -d gen=20 FS_Loci.slim
slim -d seed=78514 -d repl=2 -d loci=300 -d region=2033 -d mysample=100 -d h=0.1 -d SD=4 -d rang=51 -d gen=20 FS_Loci.slim
slim -d seed=11519 -d repl=3 -d loci=300 -d region=2033 -d mysample=100 -d h=0.1 -d SD=4 -d rang=51 -d gen=20 FS_Loci.slim
slim -d seed=85732 -d repl=4 -d loci=300 -d region=2033 -d mysample=100 -d h=0.1 -d SD=4 -d rang=51 -d gen=20 FS_Loci.slim
slim -d seed=65741 -d repl=5 -d loci=300 -d region=2033 -d mysample=100 -d h=0.1 -d SD=4 -d rang=51 -d gen=20 FS_Loci.slim

slim -d seed=2515 -d repl=1 -d loci=300 -d region=2033 -d mysample=100 -d h=0.8 -d SD=4 -d rang=51 -d gen=20 FS_Loci.slim
slim -d seed=78514 -d repl=2 -d loci=300 -d region=2033 -d mysample=100 -d h=0.8 -d SD=4 -d rang=51 -d gen=20 FS_Loci.slim
slim -d seed=11519 -d repl=3 -d loci=300 -d region=2033 -d mysample=100 -d h=0.8 -d SD=4 -d rang=51 -d gen=20 FS_Loci.slim
slim -d seed=85732 -d repl=4 -d loci=300 -d region=2033 -d mysample=100 -d h=0.8 -d SD=4 -d rang=51 -d gen=20 FS_Loci.slim
slim -d seed=65741 -d repl=5 -d loci=300 -d region=2033 -d mysample=100 -d h=0.8 -d SD=4 -d rang=51 -d gen=20 FS_Loci.slim



slim -d seed=2515 -d repl=1 -d loci=300 -d region=2033 -d mysample=300 -d h=0.5 -d SD=4 -d rang=51 -d gen=20 FS_Loci.slim      
slim -d seed=78514 -d repl=2 -d loci=300 -d region=2033 -d mysample=300 -d h=0.5 -d SD=4 -d rang=51 -d gen=20 FS_Loci.slim
slim -d seed=11519 -d repl=3 -d loci=300 -d region=2033 -d mysample=300 -d h=0.5 -d SD=4 -d rang=51 -d gen=20 FS_Loci.slim
slim -d seed=85732 -d repl=4 -d loci=300 -d region=2033 -d mysample=300 -d h=0.5 -d SD=4 -d rang=51 -d gen=20 FS_Loci.slim      
slim -d seed=65741 -d repl=5 -d loci=300 -d region=2033 -d mysample=300 -d h=0.5 -d SD=4 -d rang=51 -d gen=20 FS_Loci.slim

slim -d seed=2515 -d repl=1 -d loci=300 -d region=2033 -d mysample=300 -d h=0.1 -d SD=4 -d rang=51 -d gen=20 FS_Loci.slim
slim -d seed=78514 -d repl=2 -d loci=300 -d region=2033 -d mysample=300 -d h=0.1 -d SD=4 -d rang=51 -d gen=20 FS_Loci.slim
slim -d seed=11519 -d repl=3 -d loci=300 -d region=2033 -d mysample=300 -d h=0.1 -d SD=4 -d rang=51 -d gen=20 FS_Loci.slim
slim -d seed=85732 -d repl=4 -d loci=300 -d region=2033 -d mysample=300 -d h=0.1 -d SD=4 -d rang=51 -d gen=20 FS_Loci.slim
slim -d seed=65741 -d repl=5 -d loci=300 -d region=2033 -d mysample=300 -d h=0.1 -d SD=4 -d rang=51 -d gen=20 FS_Loci.slim

slim -d seed=2515 -d repl=1 -d loci=300 -d region=2033 -d mysample=300 -d h=0.8 -d SD=4 -d rang=51 -d gen=20 FS_Loci.slim
slim -d seed=78514 -d repl=2 -d loci=300 -d region=2033 -d mysample=300 -d h=0.8 -d SD=4 -d rang=51 -d gen=20 FS_Loci.slim
slim -d seed=11519 -d repl=3 -d loci=300 -d region=2033 -d mysample=300 -d h=0.8 -d SD=4 -d rang=51 -d gen=20 FS_Loci.slim
slim -d seed=85732 -d repl=4 -d loci=300 -d region=2033 -d mysample=300 -d h=0.8 -d SD=4 -d rang=51 -d gen=20 FS_Loci.slim
slim -d seed=65741 -d repl=5 -d loci=300 -d region=2033 -d mysample=300 -d h=0.8 -d SD=4 -d rang=51 -d gen=20 FS_Loci.slim

else

  echo "All is well, Boss.  The ${output} file is there."
fi
echo "=== Fini finito! End of SLiM QTLs Constant Selection run >" $(date)
