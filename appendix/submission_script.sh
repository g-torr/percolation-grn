#!/bin/bash -l

#SBATCH -p nms_research

# How much memory do you need **per job**?
#SBATCH --mem=32G
#SBATCH --nodes=1
# Number of cores you need 
#SBATCH --cpus-per-task=20
#Work in the current directory. 
#$ -cwd
  
# load any modules you need 
module load devtools/anaconda/2019.3-python3.7.3 

#set email address
#$ -M giuseppe.torrisi@kcl.ac.uk

#send email when job finishes
#$ -m e 

# Run thing that uses multiple cores 
#python /users/k1762355/brc_scratch/dynamics/regular\ graph/spin-glass/multiple_states/cycle.py  --threads 40
# Time reserved for execution 
#$ -l h_rt=50:00:00

python knockout_cascade.py --threads 20
