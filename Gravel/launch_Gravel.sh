#!/bin/bash
#SBATCH -n 1
#SBATCH -t 30:00:00
#SBATCH --mem=4GB

module purge
module load python/3.6.4
module load slim/3.3.1

short=$1

python -u  SLiM_Gravel.py \
-L 80000000 \
-r Human_sims/Gravel_2011_frame_sample.slim \
-N 40 \
-c GravelSimple \
-s $short \

