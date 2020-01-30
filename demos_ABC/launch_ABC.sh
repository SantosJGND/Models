#!/bin/bash
#SBATCH -t 15 -p debug -q wildfire

module purge
module load python/3.6.4
module load slim/3.3.1

short=$1

python -u  SLiM_demosABC.py \
-d "demos/PM2013_M3.txt" \
-L 100000  \
-r demos_mat/template_matVar.slim \
-N 4 \
-c ABC \
-s $short \

