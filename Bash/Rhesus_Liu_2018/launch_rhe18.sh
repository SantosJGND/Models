#!/bin/bash
#SBATCH -n 5
#SBATCH -t 50:00:00
#SBATCH --mem=15GB


module purge
module load python/3.6.4
module load slim/3.3.1

short=$1

python -u  SLiM_demosABC.py \
-d "demos/rhesus_liu18.txt" \
-L 10000000 \
-R demos_mat/template_matVar.slim \
-N 5 \
-c ABC \
-r 1 \
--mem '25GB' \
-t '25:00:00' \
--nodes 6 \
-s $short \

