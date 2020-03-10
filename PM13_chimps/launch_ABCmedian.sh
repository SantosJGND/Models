#!/bin/bash
#SBATCH -n 1
#SBATCH -t 30:00:00
#SBATCH --mem=4GB

module purge
module load python/3.6.4
module load slim/3.3.1

short=$1

python -u SLiM_ABCmedian.py \
-L 1000000 \
-r M4A_median.slim \
-N 40 \
-c PMchimpMedian \
-s $short \

