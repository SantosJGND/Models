#!/bin/bash
#SBATCH -n 10
#SBATCH -t 30:00:00
#SBATCH --mem=44GB

module purge
module load python/3.6.4

echo mcount_plot_VarIncrement.py

python -u mcount_plot_VarIncrement.py \
--species human \
--samp 800


