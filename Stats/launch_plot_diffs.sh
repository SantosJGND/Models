#!/bin/bash
#SBATCH -n 10
#SBATCH -t 25:00:00
#SBATCH --mem=40GB

module purge
module load python/3.6.4

echo mcount_plot_diffs.py

python -u mcount_plot_diffs.py \
--species human \
--stepup increment \
--samp 800 \
--steps 50 \
--reps 5 \
