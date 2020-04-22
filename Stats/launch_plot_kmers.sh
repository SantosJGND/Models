#!/bin/bash
#SBATCH -n 10
#SBATCH -t 80:00:00
#SBATCH --mem=40GB

module purge
module load python/3.6.4

echo mcount_plot_kmers.py

python -u mcount_plot_kmers.py \
--species chimp \
--stepup increment \
--samp 30 \
--steps 30 \
--reps 1 \
