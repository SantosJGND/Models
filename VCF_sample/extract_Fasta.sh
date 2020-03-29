#!/bin/bash
#SBATCH -n 1
#SBATCH -t 2:00:00
#SBATCH --mem=3GB

module purge

module load python/3.6.4

chr_request="1"
fasta="/home/jgarc235/Fastas/panTro2.fa.gz"

python -u /home/jgarc235/Rhesus/bash_commands/Extract/Extract_chr_Fasta.py -c $chr_request -f $fasta
