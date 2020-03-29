#!/bin/bash
#SBATCH -n 2
#SBATCH -t 12:00:00
#SBATCH --mem=9GB

module purge
module load vcftools/0.1.12a
module load tabix/0.2.6
module load python/3.6.4

#vcf_file="/home/jgarc235/Human/data_vcf/phase1_chr1_filtered.vcf.gz.recode.vcf.gz"
vcf_file="/scratch/jgarc235/chr1_pt.vcf.gz"
N=50
L=1000000

assembly="chr1_panTro2"
batch="GAGP_Pt"
out_dir="/home/jgarc235/Chimp/vcf_data/chimp_1MB/"
diffs=""
ids="/home/jgarc235/Chimp/chimp_ind_assignments.txt"

chrom="1"

echo $vcf_file

python -u ordered_extractions.py -v $vcf_file -c $chrom -n $N -l $L -a $assembly \
-b $batch -i $ids -o $out_dir

