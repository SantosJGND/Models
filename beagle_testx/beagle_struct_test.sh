#!/bin/bash
#SBATCH --nodes=2
#SBATCH --time=6:00:00
#SBATCH --job-name=BeagleStrucTest
#SBATCH --mem=7GB

module purge
module load beagle/5.1
module load vcftools/0.1.17
module load python/3.6.2

cd /home/x_garciaj/beagle

source ~/miniconda3/etc/profile.d/conda.sh
conda activate ./env

HOME=/home/x_garciaj/beagle/

slim_dir=$HOME$1
indfile=$HOME"inds.txt"
ind_assignments=$HOME$2
outfile=$HOME"diffs_structstats_db.txt"

nsamp=100

tkey=$RANDOM
temp_dir=$HOME"temp"$tkey

mkdir $temp_dir

temp_dir=$temp_dir"/"

tempfile=$temp_dir"tempfile.txt"
temptwo=$temp_dir"temptwo.txt"
tempsamp=$temp_dir"tempsamp.vcf"
tempout=$temp_dir"OUTPUT.vcf"
temphase=$temp_dir"tempphase"

echo $temp_dir
echo $tempfile
echo $tempout
echo $temphase

cd $slim_dir

for file in ./*/*.vcf.gz; do
echo $file

for rep in {1..10}; do

shuf -n $nsamp $ind_assignments > $temptwo

cut -f1 $temptwo > $tempfile

vcf-subset -c $tempfile $file > $tempsamp

sed '/^##/! s/|/\//g' $tempsamp > $tempout

java -jar $BEAGLE/beagle.jar gt=$tempout out=$temphase burnin=18 iterations=24

python -u /home/x_garciaj/beagle/compare_vcfs.py \
$tempout $temphase".vcf.gz" \
--name $file \
--out $outfile

for pop in `cut -f2 $temptwo | awk '!_[$1]++'`; do

echo $pop
grep $pop $temptwo | cut -f1 > $tempfile

vcf-subset -c $tempfile $file > $tempsamp

sed '/^##/! s/|/\//g' $tempsamp > $tempout

java -jar $BEAGLE/beagle.jar gt=$tempout out=$temphase burnin=18 iterations=24

python -u /home/x_garciaj/beagle/compare_vcfs.py \
$tempout $temphase".vcf.gz" \
--name $file"."$pop \
--out $outfile
done

done

rm $temp_dir*
done

rm -r $temp_dir

conda deactivate
