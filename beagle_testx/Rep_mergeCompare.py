
import os
from os import path
from tools.input_utilities import read_vcf_allel
import numpy as np
import itertools as it

from tools.pop_gen_stats import (
	PI_Taj, Watt_est, TajD
	)

import argparse

parser = argparse.ArgumentParser()

parser.add_argument("--vcf",type= str,default= '../temp/OUTPUT.vcf.gz',
	help = "vcf file name")

parser.add_argument("--dir",type= str,default= '../temp/',
	help = "temp_dir")

parser.add_argument("--pref",type= str,default= 'temp',
	help = "prefix of phased files.")

parser.add_argument("--nrep",type= int,default= 10,
	help = "number of repeats")

parser.add_argument("--name",type= str,default= 'filename',
	help = "vcf file name")

parser.add_argument("--out",type= str,default= 'table.txt',
	help = "output_table")

args = parser.parse_args()


######
######
geno_ref, summary, Names= read_vcf_allel(args.vcf,haps_extract= True)

######
######
temp_dir= args.dir
temp_pref= args.pref
temp_pref= temp_pref.split('/')[-1]
temp_range= args.nrep

vcf_arrays= []

print('###')
print(args.pref)
print(temp_dir)
print(os.listdir(temp_dir))

for file in os.listdir(temp_dir):
	if os.path.isfile(os.path.join(temp_dir,file)) and temp_pref in file and 'vcf' in file:
		vcf_file= os.path.join(temp_dir,file)

		genotype, summary, Names= read_vcf_allel(vcf_file,haps_extract= True)
		vcf_arrays.append(genotype)


print('{} temp phased files found.'.format(len(vcf_arrays)))

######
######

header=['file','nsamp','ncells','nphas','ndiffs','std']
nmax= 20
proxy_range= list(range(len(vcf_arrays)))

for idx in range(len(vcf_arrays)):

	xrang= idx +1
	poss_combs= it.combinations(proxy_range,xrang)
	poss_combs= list(poss_combs)

	if len(poss_combs) > nmax:
		prox_idx= list(range(len(poss_combs)))
		prox_idx= np.random.choice(prox_idx,10)
		poss_combs= [poss_combs[x] for x in prox_idx]

	diffs_vector= []

	for combs in poss_combs:

		phase_select= [vcf_arrays[x] for x in combs]

		phase_select= np.array(phase_select)

		phase= np.median(phase_select,axis= 0)
		final= np.random.binomial(1, phase) 

		vcf_comp= geno_ref != final

		vcf_comp= np.array(vcf_comp,dtype= int)
		ndiff= np.sum(vcf_comp)

		diffs_vector.append(ndiff)

	ncells= np.prod(vcf_comp.shape)

	filename= args.out

	if not path.exists(filename):
		os.makedirs(os.path.dirname(filename), exist_ok=True)
		with open(filename,'w') as fp:
			fp.write('\t'.join(header) + '\n')

	with open(filename,"a") as fp:
		line= [args.name,int(vcf_comp.shape[0] / 2), ncells, xrang, np.mean(diffs_vector), np.std(diffs_vector)]
		line= np.array(line,dtype=str)

		fp.write('\t'.join(list(line)) + '\n')







