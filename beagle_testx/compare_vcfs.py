
import os
from os import path
from tools.input_utilities import read_vcf_allel
import numpy as np

from tools.pop_gen_stats import (
	PI_Taj, Watt_est, TajD
	)

import argparse

parser = argparse.ArgumentParser()

parser.add_argument("input",type=str,metavar= 'N',nargs= '+',
                    help = "Reference files to read. Any number can be given.")

parser.add_argument("--name",type= str,default= 'filename',
	help = "vcf file name")

parser.add_argument("--out",type= str,default= 'table.txt',
	help = "output_table")

args = parser.parse_args()


vcf_arrays= []

for vcf_file in args.input:
	genotype, summary, Names= read_vcf_allel(vcf_file,haps_extract= True)
	vcf_arrays.append(genotype)

gen_proxy= np.sum(vcf_arrays[0],axis= 0)
gen_proxy= gen_proxy != 0
gen_proxy= vcf_arrays[0][:,gen_proxy]
print(gen_proxy.shape)
stats_popg= [PI_Taj(gen_proxy), Watt_est(gen_proxy), TajD(gen_proxy)]

vcf_comp= vcf_arrays[0] != vcf_arrays[1]
vcf_comp= np.array(vcf_comp,dtype= int)
ndiff= np.sum(vcf_comp)
ncells= np.prod(vcf_comp.shape)

filename= args.out

if not path.exists(filename):
	os.makedirs(os.path.dirname(filename), exist_ok=True)
	with open(filename,'w') as fp:
		fp.write('\t'.join(['file','nsamp','ndiffs','ncells','piT','piW','TajD']) + '\n')

with open(filename,"a") as fp:
	line= [args.name,int(vcf_comp.shape[0] / 2),ndiff,ncells,*stats_popg]
	line= np.array(line,dtype=str)

	fp.write('\t'.join(list(line)) + '\n')







