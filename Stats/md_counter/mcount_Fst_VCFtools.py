from tools.input_utilities import (
    MC_sample_matrix_simple, MC_sample_matrix_dict
)

from tools.compare_utilities import (
    gzip_request
)

from tools.fasta_utilities import (
    get_mutations, kmer_comp_index, kmer_mut_index
    )


from tools.mcounter_tools import (
    Fst_tools, get_pops
    )

import re
import pandas as pd
import os
import numpy as np
import itertools as it
import collections
import time

def recursively_default_dict():
    return collections.defaultdict(recursively_default_dict)


import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt


############## data set info
from INFO_db import INFO_dict


##########
##########
##
##
fig_dir= 'Figures/FST'
os.makedirs(fig_dir, exist_ok=True)
fig_dir= fig_dir + '/'

##########
##########

species= 'rhesus'

pop_seg_dict= INFO_dict[species]['pop_dict']

## Reql Data
main_dir= os.getcwd() + '/'
count_dir= main_dir + 'mutation_counter/count/'
dir_launch= main_dir + 'mutation_counter'
sims_dir= INFO_dict[species]['dirs']['rd']
indfile= INFO_dict[species]['w_indfile']['rd']

sample_sim= 0
sim_del= 'C'

chrom_idx= 0
request_processed= gzip_request(dir_check= sims_dir,str_format= '',requested= ['.vcf','fa'], func= 'gzip')

vcft_func= 'weir-fst-pop'

fst_dict= Fst_tools({},sim_dir= sims_dir, indfile= indfile,fig_dir= fig_dir,vcft_func= vcft_func,sample_sim=sample_sim)


####### Sort
combs= list(fst_dict[0].keys())
combs= [x for x in combs if sum([int(y in pop_seg_dict.keys()) for y in x]) == 2]
stats= {
    x: [] for x in combs
}

for idx in fst_dict.keys():
    for comb,g in fst_dict[idx].items():
        if comb in combs:
            if g:
                mean_fst= np.mean(g)
                stats[comb].append(mean_fst)

####### PLOT

tag= sims_dir.split('/')[-2]

fig,axis = plt.subplots(1,1, figsize=(10, 10))

axis.set_xlabel('pop combinations')

axis.set_xticks(list(range(1,len(combs)+1)))
axis.set_xticklabels(['{}'.format(x) for x in combs])
plt.xticks(rotation=45,fontsize= 20)
plt.yticks(fontsize= 20)

# Create an axes instance
#ax = fig.add_axes(pops)

violin_data= [stats[z] for z in combs]
# Create the boxplot
#bp = ax.violinplot(violin_data)
plt.title('segregating count')
axis.violinplot(violin_data,showmeans=False, showmedians=True,
        showextrema=False)

plt.xlabel('pop pairs',fontsize=20)

plt.ylabel('FST',fontsize= 20)
plt.ylim(0,.2)

#plt.ylim(0,.5)
plt.savefig(fig_dir + 'violinFST_{}_rd.png'.format(tag),bbox_inches='tight')
plt.close()


#############################################
#############################################
## Sims
main_dir= os.getcwd() + '/'
count_dir= main_dir + 'mutation_counter/count/'
dir_launch= main_dir + 'mutation_counter'
sims_dir= INFO_dict[species]['dirs']['sims']
indfile= INFO_dict[species]['w_indfile']['sims']

sample_sim= 0
sim_del= 'C'

chrom_idx= 0
request_processed= gzip_request(dir_check= sims_dir,str_format= '',requested= ['.vcf','fa'], func= 'gzip')

vcft_func= 'weir-fst-pop'

fst_dict= Fst_tools(pop_seg_dict,sim_dir= sims_dir, indfile= indfile,fig_dir= fig_dir,vcft_func= vcft_func,sample_sim= sample_sim)

####### Sort
combs= list(fst_dict[0].keys())

stats= {
    x: [] for x in combs
}

for idx in fst_dict.keys():
    for comb,g in fst_dict[idx].items():
        if g:
            mean_fst= np.mean(g)
            stats[comb].append(mean_fst)

####### PLOT

tag= sims_dir.split('/')[-2]

fig,axis = plt.subplots(1,1, figsize=(10, 10))

axis.set_xlabel('pop combinations')

axis.set_xticks(list(range(1,len(combs)+1)))
axis.set_xticklabels([' / '.join(x) for x in combs])
plt.xticks(rotation=45,fontsize= 20)
plt.yticks(fontsize= 20)

# Create an axes instance
#ax = fig.add_axes(pops)

violin_data= [stats[z] for z in combs]
# Create the boxplot
#bp = ax.violinplot(violin_data)
plt.title('segregating count')
axis.violinplot(violin_data,showmeans=False, showmedians=True,
        showextrema=False)

plt.xlabel('pop pairs',fontsize=20)

plt.ylabel('FST',fontsize= 20)
plt.ylim(0,.2)
#plt.ylim(0,.5)
plt.savefig(fig_dir + 'violinFST_{}_sims.png'.format(tag),bbox_inches='tight')
plt.close()

