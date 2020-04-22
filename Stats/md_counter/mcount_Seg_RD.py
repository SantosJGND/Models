from tools.input_utilities import (
    MC_sample_matrix_simple, MC_sample_matrix_dict
)

from tools.compare_utilities import (
    gzip_request
)

from tools.fasta_utilities import (
    get_mutations, kmer_comp_index, kmer_mut_index
    )



from tools.mcounter_tools_new import (
    Fst_tools, get_pop_dict
    )

#from tools.SLiM_pipe_tools import mutation_counter_launch
import re
import pandas as pd
import os
import numpy as np
import itertools as it
import collections

def recursively_default_dict():
    return collections.defaultdict(recursively_default_dict)


import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt

##
fig_dir= 'Figures/SegCount'
os.makedirs(fig_dir, exist_ok=True)
fig_dir= fig_dir + '/'


############## data set info
from INFO_db import INFO_dict


species= 'chimp'

pop_seg_dict= INFO_dict[species]['pop_dict']
ind_file= INFO_dict[species]['g_indfile']

## directories
main_dir= os.getcwd() + '/'
count_dir= main_dir + 'mutation_counter/count/'
dir_launch= main_dir + 'mutation_counter'
muted_dir= main_dir + 'mutation_counter/data/mutation_count/'
sims_dir= INFO_dict[species]['dirs']['rd']
indfile= INFO_dict[species]['w_indfile']['rd']
sims_target= sims_dir.split('/')[-2]

diffs= False
mutlog= 'toMut.log'
min_size= 5
sampling= [2,2,5]
stepup= 'increment'
sample_sim= 0
collapsed= False

row= [64,32][int(collapsed)]
col= 3
freq_extract= False
segregating= True
sim_del= 'C'
chrom_idx= 0
request_processed= gzip_request(dir_check= sims_dir,str_format= '',requested= ['.vcf','fa'], func= 'gzip')
return_private= True
distances= 'PCA'


data, data_freqs = MC_sample_matrix_simple(min_size= min_size, samp= sampling, stepup= stepup, count_dir= count_dir, 
                        dir_launch= dir_launch,main_dir= main_dir,sim_dir= sims_dir,indfile= indfile,
                          muted_dir= muted_dir, diffs= diffs,collapsed= collapsed,row= row, col= col, 
                       exclude= False,sample_sim= sample_sim,freq_extract= freq_extract,segregating= segregating,
                       return_private= return_private, distances= distances)





pops_dict= get_pop_dict(indfile= ind_file)
pops_dict= {z:g for z,g in pops_dict.items() if z in pop_seg_dict.keys()}
data_avail= list(data.keys())
#pops= data[data_avail[0]]['counts'].keys()
pops= pops_dict.keys()
pops= sorted(list(pops))

pop_counts= {
    pop: [np.sum(data[z]['counts'][pop]) for z in data_avail] for pop in pops
}



########
########
tag= sims_dir.split('/')[-2]

fig,axis = plt.subplots(1,1)

axis.set_xlabel('populations')
axis.set_xticks(list(range(1,len(pops)+1)))
axis.set_xticklabels(['{} N:{}'.format(x,len(pops_dict[x])) for x in pops])
# Create an axes instance
#ax = fig.add_axes(pops)

violin_data= [pop_counts[z] for z in pops]
# Create the boxplot
#bp = ax.violinplot(violin_data)
plt.title('segregating count')
axis.violinplot(violin_data,showmeans=True)

plt.savefig(fig_dir + 'violin_{}_rdata.png'.format(tag),bbox_inches='tight')
plt.close()



#################### Distances
####################

pop_combs= list(it.combinations(pops,2))
comb_dict= {
    z: [] for z in pop_combs
}

for av in data_avail:
    for comb in pop_combs:
        comb_dict[comb].append(data[av]['pairDist'][comb][0])

comb_dict= {
    z: np.array(g).reshape(1,-1) for z,g in comb_dict.items()
}


tag= sims_dir.split('/')[-2]

fig,axis = plt.subplots(1,1)

axis.set_xlabel('pairs')
axis.set_xticks(list(range(1,len(pop_combs)+1)))
axis.set_xticklabels(['-'.join(x) for x in pop_combs])
plt.xticks(rotation=45)
# Create an axes instance
#ax = fig.add_axes(pops)

violin_data= [comb_dict[z][0] for z in pop_combs]
# Create the boxplot
#bp = ax.violinplot(violin_data)
plt.title('Pairwise distances')
axis.violinplot(violin_data,showmeans=True)

plt.savefig(fig_dir + 'violin_{}_distances_rdata.png'.format(tag),bbox_inches='tight')
plt.close()


#############
############# SIMS
#############

pop_seg_dict= {
    g:v for v,g in pop_seg_dict.items()
}

pop_sizes= {
    v: len(g) for v,g in pops_dict.items() if v in pop_seg_dict.values()
}

sim_pops= sorted(list(pop_sizes.keys()))
print(sim_pops)
pop_combs= list(it.combinations(sim_pops,2))

sims_dir= INFO_dict[species]['dirs']['sims']
indfile= INFO_dict[species]['w_indfile']['sims']

request_processed= gzip_request(dir_check= sims_dir,str_format= '',requested= ['.vcf','fa'], func= 'gzip')

sample_sim= 0

data_sims, data_freqs_sims = MC_sample_matrix_dict(pop_seg_dict,pop_sizes,min_size= min_size, samp= 20, stepup= stepup, count_dir= count_dir, 
                        dir_launch= dir_launch,main_dir= main_dir,sim_dir= sims_dir,indfile= indfile,
                          muted_dir= muted_dir, diffs= diffs,collapsed= collapsed,row= row, col= col,
                       exclude= False,sample_sim= sample_sim,freq_extract= freq_extract,segregating= segregating,
                       return_private= return_private)


sims_avail= list(data_sims.keys())
##pop_vector= [data_sims[av][] for av in sims_avail]

sim_pop_counts= {
    z: [data_sims[x]['Nvars'][z] for x in sims_avail] for z in sim_pops
}



#################
#################
tag= sims_dir.split('/')[-2]

fig,axis = plt.subplots(1,1)

axis.set_xlabel('populations')
axis.set_xticks(list(range(1,len(sim_pops)+1)))
axis.set_xticklabels(['{} N:{}'.format(x,len(pops_dict[x])) for x in sim_pops])
# Create an axes instance
#ax = fig.add_axes(pops)

violin_data= [sim_pop_counts[z] for z in sim_pops]
# Create the boxplot
#bp = ax.violinplot(violin_data)
plt.title('segregating count')
axis.violinplot(violin_data,showmeans=True)

plt.savefig(fig_dir + sims_target + 'violin_{}_sims.png'.format(tag),bbox_inches='tight')
plt.close()


####################
####################

comb_dict= {
    z: [] for z in pop_combs
}


for av in sims_avail:
    for comb in pop_combs:
        comb_dict[comb].append(data_sims[av]['pairDist'][comb][0])

comb_dict= {
    z: np.array(g).reshape(1,-1) for z,g in comb_dict.items()
}


tag= sims_dir.split('/')[-2]

fig,axis = plt.subplots(1,1)

axis.set_xlabel('pairs')
axis.set_xticks(list(range(1,len(pop_combs)+1)))
axis.set_xticklabels(['-'.join(x) for x in pop_combs])
plt.xticks(rotation=45)
# Create an axes instance
#ax = fig.add_axes(pops)

violin_data= [comb_dict[z][0] for z in pop_combs]
# Create the boxplot
#bp = ax.violinplot(violin_data)
plt.title('Pairwise distances')
axis.violinplot(violin_data,showmeans=True)

plt.savefig(fig_dir + sims_target + 'violin_{}_distances_sims.png'.format(tag),bbox_inches='tight')
plt.close()
