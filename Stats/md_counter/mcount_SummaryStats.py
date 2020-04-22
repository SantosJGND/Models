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
fig_dir= 'Figures/StatsComp'
os.makedirs(fig_dir, exist_ok=True)
fig_dir= fig_dir + '/'


########## species data
from INFO_db import INFO_dict

species= 'chimp'

pop_seg_dict= INFO_dict[species]['pop_dict']
ind_file= INFO_dict[species]['g_indfile']

####
#### stats_file
stats_file= {}

target_dirs= {}
##############

target= 'rd'
stats_file[target]= {}
## directories
main_dir= os.getcwd() + '/'
count_dir= main_dir + 'mutation_counter/count/'
dir_launch= main_dir + 'mutation_counter'
muted_dir= main_dir + 'mutation_counter/data/mutation_count/'
sims_dir= INFO_dict[species]['dirs'][target]
indfile= INFO_dict[species]['w_indfile'][target]

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


data, data_freqs = MC_sample_matrix_simple(min_size= min_size, samp= sampling, stepup= stepup, count_dir= count_dir, 
                        dir_launch= dir_launch,main_dir= main_dir,sim_dir= sims_dir,indfile= indfile,
                          muted_dir= muted_dir, diffs= diffs,collapsed= collapsed,row= row, col= col,
                       exclude= False,sample_sim= sample_sim,freq_extract= freq_extract,segregating= segregating)



vcft_func= 'weir-fst-pop'

fst_dict= Fst_tools({},sim_dir= sims_dir, indfile= indfile,fig_dir= fig_dir,vcft_func= vcft_func,sample_sim=sample_sim)

tag= sims_dir.split('/')[-2]
target_dirs[target]= tag
################ COUNTS
################
pops_dict= get_pops(indfile= ind_file)
pops_dict= {z:g for z,g in pops_dict.items() if z in pop_seg_dict.keys()}
data_avail= list(data.keys())
#pops= data[data_avail[0]]['counts'].keys()
pops= pops_dict.keys()
pops= sorted(list(pops))

pop_counts= {
    pop: [np.sum(data[z]['counts'][pop]) for z in data_avail] for pop in pops
}

stats_file[target]['counts']= pop_counts

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
    z: [x[0] for x in g] for z,g in comb_dict.items()
}

stats_file[target]['distances']= comb_dict

####################### FST
#######################

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

stats_file[target]["FST"]= stats

###########################################################
###########################################################
############# SIMS          ###############################
###########################################################

target= 'sims'
stats_file[target]= {}
pop_seg_dict= {
    g:v for v,g in pop_seg_dict.items()
}

pop_sizes= {
    v: len(g) for v,g in pops_dict.items() if v in pop_seg_dict.values()
}


sim_pops= sorted(list(pop_sizes.keys()))
print(sim_pops)
pop_combs= list(it.combinations(sim_pops,2))

sims_dir= INFO_dict[species]['dirs'][target]
indfile= INFO_dict[species]['w_indfile'][target]

request_processed= gzip_request(dir_check= sims_dir,str_format= '',requested= ['.vcf','fa'], func= 'gzip')

sample_sim= 0

data_sims, data_freqs_sims = MC_sample_matrix_dict(pop_seg_dict,pop_sizes,min_size= min_size, samp= 20, stepup= stepup, count_dir= count_dir, 
                        dir_launch= dir_launch,main_dir= main_dir,sim_dir= sims_dir,indfile= indfile,
                          muted_dir= muted_dir, diffs= diffs,collapsed= collapsed,row= row, col= col,
                       exclude= False,sample_sim= sample_sim,freq_extract= freq_extract,segregating= segregating)


pop_seg_dict= {
    g:v for v,g in pop_seg_dict.items()
}

fst_dict= Fst_tools(pop_seg_dict,sim_dir= sims_dir, indfile= indfile,fig_dir= fig_dir,vcft_func= vcft_func,sample_sim= sample_sim)

tag= sims_dir.split('/')[-2]
target_dirs[target]= tag

################
################ COUNTS
sims_avail= list(data_sims.keys())

sim_pop_counts= {
    z: [data_sims[x]['Nvars'][z] for x in sims_avail] for z in sim_pops
}

stats_file[target]['counts']= sim_pop_counts

####################
#################### Distances

comb_dict= {
    z: [] for z in pop_combs
}


for av in sims_avail:
    for comb in pop_combs:
        comb_dict[comb].append(data_sims[av]['pairDist'][comb][0])

comb_dict= {
    z: [x[0] for x in g] for z,g in comb_dict.items()
}

stats_file[target]['distances']= comb_dict

###################### 
###################### FST

combs= list(fst_dict[0].keys())

stats= {
    x: [] for x in combs
}

for idx in fst_dict.keys():
    for comb,g in fst_dict[idx].items():
        if g:
            mean_fst= np.mean(g)

            stats[comb].append(mean_fst)


stats_file[target]['FST']= stats


########################################################
########################################################
from scipy.stats import pearsonr

Nreps= 100
target_vector= stats_file.keys()
tag_names= '_'.join(target_vector)
reps_dict= {}

for target in target_vector: 
    for stat_av in stats_file[target].keys():
        if stat_av not in reps_dict.keys():
            reps_dict[stat_av]= {}
        print(stat_av)
        for op in stats_file[target][stat_av].keys():
            avail_samples= stats_file[target][stat_av][op]
            avail_samples= list(avail_samples)
            #print(avail_samples[:10])
            samp_here= np.random.choice(avail_samples,Nreps,replace= True)

            if op not in reps_dict[stat_av].keys():
                reps_dict[stat_av][op]= {}

            reps_dict[stat_av][op][target]= samp_here


combined_dict= {x:[] for x in stats_file.keys()}

for stat_av in reps_dict.keys():
    sd_dict= {z:g for z,g in reps_dict[stat_av].items() if sum([int(y in stats_file.keys()) for y in g.keys()]) == 2}

    sg_keys= sorted(sd_dict.keys())
    xvec= []
    yvec= []

    for op in sg_keys:
        xvec.extend(sd_dict[op]['rd'])
        yvec.extend(sd_dict[op]['sims'])

    combined_dict['rd'].extend(xvec)
    combined_dict['sims'].extend(yvec)
    
    coeff= pearsonr(xvec,yvec)[0]
    ########## plot

    fig,axis = plt.subplots(1,1, figsize=(10, 10))

    axis.set_xlabel(stat_av)
    plt.title(str(round(coeff,3)))
    plt.xticks(fontsize= 20)
    plt.yticks(fontsize= 20)

    # Create an axes instance
    #ax = fig.add_axes(pops)
    # Create the boxplot
    #bp = ax.violinplot(violin_data)
    plt.scatter(xvec,yvec)

    plt.xlabel('real data',fontsize=20)

    plt.ylabel('Simulations',fontsize= 20)

    #plt.ylim(0,.5)
    plt.savefig(fig_dir + 'statsComp_{}_{}_sims.png'.format(stat_av,tag_names),bbox_inches='tight')
    plt.close()


fig,axis = plt.subplots(1,1, figsize=(10, 10))

plt.xticks(fontsize= 20)
plt.yticks(fontsize= 20)

# Create an axes instance
#ax = fig.add_axes(pops)
# Create the boxplot
#bp = ax.violinplot(violin_data)
xvec= combined_dict['rd']
yvec= combined_dict['sims']
coeff= pearsonr(xvec,yvec)[0]
plt.scatter(xvec,yvec)

plt.xlabel('real data',fontsize=20)

plt.ylabel('Simulations',fontsize= 20)

#plt.ylim(0,.5)
plt.savefig(fig_dir + 'statsComp_ALL_{}_{}_sims.png'.format(stat_av,tag_names),bbox_inches='tight')
plt.close()


