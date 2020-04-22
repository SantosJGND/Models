from tools.input_utilities import (
    MC_sample_matrix_v1
)

from tools.mcounter_tools import (
    md_SubSampVar
)

from tools.compare_utilities import (
    gzip_request
)


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



import argparse

parser = argparse.ArgumentParser()

parser.add_argument('--species', type=str,
                    default='chimp')

parser.add_argument('--minS', type=int, # ignore pops if Nsamp < minS
                    default=5)

parser.add_argument('--samp', type=int, # if stepup= increment = max Nsamp, else = min (Nsamp).
                    default=200)

parser.add_argument('--steps', type=int, # Nber of Nsamp steps
                    default=50)

parser.add_argument('--reps', type=int, # replicates per sample size.
                    default=1)

parser.add_argument('--stepup', type=str, # type of analysis - increment / other.
                    default='increment')

parser.add_argument('--simsN', type=int, # Numberr of simulations to read from available. random. 
                    default=0)           # if 0 uses all sims in sims_dir. 

parser.add_argument('--collapsed', type=bool, # collapse mutation counts.
                    default=False)

parser.add_argument('--freqs', type=bool, # return allele frequencies.
                    default=False)

parser.add_argument('--ploidy', type=int, # 
                    default=2)

parser.add_argument('--haps', type=bool, # analyse by haplotype.
                    default=True)

parser.add_argument('-d', '--diffs', type=bool, # read file of ancestral SNPs and polarise.
                    default=False)

args = parser.parse_args()


from INFO_db import INFO_dict

species= args.species

####################################### 
#######################################
######## I. DIRECTORIES

fig_dir= 'Figures/Var_Samp'
os.makedirs(fig_dir, exist_ok=True)
fig_dir= fig_dir + '/'


main_dir= os.getcwd() + '/'
count_dir= main_dir + 'mutation_counter/count/'
dir_launch= main_dir + 'mutation_counter'
muted_dir= main_dir + 'mutation_counter/data/mutation_count/'
sims_dir= INFO_dict[species]['dirs']['sims']
sims_target= sims_dir.split('/')[-2]
mutlog= 'toMut.log'

indfile= INFO_dict[species]['w_indfile']['sims']


######################################
######################################
##### II. READ DATA

sampling= [args.samp,args.steps,args.reps]
print(sampling)

row= [64,32][int(args.collapsed)]
col= 3

request_processed= gzip_request(dir_check= sims_dir,str_format= '',requested= ['.vcf','fa'], func= 'gzip')


data, data_freqs = MC_sample_matrix_v1(min_size= args.minS, samp= sampling, stepup= args.stepup, count_dir= count_dir, 
                        dir_launch= dir_launch,main_dir= main_dir,sim_dir= sims_dir,indfile= indfile,
                          muted_dir= muted_dir, diffs= args.diffs,collapsed= args.collapsed,row= row, col= col, ploidy= args.ploidy,
                       sample_sim= args.simsN,freq_extract= args.freqs,haps_extract=args.haps)


#####################################
#####################################
########## III. COUNT COMPARISONS
##########
bsep= 'C'

stats_dict, counts_dict, ref_dict= md_SubSampVar(data,row= row,col= col,bsep= bsep)


#####################################
####################################
######## IV. Extract and Sort.
######## std distribution.
######## Across sub_samples:

from scipy.stats import norm
threshold= 2 # std

pop_diff_stats= {
    z: [stats_dict[z][x] for x in sorted([x for x in stats_dict[z].keys() if x not in ['sizes', 'counts']])] for z in stats_dict.keys()
}

pop_diff_dict= {
    z: {
        'mean': [np.mean(x) for x in g],
        'std': [np.std(x) for x in g],
    } for z,g in pop_diff_stats.items()
}

######### Across references.
ref_stats= {
    z: {
        'mean': np.mean(g),
        'ucl': np.mean(g) + threshold*np.std(g),
        'lcl': np.mean(g) - threshold*np.std(g),
        'std': np.std(g)
    } for z,g in ref_dict.items()
}

average_cdf= {
    z: [norm.cdf(x,loc= ref_stats[z]['mean'],scale= ref_stats[z]['std']) for x in pop_diff_dict[z]['mean']] for z in pop_diff_dict.keys()
}


####################################
####################################
########## IV. PLOT.
##########
########## i. ref_stats distribution - STD acrss sample sizes and for full-pops. 

pop_colors= INFO_dict[species]['pop_colors']
pop_names_dict= {g:z for z,g in INFO_dict[species]['pop_dict'].items()}

std_margin= 1
threshold= 0.01
from scipy.stats import norm

ylab= 'sum of squared differences'
xlab= 'sample size'
for pop in pop_diff_dict.keys():
    pop_name= pop_names_dict[pop]
    for plot_var in [0,1]:
        plt.figure(figsize=(10,10))
        
        surface= counts_dict[pop]['sizes']
        y= pop_diff_dict[pop]['mean']
        error= pop_diff_dict[pop]['std']
        if plot_var:
            plt.errorbar(surface,y,yerr=error,label= 'var')   
        else:
            plt.plot(surface,y,label= 'mean var')

        line_surf= [min(counts_dict[pop]['sizes']),max(counts_dict[pop]['sizes'])]

        for z in ['mean','ucl','lcl']:
            y= [ref_stats[pop][z]] * 2
            plt.plot(line_surf,y,label= z)

        plt.xlabel(xlab)

        plt.ylabel(ylab)
        plt.title('grid SSD / sample proportion - control - {}'.format(pop_name), fontsize= 20)

        plt.legend()
        plt.savefig(fig_dir + sims_target + '_' + 'gridSSD_{}_{}_sbsD.png'.format(pop_name,["",'wSD'][plot_var]),bbox_inches='tight')
        plt.close()


###########################################
###########################################
#### ii. for pop distribution - histograms
nbins= 50
opacity= .7

rates_show= 3
steps= np.linspace(0,len(surface)-1,rates_show,dtype= int)


for pop in pop_diff_dict.keys():
    pop_name= pop_names_dict[pop]
    #
    dtup= [ref_dict[pop]] + [pop_diff_stats[pop][x] for x in steps]
    dtup_len= [len(x) for x in dtup]

    min_len= min(dtup_len)
    sample_dtup= [np.random.choice(x, min_len,replace= False) for x in dtup]
    #

    plt.figure(figsize=(10,10))


    plt.hist(sample_dtup[0],
         color = 'blue', label= 'full population',bins= nbins, alpha= .8)

    for idx in range(rates_show):
        plt.hist(sample_dtup[idx+1], alpha= opacity, bins= nbins,
            label= str(surface[steps[idx]]))

    plt.xlabel("MS sd - nbins: {}".format(nbins),fontsize= 20)

    plt.ylabel("frequency",fontsize= 20)
    plt.title(str(pop),fontsize= 20)

    plt.legend()
    plt.savefig(fig_dir + sims_target + '_MSsd_hist_{}_sbsD.png'.format(pop_name),bbox_inches='tight')
    plt.close()


#############################################
############### iii. convergence Across populations. 
###############
ylab= 'mean p-val'
plt.figure(figsize=(10,10))

for pop in pop_diff_dict.keys():
    pop_name= pop_names_dict[pop]
    pop_col= pop_colors[pop_name]
    x= counts_dict[pop]['sizes']
    y= [np.log(x) for x in average_cdf[pop]]
    
    plt.plot(x,y,label= pop_name,color= pop_col)

plt.plot([0,sampling[0]],[np.log(threshold)]*2,label= str(threshold))

plt.xlabel(xlab,fontsize=20)

plt.ylabel(ylab,fontsize= 20)
plt.ylim(np.log(1e-4),0.1)
plt.title('grid SSD / sample proportion - control')

plt.legend()
plt.savefig(fig_dir + sims_target + '_' + 'Variance_Convergence_sbsD.png',bbox_inches='tight')
plt.close()

