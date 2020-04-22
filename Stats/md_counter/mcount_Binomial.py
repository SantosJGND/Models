from tools.input_utilities import (
    MC_sample_matrix_v1
)

from tools.mcounter_tools_new import (
    seg_comp_v2, get_chrom_sizes
)

from tools.fasta_utilities import (
    get_mutations, kmer_comp_index, kmer_mut_index
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



##########
########## directories
chrom_size_dir= '/home/jgarc235/Rhesus/chrom_sizes/'

main_dir= os.getcwd() + '/'
count_dir= main_dir + 'mutation_counter/count/'
dir_launch= main_dir + 'mutation_counter'
muted_dir= main_dir + 'mutation_counter/data/mutation_count/'

fig_dir= 'Figures/Binomial'
os.makedirs(fig_dir, exist_ok=True)
fig_dir= fig_dir + '/'

single_figs= True

#########
### species data
from INFO_db import INFO_dict


species= 'chimp'

##
assembly= INFO_dict[species]['assembly']
chrom_sizes= get_chrom_sizes(assembly,chrom_size_dir= chrom_size_dir)
genome_size= sum(chrom_sizes.values())

prop_gen_used= .85
scale_genSize= False
##

import re
import pandas as pd

sims_dir= INFO_dict[species]['dirs']['sims']
sims_target= sims_dir.split('/')[-2]
diffs= False

mutlog= 'toMut.log'
min_size= 5
sampling= [100,100,5]
stepup= 'increment'
sample_sim= 3
collapsed= False

row= [64,32][int(collapsed)]
col= 3
freq_extract= False
segregating= True
sim_del= 'C'
chrom_idx= 0
request_processed= gzip_request(dir_check= sims_dir,str_format= '',requested= ['.vcf','fa'], func= 'gzip')

ploidy= 1
haps_extract= True

data, data_freqs = MC_sample_matrix_v1(min_size= min_size, samp= sampling, stepup= stepup, count_dir= count_dir, 
                        dir_launch= dir_launch,main_dir= main_dir,sim_dir= sims_dir,sim_del= sim_del,chrom_idx= chrom_idx,
                          muted_dir= muted_dir, diffs= diffs,collapsed= collapsed,row= row, col= col, scale_genSize= scale_genSize, prop_gen_used= prop_gen_used,
                       exclude= False,sample_sim= sample_sim,freq_extract= freq_extract,segregating= segregating,haps_extract= haps_extract,
                       ploidy= ploidy)


#################
#################

p_value= 1e-5
test_m= 'chi2'
individually= False
exclude= False
frequency_range= [0,1]
extract= 'pval'
Nbins= sampling[1]
chi_total= True


pop_asso, ref_combos= seg_comp_v2(data,p_value= p_value, test_m= test_m, individually= individually, Nbins= Nbins,
                                        exclude= exclude, frequency_range= frequency_range, extract= extract,stepup= stepup, sampling= sampling,
                                     muted_dir= muted_dir, data_freqs= data_freqs,row= row, col= col,chi_total= chi_total)



############
############

Nmax= 1
if stepup == 'increment':
    Nmax= sampling[0]
    Nbins= sampling[0]

bins= np.linspace(0,Nmax,Nbins)
bins= np.round(bins,4)
bins= [(bins[x-1],bins[x]) for x in range(1,len(bins))]


sim_stats= list(ref_combos.keys())
combs_poss= list(ref_combos[sim_stats[0]]['combs'].keys())
bins_poss= ref_combos[sim_stats[0]]['combs'][combs_poss[0]].keys()
bins_poss= list(bins_poss)

combo_grids= {
    comb: {
        z: {
            'PA': {
                v: [] for v in comb
            },
            'shared': {
                v: [] for v in comb
            }
        } for z in bins_poss
    } for comb in combs_poss
}

combo_ref= {x: {
    'PA': {
        v: [] for v in x
    },
    'shared': []
} for x in ref_combos[sim_stats[0]]['combs'].keys()}

for sim in sim_stats:
    for comb in ref_combos[sim]['combs'].keys():
        for pop in comb:
            tplpop= (sim,pop)
            combo_ref[comb]['PA'][pop].append(ref_combos[sim]['stats'][comb]['PA'][tplpop])
        
        combo_ref[comb]['shared'].append(ref_combos[sim]['stats'][comb]['shared'])

        for bi in ref_combos[sim]['combs'][comb].keys():
            for pop in comb:
                pop_sum= ref_combos[sim]['combs'][comb][bi]['PA'][pop]
                shared= ref_combos[sim]['combs'][comb][bi]['shared'][pop]
                
                combo_grids[comb][bi]['PA'][pop].extend(pop_sum)
                combo_grids[comb][bi]['shared'][pop].extend(shared)


####################
#################### Actual analysis.
## min needed file
from scipy.stats import binom

alpha= 0.01
window_pval= [5e-3,0.01]
p= 1 / 192

### rate_change
rate_steps= 100
#range_surface= np.linspace(1,20,rate_steps)
#rate_range= 1 + np.log(range_surface) * 1e-1

rate_range= np.linspace(1,10,5)

for comb_select in combo_grids.keys():


    rate_stats= recursively_default_dict()

    for rate in rate_range:
        new_rate= p * rate
        rate_stats[rate]= {
            x: [] for x in [*comb_select,'pval']
        }
        print(len(combo_grids[comb_select].keys()))

        for bi in combo_grids[comb_select].keys():

            pops_probs= []
            for pop_skew in comb_select:

                pop_ref= [x for x in comb_select if x != pop_skew][0]

                Nrep= len(combo_grids[comb_select][bi]['shared'])
                pops= list(combo_grids[comb_select][bi]['PA'].keys())

                Nt_rep= combo_grids[comb_select][bi]['shared'][pop_ref]
                Nt_skew= combo_grids[comb_select][bi]['shared'][pop_skew]
                print(Nt_rep,Nt_skew)
                #Nt= np.mean(Nt)
                Nref= combo_grids[comb_select][bi]['PA'][pop_ref]
                #Nref= np.mean(Nref)
                Nskew= combo_grids[comb_select][bi]['PA'][pop_skew]
                
                #print(Nt,Nref,Nskew)
                lbt= binom.ppf(alpha,Nt_skew,p)
                lbt= lbt + binom.ppf(alpha,Nskew,new_rate)

                probH2= 1 - binom.cdf(lbt,[Nref[x] + Nt_rep[x] for x in range(len(Nt_rep))], p)
                #print(probH2)
                pops_probs.append(probH2)

            pops_probs= np.array(pops_probs)
            #print(pops_probs.shape)
            pops_probs= np.nanmean(pops_probs,axis= 1)

            if len(probH2):
                for idx in range(len(comb_select)):
                    rate_stats[rate][comb_select[idx]].append(sum(bi[idx]) / 2)

                rate_stats[rate]['pval'].append(np.max(pops_probs))

    rate_pval_wind= {
        rate: [x for x in range(len(rate_stats[rate]['pval'])) if rate_stats[rate]['pval'][x] <= window_pval[1]] for rate in rate_stats.keys()
    }

    if single_figs:
        for i in sorted(rate_pval_wind.keys()):
            plt.figure(figsize=(20,10))
            surface= [rate_stats[i][comb_select[0]][x] for x in rate_pval_wind[i]]
            y= [rate_stats[i][comb_select[1]][x] for x in rate_pval_wind[i]]

            plt.scatter(surface,y,label= 'rate = ' + str(i))    

            plt.plot([0,sampling[0]],[np.log(alpha)]*2,label= alpha)

            plt.xlabel('Nsamp ' + comb_select[0])
            plt.ylabel('Nsamp ' + comb_select[1])

            #plt.title('grid SSD / sample proportion - control')
            
            plt.title('pval range: {}'.format('-'.join([str(x) for x in window_pval])))

            plt.xticks(fontsize=20)
            plt.yticks(fontsize=20)

            plt.legend()
            plt.savefig(fig_dir + sims_target + '_Binomial_{}_skewed.png'.format('-'.join(comb_select)),bbox_inches='tight')
            plt.close()

    else:
        plt.figure(figsize=(20,10))

    for i in sorted(rate_pval_wind.keys()):
        surface= [rate_stats[i][comb_select[0]][x] for x in rate_pval_wind[i]]
        y= [rate_stats[i][comb_select[1]][x] for x in rate_pval_wind[i]]

        plt.scatter(surface,y,label= 'rate = ' + str(i))    

    plt.plot([0,sampling[0]],[np.log(alpha)]*2,label= alpha)

    plt.xlabel('Nsamp ' + comb_select[0])
    plt.ylabel('Nsamp ' + comb_select[1])

    #plt.title('grid SSD / sample proportion - control')
    
    plt.title('pval range: {}; gen_size: {}'.format('-'.join([str(x) for x in window_pval]),genome_size))

    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)

    plt.legend()
    plt.savefig(fig_dir + sims_target + '_Binomial_{}_skewed.png'.format('-'.join(comb_select)),bbox_inches='tight')
    plt.close()

