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
from tools.mcounter_tools import mcounter_deploy
bases= 'ACGT'
ksize= 3

p_value= 1e-5
test_m= 'fisher'
individually= False
exclude= False
frequency_range= [0,1]
extract= 'pval'

pop_asso, count_data, ref_batch_dict, reference_freqs= mcounter_deploy_v2(data,p_value= p_value, test_m= test_m, individually= individually,
                                        exclude= exclude, frequency_range= frequency_range, extract= extract,
                                     muted_dir= muted_dir, sims_dir= sims_dir,data_freqs= data_freqs,collapsed= collapsed)


############## mutation grids
############## mutation grids

from functools import reduce  # forward compatibility for Python 3
import operator

from tools.fasta_utilities import (
    get_mutations, kmer_comp_index, kmer_mut_index
)


mutations= get_mutations(bases= bases,ksize= ksize)
kmers, kmer_idx= kmer_comp_index(mutations)

mut_lib= kmer_mut_index(mutations)
if collapsed:
    labels= [kmer_idx[x][0] for x in sorted(kmer_idx.keys())]
else:
    labels= ['_'.join(x) for x in mutations]

grid_labels= np.array(labels).reshape(row,4)
list_labels= grid_labels.reshape(1,np.prod(grid_labels.shape))[0]
############## process grids


available= list(count_data.keys())
subsamp= len(count_data)
avail_sub= np.random.choice(available,subsamp,replace= False)

grid_diffs= [count_data[s]['diffs'] for s in avail_sub]

### 2. calculate proportions across smulations
pop_proportions= [count_data[s]['sizes'][1] for s in avail_sub]
#pop_proportions= [round(x,3) for x in pop_proportions]
### 3. batch names
batch_names= [count_data[s]['batch'] for s in avail_sub]
batch_dict= {
    z:[x for x in range(len(avail_sub)) if batch_names[x] == z] for z in list(set(batch_names))
}



####################################################
#################################################### grid SSD

metric= 'euc_norm'
Nbins= 30
bins= np.linspace(0,1,Nbins)
bins= np.round(bins,4)
bins= [(bins[x-1],bins[x]) for x in range(1,len(bins))]

other_diffs= [count_data[s]['other'] for s in avail_sub]
pop_vector= [count_data[s]['pop'] for s in avail_sub]
pop_set= list(set(pop_vector))

#kl_diffs=[count_data[s]['kl'] for s in avail_sub]

pop_batch_dict= {
    ba: {
        pop: [x for x in batch_dict[ba] if pop_vector[x] == pop] for pop in pop_set
    } for ba in batch_dict.keys()
}

xlab= 'relative sampling'
ylab= 'mean matrix p-val'

view_sets= ['ref','anti']
grid_whole= {
   pop: {
        ba : {    
            view:{} for view in view_sets
        } for ba in batch_dict.keys()
    } for pop in pop_set
}

for i in batch_dict.keys():
    for pop in pop_batch_dict[i].keys():
        
        xprep= [pop_proportions[x] for x in pop_batch_dict[i][pop]]
        
        xprep= {
             z: [x for x in range(len(xprep)) if xprep[x] == z] for z in list(set(xprep))
        }

        ### grids
        batch_grids= [grid_diffs[x] for x in pop_batch_dict[i][pop]]
        y_prep= {
            z: [batch_grids[x] for x in xprep[z]] for z in xprep.keys()
        }

        y_prep= {
            z: [1 - np.sqrt(np.sum(x**2)) for x in y_prep[z]] for z in y_prep.keys()
        }

        surface= sorted(xprep.keys())
        y= [np.mean(y_prep[x]) for x in surface]
        error= [np.std(y_prep[x])*2 for x in surface]

        grid_whole[pop][i]['ref']= [surface,y,error]

        ###:
        
        batch_grids= [other_diffs[x] for x in pop_batch_dict[i][pop]]
        xprep= [pop_proportions[x] for x in pop_batch_dict[i][pop]]
        xprep= np.repeat(xprep,[len(x) for x in batch_grids])
        batch_grids= list(it.chain(*batch_grids))

        xprep= {
             z: [x for x in range(len(xprep)) if xprep[x] == z] for z in list(set(xprep))
        }
        y_prep= {
            z: [batch_grids[x] for x in xprep[z]] for z in xprep.keys()
        }

        y_prep= {
            z: [1-np.sqrt(np.sum(x**2)) for x in y_prep[z]] for z in y_prep.keys()
        }

        surface= sorted(xprep.keys())
        y= [np.mean(y_prep[x]) for x in surface]
        error= [np.std(y_prep[x])*2 for x in surface]
        
        grid_whole[pop][i]['anti']= [surface,y,error]


#############################################################
############################################################# STRATA 
###########

def set_SSD(set1,set2,same= False):
    '''
    return sum of squared differences between every pair of vectors across two sets.
    '''
    dists= []
    
    for ind in set1:
        
        dist_vec= [(x - ind) for x in set2] #/ np.sum(indian + x)
        dist_vec= [z**2 for z in dist_vec]
        dist_vec= [1 - np.sqrt(np.sum(x)) for x in dist_vec]
        dists.extend(dist_vec)
    
    if same:
        dists= np.array(dists)
        dists= dists.reshape(len(set1),len(set2))
        indx= np.triu_indices(len(set1),k=1)
        dists= dists[indx]
    
    return dists

    
###########
########### Pre-process reference population distances among themselves and across.

from itertools import combinations 

batch_data= {
    ba: {
        'pop': [x[0] for x in ref_batch_dict[ba]],
        'ref': [x[1] for x in ref_batch_dict[ba]]
    } for ba in ref_batch_dict.keys()
}

pop_ref_dict= {
    ba: {
        pop: [z for z in ref_batch_dict[ba] if z[0] == pop] for pop in list(batch_data[ba]['pop'])
    } for ba in ref_batch_dict.keys()
}


pop_ref_dict= {
    ba: {
        pop: {
                ref: [z[2] for z in pop_ref_dict[ba][pop] if z[1] == ref][0] for ref in list(set(batch_data[ba]['ref']))
             } for pop in list(batch_data[ba]['pop'])
    } for ba in ref_batch_dict.keys()
}


pop_anti_dict= {
    ba: {
        pop: {
            popi: ([pop_ref_dict[ba][pop][g[0]]],[pop_ref_dict[ba][popi][g[1]]]) for g in combinations(list(pop_ref_dict[ba][pop].keys()),2) for popi in pop_ref_dict[ba].keys() if pop != popi
        } for pop in pop_ref_dict[ba].keys()
        ##
    } for ba in pop_ref_dict.keys()
}



pop_ref_dict= {
    ba: {
        pop: [pop_ref_dict[ba][pop][ref] for ref in pop_ref_dict[ba][pop].keys()] for pop in pop_ref_dict[ba].keys()
        #
    } for ba in pop_ref_dict.keys()
}


pop_ref_dict= {
    ba: {
        pop: set_SSD(g,g,same= True) for pop,g in pop_ref_dict[ba].items()
        #
    } for ba in pop_ref_dict.keys()
}


###########

xlab= 'relative sampling'
ylab= 'sqrt SSD'

for pop in grid_whole.keys():
        for batch in grid_whole[pop].keys():
            for d in range(2):
                plt.figure(figsize=(10,10))
                for ep in ['ref']:
                    surface= grid_whole[pop][batch][ep][0]
                    y= grid_whole[pop][batch][ep][1]
                    error= grid_whole[pop][batch][ep][2]

                    plt.errorbar(surface,y,yerr=error,label= ep)    
            
                plt.plot([1,sampling[0]],[np.mean(pop_ref_dict[batch][pop]) - np.std(pop_ref_dict[batch][pop])]*2,label= 'up ref. diff')
                plt.plot([1,sampling[0]],[np.mean(pop_ref_dict[batch][pop])]*2,label= 'mean ref. diff')
                plt.plot([1,sampling[0]],[np.mean(pop_ref_dict[batch][pop]) + np.std(pop_ref_dict[batch][pop])]*2,label= 'lb ref. diff')

                if d>0:
                    for refp in pop_anti_dict[batch][pop].keys():
                        sts= pop_anti_dict[batch][pop][refp]
                        print(sts)
                        plt.plot([1,sampling[0]],[np.mean(sts)]*2,label= 'mean anti. diff {}'.format(refp))
                        plt.plot([1,sampling[0]],[np.mean(sts) - np.std(sts)*2]*2,label= 'lb anti. diff {}'.format(refp))

                plt.xlabel(xlab,fontsize= 20)
                plt.ylabel(ylab,fontsize=20)
                plt.xticks(fontsize= 15)
                plt.yticks(fontsize= 15)
                plt.title('grid SSD / sample proportion - control')

                plt.legend()
                plt.savefig(fig_dir + 'gridSSD_{}_control{}.png'.format(pop,['','_comparison'][d]),bbox_inches='tight')
                plt.close()


############################################################# 

for pop in grid_whole.keys():
        for batch in grid_whole[pop].keys():
            plt.figure(figsize=(10,10))
            for ep in ['ref','anti']:
                surface= grid_whole[pop][batch][ep][0]
                y= grid_whole[pop][batch][ep][1]
                error= grid_whole[pop][batch][ep][2]

                plt.errorbar(surface,y,yerr=error,label= ep)    
            
            plt.xlabel(xlab,fontsize= 20)
            plt.ylabel(ylab,fontsize=20)
            plt.xticks(fontsize= 15)
            plt.yticks(fontsize= 15)

            plt.title('grid SSD / sample proportion - control')

            plt.legend()
            plt.savefig(fig_dir + 'gridSSD_{}_control_sbsD.png'.format(pop),bbox_inches='tight')
            plt.close()

