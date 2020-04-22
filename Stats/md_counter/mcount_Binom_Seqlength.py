from tools.input_utilities import (
    MC_sample_matrix_simple, MC_sample_matrix_dict
)

from tools.compare_utilities import (
    gzip_request
)

from tools.fasta_utilities import (
    get_mutations, kmer_comp_index, kmer_mut_index
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


def get_pop_dict(indfile= 'ind_assignments.txt'):

    #ind_assignments= dir_sim + reference + '/' + indfile
    
    with open(indfile,'r') as f:
        inds= f.readlines()
    
    inds= [x.split() for x in inds]
    inds= [x for x in inds if x]
    pops= np.array(inds)[:,1]
    pop_dict= {
        z: [x for x in range(len(pops)) if pops[x] == z] for z in list(set(pops))
    }

    return pop_dict



def get_chrom_sizes(assembly= 'hg38', chrom_size_dir= '/home/jgarc235/Rhesus/chrom_sizes/', avoid= ['_','X','Y']):
    '''
    return dictionary of chrom sizes for given assembly
    '''

    filename= chrom_size_dir + assembly + '.chrom.sizes'

    with open(filename,'r') as fp: 
        lines= fp.readlines()
    
    lines= [x.split() for x in lines]
    lines= [x for x in lines if sum([int(y in x[0]) for y in avoid]) == 0]

    chrom_dict= {
        x[0].strip('chr'): int(x[1]) for x in lines
    }

    return chrom_dict


def compare_binomial(ref_combos_pval,comb_select, rate,alpha= 0.01, p= 1/96):

    for pop_skew in comb_select:
        pops_probs= []

        pop_ref= [x for x in comb_select if x != pop_skew][0]

        Nrep= len(ref_combos_pval[comb_select]['shared'][pop_ref])
        pops= list(ref_combos_pval[comb_select]['PA'].keys())
        
        Nt_rep= ref_combos_pval[comb_select]['shared'][pop_ref]
        Nt_skew= ref_combos_pval[comb_select]['shared'][pop_skew]
        #Nt= np.mean(Nt)
        Nref= ref_combos_pval[comb_select]['PA'][pop_ref]
        #Nref= np.mean(Nref)
        Nskew= ref_combos_pval[comb_select]['PA'][pop_skew]
        #Nskew= np.mean(Nskew)

        #print(Nt,Nref,Nskew)
        lbt= binom.ppf(alpha,Nt_skew,p)
        lbt= lbt + binom.ppf(alpha,Nskew,new_rate)

        probH2= 1 - binom.cdf(lbt,[Nref[x] + Nt_rep[x] for x in range(len(Nt_rep))], p)
        #print(probH2)
        pops_probs.append(probH2)
    
    pops_probs= np.array(pops_probs)
    #print(pops_probs.shape)
    pops_probs= np.max(pops_probs,axis= 0)

    return pops_probs



############## data set info
from INFO_db import INFO_dict


##############
##############
chrom_size_dir= '/home/jgarc235/Rhesus/chrom_sizes/'
main_dir= os.getcwd() + '/'
count_dir= main_dir + 'mutation_counter/count/'
dir_launch= main_dir + 'mutation_counter'
muted_dir= main_dir + 'mutation_counter/data/mutation_count/'


fig_dir= 'Figures/Binom_Seq_length'
os.makedirs(fig_dir, exist_ok=True)
fig_dir= fig_dir + '/'

##############
species= 'rhesus'

##
assembly= INFO_dict[species]['assembly']
chrom_sizes= get_chrom_sizes(assembly,chrom_size_dir= chrom_size_dir)
genome_size= sum(chrom_sizes.values())


scale_genSize= False

##

ind_file= ind_file= INFO_dict[species]['g_indfile']

pops_dict= get_pop_dict(indfile= ind_file)
pop_seg_dict= INFO_dict[species]['pop_dict']


pop_seg_dict= {
    g:v for v,g in pop_seg_dict.items()
}

pop_sizes= {
    v: len(g) for v,g in pops_dict.items()
}

######
######
from tools.input_utilities import (
    MC_sample_matrix_dict
)

import re
import pandas as pd



diffs= False
indfile= INFO_dict[species]['w_indfile']['sims']
mutlog= 'toMut.log'
min_size= 5
sampling= [100,50,5]
stepup= 'increment'
sample_sim= 5
collapsed= False

row= [64,32][int(collapsed)]
col= 3
freq_extract= False
segregating= True
sim_del= 'C'
chrom_idx= 0
return_private= True

sims_dir= INFO_dict[species]['dirs']['sims']
sims_target= sims_dir.split('/')[-2]

request_processed= gzip_request(dir_check= sims_dir,str_format= '',requested= ['.vcf','fa'], func= 'gzip')


sample_sim= 6
Lsteps= 50
samp= 5

data_sims, data_freqs_sims = MC_sample_matrix_dict(pop_seg_dict,pop_sizes,min_size= min_size, samp= samp, stepup= stepup, count_dir= count_dir, 
                        dir_launch= dir_launch,main_dir= main_dir,sim_dir= sims_dir,indfile= indfile,
                          muted_dir= muted_dir, diffs= diffs,collapsed= collapsed,row= row, col= col, scale_genSize= scale_genSize,
                       exclude= False,sample_sim= sample_sim,freq_extract= freq_extract,segregating= segregating, Lsteps= Lsteps,
                       return_private= return_private)


available= list(data_sims.keys())
len_vector= [x.split('.')[-1] for x in available]
batch_vector= [x.split('_')[0] for x in available] 
len_dict= {
    z: [x for x in range(len(len_vector))  if len_vector[x] == z] for z in list(set(len_vector))
}


ref_len= len_dict['full']
len_dict.pop('full')
len_dict= {int(z):g for z,g in len_dict.items()}

pop_combs= data_sims[available[0]]['pairDist'].keys()


#######
####### REFERENCE ANALYSIS
tag= '_ss'
ref_sims= [available[x] for x in ref_len]

ref_pvals= {}
ref_combos_pval= {}

for sim_here in ref_sims:
    sim_simple= sim_here.split(tag)[0]
    ref_pvals[sim_simple]= {}
    
    for combo in pop_combs:
        ref_pvals[sim_simple][combo]= {'PA': {}}
        pop1,pop2= combo

        ref_pair= [(sim_here, pop1),(sim_here,pop2)]

        pop_count_dict= {
            z: [] for z in combo
        }

        pop_counts= {
            g: data_sims[g[0]]['seg'][g[1]] for g in ref_pair
        }

        for idx in range(len(combo)):
            g= ref_pair[idx]

            pop_count_dict[combo[idx]].append(pop_counts[g][0])
            #ij_queue[bend].append(d)
            #d += 1

        ##
        pop_count_dict= {
            z:np.array(g) for z,g in pop_count_dict.items()
        }

        pop_count_shared= pop_count_dict[pop1] - pop_count_dict[pop2]


        PA_dict= {
            combo[z]: np.array(pop_count_shared == [1,-1][int(z)],dtype= int) for z in [0,1]
        }
        
        shared= np.array(pop_count_shared == 0,dtype= int)

        if combo not in ref_combos_pval.keys():
            ref_combos_pval[combo]= {
                'PA': {z:[] for z in combo},
                'shared': []
            }
            

        ref_combos_pval[combo]['shared'].append(shared)
        ref_pvals[sim_simple][combo]['shared']= shared
        
        for pop in combo:
            ref_combos_pval[combo]['PA'][pop].append(PA_dict)
            ref_pvals[sim_simple][combo]['PA'][pop]= PA_dict[pop]


##########
##########

PA_shared_dict= {x: {
    z: {
        'PA': {
            v: [] for v in z
        },
        'shared': {
            v: [] for v in z
        }
    } for z in pop_combs
} for x in len_dict.keys()}


for ls in len_dict.keys():
    same_len= len_dict[ls]
    
    for combo in pop_combs:
        
        for idx in same_len:
            
            sim_here= available[idx]
            sim_master= sim_here.split(tag)[0]
            
            pop1,pop2= combo
                       
            ref_pair= [(sim_here, pop1),(sim_here,pop2)]

            pop_count_dict= {
                z: [] for z in combo
            }

            pop_counts= {
                g: data_sims[g[0]]['seg'][g[1]] for g in ref_pair
            }

            for idx in range(len(combo)):
                g= ref_pair[idx]

                pop_count_dict[combo[idx]].append(pop_counts[g][0])
                #ij_queue[bend].append(d)
                #d += 1

            ##
            pop_count_dict= {
                z:np.array(g) for z,g in pop_count_dict.items()
            }
            
            PA_dict= {
                z: g + ref_pvals[sim_master][combo]['PA'][z].T[:len(g.T)].T for z,g in pop_count_dict.items()
            }
            
            PA_dict= {
                z: np.array(g == 2,dtype= int) for z,g in PA_dict.items()
            }
            
            PA_dict= {z:np.sum(g) for z,g in PA_dict.items()}
            
            shared_dict= {
                z: g + ref_pvals[sim_master][combo]['shared'].T[:len(g.T)].T for z,g in pop_count_dict.items()
            }
            
            shared_dict= {
                z: np.array(g == 2,dtype= int) for z,g in shared_dict.items()
            }
            
            shared_dict= {z:np.sum(g) for z,g in shared_dict.items()}
            
            for pop in combo:
                PA_shared_dict[ls][combo]['shared'][pop].append(shared_dict[pop])
                PA_shared_dict[ls][combo]['PA'][pop].append(PA_dict[pop])



###########
###########
alpha= 0.01
p= 1 / 96
threshold= 0.01
rate_steps= 7
rate_max= 10

####################
####################
rate_steps= 100
range_surface= np.linspace(1,20,rate_steps)
rate_range= 1 + np.log(range_surface) * 1e-1

ref_rate_stats= {}

full_pop_out_file= 'full_ratePval.txt'
#####################
### Attributing pvals

from scipy.stats import binom

rate_steps= 100
range_surface= np.linspace(1,20,rate_steps)
rate_range= 1 + np.log(range_surface) * 1e-1


for comb_select in pop_combs:
    rate_range= np.linspace(2,rate_max,rate_steps)

    rate_stats= recursively_default_dict()
    ref_rate_stats[combo]= {}

    for rate in rate_range:
        new_rate= p * rate
        rate_stats[rate]= {
            x: [] for x in ['L','pval','std']
        }

        ## REFERENCES
        if rate not in ref_rate_stats.keys():
            ref_rate_stats[rate]= {}

        pops_probs= compare_binomial(ref_combos_pval,comb_select,rate,alpha= alpha,p= p)

        if len(pops_probs):
            
            ref_rate_stats[combo][rate]= {
                'mean': np.mean(pops_probs),
                'std': np.mean(pops_probs)
            }
        ##
        # SUBSAMPLES
        for ls in sorted(PA_shared_dict.keys()):
            
            pops_probs= compare_binomial(PA_shared_dict[ls],comb_select,rate,alpha= alpha,p= p)
            
            if len(pops_probs):
                
                rate_stats[rate]['L'].append(ls)
                
                rate_stats[rate]['pval'].append(np.mean(pops_probs))
                rate_stats[rate]['std'].append(np.std(pops_probs))

    ###########
    ###########
    thresh_dict= {
        rate: [x for x in range(len(rate_stats[rate]['pval'])) if abs(rate_stats[rate]['pval'][x] - threshold) - 2 * rate_stats[rate]['std'][x] < 0] for rate in rate_stats.keys()
    } 

    thresh_dict= {
        x: g for x,g in thresh_dict.items() if g
    }

    thresh_dict= {
        x: [rate_stats[rate]['L'][x] for x in g] for x,g in thresh_dict.items()
    }
    
    ###########
    ###########
    ###########

    for show_var in [0,1]:

        for zoomin in [0,1]:

            plt.figure(figsize=(20,10))

            for rate in rate_stats.keys():
                surface= rate_stats[rate]['L']
                y= rate_stats[rate]['pval']
                error= rate_stats[rate]['std']

                if show_var:
                    plt.errorbar(surface,y,yerr=error,label= 'rate: ' + str(round(rate,3)))
                else:

                    plt.plot(surface,y,label= 'rate = ' + str(round(rate,3)))    

            if zoomin:
                plt.ylim(0,threshold + 0.005)

            if zoomin:
                plt.plot([0,max(surface)],[threshold]*2,label= alpha)

            plt.xlabel('pval')
            plt.ylabel('sequence length (bp)')

            #plt.title('grid SSD / sample proportion - control')
            
            plt.xticks(fontsize=20)
            plt.yticks(fontsize=20)

            plt.legend()
            plt.savefig(fig_dir + sims_target + '_BinomSeqL_{}_{}_{}.png'.format('-'.join(comb_select),['','std'][int(show_var)],['','zoom'][int(zoomin)]),bbox_inches='tight')
            plt.close()


    #########
    #########
    #########

    plt.figure(figsize=(20,10))

    hists= [thresh_dict[rate] for rate in sorted(thresh_dict.keys())]
    hist_labels= sorted(thresh_dict.keys())
    hist_labels= [round(x,2) for x in hist_labels]
    plt.boxplot(hists,labels= hist_labels)

    plt.ylabel('Seq Length',fontsize= 20)
    plt.xlabel('rate', fontsize= 20)

    #plt.title('grid SSD / sample proportion - control')

    plt.title('-'.join(comb_select) + ' pval thresh: {}'.format(threshold))

    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)

    plt.legend()
    plt.savefig(fig_dir + sims_target + '_BinomSeqL_hists_{}.png'.format('-'.join(comb_select)),bbox_inches='tight')
    plt.close()

    
    ##########
    ##########
    ##########
    thresh_sum_dict= {}

    for thresh in np.linspace(0.0001,0.05,5):
        thresh_dict= {
            rate: [x for x in range(len(rate_stats[rate]['pval'])) if rate_stats[rate]['pval'][x] + rate_stats[rate]['std'][x] - threshold < 0] for rate in rate_stats.keys()
        } 

        thresh_dict= {
            x: g for x,g in thresh_dict.items() if g
        }

        thresh_dict= {
            rate: [rate_stats[rate]['L'][x] for x in g] for rate,g in thresh_dict.items()
        }

        thresh_dict= {
            x: g[0] for x,g in thresh_dict.items()
        }
        thresh_sum_dict[thresh]= thresh_dict


    plt.figure(figsize=(10,10))
    for thresh in thresh_sum_dict.keys():
        surface= sorted(thresh_dict.keys())
        y= [thresh_dict[rate] for rate in sorted(thresh_dict.keys())]

        plt.plot(surface,y,label= str(thresh))

    plt.ylabel('Seq Length',fontsize= 20)
    plt.xlabel('rate', fontsize= 20)

    plt.title('-'.join(comb_select) + 'seq length for pval < {}'.format(threshold))

    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)

    plt.legend()
    plt.savefig(fig_dir + sims_target + '_BinomSeqL_threshold_{}.png'.format('-'.join(comb_select)),bbox_inches='tight')
    plt.close()


with open(fig_dir + full_pop_out_file,'w') as fp:
    header= ['comb','mean','std']
    fp.write('\t'.join(header) + '\n')

    for cb in ref_rate_stats.keys():
        comb_name= '-'.join(cb)
        for rate in sorted(ref_rate_stats[cb].keys()):
            mean_stat= ref_rate_stats[cb][rate]['mean']
            sd_stat= ref_rate_stats[cb][rate]['std']
            fp.write('\t'.join([comb_name,mean_stat,sd_stat]) + '\n')

