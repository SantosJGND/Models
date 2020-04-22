from tools.input_utilities import (
    MC_sample_matrix_v1
)

from tools.mcounter_tools import (
    mcounter_deploy_v2
)

from tools.fasta_utilities import (
    get_mutations, kmer_comp_index, kmer_mut_index
    )

from tools.compare_utilities import (
    gzip_request
)


#from tools.SLiM_pipe_tools import mutation_counter_launch
import re
import pandas as pd
import os
import numpy as np
import itertools as it
import collections
from scipy.stats import norm


def recursively_default_dict():
    return collections.defaultdict(recursively_default_dict)



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


###############
############### plots

import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt

from INFO_db import INFO_dict



########
########
fig_dir= 'Figures/kmers'
os.makedirs(fig_dir, exist_ok=True)
fig_dir= fig_dir + '/'


## directories
main_dir= os.getcwd() + '/'
count_dir= main_dir + 'mutation_counter/count/'
dir_launch= main_dir + 'mutation_counter'
muted_dir= main_dir + 'mutation_counter/data/mutation_count/'
sims_dir= INFO_dict[args.species]['dirs']['sims']
sims_target= sims_dir.split('/')[-2]
mutlog= 'toMut.log'

indfile= INFO_dict[args.species]['w_indfile']['sims']


######################################
######################################
##### READ DATA


sampling= [args.samp,args.steps,args.reps]

row= [64,32][int(args.collapsed)]
col= 3
ksize= 3
bases= 'ACGT'

request_processed= gzip_request(dir_check= sims_dir,str_format= '',requested= ['.vcf','fa'], func= 'gzip')


data, data_freqs = MC_sample_matrix_v1(min_size= args.minS, samp= sampling, stepup= args.stepup, count_dir= count_dir, 
                        dir_launch= dir_launch,main_dir= main_dir,sim_dir= sims_dir,indfile= indfile, ploidy= args.ploidy,
                          muted_dir= muted_dir, diffs= args.diffs,collapsed= args.collapsed,row= row, col= col,
                       sample_sim= args.simsN,freq_extract= args.freqs,haps_extract=args.haps, ksize= ksize, bases= bases)


#######################################
#######################################
###### Count and compare. 

p_value= 1e-5
test_m= 'chi2'
exclude= False
frequency_range= [0,1]
extract= 'pval'

pop_asso, count_data, ref_batch_dict, reference_freqs= mcounter_deploy_v2(data,p_value= p_value, test_m= test_m, bases= bases, ksize= ksize,
                                        frequency_range= frequency_range, extract= extract, muted_dir= muted_dir, row= row, col= col,
                                     sims_dir= sims_dir,data_freqs= data_freqs,collapsed= args.collapsed)


############## mutation grids

from functools import reduce  # forward compatibility for Python 3
import operator

from tools.fasta_utilities import (
    get_mutations, kmer_comp_index, kmer_mut_index
)

mutations= get_mutations(bases= bases,ksize= ksize)
kmers, kmer_idx= kmer_comp_index(mutations)

mut_lib= kmer_mut_index(mutations)
if args.collapsed:
    labels= [kmer_idx[x][0] for x in sorted(kmer_idx.keys())]

else:
    labels= ['_'.join(x) for x in mutations]

grid_labels= np.array(labels).reshape(row,col)
list_labels= grid_labels.reshape(1,np.prod(grid_labels.shape))[0]
##############
############## process grids

available= list(count_data.keys())
subsamp= len(count_data)
avail_sub= np.random.choice(available,subsamp,replace= False)

### 1. extract grids
grids= [count_data[s]['grids'] for s in avail_sub]
grid_shape= grids[0].shape
grid_total= np.prod(grid_shape)

grid_props= [count_data[s]['props'] for s in avail_sub]
grid_diffs= [count_data[s]['diffs'] for s in avail_sub]
## extract statistics per mutation.
mut_grid= {}
mut_diffs= {}
mut_props= {}

for r in range(grid_shape[0]):
    for c in range(grid_shape[1]):
        mut= grid_labels[r,c]
        mut_grid[mut]= []
        mut_diffs[mut]= []
        mut_props[mut]= []

        for idx in range(len(avail_sub)):
            mut_grid[mut].append(grids[idx][r,c])
            mut_diffs[mut].append(grid_diffs[idx][r,c])
            mut_props[mut].append(grid_props[idx][r,c])

###
### 2. calculate proportions across smulations

pop_sizes= [count_data[s]['sizes'][1]  for s in avail_sub]

pop_proportions= pop_sizes
#pop_proportions= [count_data[s]['sizes'][1] / count_data[s]['sizes'][0] for s in avail_sub]

pop_proportions= [round(x,3) for x in pop_proportions]
### 3. batch names
batch_names= [count_data[s]['pop'] for s in avail_sub]

batch_dict= {
    z:[x for x in range(len(avail_sub)) if batch_names[x] == z] for z in list(set(batch_names))
}
##
pop_vector= [count_data[s]['pop'] for s in avail_sub]
batch_pop_dict= {
    z: {
        p: [x for x in batch_dict[z] if pop_vector[x] == p] for p in list(set(pop_vector))
    } for z in batch_dict.keys()
}



################################################
################################################
######## PLOT

pop_colors= INFO_dict[args.species]['pop_colors']
pop_names_dict= {g:z for z,g in INFO_dict[args.species]['pop_dict'].items()}

##
sampling_str= str(sampling[0])
## plot first for every mutation context.
for kmer in mut_grid.keys():
    fig_kmer= fig_dir + '/' + kmer
    os.makedirs(fig_kmer, exist_ok=True)
    fig_kmer= fig_kmer + '/'
    
    plot_data= {
        'pval': mut_grid[kmer],
        'diffs': mut_diffs[kmer],
        'props': mut_props[kmer]
    }
    
    for strata in plot_data.keys():
        batch_hold= {}
        ydata= plot_data[strata]
        
        xlab= 'relative sampling'
        ylab= 'mean matrix {}'.format(strata)

        d=0

        colors= ['ro', 'bs', 'g^']

        append_all= []
        for i in batch_dict.keys():
            
            xprep= [pop_proportions[x] for x in batch_dict[i]]
            #print(xprep)
            xprep= {
                z: [x for x in range(len(xprep)) if xprep[x] == z] for z in list(set(xprep))
            }
            
            yprep= [ydata[x] for x in batch_dict[i]]
            yprep= {
                z: [yprep[x] for x in xprep[z]] for z in xprep.keys()
            }
            
            x= sorted(xprep.keys())
            y= [np.nanmean(yprep[z]) for z in x]
            error= [np.nanstd(yprep[z]) for z in x]

            batch_hold[i]= {
                'x': x,
                'y': y,
                'error': error
            }

            plt.figure(figsize=(10, 10))

            plt.errorbar(x,y,yerr=error)
            plt.xlabel(xlab + ' {} comparisons'.format(len(batch_dict[i])),fontsize=20)
            #plt.ylim(0,1.5)
            plt.ylabel(ylab,fontsize=20)
            plt.xticks(fontsize= 15)
            plt.yticks(fontsize= 15)
            plt.title(i,fontsize= 20)

            plt.savefig(fig_kmer + sims_target + '_' + '{}_{}_{}.png'.format(kmer,i, strata),bbox_inches='tight')
            plt.close()
            
            d += 1

        plt.figure(figsize=(10, 10))

        for i in batch_hold.keys():
            
            plt.errorbar(batch_hold[i]['x'],batch_hold[i]['y'],yerr=batch_hold[i]['error'],label= i)

        plt.xlabel(xlab,fontsize=20)
        #plt.ylim(0,1.5)
        plt.ylabel(ylab,fontsize=20)
        plt.xticks(fontsize= 15)
        plt.yticks(fontsize= 15)
        plt.title('combined stats')
        plt.legend()

        plt.savefig(fig_kmer + sims_target + '_' + '{}_{}_combined_{}_{}.png'.format(args.stepup,sampling_str,kmer,strata),bbox_inches='tight')
        plt.close()



########################################## GEN

xlab= 'relative sampling'
ylab_temp= 'mean matrix {}'

keys_get= ['pval','props','diffs']

compound_kmer= {
        y: {
                g: {
                    'mean':[],
                    'std':[]
                } for g in batch_dict.keys()
            } for y in ['pval','diffs','props']
        }



plot_data= {
    batch:{
        kmer: {
                'pval': mut_grid[kmer],
                'props': mut_props[kmer],
                'diffs': mut_diffs[kmer]
        } for kmer in mut_diffs.keys() if 'CG' not in kmer
    } for batch in batch_dict.keys()
}



for i in batch_dict.keys(): 
    
    xprep= [pop_sizes[x] for x in batch_dict[i]]
    #print(xprep)
    xprep= {
        z: [x for x in range(len(xprep)) if xprep[x] == z] for z in list(set(xprep))
    }
    
    for strata in keys_get:
        for plt_type in ['lines','markers']:
            plt.figure(figsize=(10, 10))
            ylab= ylab_temp.format(strata)
            for kmer in plot_data[i].keys():
                
                ydata= plot_data[i][kmer][strata]
                

                d=0

                colors= ['ro', 'bs', 'g^']

                append_all= []


                yprep= [ydata[x] for x in batch_dict[i]]
                yprep= {
                    z: [yprep[x] for x in xprep[z]] for z in xprep.keys()
                }
                
                if plt_type == 'lines':
                    xvec= sorted(xprep.keys())
                    ydat= yprep[z]
                    ydat= [x**2 for x in ydat]
                    y= [np.nanmean(ydat) for z in xvec]
                    error= [np.nanstd(ydat) for z in xvec]
                    
                    compound_kmer[strata][i]['mean'].append(y)
                    compound_kmer[strata][i]['std'].append(error)
                    
                    plt.errorbar(xvec,y,yerr=error)

                else:
                    xvec= sorted(xprep.keys())
                    y= [yprep[z] for z in xvec]

                    xvec= [[xvec[x]] * len(y[x]) for x in range(len(xvec))]
                    xvec= list(it.chain(*xvec))
                    y= list(it.chain(*y))

                    plt.plot(xvec,y,'.')



            plt.xlabel(xlab,fontsize=20)
            #plt.ylim(0,1.5)
            plt.ylabel(ylab,fontsize=20)
            plt.xticks(fontsize= 15)
            plt.yticks(fontsize= 15)
            plt.title('combined stats')

            plt.savefig(fig_dir + sims_target + '_' + '{}_{}_combined_{}_{}_{}.png'.format(args.stepup,sampling_str,plt_type,i,strata),bbox_inches='tight')
            plt.close()



######################
######################



for strata in ['pval','diffs','props']:

    plt.figure(figsize=(10,10))
    ylab= ylab_temp.format(strata)
    for batch in batch_dict.keys():
        global_x= sorted(set([pop_proportions[x] for x in batch_dict[batch]]))
        
        global_y= compound_kmer[strata][batch]['mean']
        global_y= np.array(global_y)
        global_y= np.ma.masked_where(global_y == np.inf, global_y)
        global_y= list(np.mean(global_y,axis= 0))
        
        global_error= compound_kmer[strata][batch]['std']
        global_error= np.array(global_error)
        global_error= np.ma.masked_where(global_error == np.inf, global_error)
        global_error= list(np.std(global_error,axis= 0))
        
        plt.errorbar(global_x,global_y,yerr=global_error,label= batch)

    plt.xlabel(xlab,fontsize=0)
    #plt.ylim(0,1.5)
    plt.ylabel(ylab,fontsize=20)
    plt.xticks(fontsize= 15)
    plt.yticks(fontsize= 15)
    plt.title('combined stats',fontsize=20)

    plt.legend()
    plt.savefig(fig_dir + sims_target + '_' + '{}_{}_combined_{}_{}.png'.format(args.stepup,sampling_str,'kmers',strata),bbox_inches='tight')
    plt.close()





for strata in ['pval','diffs','props']:


    plt.figure(figsize=(10,10))
    
    for batch in batch_dict.keys():
        global_x= sorted(set([pop_proportions[x] for x in batch_dict[batch]]))
        
        global_error= compound_kmer[strata][batch]['std']
        global_error= np.array(global_error)
        global_error= np.ma.masked_where(global_error == np.inf, global_error)
        global_error= list(np.nanmean(global_error,axis= 0))
        
        plt.plot(global_x,global_error,label= batch)
    
    plt.xlabel(xlab,fontsize=20)
    #plt.ylim(0,1.5)
    plt.ylabel('variance',fontsize=20)
    plt.xticks(fontsize= 15)
    plt.yticks(fontsize= 15)
    plt.title('combined stats',fontsize=20)
    
    plt.legend()
    plt.savefig(fig_dir + sims_target + '_' + '{}_{}_combined_{}_{}.png'.format(args.stepup,sampling_str,'variance',strata),bbox_inches='tight')
    plt.close()


############
############ p-value of 0 averaged across mutation types.


strata= 'diffs'

plt.figure(figsize=(10,10))

for batch in batch_dict.keys():
    global_x= sorted(set([pop_proportions[x] for x in batch_dict[batch]]))
    global_means= compound_kmer[strata][batch]['mean']
    global_means= np.array(global_means)
    print(global_means.shape)
    #global_means= np.ma.masked_where(global_means == np.inf, global_means)
    #global_means= list(np.nanmean(global_means,axis= 0))
    
    global_error= compound_kmer[strata][batch]['std']
    global_error= np.array(global_error)
    print(global_error.shape)
    #global_error= np.ma.masked_where(global_error == np.inf, global_error)
    #global_error= list(np.nanmean(global_error,axis= 0))
    
    print(len(global_error))
    print([len(x) for x in global_error])
    print(global_means.shape)
    #print(global_means[:20])
    #print(global_error[:20])
    pval_array= []
    global_sd= []
    surface= []
    for idx in range(len(global_x)):
        pval= []

        for rep in range(len(global_means)):
            #print(global_means[rep,idx])
            #print(global_error[rep,idx])
            pval_i= norm.cdf(0,loc= global_means[rep,idx],scale= global_error[rep,idx])
            pval.append(pval_i)

        #pval= [np.log(x) for x in pval]
        if pval:
            pval_array.append(np.nanmean(pval))
            global_sd.append(np.nanstd(pval))
            surface.append(global_x[idx])

    dol= [x for x in range(len(pval_array))]
    #
    print(pval_array[:20])
    #pval_array= [np.log(x) for x in pval_array]

    global_x= [global_x[x] for x in dol]
    
    plt.errorbar(surface,pval_array,yerr= global_sd,label= batch)

#plt.plot([0,sampling[0]],[np.log(threshold)]*2,label= str(threshold))

plt.xlabel(xlab,fontsize=20)
#plt.ylim(0,1.5)
plt.ylabel('variance',fontsize=20)
plt.xticks(fontsize= 15)
plt.yticks(fontsize= 15)

#plt.ylim(np.log(1e-4),0.1)

plt.title('True Value average p-val',fontsize=20)

plt.legend()
plt.savefig(fig_dir + sims_target + '_' + '{}_{}_combined_{}_{}.png'.format(args.stepup,sampling_str,'pvalue',strata),bbox_inches='tight')
plt.close()

