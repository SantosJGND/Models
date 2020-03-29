
import os

from tools.SLiM_pipe_tools import (
    read_chrom_sizes, region_samplev2,
    fasta_RextractUnif, write_fastaEx, process_recipe,
    SLiM_dispenserv3, 
)



###########

if __name__ == '__main__':

	import argparse

	parser = argparse.ArgumentParser()

	parser.add_argument('-r', '--recipe', type=str,
	                    default='Human_sims/Gravel_2011_frame_sample.slim')

	parser.add_argument('-c', '--cookbook', type=str,
	                    default='simple2split')

	parser.add_argument('-s', '--short', type=str,
	                    default='Gravel')

	parser.add_argument('-a', '--assembly', type=str,
	                    default='hg38')

	parser.add_argument('--mem', type=str,
	                    default='15GB')

	parser.add_argument('-t', type=str,
	                    default='30:00:00')

	parser.add_argument('--nodes', type=int,
	                    default= 4)

	parser.add_argument('--rec', type=float,
	                    default= 1e-8)

	parser.add_argument('-N', type=int,
	                    default= 40)

	parser.add_argument('-L', type=int,
	                    default= 1000000)

	#parser.add_argument('-a', '--annotations', type=str, nargs='+', default= [])

	args = parser.parse_args()

	## directories
	main_dir= os.getcwd() + '/'
	slim_dir= ''
	fastas_dir= '/home/jgarc235/Fastas/'
	##

	## sub-directories.

	dir_data= main_dir + 'mutation_counter/data/gravel_10MB/'
	count_dir= main_dir + 'mutation_counter/count/'
	dir_launch= main_dir + 'mutation_counter'
	slim_soft= slim_dir + 'sim*'

	summary_file= 'sims.log'
	mutlog= 'toMut.log'

	#
	##
	## SLiM recipe.
	#sim_dir= main_dir + 'Recipes/Human_sims/'
	sim_recipe= main_dir + 'Recipes/' + args.recipe
	##
	##
	#

	###########   ##############################   #############
	############################################################

	## Simulation tag names, assembly to select from.
	batch_name= args.short
	assembly= args.assembly

	## files & variables
	## fasta segment lengths; number of segments / sims.
	L= args.L
	N= args.N


	############################################################
	########################      ##############################


	## Read chrom_sizes file to decide where to sample files from. 
	chrom_sizes= read_chrom_sizes(assembly)

	## Sample fasta.
	##
	fasta= fastas_dir + assembly + '.fa.gz'
	rseqs= region_samplev2(L, chrom_sizes,N, fasta)

	from tools.SLiM_pipe_tools import SLiM_dispenserv1
	from tools.cookbook import cook_constants_Gravel2sampleRange

	selected_book= 'cook_constants_' + args.cookbook

	import tools.cookbook
	book= getattr(tools.cookbook, selected_book)

	
	cookargs= {
		"s1": 500,
		"s2": 500,
		"s3": 500,
		'rec': args.rec
		}

	sim_store, cookID= book(rseqs,dir_data= dir_data,
	               slim_dir= slim_dir, batch_name= batch_name,**cookargs)

	print('launch SLiM jobs.')
	SLiM_dispenserv3(sim_store, sim_recipe, cookID= cookID, slim_dir= slim_dir, batch_name= batch_name,
	                    ID= cookID, L= L, logSims= summary_file, mutlog= mutlog)


	#########                                      ##############
	#############################################################

