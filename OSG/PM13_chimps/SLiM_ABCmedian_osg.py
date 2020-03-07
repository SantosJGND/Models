
from tools.SLiM_pipe_tools import (
    read_chrom_sizes, region_samplev2,
    fasta_RextractUnif, write_fastaEx, process_recipe,
    SLiM_osg_dispenser, 
)

from tools.ABC_utilities import demo_to_recipe 

###########

if __name__ == '__main__':

        import argparse
        
        parser = argparse.ArgumentParser()
        
        parser.add_argument('-R', '--recipe', type=str,
                            default='Human_sims/Gravel_2011_frame_sample.slim')

        parser.add_argument('-d', '--demo', type=str,
                            default='demos/P&M2013_M3.txt')
        
        parser.add_argument('-c', '--cookbook', type=str,
                            default='SizeChange')
        
        parser.add_argument('-s', '--short', type=str,
                            default='Gravel')
        
        parser.add_argument('-a', '--assembly', type=str,
                            default='panTro5')

        parser.add_argument('-r','--rescale', type=float,
                            default= 1)
                
        parser.add_argument('-N', type=int,
                            default= 40)
        
        parser.add_argument('-L', type=int,
                            default= 1000000)
        
        parser.add_argument('--NeC', type=int,
                            default= 20000)
        
        parser.add_argument('--Nef', type=int,
                            default= 400000)
        
        parser.add_argument('--rate', type=float,
                            default= 1.03)
        
        parser.add_argument('--mu', type=float,
                            default= 1e-8)

        parser.add_argument('--rec', type=float,
                            default= 1e-8)
                
        parser.add_argument('--mfile', type=str,
                            default= 'mut_matrix_v0.txt')
        
        
        #parser.add_argument('-a', '--annotations', type=str, nargs='+', default= [])

        args = parser.parse_args()

        ## directories
        import os
        main_dir= os.getcwd() + '/'
        slim_dir= '/home/douradojns/SLiM/'
        fastas_dir= '/home/douradojns/SLiM/Fastas/'
        chrom_sizes_dir= '/home/douradojns/SLiM/chrom_sizes/'
        ##

        ## sub-directories.
        
        dir_data= main_dir + 'mutation_counter/data/sims/'
        count_dir= main_dir + 'mutation_counter/count/'
        dir_launch= main_dir + 'mutation_counter'
        slim_soft= slim_dir + 'sim*'
        matrix_dir= '' #'/' + main_dir + 'mut_matices/'
        log_dir= 'log'

        summary_file= 'sims.log'
        mutlog= 'toMut.log'


        os.makedirs(log_dir, exist_ok=True)

        #
        ##
        ## SLiM recipe.
        #sim_dir= main_dir + 'Recipes/Human_sims/'
        sim_template= main_dir + 'Recipes/' + args.recipe
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
        chrom_sizes= read_chrom_sizes(assembly, size_dir= chrom_sizes_dir)

        ## Sample fasta.
        ##
        fasta= fastas_dir + assembly + '.fa.gz'
        rseqs= region_samplev2(L, chrom_sizes,N, fasta)
        
        from tools.SLiM_pipe_tools import SLiM_dispenserv2

        ## Perform Simulations
        ## I. get cookfunction and args:
        selected_book= 'cook_constants_' + args.cookbook

        import tools.cookbook
        book= getattr(tools.cookbook, selected_book)

        cookargs= {
            "s1": 500,
            "s2": 500,
            "s3": 500,
            "s4": 500,
            "mu": args.mu,
            "rec": args.rec
        }

        sim_store, cookID= book(rseqs,dir_data= dir_data,
                       slim_dir= slim_dir, batch_name= batch_name,**cookargs)

        print('launch SLiM jobs.')
        SLiM_osg_dispenser(sim_store, sim_recipe= sim_template,cookID= cookID, slim_dir= slim_dir, batch_name= batch_name,
                            ID= cookID, L= L, logSims= summary_file, mutlog= mutlog,dir_data= dir_data,
                            cpus= 1,Nmem= 1,mem= 'GB',diskN= 1,diskS= 'GB',log_dir= log_dir)
        
        #########                                      ##############
        #############################################################
