## VCF sampling

It might be desirable when running simulations to compare them to real data, usually that data on which the models used are based. 

In order to facilitate the extraction of the same statistics across simulated and real data it is then useful to prepare these data in the same way. 

Below is a list of elements necessary for analysis and the formats used here:

**Required elements**: 

- genotype data - `VCF`;

**Optional elements**:

- reference sequence - `oneline fasta` (*);

- individual to population assignment - _e.g._ [ind_assignments.txt](chimp_ind_assignments.py);

- ancestral differences - [specific to real data] see [anc_diffs_file](diffs_example.py).


The pipelines proposed in this repository for the deployment of simulations using the software SLiM v3 (see [../Bash](../Bash/), [../OSG](../OSG/)) create replicate directories carrying the elements listed above except for ancestral differences. 

The python script `ordered_extractions.py` generates serializes the extraction of windows of genotypic data. 

**arguments**:

> `-v` or `--vcf`: Parent VCF file to extract from.

> `-l` or `--length`: length, in bp, of the windows to extract.

> `-n` or `--number`: number of windows to extract.

> `-a` or `--assembly`: reference fasta corresponding to the assembly used for `VCF`. 

> `-b` or `--batch`: instance ID. Appended to directy and file names across windows. 

> `-i` or `--ids`: individual to population assignment file.

> `-o` or `--out`: output directory where to store windows to. 

**optional arguments**:

> `-d` or `--diffs`: ancestral differences file. 

> `-c` or `--chrom`: extract from a single chromosome. 


**example deployment**


- bash deployment: [vcf_sampler.sh](vcf_sampler.sh)

##############################################################################
##############################################################################

[1] - Extractions are performed using VFCtools. 

[2] - Fasta and chrom_size file directories are currently hardcoded. 

[3] - Windows are extracted at random from across the genome (or chromosome) using the chrom_sizes.sizes file. 

[4] - Positions in VCF file and differences file are reset to the initial position of each region.

