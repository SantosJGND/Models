## Cross-window data summaries.

This section of the repository will deal with the downstream analysis of data prepared using the tools provided in sections [OSG](../OSG), [Bash](../Bash) ou [VCF_sampler](../VCF_sample). 

As described in those repositories, data is generated under the same format, with windows of genetic data stored w/ the corresponding reference sequence and individual to population assignment. File and directory names follow specific rules, see [data](/data) for an example.

The purposes of data simulation or extraction by windows will vary from project to project. Even within the same project the same data will often be analysed in different ways. It is important in this case to ensure that various analyses are coherent, i.e. that the same steps across analyses are completed by the exact same way. But it is also important to set up analyses in a way that will facilitate the repurposing of internal elements, and that updates or corrections to shared functions are carried automatically across platforms. 

The idea to identify which steps can be shared and follow a functional approach that allows these steps to be stored in shared tools directories. 

The following functions were set up with these purposes in mind. 

### I. INFO database.

The first thing is to ensure that different analyses are run on the same data. Having source directories managed from a single data base prevents us from having to re-type these every time, makeing things easier, but also avoiding mistakes.

This pipeline relies for this purpose on the INFO_db.py data base. The root keys of INFO_db are currently set to be **species** names. This is because several aspects of data will be shared between simulations and data extractions of data from the same species. Within each species' entry we find the keys:

> shared entries: 
- dirs: `dict`, this dictionary holds the directories where data windows were stored. Directories each have a key, e.g. "sims", indicating the type or batch.

- w_indfile: `dict`, the names of files holding individual to population association within each window. these can be batch specific (see [VCF_sample](../VCF_sample)). 

> species' specific data:

- g_indile: `str` , table of individual to population assignment from real data. While populations might not all be simulated (e.g. the AFR populaiton in the Human model, see [Bash](../Bash), but it is assumed that the populations are the same. 

- pop_dict: `dict`, association between population names in real and simulated data. 

- pop_colors: `dict`, population specific colors to be used in plots. 

- assembly: `str`, assembly name. 


### II. Input_utilities.

The first large axis of convergence among applications is reading the input. To be specific, different statistics might be extracted from the same data. Or the same statistics might be extracted from the same data at different levels (e.g. by population and by SNP).

The functions held in [input_utilities.py](tools/input_utilities.py) serve this purpose. Each takes as main input the directory where data is to be read. other important arguments pertain to how information is extracted. 

> **general arguments**: 
- frequency_range: `list`, alleles w/ frequencies outside this range are excluded. defaults to [0,1], open brackets.
- haps_extract: `bool`, use for phased data ploidy == 2. If True ALT allele counts are summed by ind and position, False haps are used as observations.  
- return_private: `bool`, excludes shared alleles if true.

> **mutation-type arguments**:

- ksize & bases: `int` & `str`, control mutation-type to count. 
- collapse: `bool`, collapse mutation counts by complementarity if true. 
- row & col: `int`, shape of mutation count matrices returned, numpy. 
- single: `bool`, if true, returns mutation type counts by individual as well as population. 

> **sub-sampling**: 
- stepup: `str`, if "increment", ranges between 2 and `args.samp`. Otherwise ranges between `args.samp` and Npop = population size. 
- steps: number of steps in sample size between the range determined by `args.stepup`.
- reps: `int`, number of replicates by sample size by population.

#### II.a. input_cofactors.py.

Functions in input_utilities set the platform for reading and storing data. Functions in input_cofactors perform the work.


### III. Mcounter tools.

Functions to extract statistics. The same structure is applied as for input utilities. [mcounter_utilities](tools/mcounter_utilities.py) functions set the platform parse and structure data, calculations are handed to [mcounter_cofactor](tools/mcounter_cofactors.py) functions. 

Mcounter functions have specific purposes that will be specific to the scientific question at hand and we do not go here into detail. 