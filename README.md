## SLiM simulations.

This repository contains scripts to faclitate simulations using SLiM, as well as scripts and notebooks to analyse the output.

Each directory is concerned with some particular feature of simulations. These problems were tackled one at a time and each benefitted from lessons already learned, so that the functions overlap, and simulations become increasingly complete.

The directories are listed below, in order of chronological development, roughly in increasing order of complexity.

### ABC simulations - .
This repository was motivated by the need to accurately translate the demographic inferences of a published deomographic model into SLiM simulations.

#### I. `demos`
Prado-Martinez et al. (2013) studied the evolution of chimpanzee subspecies. Using genetic data and simulations, the authors assessed the likelihood of several demographic models. They concluded on two most likely models, for which they provide summary statistics on the parameters.

We provide a set of functions to process this information into SLiM recipes. The input format required `demos` is described. 
    
> directory notebook: [link](https://github.com/SantosJGND/SLiM/tree/master/demos_ABC)

> Chimp demos: [PM2013_M4A](demos_ABC/demos/PM2013_M4A.txt); [PM2013_M3](demos_ABC/demos/PM2013_M3.txt)

#### II. Converting to demos format
Liu _et al._ (2018) studied the evolution of the rhesus monkey (_Macaca mulatta_). The authors performed demographic inference on 88 sequences from five subspecies and published samples of the best fitting parameters for the model tested.

We convert this data to the format `demos`.

notebook: [link](https://nbviewer.jupyter.org/github/SantosJGND/SLiM/blob/f6320157e0e6ea6a64503690e76519efda73b5e0/Rhesus_Liu_2018/Parameter_input.ipynb)

and generate SLiM recipes.

notebook: [link](https://nbviewer.jupyter.org/github/SantosJGND/SLiM/blob/f6320157e0e6ea6a64503690e76519efda73b5e0/Rhesus_Liu_2018/ABC_demo.ipynb)

> Rhesus demos: [rhesus_liu18](Rhesus_Liu_2018/rhesus_liu18.txt)


### OSG

Deployment to the Open Science Grid.

> [README.md](OSG/)