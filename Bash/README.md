### I. `demos`

Prado-Martinez et al. (2013) studied the evolution of chimpanzee subspecies. Using genetic data and simulations, the authors assessed the likelihood of several demographic models. They concluded on two most likely models, for which they provide summary statistics on the parameters.

We provide a set of functions to process this information into SLiM recipes. The input format required `demos` is described. 
    
> directory README: [link](https://github.com/SantosJGND/SLiM/tree/master/Bash/demos_ABC)

> Chimp demos: [PM2013_M4A](Bash/demos_ABC/demos/PM2013_M4A.txt); [PM2013_M3](Bash/demos_ABC/demos/PM2013_M3.txt)

### II. Converting to demos format
Liu _et al._ (2018) studied the evolution of the rhesus monkey (_Macaca mulatta_). The authors performed demographic inference on 88 sequences from five subspecies and published samples of the best fitting parameters for the model tested.

We convert this data to the format `demos`.

notebook: [link](https://nbviewer.jupyter.org/github/SantosJGND/SLiM/blob/master/Bash/Rhesus_Liu_2018/Parameter_input.ipynb)

and generate SLiM recipes.

notebook: [link](https://nbviewer.jupyter.org/github/SantosJGND/SLiM/blob/master/Bash/Rhesus_Liu_2018/ABC_demo.ipynb)

> Rhesus demos: [rhesus_liu18](Bash/Rhesus_Liu_2018/rhesus_liu18.txt)
