
### OSG

The Open Science Grid is a distributed computing partnership for data-intensive research.

OSG is a useful tool for when many jobs must be deployed, such as in simulation studies. 

deploying to the OSG network requires only an additional layer compared to submitting to a regular computing cluster. 

This layer consists of a `submit` file describing the executable to run, arguments, necessary input files and memory requirements. 

- visit the OSG [website](https://support.opensciencegrid.org/support/home) for more information.

## This repertoire.

OSG submission allows jobs to be replicated. However in certain settings , for example for ranging parameter variables, it can still be productive to have the writting of submit files be automated.

The function `osg_template` performs this task. 

See [osg_submit.py](osg_submit.py). 

## Application. 

We used the osg_template function to perform simulations of primate species' evolution. 

- Human deployment: [launch_Gravel_osg.sh](Gravel/launch_Gravel_osg.sh);
- Chimp deployment: [launch_ABC.sh](demos_ABC/launch_ABC.sh);
- Rhesus deployment: [launch_Rhesus_Liu18.sh](Rhesus_Liu18/launch_Rhesus_Liu18.sh);

