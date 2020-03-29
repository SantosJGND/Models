#!/bin/bash

#module purge
module load python/3.7.0
module load py-numpy/1.15.2-py3.7
module load py-scipy/1.1.0-py3.7

#module load slim/3.3.1

short=$1

python -u  SLiM_Rhesus_Liu18_osg.py \
-d "rhesus_liu18.txt" \
-L 1000000  \
-R demos_mat/template_matVar2.slim \
-N 10 \
-c ABC \
-r 1 \
-s $short \

