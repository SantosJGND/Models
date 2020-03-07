#!/bin/bash

#module purge
module load python/3.7.0
module load py-numpy/1.15.2-py3.7
module load py-scipy/1.1.0-py3.7

#module load slim/3.3.1

short=$1

python -u  SLiM_ABCmedian_osg.py \
-L 1000000  \
-R M4A_median.slim \
-N 10 \
-c PMchimpMedian \
-s $short \

