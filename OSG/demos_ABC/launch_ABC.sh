#!/bin/bash

#module purge
module load python/3.7.0
module load py-numpy/1.15.2-py3.7
module load py-scipy/1.1.0-py3.7

#module load slim/3.3.1

short=$1

python -u  SLiM_demosABC_osg.py \
-L 1000000  \
-R M4Amedian.slim \
-N 2 \
-c PMchimpMedian \
-s $short \
--cpus 3 \
--nmem 4 \
