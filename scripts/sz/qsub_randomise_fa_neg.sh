#!/bin/bash
# 
#$ -V
#$ -N randomise-fa-neg
#$ -S /bin/bash
#$ -j y
#$ -o /share/work1/szhang/projects/492/tbss/stats/082014/gamma-cognglobal/randomise-fa-neg.log
#$ -wd /share/work1/szhang/projects/492/tbss/stats/082014/gamma-cognglobal
#
randomise -i all_FA_skeletonised -o fa-gamma-cognglobal-neg -m mean_FA_skeleton_mask -d gamma-cognglobal-neg.mat -t gamma-cognglobal-neg.con -n 5000 --T2
