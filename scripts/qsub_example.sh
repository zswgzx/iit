#!/bin/bash
# 
#$ -V
#$ -N randomise-fa-neg
#$ -S /bin/bash
#$ -j y
#$ -o /share/work1/szhang/projects/492/tbss/stats/082014/largeAlpha-cognglobal/randomise-fa-neg.log
#$ -wd /share/work1/szhang/projects/492/tbss/stats/082014/largeAlpha-cognglobal
#$ -l h=compute-0-6
#
randomise -i all_FA_skeletonised -o fa-largeAlpha-cognglobal-neg -m mean_FA_skeleton_mask -d largeAlpha-cognglobal-neg.mat -t largeAlpha-cognglobal-neg.con -n 5000 --T2
