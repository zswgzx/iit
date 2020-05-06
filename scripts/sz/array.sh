#!/bin/bash

NUMSUB=52
NSPLIT=2

# open file that contains all subject names
array=(`cat subjects-all | tr '\n' ' '`)

# main loop
for (( i=0; i<$NSPLIT; i++ )) ;do
    sub=${array[i]}
    for (( j=`echo "$i+1" | bc -l`; j<$NUMSUB; j++ )) ;do
      trg=${array[j]}
      echo $sub " -> " $trg
    done
done

exit
