#!/bin/bash

temp=0
sep=$( printf '\t') ## adjust as needed

while IFS=$sep read a b
do
   printf "%s$sep%s$sep%s\n" "$temp" "$a" "$b"
   let "temp+=1"
done < AtlasLabelPairs.txt