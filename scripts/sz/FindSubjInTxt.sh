#!/bin/bash

n=195

#Declare array
declare -a NUM SUB

let count=0
while read LINE;do
  let line=$((LINE-1))
  NUM[$count]=$line
  ((count++))
done < group-$n.txt

#Open file for reading to array
let count=0
while read LINE;do
  SUB[$count]=$LINE
  ((count++))
done < subjects-204.txt

for (( i=0;i<$n;i++ ));do
  loc=${NUM[$i]}
  echo ${SUB[$loc]}
done
