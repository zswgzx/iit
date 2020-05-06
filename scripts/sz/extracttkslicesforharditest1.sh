#!/bin/bash

#This script create slices from all dw directions of all datasets, for hardi test
#data used is 45 slices smth_x/sq? image from Anna, deformation is from subject-based(ICBM_JD) registration of TK 

# check # of arguments
if [[ $# < 4 ]]; then
    echo "Input: subjects-67.txt, dm0.hdr"
    echo "Usage: $0 <starting slice> <# slices to process> <x or sq> <b0 exist? y:0, n:1>"
    echo "e.g. $0 71 90 sq 0"
    exit 1
fi

sad=$(awk '$1' subjects-67.txt)
#sub='TRT007'
subjnum=1
startslice=$1 
nslice=$2
temp=$startslice

#step1:make dir
cd ~/Desktop/dtigui_template/harditest/dwdata
for ((j=1;j<=$nslice;j++))
do
    mkdir slice$temp
    let "temp+=1"
done
cd -

#step2: subject-wise loop
for sub in $sad; do
    cd $sub/
    mv ~/Desktop/dtigui_template/harditest/{$sub*.df.nii.gz,warpfa.nii.gz,renamedwdata1} ./    

    #deal with DWs 
    for ((i=$4;i<=12;i++ ))	#start from 0 or 1,depending on presence of b0
    do
        #add hdr 
	cp ../dm0.hdr smth2_$3_d${i}.hdr
        #apply deformation (amazing, the hdr file doesn't need to be ziped!)
	deformationScalarVolume -in smth2_$3_d${i}.img.gz -trans $sub-combined.df.nii.gz -out dmwarp${i}.nii.gz -target warpfa.nii.gz
        #extract 1 slice
	fslslice dmwarp${i}.nii.gz dmwarp${i}
        renamedwdata1 $startslice $nslice $i $subjnum
        rm dmwarp*.gz
    done

    mv $sub*.df.nii.gz warpfa.nii.gz renamedwdata1 ~/Desktop/dtigui_template/harditest/
    let "subjnum += 1"
    rm *hdr
    cd ../    
    printf "%s done\n" $sub

done
