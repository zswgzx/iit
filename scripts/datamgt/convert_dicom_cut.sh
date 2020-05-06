#!/bin/bash
FILE=mg_zips.dat
#working_folder=/mri/invivo/raw/mg
while read line; do
        parent_folder=$(dirname "$line")
 	date_visit_projid=$(echo $parent_folder | cut -d/ -f2 )
        datestr=$(echo $date_visit_projid | cut -d_ -f1 )
        visit=$(echo $date_visit_projid | cut -d_ -f2 )
        projid=$(echo $date_visit_projid | cut -d_ -f3 )
        projid_visit=$(echo $projid"_"$visit)
        echo "Converting $parent_folder"
        unzip -q $line -d $parent_folder/dicom
        mkdir $parent_folder/$projid_visit
        ./nii_convert.sh $parent_folder $projid_visit
        cd $parent_folder
        zipname=$(echo $projid_visit"_nii.zip")
        zip -r $zipname ./$projid_visit -i \*.nii
        cd ..
        rm -rf $parent_folder/dicom
        rm -rf $parent_folder/$projid_visit        
	echo ""
        echo ""
done < $FILE
