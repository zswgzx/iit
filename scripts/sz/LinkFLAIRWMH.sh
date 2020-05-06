for sub in `cat ../subjects`;do
	ln -s ~/work/dti-t1-wmh/wmh-voxelwise/FLAIR-T1-WMH-raw/$sub/$sub-FL-orig-masked.nii.gz
	ln -s ~/work/dti-t1-wmh/wmh-voxelwise/FLAIR-T1-WMH-raw/$sub/$sub-FL-wml-masked.nii.gz
	echo $sub done
done
