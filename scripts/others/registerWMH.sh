while read directory
do 
visit=${directory:0:6}
projid=${directory:10:18}
#printf -v T1 "/data/aging_dti/%s/co*" "$directory" #"$visit" "$projid"
echo $directory
printf -v AM "/data01/aging_dti01/tbss/registration/Amplitude_maps/%s.hdr" "$projid"
#echo $AM
printf -v FL "/data01/aging_dti01/tbss/registration/FLAIR/%s/%s.FL_orig_brain.nii" "$directory" "$directory"
#echo $FL
printf -v FL_out "/data01/aging_dti01/tbss/registration/registered_FLAIR/%s.nii" "$projid"
#echo $FL_out
printf -v WMH "/data01/aging_dti01/tbss/registration/FLAIR/%s/%s.FL.wml_06.hdr" "$directory" "$directory"
printf -v WMH_out "/data01/aging_dti01/tbss/registration/registered_WMH/%s.nii" "$projid"
rm $FL_out
rm $WMH_out
flirt -in $FL -ref $AM -o $FL_out -searchcost normmi -cost normmi -dof 6 -omat reg_params.mat
flirt -in $WMH -ref $AM -out $WMH_out -applyxfm -init reg_params.mat
printf -v FL_out_gz "%s.gz" $FL_out
printf -v WMH_out_gz "%s.gz" $WMH_out
gunzip $FL_out_gz
gunzip $WMH_out_gz


done < directories.txt
