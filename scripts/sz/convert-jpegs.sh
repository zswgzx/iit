EXT_RAW='ppm'
EXT_CONV='jpeg'

for i in `ls -l| grep dr | awk '{print $9}'`;do
    for fname in `ls $i/*.${EXT_RAW}`;do
	fname1=`echo ${fname%.*}`
	~/work/software/afni/linux_gcc33_64/cjpeg -outfile $fname1.${EXT_CONV} -progressive $fname1.${EXT_RAW}
    done
    echo "$EXT_RAW files in folder $i converted to $EXT_CONV format"
    rm $i/*.${EXT_RAW}
done
