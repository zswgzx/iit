while read directory
do 
visit=${directory:0:6}
projid=${directory:10:18}

# for /data/aging_dti/ in subjects0.txt, see list0.txt for details
# 24 subjects below
<< EOF
cp /data/aging_dti/$directory/DICOM/PA000000/ST000000/SE000003_proc/SE000003_DMC_R1_SAVE/SE000003_DMC_R1_DT.hdr $projid.hdr
cp /data/aging_dti/$directory/DICOM/PA000000/ST000000/SE000003_proc/SE000003_DMC_R1_SAVE/SE000003_DMC_R1_DT.img $projid.img
EOF

# 1 subject below
<<EOF
cp /data/aging_dti/$directory/DICOM/PA000000/ST000000/SE000004_proc/SE000004_DMC_R1_SAVE/SE000004_DMC_R1_DT.hdr $projid.hdr
cp /data/aging_dti/$directory/DICOM/PA000000/ST000000/SE000004_proc/SE000004_DMC_R1_SAVE/SE000004_DMC_R1_DT.img $projid.img
EOF

# 133 subjects below
<<EOF
cp /data/aging_dti/$directory/${visit}_${projid}/DICOM/PA000000/ST000000/SE000003_proc/SE000003_DMC_R1_SAVE/SE000003_DMC_R1_DT.hdr $projid.hdr
cp /data/aging_dti/$directory/${visit}_${projid}/DICOM/PA000000/ST000000/SE000003_proc/SE000003_DMC_R1_SAVE/SE000003_DMC_R1_DT.img $projid.img
EOF

# 3 subjects below
<<EOF
cp /data/aging_dti/$directory/${projid}/DICOM/PA000000/ST000000/SE000003_proc/SE000003_DMC_R1_SAVE/SE000003_DMC_R1_DT.hdr $projid.hdr
cp /data/aging_dti/$directory/${projid}/DICOM/PA000000/ST000000/SE000003_proc/SE000003_DMC_R1_SAVE/SE000003_DMC_R1_DT.img $projid.img
EOF

# 1 subject below
<<EOF
cp /data/aging_dti/$directory/SE000003_proc/SE000003_DMC_R1_SAVE/SE000003_DMC_R1_DT.hdr $projid.hdr
cp /data/aging_dti/$directory/SE000003_proc/SE000003_DMC_R1_SAVE/SE000003_DMC_R1_DT.img $projid.img
EOF

# 1 subject below
<<EOF
cp /data/aging_dti/$directory/${visit}_${projid}/DICOM/PA000000/ST000000/SE000005_proc/SE000005_DMC_R1_SAVE/SE000005_DMC_R1_DT.hdr $projid.hdr
cp /data/aging_dti/$directory/${visit}_${projid}/DICOM/PA000000/ST000000/SE000005_proc/SE000005_DMC_R1_SAVE/SE000005_DMC_R1_DT.img $projid.img
EOF

# 1 subject below
<<EOF
cp /data/aging_dti/$directory/${projid}.MAP/DICOM/PA000000/ST000000/SE000003_proc/SE000003_DMC_R1_SAVE/SE000003_DMC_R1_DT.hdr $projid.hdr
cp /data/aging_dti/$directory/${projid}.MAP/DICOM/PA000000/ST000000/SE000003_proc/SE000003_DMC_R1_SAVE/SE000003_DMC_R1_DT.img $projid.img
EOF

# 1 subject below
<<EOF
cp /data/aging_dti/$directory/${visit}_${projid}/DICOM/PA000000/ST000000/SE000004_proc/SE000004_DMC_R1_SAVE/SE000004_DMC_R1_DT.hdr $projid.hdr
cp /data/aging_dti/$directory/${visit}_${projid}/DICOM/PA000000/ST000000/SE000004_proc/SE000004_DMC_R1_SAVE/SE000004_DMC_R1_DT.img $projid.img
EOF

# for /data01/aging_dti01/ in subjects1.txt, see list1.txt for details
# 29 subjects below
<<EOF
cp /data01/aging_dti01/$directory/DICOM/PA000000/ST000000/SE000003_proc/SE000003_DMC_R1_SAVE/SE000003_DMC_R1_DT.hdr $projid.hdr
cp /data01/aging_dti01/$directory/DICOM/PA000000/ST000000/SE000003_proc/SE000003_DMC_R1_SAVE/SE000003_DMC_R1_DT.img $projid.img
EOF

# 1 subject below
<<EOF
cp /data01/aging_dti01/$directory/DICOM/PA000000/ST000000/SE000005_proc/SE000005_DMC_R1_SAVE/SE000005_DMC_R1_DT.hdr $projid.hdr
cp /data01/aging_dti01/$directory/DICOM/PA000000/ST000000/SE000005_proc/SE000005_DMC_R1_SAVE/SE000005_DMC_R1_DT.img $projid.img
EOF

# 23 subjects below
<<EOF
cp /data01/aging_dti01/$directory/PA000000/ST000000/SE000003_proc/SE000003_DMC_R1_SAVE/SE000003_DMC_R1_DT.hdr $projid.hdr
cp /data01/aging_dti01/$directory/PA000000/ST000000/SE000003_proc/SE000003_DMC_R1_SAVE/SE000003_DMC_R1_DT.img $projid.img
EOF

# for /data02/aging_dti02/ in subjects2.txt, see list2.txt for details
# 21 subjects below
<<EOF
cp /data02/aging_dti02/$directory/DICOM/PA000000/ST000000/SE000003_proc/SE000003_DMC_R1_SAVE/SE000003_DMC_R1_DT.hdr $projid.hdr
cp /data02/aging_dti02/$directory/DICOM/PA000000/ST000000/SE000003_proc/SE000003_DMC_R1_SAVE/SE000003_DMC_R1_DT.img $projid.img
EOF

# 22 subjects below
<<EOF
cp /data02/aging_dti02/$directory/PA000000/ST000000/SE000003_proc/SE000003_DMC_R1_SAVE/SE000003_DMC_R1_DT.hdr $projid.hdr
cp /data02/aging_dti02/$directory/PA000000/ST000000/SE000003_proc/SE000003_DMC_R1_SAVE/SE000003_DMC_R1_DT.img $projid.img
EOF

# 1 subject below
<<EOF
cp /data02/aging_dti02/$directory/SE000003_proc/SE000003_DMC_R1_SAVE/SE000003_DMC_R1_DT.hdr $projid.hdr
cp /data02/aging_dti02/$directory/SE000003_proc/SE000003_DMC_R1_SAVE/SE000003_DMC_R1_DT.img $projid.img
EOF

done < subjects0.txt
