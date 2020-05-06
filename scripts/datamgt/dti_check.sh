#!/bin/bash

# script: dti_check.sh
# purpose: check the data stored on the IIT Buffalo is "DTI Complete"
# For each MRI/scan/projid, the following must be present 
# 1. the co* files in the ROOT  [co(projid).nii & o(projid).nii]
# 2. in _proc, the files named SE0000?_DMC{template,structural}.nii
# 3. in _proc, the files in dir SE0000?_DMC_R1_SAVE/*

# Optional:
# 4. transfer into _proc_hold, the files in dir SE0000?_DMC_R1_FSL_RAW_UNSORTED/*


LOG=/work1/dti/dti_check.log

DTI_ROOT=/work1/dti/

#foreach folder in srce, copy to dest
#DIR_NAME="090211_01_54519335"

while read DIR_NAME; do

FULL_PATH="$DTI_ROOT$DIR_NAME"

#check if directory exists	
if [ ! -d "$FULL_PATH" ]; then
  echo "$DIR_NAME is missing root directory; not at $FULL_PATH"
else
  #now check 1-4
  # 1. co files
  CO_COUNT=`ls -1 $FULL_PATH/co*.nii | wc -l `

  # 2 structural
  STRU_COUNT=`ls -1 $FULL_PATH/*structural*.nii | wc -l `
  TEMP_COUNT=`ls -1 $FULL_PATH/*template*.nii | wc -l `

  # 3 _DMC_R1_SAVE
  #get the full path to DMC_R1_SAVE  (dir)
  DMC_R1_SAVE_PATH=`find $FULL_PATH -type d -name "*_DMC_R1_SAVE"`   

  if [ "$DMC_R1_SAVE_PATH" ]; then
    DMC_R1_SIZE=`du -sk $DMC_R1_SAVE_PATH | awk '{print $1}'`
  else 
    DMC_R1_SIZE=0
  fi

  #output whether each directory is correct
  if [[ $CO_COUNT -ne 0 && $STRU_COUNT -ne 0 && $TEMP_COUNT -ne 0 && $DMC_R1_SIZE -gt 50000 ]]; then
    echo "$DIR_NAME [CORRECT] " >> $LOG
  else
    echo "$DIR_NAME [HAS ERROR]   " >> $LOG

    if [[ $CO_COUNT -eq 0 ]]; then
      echo "$DIR_NAME [ERROR] no co*.nii files" >> $LOG  
    fi

    if [[ $STRU_COUNT -eq 0 ]]; then
      echo "$DIR_NAME [ERROR] no *DMCstructural.nii file" >> $LOG  
    fi

    if [[ $TEMP_COUNT -eq 0 ]]; then
      echo "$DIR_NAME [ERROR] no *DMCtemplate.nii file" >> $LOG  
    fi

    if [[ $DMC_R1_SIZE -le 50000 ]]; then
      echo "$DIR_NAME [ERROR] DMC_R1_SAVE is $DMC_R1_SIZE kilobytes" >> $LOG  
    fi  
  fi
fi

done < dirs_check.txt
