/*********************************************************************
 * SZ modified from code of Kate Fissell, Univ of Pittsburgh, May 2005.
 *********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include "nifti1.h"

#define MIN_HEADER_SIZE 348

int main(int argc,char *argv[]) 
{
  nifti_1_header hdr;
  nifti1_extender pad;
  FILE *fpin;
  int ret,i;

  // process commandline parameters
  if (argc != 2) {
    fprintf(stderr, "\n%s - output all header info from nii file",argv[0]);
    fprintf(stderr, "\nUsage: %s <raw file>\n",argv[0]);
    fprintf(stderr, "\ne.g. %s raw.nii\n",argv[0]);
    exit(1);
  }

  // open and read header
  fpin = fopen(argv[1],"rb");
  if (fpin == NULL) {
    fprintf(stderr, "\nError opening header file %s\n",argv[1]);
    exit(1);
  }

  ret = fread(&hdr, MIN_HEADER_SIZE, 1, fpin);
  if (ret != 1) {
    fprintf(stderr, "\nError reading header file %s\n",argv[1]);
    exit(1);
  }

  //print hdr info
  fprintf(stderr,"sizeof_hdr\t%d\n",hdr.sizeof_hdr);
  //fprintf(stderr,"data_type\t%s\n",hdr.data_type);
  //fprintf(stderr,"db_name\t%s\n",hdr.db_name);
  //fprintf(stderr,"extents\t%d\n",hdr.extents);
  //fprintf(stderr,"session_error\t%d\n",hdr.session_error);
  //fprintf(stderr,"regular\t%c\n",hdr.regular);
  fprintf(stderr,"dim_info(value)\t%d\n\n",hdr.dim_info);

  for (i=0;i<8;i++) fprintf(stderr,"dim[%d]\t%d\n",i,hdr.dim[i]);
  fprintf(stderr,"intent_p1\t%f\n",hdr.intent_p1);
  fprintf(stderr,"intent_p2\t%f\n",hdr.intent_p2);
  fprintf(stderr,"intent_p3\t%f\n",hdr.intent_p3);
  fprintf(stderr,"intent_code\t%d\n",hdr.intent_code);
  fprintf(stderr,"datatype\t%d\n",hdr.datatype);
  fprintf(stderr,"bitpix\t%d\n",hdr.bitpix);
  fprintf(stderr,"slice_start\t%d\n",hdr.slice_start);
  for (i=0;i<8;i++) fprintf(stderr,"pixdim[%d]\t%f\n",i,hdr.pixdim[i]);
  fprintf(stderr,"vox_offset\t%f\n",hdr.vox_offset);
  fprintf(stderr,"scl_slope\t%f\n",hdr.scl_slope);
  fprintf(stderr,"scl_inter\t%f\n",hdr.scl_inter);
  fprintf(stderr,"slice_end\t%d\n",hdr.slice_end);
  fprintf(stderr,"slice_code\t%d\n",hdr.slice_code);
  fprintf(stderr,"xyzt_units\t%d\n",hdr.xyzt_units);
  fprintf(stderr,"(xyz_units: 1=METER, 2=MM, 3=MICRON\nt_units: 8=SEC, 16=MSEC, 24=USEC, 32=HZ, 40=PPM, 48=RADS)\n");
  fprintf(stderr,"cal_max\t%f\n",hdr.cal_max);
  fprintf(stderr,"cal_min\t%f\n",hdr.cal_min);
  fprintf(stderr,"slice_duration\t%f\n",hdr.slice_duration);
  fprintf(stderr,"toffset\t%f\n\n",hdr.toffset);
  //fprintf(stderr,"glmax\t%d\n",hdr.glmax);
  //fprintf(stderr,"glmin\t%d\n",hdr.glmin);

  fprintf(stderr,"descrip\t%s\n",hdr.descrip);
  fprintf(stderr,"aux_file\t%s\n\n",hdr.aux_file);

  fprintf(stderr,"qform_code\t%d\n",hdr.qform_code);
  fprintf(stderr,"sform_code\t%d\n",hdr.sform_code);
  fprintf(stderr,"(q/sform_code: 0=UNKNOWN, 1=SCANNER_ANAT, 2=ALIGNED_ANAT, 3=TALAIRACH, 4=MNI_152)\n");
  fprintf(stderr,"qfac(if qform_code)\t%f\n",hdr.pixdim[0]);
  fprintf(stderr,"quatern_b\t%f\n",hdr.quatern_b);
  fprintf(stderr,"quatern_c\t%f\n",hdr.quatern_c);
  fprintf(stderr,"quatern_d\t%f\n",hdr.quatern_d);
  fprintf(stderr,"qoffset_x\t%f\n",hdr.qoffset_x);
  fprintf(stderr,"qoffset_y\t%f\n",hdr.qoffset_y);
  fprintf(stderr,"qoffset_z\t%f\n",hdr.qoffset_z);
  fprintf(stderr,"srow_x[1-4]\t%f  %f  %f  %f\n",hdr.srow_x[0],hdr.srow_x[1],hdr.srow_x[2],hdr.srow_x[3]);
  fprintf(stderr,"srow_y[1-4]\t%f  %f  %f  %f\n",hdr.srow_y[0],hdr.srow_y[1],hdr.srow_y[2],hdr.srow_y[3]);
  fprintf(stderr,"srow_z[1-4]\t%f  %f  %f  %f\n\n",hdr.srow_z[0],hdr.srow_z[1],hdr.srow_z[2],hdr.srow_z[3]);

  fprintf(stderr,"intent_name\t%s\n",hdr.intent_name);
  fprintf(stderr,"magic\t%s\n",hdr.magic);

  //print extender
  ret = fread(&pad, 4, 1, fpin);
  if (ret != 1) {
    fprintf(stderr, "\nError reading header file %s\n",argv[1]);
    exit(1);
  }
  fclose(fpin);
  fprintf(stderr,"extender\t%d  %d  %d  %d\n",pad.extension[0],pad.extension[1],pad.extension[2],pad.extension[3]);

  return 0;
}

