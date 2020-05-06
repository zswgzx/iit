/*********************************************************************
 * SZ modified from code of Kate Fissell, Univ of Pittsburgh, May 2005.
 *********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "nifti1.h"

#define MIN_HEADER_SIZE 348
#define NII_HEADER_SIZE 352
#define FG fgets(c120,120,fpin)

int main(int argc,char *argv[]) 
{
  nifti_1_header hdr;
  nifti1_extender pad;
  FILE *fpin,*fpout;
  int ret,i;
  char c120[120],descript[80]="SZhang";

  // process commandline parameters
  if (argc != 3) {
    fprintf(stderr, "\n%s - set all header info to nii file",argv[0]);
    fprintf(stderr, "\nUsage: %s <setting file> <output file>\n",argv[0]);
    fprintf(stderr, "\ne.g. %s input.txt raw.nii\n",argv[0]);
    exit(1);
  }

  // open and read header
  fpin = fopen(argv[1],"rb");
  if (fpin == NULL) {
    fprintf(stderr, "\nError opening setting file %s\n",argv[1]);
    exit(1);
  }

 fscanf(fpin,"%d",&hdr.sizeof_hdr);FG;
 FG;FG;FG;FG;FG;
 hdr.dim_info='\0';FG;
 for (i=0;i<8;i++) {fscanf(fpin,"%hd",&hdr.dim[i]);FG;}
 fscanf(fpin,"%f",&hdr.intent_p1);FG;
 fscanf(fpin,"%f",&hdr.intent_p2);FG;
 fscanf(fpin,"%f",&hdr.intent_p3);FG;
 fscanf(fpin,"%hd",&hdr.intent_code);FG;
 fscanf(fpin,"%hd",&hdr.datatype);FG;
 fscanf(fpin,"%hd",&hdr.bitpix);FG;
 fscanf(fpin,"%hd",&hdr.slice_start);FG;
 for (i=0;i<8;i++) {fscanf(fpin,"%f",&hdr.pixdim[i]);FG;}
 fscanf(fpin,"%f",&hdr.vox_offset);FG;
 fscanf(fpin,"%f",&hdr.scl_slope);FG;
 fscanf(fpin,"%f",&hdr.scl_inter);FG;
 fscanf(fpin,"%hd",&hdr.slice_end);FG;
 fscanf(fpin,"%d",&hdr.slice_code);FG;
 fscanf(fpin,"%d",&hdr.xyzt_units);FG;
 fscanf(fpin,"%f",&hdr.cal_max);FG;
 fscanf(fpin,"%f",&hdr.cal_min);FG;
 fscanf(fpin,"%f",&hdr.slice_duration);FG;
 fscanf(fpin,"%f",&hdr.toffset);FG;
 FG;FG;

 fscanf(fpin,"%79s",hdr.descrip);
 if (!strncmp(hdr.descrip,"n/a",3)) {
   strncpy(hdr.descrip,descript,6);
   hdr.descrip[6]='\0';
 }
 hdr.descrip[79]='\0';
 //hdr.descrip[0]='\0';
 FG;
 fscanf(fpin,"%23s",hdr.aux_file);
 if (!strncmp(hdr.aux_file,"n/a",3)) hdr.aux_file[0]='\0';
 hdr.aux_file[23]='\0';
 //hdr.aux_file[0]='\0';
 FG;

 fscanf(fpin,"%hd",&hdr.qform_code);FG;
 fscanf(fpin,"%hd",&hdr.sform_code);FG;
 fscanf(fpin,"%f",&hdr.quatern_b);FG;
 fscanf(fpin,"%f",&hdr.quatern_c);FG;
 fscanf(fpin,"%f",&hdr.quatern_d);FG;
 fscanf(fpin,"%f",&hdr.qoffset_x);FG;
 fscanf(fpin,"%f",&hdr.qoffset_y);FG;
 fscanf(fpin,"%f",&hdr.qoffset_z);FG;
 fscanf(fpin,"%f %f %f %f",&hdr.srow_x[0],&hdr.srow_x[1],&hdr.srow_x[2],&hdr.srow_x[3]);FG;
 fscanf(fpin,"%f %f %f %f",&hdr.srow_y[0],&hdr.srow_y[1],&hdr.srow_y[2],&hdr.srow_y[3]);FG;
 fscanf(fpin,"%f %f %f %f",&hdr.srow_z[0],&hdr.srow_z[1],&hdr.srow_z[2],&hdr.srow_z[3]);FG;

 fscanf(fpin,"%15s",hdr.intent_name);
 if (!strncmp(hdr.intent_name,"n/a",3)) hdr.intent_name[0]='\0';
 hdr.intent_name[15]='\0';
 //hdr.intent_name[0]='\0';
 FG;
 fscanf(fpin,"%3s",hdr.magic);FG;hdr.magic[3]='\0';
 fscanf(fpin,"%d",&pad.extension[0]);FG;
 for (i=1;i<=3;i++) pad.extension[1]=pad.extension[2]=pad.extension[3]=0;
 fclose(fpin);

 //output and write header
 fpout = fopen(argv[2],"wb");
 if (fpout == NULL) {
   fprintf(stderr, "\nError opening header file %s\n",argv[2]);
   exit(1);
 }

 ret = fwrite(&hdr, MIN_HEADER_SIZE, 1, fpout);
 if (ret != 1) {
   fprintf(stderr, "\nError writing header file %s\n",argv[2]);
   exit(1);
 }
 ret = fwrite(&pad, 4, 1, fpout);
 if (ret != 1) {
   fprintf(stderr, "\nError writing header file extension to %s\n",argv[2]);
   exit(1);
 }
 fclose(fpout);
  /*
  ret = fread(&hdr, MIN_HEADER_SIZE, 1, fpin);
  if (ret != 1) {
    fprintf(stderr, "\nError reading header file %s\n",argv[1]);
    exit(1);
  }
  */
  return 0;
}

