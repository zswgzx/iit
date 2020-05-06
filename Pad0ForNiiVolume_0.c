/*********************************************************************
 * gcc Pad0ForNiiVolume_x.c -o Pad0ForNiiVolume_x (Mac) 
 * To change to your datatype, change the lines:
 * typedef float MY_DATATYPE;
 * hdr.datatype = NIFTI_TYPE_FLOAT32;
 * hdr.bitpix = 32;
 *********************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "nifti1.h"

typedef float MY_DATATYPE;
#define MIN_HEADER_SIZE 348
#define NII_HEADER_SIZE 352

void pad_zeros(char *infile,char *outfile,int left,int right,int top,int down,int inf,int sup)
{
  nifti_1_header hdr;
  nifti1_extender pad={0,0,0,0};
  FILE *fpin,*fpout;
  int ret,i,j,k,l,rawx,rawy,count;
  MY_DATATYPE zero=0,data,*buf;

/********** open and read header */
  fpin = fopen(infile,"rb");
  if (fpin == NULL) 
    {
      fprintf(stderr, "\nError opening file %s\n",infile);
      exit(1);
    }
  ret = fread(&hdr, MIN_HEADER_SIZE, 1, fpin);
  if (ret != 1) 
    {
      fprintf(stderr, "\nError reading file %s\n",infile);
      exit(1);
    }

 /*temporary settings*/
  rawx=hdr.dim[1]+left;
  rawx+=right;
  rawy=hdr.dim[2]+top;
  rawy+=down;
  ret=hdr.dim[3]+inf;
  ret+=sup;

/********** print a little header information */
  fprintf(stderr, "%s header info:",infile);
  fprintf(stderr, "\nXYZT dimensions: %d %d %d %d",hdr.dim[1],hdr.dim[2],hdr.dim[3],hdr.dim[4]);
  fprintf(stderr, "\nDatatype code and bits/pixel: %d %d",hdr.datatype,hdr.bitpix);
  fprintf(stderr, "\nPad to dim: %d by %d by %d\n",rawx,rawy,ret);

  rawx=hdr.dim[1];
  rawy=hdr.dim[2];

/* modify dimension parameters*/
  hdr.dim[1] += ((short)left+(short)right);
  hdr.dim[2] += ((short)top+(short)down);
  hdr.dim[3] += ((short)inf+(short)sup);
  hdr.xyzt_units=NIFTI_UNITS_MM;

/********** write first 348 bytes of header   */
  fpout = fopen(outfile,"wb");
  if (fpout == NULL) 
    {
      fprintf(stderr, "\nError creating file %s\n",outfile);
      exit(1);
    }

  ret = fwrite(&hdr, MIN_HEADER_SIZE, 1, fpout);
  if (ret != 1) 
    {
      fprintf(stderr, "\nError writing header file %s\n",outfile);
      exit(1);
    }
  /* write extender*/
  ret = fwrite(&pad, 4, 1, fpout);
  if (ret != 1) 
    {
      fprintf(stderr, "\nError writing header file extension pad %s\n",outfile);
      exit(1);
    }
/********** open the datafile, jump to data offset */
  ret = fseek(fpin, (long)(hdr.vox_offset), SEEK_SET);
  if (ret != 0) 
    {
      fprintf(stderr, "\nError doing fseek() to %ld in data file %s\n",(long)(hdr.vox_offset), infile);
      exit(1);
    }
/*allocate buffer for one slice*/
  buf = (MY_DATATYPE *) malloc(sizeof(MY_DATATYPE)*hdr.dim[1]*hdr.dim[2]);
  if (buf == NULL) 
    {
      printf("*** Failed to allocate buffer\n");
      exit(1);
    }
/*copy data*/
  for (i=0;i<hdr.dim[4];i++)    //iter for volume
    {
      for (j=0;j<hdr.dim[3];j++)    //slice
	{
	  if (j<inf || j>hdr.dim[3]-1-sup) //inf and sup padding
	    {
	      for (k=0;k<hdr.dim[1]*hdr.dim[2];k++) buf[k]=zero;
	    }
	  else
	    {
	      for (k=0;k<top;k++) //top padding
		for (l=0;l<hdr.dim[1];l++)
		  buf[k*hdr.dim[2]+l]=zero;

	      for (k=top;k<rawy+top;k++)
		{
		  for (l=0;l<left;l++) buf[k*hdr.dim[2]+l]=zero;

		  for (l=left;l<left+rawx;l++)
		    {
		      ret=fread(&data,sizeof(MY_DATATYPE),1,fpin);
		      if (!ret)
			{
			  fprintf(stderr, "\nError reading volume %d at (%d,%d,%d) (start from 1, not 0)\n",i+1,l-left+1,k-top+1,j-inf+1);
			  exit(1);
			}
		      buf[k*hdr.dim[2]+l]=data;
		    }

		  for (l=left+rawx;l<hdr.dim[1];l++) buf[k*hdr.dim[2]+l]=zero;
		}

	      for (k=rawy+top;k<hdr.dim[2];k++)//down padding
		for (l=0;l<hdr.dim[1];l++)
		  buf[k*hdr.dim[2]+l]=zero;
	    }

	  ret=fwrite(buf,sizeof(MY_DATATYPE),hdr.dim[1]*hdr.dim[2],fpout);
	  if (ret != hdr.dim[1]*hdr.dim[2])
	    {
	      fprintf(stderr, "\nError writing volume %d at slice %d (start from 1, not 0)\n",i+1,j+1);
	      exit(1);
	    }
	}
    }

  free(buf);
  fclose(fpin);
  fclose(fpout);
}
/**********************************************************************
 *
 **********************************************************************/
int main(int argc,char *argv[]) 
{
  int ax_left,ax_right,ax_top,ax_down,ax_inf,ax_sup;
/********** process commandline parameters */
  if (argc != 9) 
    {
      fprintf(stderr, "%s - pad zeros around raw nii volume file",argv[0]);
      fprintf(stderr,"\nInput: nifti1.h and rawfile.nii, data type float(need to modify code for other types)\nOutput: user specified name, same data type");
      fprintf(stderr, "\n\nUsage: %s <raw file> <output file> <axial_left> <axial_right> <axial_top> <axial_down> <inferior> <superior>\n",argv[0]);
      fprintf(stderr, "\ne.g. %s raw_DMC_R1_tensor.nii(240*240*108->256*256*128) new.nii 8 8 8 8 10 10\n",argv[0]);
      exit(1);
    }

 ax_left=atoi(argv[3]);
 ax_right=atoi(argv[4]); 
 ax_top=atoi(argv[5]); 
 ax_down=atoi(argv[6]); 
 ax_inf=atoi(argv[7]);
 ax_sup=atoi(argv[8]);

 pad_zeros(argv[1],argv[2],ax_left,ax_right,ax_top,ax_down,ax_inf,ax_sup);
 return 0;
}
