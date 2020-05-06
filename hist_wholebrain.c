/*sz last updated: Dec 2, 09
Func: write  data to a txt file for matlab histogram process
Note: 
*/

#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>

#define XRES 256 
#define YRES 256
//#define ZRES 256
//#define NII_HEADER 352

int main(int argc, char **argv)
{
  FILE *fp_std,*fp_trstd,*fp_txt;
  float std,THR;
  char filename[200],type_flag;
  int  s, r, c, i, j,seek_flag,ZRES;

  // assert number of arguments
  if (argc != 4) 
    {
      printf("%s - write non-zeros intensity voxels to a txt file for histogram in matlab, it's for whole volume\n",argv[0]);
      fprintf(stderr, "Usage: %s <image> <txt> <#slices>\n",argv[0]);
      fprintf(stderr, "i.e. %s hp_fa_std.img hp_fa_std.txt 256\n",argv[0]);
      exit(1);
    }
  ZRES=atoi(argv[3]);

  // open std
  sprintf(filename,argv[1]);
  if((fp_std=fopen(filename,"rb"))==NULL)
    {
      fprintf(stderr,"Error opening file %s.\n",filename);
      exit(1);
    }

  // open txt
  sprintf(filename,argv[2]);
  if((fp_txt=fopen(filename,"w"))==NULL)
    {
      fprintf(stderr,"Error opening file %s.\n",filename);
      exit(1);
    }

  //-----------------------------------------
  //vox-wise comparison
  for(s=1;s<=ZRES;s++)
    {
      for(r=1;r<=YRES;r++)
	{
	  for(c=1;c<=XRES;c++)
	    {
	      //read intensities
	      fread(&std,sizeof(float),1,fp_std);

	      //check value 
	      if (std>0)	  fprintf(fp_txt,"%.8f\n",std);
	    }	
	}
    }

  fclose(fp_std);
  fclose(fp_txt);

  return 0;
}
