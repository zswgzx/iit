/*sz last updated: Dec 2, 09
Func: write  data to a txt file for matlab histogram process
*/

#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>

int main(int argc, char **argv)
{
  FILE *fp_img,*fp_mask,*fp_txt;
  float std,THR;
  char filename[200],type_flag,mask;
  int  s, r, c, i, j,seek_flag,ZRES;

  // assert number of arguments
  if (argc != 3) {
    printf("%s - write non-zeros intensity voxels to a txt file for histogram in matlab, it's for whole volume\n",argv[0]);
    fprintf(stderr, "Usage: %s <image> <txt>\n",argv[0]);
    fprintf(stderr, "i.e. %s hp_fa_std.img hp_fa_std.txt\n",argv[0]);
    exit(1);
  }
 
  // open std
  sprintf(filename,argv[1]);
  if((fp_img=fopen(filename,"rb"))==NULL) {
    fprintf(stderr,"Error opening file %s.\n",filename);
    exit(1);
  }

  // open txt
  sprintf(filename,argv[2]);
  if((fp_txt=fopen(filename,"w"))==NULL) {
    fprintf(stderr,"Error opening file %s.\n",filename);
    exit(1);
  }

  //-----------------------------------------
  //vox-wise comparison
  while(feof(fp_img)==0) {
    //read intensities
    fread(&std,sizeof(float),1,fp_img);
  
    //check value 
    if (std>0)	  fprintf(fp_txt,"%.8f\n",std);
    
  }

  fclose(fp_img);
  fclose(fp_txt);

  return 0;
}
