/* SZ 
revised from xcorr4img-v1.c, should replace it
*/

#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include <math.h>

float pair_matching(char a[],char b[], char arg[], int XRES,int YRES,int ZRES,int numsub)
{
  FILE *fp_org[2],*fp_mask;
  char filename[130],mask;
  int i,s,r,c;
  float raw[2],sum_denom1=0,sum_denom2=0,sum_numer=0,xcorr;
  const char *path_sub="reg%s_%s.img";

  //open data files from subject a
  sprintf(filename,path_sub,a,arg);
  if((fp_org[0]=fopen(filename,"rb"))==NULL) {
    fprintf(stderr,"Error opening file %s\n",filename);
    exit(1);
  }
  
  //open data files from subject b
  sprintf(filename,path_sub,b,arg);
  if((fp_org[1]=fopen(filename,"rb"))==NULL) {
    fprintf(stderr,"Error opening file %s\n",filename);
    exit(1);
  }

  //open mask file
  sprintf(filename,"commonmask%d.img",numsub);
  if((fp_mask=fopen(filename,"rb"))==NULL) {
    fprintf(stderr,"Error opening file %s\n",filename);
    exit(1);
  }
  
  // voxel based loop
  for(s=1;s<=ZRES;s++)
    for(r=1;r<=YRES;r++)
      for(c=1;c<=XRES;c++) {
 	  //read data from both subjects and mask
	  if ((fread(&raw[0],sizeof(float),1,fp_org[0]))!=1) {
	    fprintf(stderr,"Error reading raw file 1(s%d r%d c%d).\n",s,r,c);
	    exit(1);
	  }
	  if ((fread(&raw[1],sizeof(float),1,fp_org[1]))!=1) {
	    fprintf(stderr,"Error reading raw file 2(s%d r%d c%d).\n",s,r,c);
	    exit(1);
	  }
	  if ((fread(&mask,sizeof(char),1,fp_mask))!=1) {
	    fprintf(stderr,"Error reading mask(s%d r%d c%d).\n",s,r,c);
	    exit(1);
	  }

	  // main comparison
	  if (mask != 0) {
	    //xcorr of numerator and denominator
	    sum_numer += raw[0]*raw[1];
	    sum_denom1 += pow(raw[0],2);
	    sum_denom2 += pow(raw[1],2);
	  }
      }

  xcorr=sum_numer/sqrt(sum_denom1*sum_denom2);	 

  // clean up file pointers
  for (i=0;i<2;i++) fclose(fp_org[i]);
  fclose(fp_mask);

  return xcorr;
}
//**************************************************************************************//

int main(int argc, char **argv)
{ 
  char filename[180],subname[20];
  FILE  *fp_sub,*fp_results;
  int i,j,subc=0,xres,yres,zres,numsub;
  float temp;

  if (argc!=6) {
    fprintf(stderr,"%s - output cross correlation values of all pairs of scalar maps within common mask to a text file\n",argv[0]);
    fprintf(stderr, "\nInput: subjects-??.txt,reg[subj]_[fa etc.].img (float), commonmask??.img (char)\nOutput: results_[fa etc.].txt\n");
    fprintf(stderr,"\nUsage: %s <fa or trace> <xres> <yres> <zres> <subjects #>", argv[0]);
    fprintf(stderr,"\ne.g. %s fa 256 256 181 22\n", argv[0]);
    exit(1);
  }

  xres=atoi(argv[2]);
  yres=atoi(argv[3]);
  zres=atoi(argv[4]);
  numsub=atoi(argv[5]);
  char sub[numsub][20];

  //////////////////////////////////////////////////////////
  //open file that contains all subject names
  sprintf(filename,"subjects-%d.txt",numsub);
  if((fp_sub=fopen(filename,"r"))==NULL) {
    fprintf(stderr,"Error opening file subjects.txt\n");
    exit(1);
  }

  // store subject names to arrays
  while (fgets(subname, sizeof(subname), fp_sub)) {
    // cut newline from the end of subname
    char *tmp = subname;
    while (*tmp) {
      if (*tmp == '\n') *tmp = '\0';
      tmp++;
    }
    strcpy(sub[subc],subname);
    if (++subc == numsub) break;
   }

  fclose(fp_sub);
  ////////////////////////////////////////////////////////
  // open file to write results
  sprintf(filename,"results_%s.txt",argv[1]);
  if((fp_results=fopen(filename,"w"))==NULL) {
    fprintf(stderr,"Error opening file %s\n",filename);
    exit(1);
  }
  fprintf(fp_results,"xcorr_%s\n",argv[1]);

  //main loop of comparison
  for (i=0;i<numsub-1;i++) {
    for (j=i+1;j<numsub;j++) {

      printf("Comparing %s and %s:\n",sub[i],sub[j]);

      temp = pair_matching(sub[i],sub[j],argv[1],xres,yres,zres,numsub);
      fprintf(fp_results,"%f\n",temp);

      //show info
      printf("xcorr_%s=%f\n",argv[1],temp);
      printf("---------------------------------------------\n");	
      
    }
    printf("\n************Subject %s done*************\n\n",sub[i]);
  }

  fclose(fp_results);
  return 0;
}



	
