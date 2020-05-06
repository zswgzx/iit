/* SZ last updated: 150810
revised from scalarmap-std-v1.c, should replace it
*/

#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include <math.h>

#define HEADER 352

float stdev(float *xc, int num)
{
  //find the standard deviation of a 1D array, ignore 0 in it!
  float avg,stdv,sum=0.0;
  int i,count=0;

  //find average
  for (i=0;i<num;i++) {
    if (xc[i]!=0) {
      count+=1;
      sum += xc[i];
    }
  }
  if (count!=0) {
    if (count==1) stdv=0.0; //exclude the situation with only 1 nonzero value
    else {
      avg=sum/count;
      sum=0;
      //sum of square
      for (i=0;i<num;i++) {
	if (xc[i]!=0)      sum += pow(xc[i]-avg,2);
      }
      stdv = sqrt(sum/(count-1));
    }
  }
  else stdv=0.0;

  //printf("avg= %f,sum=%f,count=%d ",avg,sum,count); //for debug use
  return stdv;
}

//******************************************************************************//

int main(int argc, char **argv)
{ 
  char filename[180],subname[20];
  FILE  *fp_sub,*fp_results;
  int xres,yres,zres,numsub,hdr_flag,s,r,c,i,subc=0,checkseek;
  float stdv;

  // assert number of arguments
  if (argc != 8) {
    fprintf(stderr, "%s - standard dev. of scalar maps(fa, trace, etc.)\n\nInput: subjects-??.txt,scalar maps (float)\nOutput: std.img (float)\n\n",argv[0]);
    fprintf(stderr, "Usage: %s <xdim> <ydim> <zdim> <subjects#> <header_flag,1 or 0, exist or not> <filename_prefix> <filename_suffix>\n",argv[0]);
    fprintf(stderr, "e.g. %s 181 217 181 22 1 reg -fa.nii \n",argv[0]);
    fprintf(stderr, "  or %s 181 217 181 22 0 reg -fa.img \n",argv[0]);
    exit(1);
  }

  xres=atoi(argv[1]);
  yres=atoi(argv[2]);
  zres=atoi(argv[3]);
  numsub=atoi(argv[4]);
  hdr_flag=atoi(argv[5]);

  FILE *fp[numsub];
  float data[numsub];
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
  if((fp_results=fopen("std.img","wb"))==NULL) {
    fprintf(stderr,"Error opening file std.img\n");
    exit(1);
  }
   
  for (i=0;i<numsub;i++) {
    sprintf(filename,"%s%s%s",argv[6],sub[i],argv[7]);
    if((fp[i]=fopen(filename,"rb"))==NULL) {
      fprintf(stderr,"Error opening file %s\n",filename);
      exit(1);
    }
  }

  if (hdr_flag == 1) {
    for (i=0;i<numsub;i++) {
      checkseek=fseek(fp[i],HEADER,SEEK_SET);
      if(checkseek!=0) {
	fprintf(stderr,"fseek error for %d.\n",i+1);
	exit(1);
      }
    }
  }

  //vox-wise main loop
  for(s=1;s<=zres;s++)
    for(r=1;r<=yres;r++)
      for(c=1;c<=xres;c++) {
	for (i=0;i<numsub;i++) fread(&data[i],sizeof(float),1,fp[i]);
	stdv=stdev(data,numsub);
	fwrite(&stdv,1,sizeof(float),fp_results);
      }

  for (i=0;i<numsub;i++) fclose(fp[i]);
  fclose(fp_results);

  return 0;
}



	
