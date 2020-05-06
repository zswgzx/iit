/*SZ  
revised from subjcou-v1.c, should replace it

structure:
read names of all subjects
allocate memory for 1 slice of data and result
open necessary files
for (each slice)
    read 1 slice of data from all subjects and meandyad_v1
    for (each voxel in the slice)
       calculate dot product between meandyad_v1 and v1 of each subject
       enforce cos boundary cond., find abs of cos and convert it to angle
       sort the angles and find % of cone of uncertainty in the voxel
    write 1 slice of result
clean up

Notice: 
1. cou should be between 0-90 degrees, so the absolute value of the inner product is of interest; also, the range of acos() is 0 to pi
2. use fabs instead of abs for getting float point absolute value;
3. when using qsort for float point values comparison, scale the returning results to significant level so it won't fail
*/ 

#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include <math.h>

#define PI 3.1415926535 
#define PCTG 95 //percentage for cone of uncertainty
#define FSCALE 1000 //float point scaling factor for compare function in qsort
#define FREE_ARG char*
#define NR_END 1
//#define HEAD 352

void nrerror(char error_text[])
/* Numerical Recipes standard error handler */
{
  fprintf(stderr,"Numerical Recipes run-time error...\n");
  fprintf(stderr,"%s\n",error_text);
  fprintf(stderr,"...now exiting to system...\n");
  exit(1);
}

float *fvector(long nl, long nh)
/* allocate a float vector with subscript range v[nl..nh] */
{
  float *v;

  v=(float *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(float)));
  if (!v) nrerror("allocation failure in fvector()");
  return v-nl+NR_END;
}

void free_fvector(float *v, long nl, long nh)
/* free a float vector allocated with fvector() */
{
  free((FREE_ARG) (v+nl-NR_END));
}

int compare (const void * a, const void * b)
{
  return (int)(FSCALE*( *(float *)a - *(float *)b ));
}
//////////////////////////////////////////////////////////////////////////
int main(int argc, char **argv)
{
  FILE *fp_sub,*fp_std;
  float *v1,*v2,*v3,*cu;
  char filename[200], subname[25];
  unsigned int  i, j, z,xres,yres,zres,subc = 0,numsub,pos;
  unsigned long size;
  const char *v1file = "reg%s_v1.img";

  if (argc!=5) {
    fprintf(stderr,"%s - calculate cone of uncertainty voxel by voxel from all subjects' primary eigenvectors\n\n",argv[0]);
    fprintf(stderr,"Input:meandyad_v1.img, reg[subj]_v1.img (volume by volume, float), subjects-??.txt\nOutput:cou.img (float)\n");
    fprintf(stderr,"\nUsage: %s <xdim> <ydim> <zdim> <# subjects>", argv[0]);
    fprintf(stderr,"\ne.g. %s 181 217 181 18\n", argv[0]);
    exit(1);
  }
 
  xres = atoi(argv[1]);
  yres = atoi(argv[2]);
  zres = atoi(argv[3]);
  numsub = atoi(argv[4]);
  size=xres*yres;
  pos=numsub*PCTG/100;

  FILE *fp_e1[numsub];
  float *e11[numsub],*e12[numsub],*e13[numsub],cos[numsub],data[numsub];
  char sub[numsub][25];

  //////////////////////////////////////////////////////////
  //open file that contains all subject names
  sprintf(filename,"subjects-%d.txt",numsub);
  if((fp_sub=fopen(filename,"r"))==NULL) {
    fprintf(stderr,"Error opening file %s\n",filename);
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
  //////////////////////////////////////////////////////////
  //allocation for vector
  for(i=0;i<numsub;i++) {
    e11[i]=fvector(1,size);
    e12[i]=fvector(1,size);
    e13[i]=fvector(1,size);
  }

  v1=fvector(1,size);
  v2=fvector(1,size);
  v3=fvector(1,size);
  cu=fvector(1,size);
 
  //open necessary files
  if((fp_std=fopen("cou.img","wb"))==NULL) {
    fprintf(stderr,"Error opening cou.img\n");
    exit(1);
  }

  if((fp_sub=fopen("meandyad_v1.img","rb"))==NULL) {
    fprintf(stderr,"Error opening meandyad_v1.img\n");
    exit(1);
  }

  for (j=0;j<numsub;j++) {
    sprintf(filename,v1file,sub[j]);
    if((fp_e1[j]=fopen(filename,"rb"))==NULL) {
      fprintf(stderr,"Error opening file %s\n",filename);
      exit(1);
    }
  }
  //////////////////////////////////////////////////////////
  //slice-wise operation
  for(z=0;z<zres;z++) {
    //read all samples' data of 1 slice
    for(i=0;i<numsub;i++) {
     
      fseek(fp_e1[i],z*size*sizeof(float),SEEK_SET);
      if(fread(e11[i],sizeof(float),size,fp_e1[i])!=size) {
	fprintf(stderr,"Error reading 1 from e1 file of subj%d, slice%d (1-offset).\n",i+1,z);
	exit(1);
      }

      fseek(fp_e1[i],sizeof(float)*((z+zres)*size),SEEK_SET);
      if(fread(e12[i],sizeof(float),size,fp_e1[i])!=size) {
	fprintf(stderr,"Error reading 2 from e1 file of subj%d, slice%d (1-offset).\n",i+1,z);
	exit(1);
      }

      fseek(fp_e1[i],sizeof(float)*((z+2*zres)*size),SEEK_SET);
      if(fread(e13[i],sizeof(float),size,fp_e1[i])!=size) {
	fprintf(stderr,"Error reading 3 from e1 file of subj%d, slice%d (1-offset).\n",i+1,z);
	exit(1);
      }
    }

    fseek(fp_sub,z*size*sizeof(float),SEEK_SET);
    if(fread(v1,sizeof(float),size,fp_sub)!=size) {
      fprintf(stderr,"Error reading 1 from dyad file slice%d (1-offset).\n",z);
      exit(1);
    }

    fseek(fp_sub,sizeof(float)*((z+zres)*size),SEEK_SET);
    if(fread(v2,sizeof(float),size,fp_sub)!=size) {
      fprintf(stderr,"Error reading 2 from dyad file slice%d (1-offset).\n",z);
      exit(1);
    }

    fseek(fp_sub,sizeof(float)*((z+2*zres)*size),SEEK_SET);
    if(fread(v3,sizeof(float),size,fp_sub)!=size) {
      fprintf(stderr,"Error reading 3 from dyad file slice%d (1-offset).\n",z);
      exit(1);
    }

    for(j=0;j<size;j++) {
      for(i=0;i<numsub;i++) {
	//a(.)b=||a||*||b||*cos(theta) - dot product
	cos[i]=e11[i][j]*v1[j]+e12[i][j]*v2[j]+e13[i][j]*v3[j];
	if (cos[i]) 
	  cos[i]/=sqrt(pow(e11[i][j],2)+pow(e12[i][j],2)+pow(e13[i][j],2))*sqrt(pow(v1[j],2)+pow(v2[j],2)+pow(v3[j],2));
	if (cos[i]>1 || cos[i]<-1)
	  cos[i]=((cos[i]>0)?1:-1); // enforcing cosine boundary condition
	//printf("impossible cosine value cos[%d]=%f at [%d,%d,%d] (0-offset)\n",i,cos[i],(j+1)%xres-1,1+j-xres*(int)((j+1)/xres),z-1);

	if (cos[i]!=0) {
	  cos[i]=fabs(cos[i]);
	  data[i]=acos(cos[i])*180.0/PI;
	  //if (cos[i]<0.03) data[i]=0;
	} else
	  data[i]=0.;
      }

      qsort(data,numsub,sizeof(float),compare);
      cu[j]=data[pos-1];
      /*
      //debug
      if (z==128 && j==127+xres*128) {
	printf("-------------------------------\n");
	for (i=0;i<numsub;i++) {
	  printf("e%d:%f %f %f\n",i+1,e11[i][j],e12[i][j],e13[i][j]);
	}
	printf("v=%f %f %f\n",v1[j],v2[j],v3[j]);
	for (i=0;i<numsub;i++) {
	  printf("cos[%d]=%f\tdata[%d]=%f\n",i,cos[i],i,data[i]);
	}
	printf("pos=%d\tcu[pos]=%f\n",pos,cu[j]);
      }
      */
    }

    if(fwrite(cu,sizeof(float),size,fp_std)!=size) {
      fprintf(stderr,"Error writing cou file %hd.\n",i);
      exit(1);
    }
    //printf("%d slice done\n", z);
  }

  for(i=0;i<numsub;i++) {
    fclose(fp_e1[i]);
      
    free_fvector(e11[i],1,size);
    free_fvector(e12[i],1,size);
    free_fvector(e13[i],1,size);  
  }

  free_fvector(v1,1,size);
  free_fvector(v2,1,size);
  free_fvector(v3,1,size);
  free_fvector(cu,1,size);
  fclose(fp_sub);
  fclose(fp_std);

  return 0;
}
