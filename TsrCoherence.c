/*SZ
revised from HP's coherence.c, should replace coh4img_v1.c
Notice: (NaN is observed sometimes at the edge of the output image, if beta2+beta3<0) k=1 voxels should be masked out because it means either all images allign perfect(rare) or just 1 of them is picked and all others are ignored.

structure:
get subjects' names and open corresponding v1 files
open output files and allocate memories
for (each slice)
    read vectors from the slice of all subjects
    for (each voxel in the slice)
        calculate dyadic tensor and diagonize it to find coherence (Jones, 02)
	count # of subjects that are involved in calculation
    write coherence and count
clean up
*/ 

#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include <math.h>

#define FREE_ARG char*
#define NR_END 1
#define TINY 1e-7
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

double *vector(long nl, long nh)
/* allocate a double vector with subscript range v[nl..nh] */
{
  double *v;

  v=(double *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(double)));
  if (!v) nrerror("allocation failure in vector()");
  return v-nl+NR_END;
}

void free_vector(double *v, long nl, long nh)
/* free a double vector allocated with vector() */
{
  free((FREE_ARG) (v+nl-NR_END));
}

double **matrix(long nrl, long nrh, long ncl, long nch)
/* allocate a double matrix with subscript range m[nrl..nrh][ncl..nch] */
{
  long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
  double **m;

  /* allocate pointers to rows */
  m=(double **) malloc((size_t)((nrow+NR_END)*sizeof(double*)));
  if (!m) nrerror("allocation failure 1 in matrix()");
  m += NR_END;
  m -= nrl;
  
  /* allocate rows and set pointers to them */
  m[nrl]=(double *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(double)));
  if (!m[nrl]) nrerror("allocation failure 2 in matrix()");
  m[nrl] += NR_END;
  m[nrl] -= ncl;
  
  for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

  /* return pointer to array of pointers to rows */
  return m;
}

void free_matrix(double **m, long nrl, long nrh, long ncl, long nch)
/* free a double matrix allocated by matrix() */
{
  free((FREE_ARG) (m[nrl]+ncl-NR_END));
  free((FREE_ARG) (m+nrl-NR_END));
}

double determine_vector(double *vector,int dim)
{
  int i;
  double det=0.0;
  for(i=1;i<=dim;i++)
    det=det+vector[i]*vector[i];
  
  return (double)(sqrt(det));
}
//eigen system calculation processes below
static double sqrarg;
#define SQR(a) ((sqrarg=(a)) == 0.0 ? 0.0 : sqrarg*sqrarg)

double pythag(double a, double b)
{
  double absa,absb;
  absa=fabs(a);
  absb=fabs(b);
  if (absa > absb) return absa*sqrt(1.0+SQR(absb/absa));
  else return (absb == 0.0 ? 0.0 : absb*sqrt(1.0+SQR(absa/absb)));
}

#undef SQR

void tred2(double **a, int n, double d[], double e[])
{
  int l,k,j,i;
  double scale,hh,h,g,f;

  for (i=n;i>=2;i--) {
    l=i-1;
    h=scale=0.0;
    if (l > 1) {
      for (k=1;k<=l;k++)
	scale += fabs(a[i][k]);
      if (scale == 0.0)
	e[i]=a[i][l];
      else {
	for (k=1;k<=l;k++) {
	  a[i][k] /= scale;
	  h += a[i][k]*a[i][k];
	}
	f=a[i][l];
	g=(f >= 0.0 ? -sqrt(h) : sqrt(h));
	e[i]=scale*g;
	h -= f*g;
	a[i][l]=f-g;
	f=0.0;
	for (j=1;j<=l;j++) {
	  a[j][i]=a[i][j]/h;
	  g=0.0;
	  for (k=1;k<=j;k++)
	    g += a[j][k]*a[i][k];
	  for (k=j+1;k<=l;k++)
	    g += a[k][j]*a[i][k];
	  e[j]=g/h;
	  f += e[j]*a[i][j];
	}
	hh=f/(h+h);
	for (j=1;j<=l;j++) {
	  f=a[i][j];
	  e[j]=g=e[j]-hh*f;
	  for (k=1;k<=j;k++)
	    a[j][k] -= (f*e[k]+g*a[i][k]);
	}
      }
    } else
      e[i]=a[i][l];
    d[i]=h;
  }
  d[1]=0.0;
  e[1]=0.0;
  /* Contents of this loop can be omitted if eigenvectors not
     wanted except for statement d[i]=a[i][i]; */
  for (i=1;i<=n;i++) {
    l=i-1;
    if (d[i]) {
      for (j=1;j<=l;j++) {
	g=0.0;
	for (k=1;k<=l;k++)
	  g += a[i][k]*a[k][j];
	for (k=1;k<=l;k++)
	  a[k][j] -= g*a[k][i];
      }
    }
    d[i]=a[i][i];
    a[i][i]=1.0;
    for (j=1;j<=l;j++) a[j][i]=a[i][j]=0.0;
  }
}

#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))
//static double maxarg1,maxarg2;

void tqli(double d[], double e[], int n, double **z)
{
  double pythag(double a, double b);
  int m,l,iter,i,k;
  double s,r,p,g,f,dd,c,b;

  for (i=2;i<=n;i++) e[i-1]=e[i];
  e[n]=0.0;
  for (l=1;l<=n;l++) {
    iter=0;
    do {
      for (m=l;m<=n-1;m++) {
	dd=fabs(d[m])+fabs(d[m+1]);
	if ((double)(fabs(e[m])+dd) == dd) break;
      }
      if (m != l) {
	if (iter++ == 30) nrerror("Too many iterations in tqli");
	g=(d[l+1]-d[l])/(2.0*e[l]);
	r=pythag(g,1.0);
	g=d[m]-d[l]+e[l]/(g+SIGN(r,g));
	s=c=1.0;
	p=0.0;
	for (i=m-1;i>=l;i--) {
	  f=s*e[i];
	  b=c*e[i];
	  e[i+1]=(r=pythag(f,g));
	  if (r == 0.0) {
	    d[i+1] -= p;
	    e[m]=0.0;
	    break;
	  }
	  s=f/r;
	  c=g/r;
	  g=d[i+1]-p;
	  r=(d[i]-g)*s+2.0*c*b;
	  d[i+1]=g+(p=s*r);
	  g=c*r-b;
	  for (k=1;k<=n;k++) {
	    f=z[k][i+1];
	    z[k][i+1]=s*z[k][i]+c*f;
	    z[k][i]=c*z[k][i]-s*f;
	  }
	}
	if (r == 0.0 && i >= l) continue;
	d[l] -= p;
	e[l]=g;
	e[m]=0.0;
      }
    } while (m != l);
  }
}

#undef SIGN

void eigsrt(double d[], double **v, int n)
{
  int k,j,i;
  double p;
  
  for (i=1;i<n;i++) {
    p=d[k=i];
    for (j=i+1;j<=n;j++)
      if (d[j] >= p) p=d[k=j];
    if (k != i) {
      d[k]=d[i];
      d[i]=p;
      for (j=1;j<=n;j++) {
	p=v[j][i];
	v[j][i]=v[j][k];
	v[j][k]=p;
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////
int main(int argc, char **argv)
{
  FILE *fp_sub,*fp_coherence,*fp_count,*fp_meandyad;
  float *k,*meandyad1,*meandyad2,*meandyad3;
  char filename[100], subname[60],*subcount,count;
  unsigned int  i, j, z, subc=0,xres,yres,zres,numsub,size;
  double **dyadic,*d1,*e1,*e_temp,beta1, beta2, beta3;
  const char *v1file = "reg%s_v1.img";
	
  // assert number of arguments
  if (argc != 5) {
    fprintf(stderr, "%s -calculate coherence from v1 maps (not masked) and primary eigenvector of mean dyadic tensor (volume by volume, float)\n",argv[0]);
    fprintf(stderr, "\nInput: subjects-??.txt,reg[subj]_v1.img\nOutput: coherence.img (float), cohcount.img (char), meandyad_v1.img\n");
    fprintf(stderr, "\nUsage: %s <xdim> <ydim> <zdim> <#subjects>\n",argv[0]);
    fprintf(stderr, "e.g. %s 256 256 256 22\n",argv[0]);
    exit(1);
  }
 
  xres=atoi(argv[1]);
  yres=atoi(argv[2]);
  zres=atoi(argv[3]);
  numsub=atoi(argv[4]);
  size=xres*yres;

  float *e11[numsub],*e12[numsub],*e13[numsub];
  FILE *fp_e1[numsub];
  double det[numsub];
  //////////////////////////
  // get names of subjects from subjects.txt
  sprintf(filename,"subjects-%d.txt",numsub);
  if((fp_sub=fopen(filename,"r"))==NULL) {
    fprintf(stderr,"Error opening file subjects.txt\n");
    exit(1);
  }

  // open subject e1 volumes
  while (fgets(subname, sizeof(subname), fp_sub)) {
    // cut newline from the end of subname
    char *tmp = subname;
    while (*tmp) {
      if (*tmp == '\n') *tmp = '\0';
      tmp++;
    }
    sprintf(filename, v1file, subname);
    if((fp_e1[subc]=fopen(filename,"rb"))==NULL) {
      fprintf(stderr,"Error opening file %s\n", filename);
      exit(1);
    }

    if (++subc == numsub) break;
  }	 
  /////////////////////////
  // open volume for writing coherence
  if((fp_coherence=fopen("coherence.img","wb"))==NULL) {
    fprintf(stderr,"Error opening file coherence.img.\n");
    exit(1);
  }

  if((fp_count=fopen("cohcount.img","wb"))==NULL) {
    fprintf(stderr,"Error opening file cohcount.img.\n");
    exit(1);
  }

  sprintf(filename,"meandyad_v1.img");
    if((fp_meandyad=fopen(filename,"wb"))==NULL) {
      fprintf(stderr,"Error opening file %s.\n",filename);
      exit(1);
    }
 
  for(i=0;i<numsub;i++) {
    e11[i]=fvector(1,size);
    e12[i]=fvector(1,size);
    e13[i]=fvector(1,size);
  }
  e_temp=vector(1,3);
  d1=vector(1,3);
  e1=vector(1,3);
  dyadic=matrix(1,3,1,3);

  subcount=(char *)malloc(sizeof(char)*size);
  if (subcount==NULL) {
    printf("subcount allocation failure\n");
    exit(1);
  }

  k=(float *)malloc(sizeof(float)*size);
  if (k==NULL) {
    printf("k allocation failure\n");
    exit(1);
  }

  // allocate buffer to store slice data
  meandyad1 = malloc(sizeof(float)*size*zres);
  if (meandyad1 == NULL) {
    printf("*** Failed to allocate meandyad buffer 1\n");
    exit(1);
  }

  meandyad2 = malloc(sizeof(float)*size*zres);
  if (meandyad2 == NULL) {
    printf("*** Failed to allocate meandyad buffer 2\n");
    exit(1);
  }

  meandyad3 = malloc(sizeof(float)*size*zres);
  if (meandyad3 == NULL) {
    printf("*** Failed to allocate meandyad buffer 3\n");
    exit(1);
  }

  //////////////////////////
  //slice-wise loop
  for(z=1;z<=zres;z++) {
    //read 1 slice of vectors from all subj
    for(i=0;i<numsub;i++) {
      fseek(fp_e1[i],(z-1)*size*sizeof(float),SEEK_SET);
      if(fread(e11[i],sizeof(float),size,fp_e1[i])!=size) {
	fprintf(stderr,"Error reading 1 from e1 file of subj%d slice%d (1-offset).\n",i,z);
	exit(1);
      }

      fseek(fp_e1[i],(z-1)*size*sizeof(float) + zres*size*sizeof(float),SEEK_SET);
      if(fread(e12[i],sizeof(float),size,fp_e1[i])!=size) {
	fprintf(stderr,"Error reading 2 from e1 file of subj%d slice%d (1-offset).\n",i,z);
	exit(1);
      }

      fseek(fp_e1[i],(z-1)*size*sizeof(float) + 2*zres*size*sizeof(float),SEEK_SET);
      if(fread(e13[i],sizeof(float),size,fp_e1[i])!=size) {
	fprintf(stderr,"Error reading 2 from e1 file of subj%d slice%d (1-offset).\n",i,z);
	exit(1);
      }
    }

    //voxel-wise in-plane processes
    for(j=0;j<size;j++) {
      dyadic[1][1]=dyadic[2][2]=dyadic[3][3]=dyadic[1][2]=dyadic[1][3]=dyadic[2][3]=dyadic[2][1]=dyadic[3][2]=dyadic[3][1]=0.0;
      count=0;
      k[j]=0.;

      for(i=0;i<numsub;i++) {
	dyadic[1][1]+=e11[i][j]*e11[i][j];
	dyadic[2][2]+=e12[i][j]*e12[i][j];
	dyadic[3][3]+=e13[i][j]*e13[i][j];
	dyadic[1][2]+=e11[i][j]*e12[i][j];
	dyadic[1][3]+=e11[i][j]*e13[i][j];
	dyadic[2][3]+=e12[i][j]*e13[i][j];

	e_temp[1]=(double)e11[i][j];
	e_temp[2]=(double)e12[i][j];
	e_temp[3]=(double)e13[i][j];
	det[i]=determine_vector(e_temp,3);
	if (det[i] != 0) count+=1;
      }

      if(count!=0) {
	dyadic[1][1]/=count;
	dyadic[2][2]/=count;
	dyadic[3][3]/=count;
	dyadic[1][2]/=count;
	dyadic[1][3]/=count;
	dyadic[2][3]/=count;
	dyadic[2][1]=dyadic[1][2];
	dyadic[3][1]=dyadic[1][3];
	dyadic[3][2]=dyadic[2][3];
	      
	tred2(dyadic,3,d1,e1);
	tqli(d1,e1,3,dyadic);
	eigsrt(d1,dyadic,3);
	      
	beta1=d1[1];
	beta2=d1[2];
	beta3=d1[3];

	if (beta2<TINY && beta2>-TINY) beta2=0;
	if (beta3<TINY && beta3>-TINY) beta3=0;

	k[j]=(float)(1.0-sqrt((beta2+beta3)/(2.0*beta1)));
	meandyad1[j+(z-1)*size]=(float)dyadic[1][1];
	meandyad2[j+(z-1)*size]=(float)dyadic[2][1];
	meandyad3[j+(z-1)*size]=(float)dyadic[3][1];
      }
     else {
       k[j]=meandyad1[j+(z-1)*size]=meandyad2[j+(z-1)*size]=meandyad3[j+(z-1)*size]=0.0;
     }
      subcount[j]=count;
      /*
      //debug
      if (z==128+1 && j==127+xres*128) {
       for (i=0;i<numsub;i++) printf("e%d=[%f\t%f\t%f]\n",i+1,e11[i][j],e12[i][j],e13[i][j]);
       printf("count=%d\n",count);
       printf("beta1,2,3=%f\t%f\t%f\n",beta1,beta2,beta3);
       printf("k=%f\n",k[j]);
       printf("meandyad_v1=[%f\t%f\t%f]\n",meandyad1[j+(z-1)*size],meandyad2[j+(z-1)*size],meandyad3[j+(z-1)*size]);
      }	
      */
    } 

    if ((fwrite(k,sizeof(float),size,fp_coherence))!=size) {
      fprintf(stderr,"Error writing slice%d (1-offset) to cohrence.img.\n",z);
      exit(1);
    }

    if ((fwrite(subcount,sizeof(char),size,fp_count))!=size) {
      fprintf(stderr,"Error writing slice%d (1-offset) to cohcount.img.\n",z);
      exit(1);
    }
  }

  if ((fwrite(meandyad1,sizeof(float),size*zres,fp_meandyad))!=size*zres) {
    fprintf(stderr,"Error writing 1st volume to meandyad_v1.img.\n");
    exit(1);
  }

  if ((fwrite(meandyad2,sizeof(float),size*zres,fp_meandyad))!=size*zres) {
    fprintf(stderr,"Error writing 2nd volume to meandyad_v1.img.\n");
    exit(1);
  }

  if ((fwrite(meandyad3,sizeof(float),size*zres,fp_meandyad))!=size*zres) {
    fprintf(stderr,"Error writing 3rd volume to meandyad_v1.img.\n");
    exit(1);
  }

  for(i=0;i<numsub;i++) {
    fclose(fp_e1[i]);
      
    free_fvector(e11[i],1,size);
    free_fvector(e12[i],1,size);
    free_fvector(e13[i],1,size);
  }

  fclose(fp_meandyad);
  free(meandyad1);
  free(meandyad2);
  free(meandyad3);

  free_vector(e_temp,1,3);
  free_matrix(dyadic,1,3,1,3);
  free_vector(d1,1,3);
  free_vector(e1,1,3);
  free(subcount);
  free(k);

  fclose(fp_coherence);
  fclose(fp_sub);
  fclose(fp_count);

  return 0;
}
