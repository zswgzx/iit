/* SZ 
revised from dted_logtsr-v1.c, should replace it and pair_match-v1.c

structure:
check # arguments
read setting file
store subjects' names
initialize chose scalar maps
compare each pair of subjects for similarity metrics (pair_matching)
average selected files

inside pair_matching:
open chosen scalar maps
open two subjects' tensor files and commonmask??.img
for (each voxel)
    if (in mask) diagonize both tensors
    read cohsen metrics
    calculate the new value of metrics in the voxel and update it in raw files
clean up

Output:  available indices are- Euclidean distance between tensors or deviatorics, dot product between tensors or deviatorics, OVL (overlap of eigval-eigvec pairs), cosine of angle between two PD weighted by geometric mean. 6 maps (all in float type, same dim) and 1 txt file are created and updated.

Notes: better split subjects in 20:47 to make the load even.
*/

#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include <math.h>

#define NR_END 1
#define FREE_ARG char*  
#define FG fgets(c120,120,fp_config)

float results[4]={0,0,0,0};

void nrerror(char error_text[])
/* Numerical Recipes standard error handler */
{
  fprintf(stderr,"Numerical Recipes run-time error...\n");
  fprintf(stderr,"%s\n",error_text);
  fprintf(stderr,"...now exiting to system...\n");
  exit(1);
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

float euclidean_dist(float (*xc1)[6], float (*xc2)[6])
{
  float sum=0.0,*p1,*p2;
  char i=0;

  p1=(float *)xc1;
  p2=(float *)xc2;
  for (;i<6;i++) {
    sum += (*p1-*p2)*(*p1-*p2)*((i<3)?1:2);
    ++p1;
    ++p2;
  }

  return sqrt(sum);
}
  
float tensor_dotproduct(float (*xc1)[6], float (*xc2)[6])
{
  float sum=0.0,*p1,*p2;
  char i=0;

  p1=(float *)xc1;
  p2=(float *)xc2;
  for (;i<6;i++) {
    sum += (*p1)*(*p2)*((i<3)?1:2);
    ++p1;
    ++p2;
  }

  return sum;
}

FILE *open_out_file(const char *file_prefix)
{
  FILE *fp;
  char filename[130];
  const char *format="%s.img";

  sprintf(filename,format,file_prefix);
  if((fp=fopen(filename,"rb+"))==NULL) {
    fprintf(stderr,"Error reading file %s for update\n",filename);
    exit(1);
  }
  return fp;
}

FILE *write_out_file(const char *file_prefix)
{
  FILE *fp;
  char filename[130];
  const char *format="%s.img";

  sprintf(filename,format,file_prefix);
  if((fp=fopen(filename,"wb"))==NULL) {
    fprintf(stderr,"Error creating file %s\n",filename);
    exit(1);
  }
  return fp;
}

FILE *open_in_file(const char *subname,const char *suffix)
{
  FILE *fp;
  char filename[130];
  const char *format="%s%s";

  sprintf(filename,format,subname,suffix);
  if((fp=fopen(filename,"rb"))==NULL) {
    fprintf(stderr,"Error opening file %s\n",filename);
    exit(1);
  }
  return fp;
}
/*--------matrix manipulations below---------*/
void matrix_mult_matrix(double **R,double **matrix1,double **matrix2,int dim)
{
  int i,j,k;
  for(i=1;i<=dim;i++)
    for(j=1;j<=dim;j++) {
      R[i][j]=0.0;
      for(k=1;k<=dim;k++)
	R[i][j]+=matrix1[i][k]*matrix2[k][j];
    }
}

void transpose_matrix(double **R,double **matrix,int dim)
{
  int i,j;
  for(i=1;i<=dim;i++)
    for(j=1;j<=dim;j++)
      R[i][j]=matrix[j][i];
}

void gettensor(double V1[],double V2[],double V3[],double L1, double L2,double L3,double **D)
{
  double **E,**L,**E_transpose,**D_temp;
  short i;

  E=matrix(1,3,1,3);
  L=matrix(1,3,1,3);
  E_transpose=matrix(1,3,1,3);
  D_temp=matrix(1,3,1,3);

  for(i=1;i<=3;i++) {
    E[i][1]=V1[i];
    E[i][2]=V2[i];
    E[i][3]=V3[i];
    L[i][1]=L[i][2]=L[i][3]=0.0;
  }
  L[1][1]=L1;
  L[2][2]=L2;
  L[3][3]=L3;

  transpose_matrix(E_transpose,E,3);
  matrix_mult_matrix(D_temp,E,L,3);
  matrix_mult_matrix(D,D_temp,E_transpose,3);

  free_matrix(E,1,3,1,3);
  free_matrix(L,1,3,1,3);
  free_matrix(E_transpose,1,3,1,3);
  free_matrix(D_temp,1,3,1,3);
}
/* for eigensystem calculation below */
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

float *pair_matching(char a[],char b[],int xcorr_flag,int dted_flag,int dtdp_flag,int dved_flag,int dvdp_flag,int ovl_flag,int cos_flag,int xres,int yres, int zres,int numsub,int logtsr_flag)
{
  FILE *fp_tensor[2],*fp_mask, *fp_dted, *fp_dved, *fp_ovl,*fp_count;//*fp_dtdp, *fp_dvdp,*fp_cos
  char filename[130],mask;//comm[100]
  unsigned int i,j,s,r,c;
  unsigned short scount;
  unsigned long count=0,pos;
  float tensor[2][6],l1[2],l2[2],l3[2],v1[2][3],v2[2][3],v3[2][3],deviatoric[2][6],trace[2],dted,dved,povl,sum=0,temp[6]={0,0,0,0,0,0},logtsr[2][6];//cos,dvdp,dtdp,logDeviaTsr[2][6],fa[2]
  extern float results[4];
  double d1[4],e1[4],**a1,vec1[4],vec2[4],vec3[4];

  // open chosen scalar maps
  if (dted_flag == 1) fp_dted=open_out_file("TensorDist");
  if (dved_flag == 1) fp_dved=open_out_file("DeviatoricDist");
  if (ovl_flag == 1) fp_ovl=open_out_file("OVL");
  fp_count=open_out_file("count");

  //open tensor files from subject a and b
  fp_tensor[0]=open_in_file(a,".img");
  fp_tensor[1]=open_in_file(b,".img");

  //open mask file
  sprintf(filename,"commonmask%d.img",numsub);
  if((fp_mask=fopen(filename,"rb"))==NULL) {
    fprintf(stderr,"Error opening file %s\n",filename);
    exit(1);
  }

  a1=matrix(1,3,1,3);
  // voxel based loop
  for(s=0;s<zres;s++)
    for(r=0;r<yres;r++)
      for(c=0;c<xres;c++) {
	pos=s*xres*yres+r*xres+c;
	//read mask
	//fseek(fp_mask,sizeof(char)*pos,SEEK_SET);
	if ((fread(&mask,sizeof(char),1,fp_mask))!=1) {
	  fprintf(stderr,"Error reading mask at [%d,%d,%d].\n",c,r,s);
	  exit(1);
	}

	// main comparison
	if (mask != 0) {
	  //read tensor
	  for (i=0;i<2;i++) {
	    fseek(fp_tensor[i],6*sizeof(float)*pos,SEEK_SET);
	    if((fread(&tensor[i][0],sizeof(float),6,fp_tensor[i]))!=6) {
	      fprintf(stderr,"Error reading tensorfile %d.\n",i+1);
	      exit(1);
	    }
	    
	    a1[1][1]=tensor[i][0];
	    a1[2][2]=tensor[i][1]; 
	    a1[3][3]=tensor[i][2];

	    a1[1][2]=a1[2][1]=tensor[i][3];
	    a1[1][3]=a1[3][1]=tensor[i][4];
	    a1[2][3]=a1[3][2]=tensor[i][5];

	    tred2(a1,3,d1,e1);
	    tqli(d1,e1,3,a1);  
	    eigsrt(d1,a1,3);

	    l1[i]=d1[1];
	    l2[i]=d1[2];
	    l3[i]=d1[3];

	    v1[i][0]=vec1[1]=a1[1][1];
	    v1[i][1]=vec1[2]=a1[2][1];
	    v1[i][2]=vec1[3]=a1[3][1];

	    v2[i][0]=vec2[1]=a1[1][2];
	    v2[i][1]=vec2[2]=a1[2][2];
	    v2[i][2]=vec2[3]=a1[3][2];

	    v3[i][0]=vec3[1]=a1[1][3];
	    v3[i][1]=vec3[2]=a1[2][3];
	    v3[i][2]=vec3[3]=a1[3][3];

	    if (l1[i]<=0 || l2[i]<=0 || l3[i]<=0)
	      printf("nonpositive diffusivity detected in subj %d [%3d,%3d,%3d] (0-offset)\n",i+1,c,r,s);
	    e1[1]=log(d1[1]);e1[2]=log(d1[2]);e1[3]=log(d1[3]);
	    gettensor(vec1,vec2,vec3,e1[1],e1[2],e1[3],a1);

	    logtsr[i][0]=a1[1][1];
	    logtsr[i][1]=a1[2][2];
	    logtsr[i][2]=a1[3][3];
	    logtsr[i][3]=a1[1][2];
	    logtsr[i][4]=a1[1][3];
	    logtsr[i][5]=a1[2][3];
	    /*
	    //debug
	    if (c==129 && r==109 && s==75) {
	      printf("tensor %d = %f %f %f %f %f %f\n",i+1, tensor[i][0],tensor[i][1],tensor[i][2],tensor[i][3],tensor[i][4],tensor[i][5]);
	      printf("log tensor %d = %f %f %f %f %f %f\n",i+1, logtsr[i][0],logtsr[i][1],logtsr[i][2],logtsr[i][3],logtsr[i][4],logtsr[i][5]);
	      printf("l1 = %f l2 = %f l3 = %f\n", l1[i],l2[i],l3[i]);
	      printf("v1 = %f %f %f \n", v1[i][0],v1[i][1],v1[i][2]);
	      printf("v2 = %f %f %f \n", v2[i][0],v2[i][1],v2[i][2]);
	      printf("v3 = %f %f %f \n", v3[i][0],v3[i][1],v3[i][2]);
	      printf("--------------------------------------------\n");
	    }
	    */
	    //calculate (log) deviatoric tensor (anisotropic part of raw tsr)
	    trace[i]=(tensor[i][0]+tensor[i][1]+tensor[i][2])/3.0;
	    for (j=0;j<6;j++)  deviatoric[i][j]=tensor[i][j]-((j<3)?trace[i]:0);
	    /*
	    //log deviatoric tsr
	    if (d1[1]-trace[i]<=0 || d1[2]-trace[i]<=0 || d1[3]-trace[i]<=0)
	      printf("nonpositive deviatoric part detected in subj %d [%3d,%3d,%3d] (0-offset)\n",i+1,c,r,s);
	    for (j=1;j<=3;j++) e1[j]=log(d1[j]-trace[i]);
	    gettensor(vec1,vec2,vec3,e1[1],e1[2],e1[3],a1);
	    logDeviaTsr[i][0]=a1[1][1];
	    logDeviaTsr[i][1]=a1[2][2];
	    logDeviaTsr[i][2]=a1[3][3];
	    logDeviaTsr[i][3]=a1[1][2];
	    logDeviaTsr[i][4]=a1[1][3];
	    logDeviaTsr[i][5]=a1[2][3];
	    */
	  }

	  //read data from scalor maps for update
	  if (dted_flag == 1) {
	    fseek(fp_dted,sizeof(float)*pos,SEEK_SET);
	    fread(&temp[0],sizeof(float),1,fp_dted);
	  }
	  if (dved_flag == 1) {
	    fseek(fp_dved,sizeof(float)*pos,SEEK_SET);
	    fread(&temp[2],sizeof(float),1,fp_dved);
	  }
	  if (ovl_flag == 1) {
	    fseek(fp_ovl,sizeof(float)*pos,SEEK_SET);
	    fread(&temp[5],sizeof(float),1,fp_ovl);
	  }
	      	      
	  //Euclidean dist between DT
	  if (!logtsr_flag) dted=euclidean_dist(tensor,(tensor+1));
	  else dted=euclidean_dist(logtsr,(logtsr+1));

	  //Euclidean dist between deviatoric
	  dved=euclidean_dist(deviatoric,(deviatoric+1));
	  //if (!logtsr_flag) dved=euclidean_dist(deviatoric,(deviatoric+1));
	  //else dved=euclidean_dist(logDeviaTsr,(logDeviaTsr+1));
	      
	  //calculate fa from eigenvalues
	  //for (i=0;i<2;i++)  fa[i]=sqrt(1.5*(pow((l1[i]-trace[i]),2)+pow((l2[i]-trace[i]),2)+pow((l3[i]-trace[i]),2))/(pow(l1[i],2)+pow(l2[i],2)+pow(l3[i],2)));

	  //nonpositive eigenvalue is NOT acceptable!
	  if ((l1[0]>0 && l2[0]>0 && l3[0]>0)&&(l1[1]>0 && l2[1]>0 && l3[1]>0)) {
	    //partial OVL
	    povl=(l1[0]*l1[1]*pow((v1[0][0]*v1[1][0]+v1[0][1]*v1[1][1]+v1[0][2]*v1[1][2]),2) +l2[0]*l2[1]*pow((v2[0][0]*v2[1][0]+v2[0][1]*v2[1][1]+v2[0][2]*v2[1][2]),2) +l3[0]*l3[1]*pow((v3[0][0]*v3[1][0]+v3[0][1]*v3[1][1]+v3[0][2]*v3[1][2]),2)) / (l1[0]*l1[1]+l2[0]*l2[1]+l3[0]*l3[1]);
	    sum += povl;
	    count += 1;
		  
	    // update count map
	    fseek(fp_count,sizeof(short)*pos,SEEK_SET);
	    if ((fread(&scount,sizeof(short),1,fp_count))!=1) {
	      fprintf(stderr,"Error reading mask.\n");
	      exit(1);
	    }
	    scount += 1;
	    fseek(fp_count,sizeof(short)*pos,SEEK_SET);
	    fwrite(&scount,sizeof(short),1,fp_count);
	  } else {
	    //printf("Fail to get cos and ovl at (row=%d,col= %d,slice= %d). (Negative diffusivity detected)\n",r,c,s);
	    //cos=0;
	    povl=0;
	  }

	  dted += temp[0];
	  dved += temp[2];
	  povl += temp[5];
	      
	  // write data in scalor map for correspond measure
	  if (dted_flag == 1) {
	    fseek(fp_dted,sizeof(float)*pos,SEEK_SET);
	    fwrite(&dted,sizeof(float),1,fp_dted);
	  }
	  if (dved_flag == 1) {
	    fseek(fp_dved,sizeof(float)*pos,SEEK_SET);
	    fwrite(&dved,sizeof(float),1,fp_dved);
	  }
	  if (ovl_flag == 1) {
	    fseek(fp_ovl,sizeof(float)*pos,SEEK_SET);
	    fwrite(&povl,sizeof(float),1,fp_ovl);
	  }
	      
	}

	for (j=0;j<6;j++) temp[j]=0;
	/*
        // for debug use
	if (c==125 && r==120 && s==75) {
	  for (i=0;i<2;i++) {
	    printf("tensor %d = %f %f %f %f %f %f\n",i+1, tensor[i][0],tensor[i][1],tensor[i][2],tensor[i][3],tensor[i][4],tensor[i][5]);
	    printf("l1 = %f l2 = %f l3 = %f\n", l1[i],l2[i],l3[i]);
	    printf("v1 = %f %f %f \n", v1[i][0],v1[i][1],v1[i][2]);
	    printf("v2 = %f %f %f \n", v2[i][0],v2[i][1],v2[i][2]);
	    printf("v3 = %f %f %f \n", v3[i][0],v3[i][1],v3[i][2]);
	    printf("fa = %f\n", fa[i]);
	    printf("deviatoric %d = %f %f %f %f %f %f\n",i+1, deviatoric[i][0],deviatoric[i][1],deviatoric[i][2],deviatoric[i][3],deviatoric[i][4],deviatoric[i][5]);
	    printf("trace = %f\n", 3*trace[i]);
	    printf("--------------------------------------------\n");
	  }

	  printf("tensor dist= %f, dot product = %f\n",dted,dtdp);
	  printf("deviatoric dist= %f, dot product = %f\n",dved,dvdp);
	  printf("cos = %f\n", cos);
	}
	*/
      }
  //printf("sum=%f, vox # = %ld, OVL = %f\n",sum,count,sum/count);

  // clean up file pointers
  for (i=0;i<2;i++) fclose(fp_tensor[i]);
  fclose(fp_mask);
  fclose(fp_count);

  if (dted_flag == 1) fclose(fp_dted);
  if (dved_flag == 1) fclose(fp_dved);
  if (ovl_flag == 1) fclose(fp_ovl);

  //return pointers to results
  if (xcorr_flag == 1) results[3] = sum/count;
  free_matrix(a1,1,3,1,3);
  return results;
}
//********************************************************************************//
int main(int argc,char **argv)
{ 
  char filename[180], c120[120],subname[20],subj1[120],subj2[120];//comm[120]
  int xcorr_flag,dted_flag,dtdp_flag,dved_flag,dvdp_flag,ovl_flag,cos_flag,logtsr_flag;
  FILE  *fp_sub,*fp_config,*fp_dted, *fp_dved, *fp_ovl,*fp_count,*fp_avgdted,*fp_avgdved,*fp_avgovl;//*fp_dtdp,*fp_dvdp,*fp_cos,*fp_results
  int s,r,c,i,j,subc=0,xres,yres,zres,numsub;
  short szero=0,scount;
  float temp[4],*pt,dtedsum,dvedsum,ovlsum,avg,*zero_vol;// zero=0.0
  //const char *path_txt="xcorr_ovl";

  // assert number of arguments
  if (argc != 7) {
    fprintf(stderr, "%s- calculate average value of similarity metrics such as DTED DVED OVL for all possible pair of subjects\n",argv[0]);
    fprintf(stderr, "\nInput:subjects-?.txt, setting, commonmask?.img (char), reg[subj].img (tensor file, float, DTIGUI format)\nOutput: user-dependent, may include any of avgdted.img, avgdved.img, avgovl.img, or avglogdted.img (all float)\n\n");
    fprintf(stderr, "Usage: %s <xdim-col> <ydim-row> <zdim-slice> <# of subjects> <log Euclidean tsr? 0:no, 1:yes> <setting txt>\n",argv[0]);
    fprintf(stderr, "e.g. %s  256 256 256 22 0 setting\n",argv[0]);
    fprintf(stderr, "  or %s  256 256 256 22 1 setting\n",argv[0]);
    exit(1);
  }

  xres=atoi(argv[1]);
  yres=atoi(argv[2]); 
  zres=atoi(argv[3]);
  numsub=atoi(argv[4]);
  logtsr_flag=atoi(argv[5]);
  if (logtsr_flag) printf("log Euclidean tensor mode activated...\n");

  char sub[numsub][20];
  ////////////////////////////////////////////////////////// 
  //read in information from the config file 
  if ((fp_config=fopen(argv[6],"r")) == NULL) {
    fprintf(stderr,"Error opening setting file\n");
    exit(1);
  }

  fscanf(fp_config,"%ud",&xcorr_flag);FG;
  fscanf(fp_config,"%ud",&dted_flag);FG;
  fscanf(fp_config,"%ud",&dtdp_flag);FG;
  fscanf(fp_config,"%ud",&dved_flag);FG;
  fscanf(fp_config,"%ud",&dvdp_flag);FG;
  fscanf(fp_config,"%ud",&ovl_flag);FG;
  fscanf(fp_config,"%ud",&cos_flag);FG;
  //fscanf(fp_config,"%hd",&anatomicaltype);FG;

  fclose(fp_config);
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
  ////////////////////////////////////////////////////////
  //memory allocation for a volume of zeros
  zero_vol=(float *)malloc(xres*yres*zres*sizeof(float));
  for(s=0;s<zres;s++)
    for(r=0;r<yres;r++)
      for(c=0;c<xres;c++)
	zero_vol[s*xres*yres+r*xres+c]=0.;
  ////////////////////////////////////////////////////////
  // similarity comparison between each two subjs
  //create chosen scalar maps
  if (dted_flag == 1) {
    fp_dted=write_out_file("TensorDist");
    if ((fwrite(zero_vol,sizeof(float),xres*yres*zres,fp_dted))!=xres*yres*zres) {
      fprintf(stderr,"Fail to initialize TensorDist.img\n");
      exit(1);
    }
    fclose(fp_dted);
  }

  if (dved_flag == 1) {
    fp_dved=write_out_file("DeviatoricDist");
    if ((fwrite(zero_vol,sizeof(float),xres*yres*zres,fp_dved))!=xres*yres*zres) {
      fprintf(stderr,"Fail to initialize DeviatoricDist.img\n");
      exit(1);
    }
    fclose(fp_dved);
  }

  if (ovl_flag == 1) {
    fp_ovl=write_out_file("OVL");
    if ((fwrite(zero_vol,sizeof(float),xres*yres*zres,fp_ovl))!=xres*yres*zres) {
      fprintf(stderr,"Fail to initialize OVL.img\n");
      exit(1);
    }
    fclose(fp_ovl);
  }

  if((fp_count=fopen("count.img","wb"))==NULL) {
    fprintf(stderr,"Error opening file count.img\n");
    exit(1);
  }

  for(s=1;s<=zres;s++)
    for(r=1;r<=yres;r++)
      for(c=1;c<=xres;c++)
	fwrite(&szero,sizeof(short),1,fp_count);
  fclose(fp_count);
  free(zero_vol);

  //main loop of comparison
  for (i=0;i<numsub-1;i++) {
    for (j=i+1;j<numsub;j++) {

      printf("Comparing %s and %s:\n",sub[i],sub[j]);
      sprintf(subj1,"reg%s",sub[i]);
      sprintf(subj2,"reg%s",sub[j]);
      pt = pair_matching(subj1,subj2,xcorr_flag,dted_flag,dtdp_flag,dved_flag,dvdp_flag,ovl_flag,cos_flag,xres,yres,zres,numsub,logtsr_flag);
	  
      if (xcorr_flag ==1) {
	temp[0] = *pt;
	temp[1] = *(pt+1);
	temp[2] = *(pt+2);
	temp[3] = *(pt+3);
	//fprintf(fp_results,"%f   %f   %f   %f\n",temp[0],temp[1],temp[2],temp[3]);
	//show info
	printf("ovl=%f\n",temp[3]);
      }
    }
  }
	  
  //if (xcorr_flag ==1)  fclose(fp_results);

  //average selected files by dividing count.img
  if (dted_flag == 1) {
    if (logtsr_flag) fp_avgdted=write_out_file("avglogdted");
    else fp_avgdted=write_out_file("avgdted");
    fp_dted=open_in_file("TensorDist",".img");
  }
 
  if (dved_flag == 1) {
    //if (logtsr_flag) fp_avgdved=write_out_file("avglogdved");
    //else fp_avgdted=write_out_file("avgdved");
    fp_avgdved=write_out_file("avgdved");
    fp_dved=open_in_file("DeviatoricDist",".img");
  }
 
  if (ovl_flag == 1) {
    fp_avgovl=write_out_file("avgovl");
    fp_ovl=open_in_file("OVL",".img");
  }

  fp_count=open_in_file("count",".img");
  printf("averaging maps...\n");

  for(s=1;s<=zres;s++)
    for(r=1;r<=yres;r++)
      for(c=1;c<=xres;c++) {
	fread(&scount,sizeof(short),1,fp_count);
	if (dted_flag == 1) {
	  fread(&dtedsum,sizeof(float),1,fp_dted);
	  avg=scount?(dtedsum/scount):0.;
	  fwrite(&avg,sizeof(float),1,fp_avgdted);
	}

	if (dved_flag == 1) {
	  fread(&dvedsum,sizeof(float),1,fp_dved);
	  avg=scount?(dvedsum/scount):0.;
	  fwrite(&avg,sizeof(float),1,fp_avgdved);
	}

	if (ovl_flag == 1) {
	  fread(&ovlsum,sizeof(float),1,fp_ovl);
	  avg=scount?(ovlsum/scount):0.;
	  fwrite(&avg,sizeof(float),1,fp_avgovl);
	}
      }
    
  system("rm count.img");
  if (dted_flag==1) {
    fclose(fp_dted);
    fclose(fp_avgdted);
    system("rm TensorDist.img");
  }
  if (dved_flag==1) {
    fclose(fp_dved);
    fclose(fp_avgdved);
    system("rm DeviatoricDist.img");
  }
  if (ovl_flag==1) {
    fclose(fp_ovl);
    fclose(fp_avgovl);
    system("rm OVL.img");
  }
  fclose(fp_count);

  return 0;
}



	
