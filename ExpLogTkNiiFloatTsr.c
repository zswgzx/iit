/* SZ */
#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include <math.h>

#define NR_END 1
#define FREE_ARG char*  
#define HDR 352

void nrerror(char error_text[])
/* Numerical Recipes standard error handler */
{
        fprintf(stderr,"Numerical Recipes run-time error...\n");
        fprintf(stderr,"%s\n",error_text);
        fprintf(stderr,"...now exiting to system...\n");
        exit(1);
}

double *vector(long nl, long nh)
/* allocate a double vector with subscript range v[nl..nh] */
{
        double *v;

        v=(double *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(double)));
        if (!v) nrerror("allocation failure in vector()");
        return v-nl+NR_END;
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

void free_vector(double *v, long nl, long nh)
/* free a double vector allocated with vector() */
{
        free((FREE_ARG) (v+nl-NR_END));
}

void free_matrix(double **m, long nrl, long nrh, long ncl, long nch)
/* free a double matrix allocated by matrix() */
{
        free((FREE_ARG) (m[nrl]+ncl-NR_END));
        free((FREE_ARG) (m+nrl-NR_END));
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
/*--------matrix manipulations below---------*/
void matrix_mult_matrix(double **R,double **matrix1,double **matrix2,int dim)
{
  int i,j,k;
  for(i=1;i<=dim;i++)
    for(j=1;j<=dim;j++)
	R[i][j]=0.0;
  for(i=1;i<=dim;i++)
    for(j=1;j<=dim;j++)
      for(k=1;k<=dim;k++)
	R[i][j]=R[i][j]+matrix1[i][k]*matrix2[k][j];
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

  for(i=1;i<=3;i++)
    {
      E[i][1]=V1[i];
      E[i][2]=V2[i];
      E[i][3]=V3[i];
      L[i][1]=0.0;
      L[i][2]=0.0;
      L[i][3]=0.0;
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
//************************************************************************************************//
int main(int argc, char **argv)
{ 
  char filename[50],ctemp;
  FILE *fp_tensor,*fp_newtsr;
  int s,r,c,i,xres,yres,zres;
  unsigned long vsize;
  double *xc,*d1,*e1,**a1,*V1,*V2,*V3,**D;
  float *data[6];
  short op;
 
  if (argc != 7) {
    fprintf(stderr, "%s-  find log or exp of tensor(float) .nii or .img format with lower-tri (xx,yx,yy,zx,zy,zz) order.\n",argv[0]);
    fprintf(stderr, "\nInput: file.{nii or img, same format}\nOutput:{log,exp}_file.nii\n");
    fprintf(stderr, "Usage: %s <prefix> <suffix> <xres> <yres> <zres> <0:log, 1:exp>\n",argv[0]);
    fprintf(stderr, "i.e. %s meanD img 256 256 256 0\n",argv[0]);
    exit(1);
  }

  xres=atoi(argv[3]);
  yres=atoi(argv[4]);
  zres=atoi(argv[5]);
  op=atoi(argv[6]);
  vsize=xres*yres*zres;

  xc=vector(1,6);
  d1=vector(1,3);
  e1=vector(1,3);
  a1=matrix(1,3,1,3);

  V1=vector(1,3);
  V2=vector(1,3);
  V3=vector(1,3);
  D=matrix(1,3,1,3);

  // open tensor file
  sprintf(filename,"%s.%s",argv[1],argv[2]);
  if((fp_tensor=fopen(filename,"rb"))==NULL) {
    fprintf(stderr,"Error opening file %s\n",filename);
    exit(1);
  }

  if (!strncmp(argv[2],"nii",3)) {
    printf("nii file, copying header to output...\n");

    if (op) sprintf(filename,"exp_%s.nii",argv[1],argv[2]);
    else sprintf(filename,"log_%s.nii",argv[1],argv[2]);
    if((fp_newtsr=fopen(filename,"wb"))==NULL) {
      fprintf(stderr,"Error opening file %s\n",filename);
      exit(1);
    }

    for (i=1;i<=HDR;i++) {
      c=fgetc(fp_tensor);
      fputc(c,fp_newtsr);
    }
  } else {
    if (op) sprintf(filename,"exp_%s.img",argv[1],argv[2]);
    else sprintf(filename,"log_%s.img",argv[1],argv[2]);
    if((fp_newtsr=fopen(filename,"wb"))==NULL) {
      fprintf(stderr,"Error opening file %s\n",filename);
      exit(1);
    }
  }

  // allocate buffers and read source data
  for (i=0;i<6;i++) {
    data[i] = malloc(sizeof(float)*vsize);
    if (data[i] == NULL) {
      printf("*** Failed to allocate buffer %d\n",i+1);
      exit(1);
    }

    if((fread(data[i],sizeof(float),vsize,fp_tensor))!=vsize) {
      fprintf(stderr,"Error reading data file in vol %d.\n",i+1);
      exit(1);
    }
  }
  fclose (fp_tensor);

  //mainloop
  for(s=1;s<=zres;s++)
    for(r=1;r<=yres;r++)
      for(c=1;c<=xres;c++) {
	for(i=1;i<=6;i++) xc[i]=data[i-1][(s-1)*xres*yres+(r-1)*xres+c-1];

	if(xc[1]!=0 || xc[2]!=0 || xc[3]!=0 || xc[4]!=0 || xc[5]!=0 || xc[6]!=0) {// make sure the tensor is not null
	  a1[1][1]=xc[1];
	  a1[2][2]=xc[3]; 
	  a1[3][3]=xc[6];
	  a1[1][2]=xc[2];
	  a1[1][3]=xc[4];
	  a1[2][3]=xc[5];
	  a1[2][1]=xc[2];
	  a1[3][1]=xc[4];
	  a1[3][2]=xc[5];

	  //eigen decomp
	  tred2(a1,3,d1,e1);
	  tqli(d1,e1,3,a1);  
	  eigsrt(d1,a1,3);

	  //exp or log
	  if (op) {
	    d1[1] = exp(d1[1]);
	    d1[2] = exp(d1[2]);
	    d1[3] = exp(d1[3]);
	  } else {
	    d1[1] = log(d1[1]);
	    d1[2] = log(d1[2]);
	    d1[3] = log(d1[3]);
	  }

	  //substitute eigenvectors
	  V1[1] = a1[1][1];
	  V1[2] = a1[2][1];
	  V1[3] = a1[3][1];

	  V2[1] = a1[1][2];
	  V2[2] = a1[2][2];
	  V2[3] = a1[3][2];

	  V3[1] = a1[1][3];
	  V3[2] = a1[2][3];
	  V3[3] = a1[3][3];
	      
	  gettensor(V1,V2,V3,d1[1],d1[2],d1[3],D); //calculate tensor from eig sys
	  // overwrite raw tensor
	  data[0][(s-1)*xres*yres+(r-1)*xres+c-1]=D[1][1];
	  data[1][(s-1)*xres*yres+(r-1)*xres+c-1]=D[2][1];
	  data[2][(s-1)*xres*yres+(r-1)*xres+c-1]=D[2][2];
	  data[3][(s-1)*xres*yres+(r-1)*xres+c-1]=D[3][1];
	  data[4][(s-1)*xres*yres+(r-1)*xres+c-1]=D[3][2];
	  data[5][(s-1)*xres*yres+(r-1)*xres+c-1]=D[3][3];
	} 
      }
  printf("Writing output files.\n");

  //write to output file
  for (i=0;i<6;i++) {
    if((fwrite(data[i],sizeof(float),vsize,fp_newtsr))!=vsize) {
      fprintf(stderr,"Error writing data to output file in vol %d.\n",i+1);
      exit(1);
    }
  }

  // clean up
  free_vector(xc,1,6);
  free_vector(d1,1,3);
  free_vector(e1,1,3);
  free_vector(V1,1,3);
  free_vector(V2,1,3);
  free_vector(V3,1,3);
  free_matrix(a1,1,3,1,3);
  free_matrix(D,1,3,1,3);
  
  for (i=0;i<6;i++) free(data[i]);
  fclose(fp_newtsr);

  return 0;
}



	
