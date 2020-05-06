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

float ***f3tensor(long nrl, long nrh, long ncl, long nch, long ndl, long ndh)
/* allocate a float 3tensor with range t[nrl..nrh][ncl..nch][ndl..ndh] */
{
  long i,j,nrow=nrh-nrl+1,ncol=nch-ncl+1,ndep=ndh-ndl+1;
  float ***t;
  
  /* allocate pointers to pointers to rows */
  t=(float ***) malloc((size_t)((nrow+NR_END)*sizeof(float**)));
  if (!t) nrerror("allocation failure 1 in d3tensor()");
  t += NR_END;
  t -= nrl;
  
  /* allocate pointers to rows and set pointers to them */
  t[nrl]=(float **) malloc((size_t)((nrow*ncol+NR_END)*sizeof(float*)));
  if (!t[nrl]) nrerror("allocation failure 2 in d3tensor()");
  t[nrl] += NR_END;
  t[nrl] -= ncl;
  
  /* allocate rows and set pointers to them */
  t[nrl][ncl]=(float *) malloc((size_t)((nrow*ncol*ndep+NR_END)*sizeof(float)));
  if (!t[nrl][ncl]) nrerror("allocation failure 3 in d3tensor()");
  t[nrl][ncl] += NR_END;
  t[nrl][ncl] -= ndl;
  
  for(j=ncl+1;j<=nch;j++) t[nrl][j]=t[nrl][j-1]+ndep;
  for(i=nrl+1;i<=nrh;i++) {
    t[i]=t[i-1]+ncol;
    t[i][ncl]=t[i-1][ncl]+ncol*ndep;
    for(j=ncl+1;j<=nch;j++) t[i][j]=t[i][j-1]+ndep;
  }
  
   /* return pointer to array of pointers to rows */
  return t;
}

void free_f3tensor(float ***t, long nrl, long nrh, long ncl, long nch, long ndl, long ndh)
/* free a float d3tensor allocated by d3tensor() */
{
  free((FREE_ARG) (t[nrl][ncl]+ndl-NR_END));
  free((FREE_ARG) (t[nrl]+ncl-NR_END));
  free((FREE_ARG) (t+nrl-NR_END));
}

double ***d3tensor(long nrl, long nrh, long ncl, long nch, long ndl, long ndh)
/* allocate a double 3tensor with range t[nrl..nrh][ncl..nch][ndl..ndh] */
{
  long i,j,nrow=nrh-nrl+1,ncol=nch-ncl+1,ndep=ndh-ndl+1;
  double ***t;
  
  /* allocate pointers to pointers to rows */
  t=(double ***) malloc((size_t)((nrow+NR_END)*sizeof(double**)));
  if (!t) nrerror("allocation failure 1 in d3tensor()");
  t += NR_END;
  t -= nrl;
  
  /* allocate pointers to rows and set pointers to them */
  t[nrl]=(double **) malloc((size_t)((nrow*ncol+NR_END)*sizeof(double*)));
  if (!t[nrl]) nrerror("allocation failure 2 in d3tensor()");
  t[nrl] += NR_END;
  t[nrl] -= ncl;
  
  /* allocate rows and set pointers to them */
  t[nrl][ncl]=(double *) malloc((size_t)((nrow*ncol*ndep+NR_END)*sizeof(double)));
  if (!t[nrl][ncl]) nrerror("allocation failure 3 in d3tensor()");
  t[nrl][ncl] += NR_END;
  t[nrl][ncl] -= ndl;
  
  for(j=ncl+1;j<=nch;j++) t[nrl][j]=t[nrl][j-1]+ndep;
  for(i=nrl+1;i<=nrh;i++) {
    t[i]=t[i-1]+ncol;
    t[i][ncl]=t[i-1][ncl]+ncol*ndep;
    for(j=ncl+1;j<=nch;j++) t[i][j]=t[i][j-1]+ndep;
  }
  
   /* return pointer to array of pointers to rows */
  return t;
}


void free_d3tensor(double ***t, long nrl, long nrh, long ncl, long nch, long ndl, long ndh)
/* free a float d3tensor allocated by d3tensor() */
{
  free((FREE_ARG) (t[nrl][ncl]+ndl-NR_END));
  free((FREE_ARG) (t[nrl]+ncl-NR_END));
  free((FREE_ARG) (t+nrl-NR_END));
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
//************************************************************************************************//
int main(int argc, char **argv)
{ 
  char filename[50];
  FILE *fp_e1,*fp_e2,*fp_e3,*fp_l1,*fp_l2,*fp_l3,*fp_tensor;
  short op;
  int s,r,c,i,xres,yres,zres;
  unsigned long vsize;
  double eigval[3],*xc,*d1,*e1,**a1,***eigvect1;
  float l1,l2,l3,*data[6],***V11,***V12,***V13, ***V21,***V22,***V23, ***V31,***V32,***V33;
 
  if (argc != 7) {
    fprintf(stderr, "%s- produce L1-3 (evalue) and/or V1-3 (evectors,volume by volume) from tensor(float) .nii or .img format with lower-tri (xx,yx,yy,zx,zy,zz) order.\n",argv[0]);
    fprintf(stderr, "\nInput: file.{nii or img, same format}\nOutput:file_{l1-3,or v1-3}.img\n");
    fprintf(stderr, "Usage: %s <prefix> <suffix> <xres> <yres> <zres> <op- 0:lv1-3, 1:l1-3, 2:v1-3, 3:v1>\n",argv[0]);
    fprintf(stderr, "i.e. %s meanD img 256 256 256 3\n",argv[0]);
    exit(1);
  }

  xres=atoi(argv[3]);
  yres=atoi(argv[4]);
  zres=atoi(argv[5]);
  op=atoi(argv[6]);
  if (op>3 || op <0) {printf("invalid op value, exiting...\n");exit(1);}
  vsize=xres*yres*zres;

  xc=vector(1,6);
  d1=vector(1,3);
  e1=vector(1,3);
  a1=matrix(1,3,1,3);

  eigvect1 = d3tensor(1,yres,1,xres,1,3);

  V11=f3tensor(1,zres,1,yres,1,xres);
  V12=f3tensor(1,zres,1,yres,1,xres);
  V13=f3tensor(1,zres,1,yres,1,xres);
  V21=f3tensor(1,zres,1,yres,1,xres);
  V22=f3tensor(1,zres,1,yres,1,xres);
  V23=f3tensor(1,zres,1,yres,1,xres);
  V31=f3tensor(1,zres,1,yres,1,xres);
  V32=f3tensor(1,zres,1,yres,1,xres);
  V33=f3tensor(1,zres,1,yres,1,xres);

  // create eigen system files based on op
  if (op==0 || op==1) {
    sprintf(filename,"%s_l1.img",argv[1]);
    if((fp_l1=fopen(filename,"wb"))==NULL) {
      fprintf(stderr,"Error openning file %s\n",filename);
      exit(1);
    }

    sprintf(filename,"%s_l2.img",argv[1]);
    if((fp_l2=fopen(filename,"wb"))==NULL) {
      fprintf(stderr,"Error openning file %s\n",filename);
      exit(1);
    }

    sprintf(filename,"%s_l3.img",argv[1]);
    if((fp_l3=fopen(filename,"wb"))==NULL) {
      fprintf(stderr,"Error openning file %s\n",filename);
      exit(1);
    }
  }

  if (op==0 || op==2 || op==3) {
    sprintf(filename,"%s_v1.img",argv[1]);
    if((fp_e1=fopen(filename,"wb"))==NULL) {
      fprintf(stderr,"Error openning file %s\n",filename);
      exit(1);
    }
  }

  if (op==0 || op==2) {
    sprintf(filename,"%s_v2.img",argv[1]);
    if((fp_e2=fopen(filename,"wb"))==NULL) {
      fprintf(stderr,"Error openning file %s\n",filename);
      exit(1);
    }

    sprintf(filename,"%s_v3.img",argv[1]);
    if((fp_e3=fopen(filename,"wb"))==NULL) {
      fprintf(stderr,"Error openning file %s\n",filename);
      exit(1);
    }
  }

  // open tensor file
  sprintf(filename,"%s.%s",argv[1],argv[2]);
  if((fp_tensor=fopen(filename,"rb"))==NULL) {
    fprintf(stderr,"Error opening file %s\n",filename);
    exit(1);
  }

  if (!strncmp(argv[2],"nii",3)) {
    printf("nii file, skip header...\n");
    if ((fseek(fp_tensor,HDR,SEEK_SET))!=0) {
      fprintf(stderr,"source file hdr fseek error\n");
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

  //mainloop
  for(s=1;s<=zres;s++)
    for(r=1;r<=yres;r++)
      for(c=1;c<=xres;c++) {
	for(i=1;i<=6;i++) xc[i]=data[i-1][(s-1)*xres*yres+(r-1)*xres+c-1];

	if(xc[2]!=0 || xc[3]!=0 || xc[4]!=0 || xc[5]!=0 || xc[6]!=0 || xc[1]!=0) {// make sure the tensor is not null
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
	    
	  eigval[0] = d1[1];
	  eigval[1] = d1[2];
	  eigval[2] = d1[3];
	  l1=(float)d1[1];
	  l2=(float)d1[2];
	  l3=(float)d1[3];
	      
	  eigvect1[r][c][1] = a1[1][1];
	  eigvect1[r][c][2] = a1[2][1];
	  eigvect1[r][c][3] = a1[3][1];
	      
	  V11[s][r][c] = (float)a1[1][1] ;
	  V12[s][r][c] = (float)a1[2][1] ;
	  V13[s][r][c] = (float)a1[3][1] ;
	  V21[s][r][c] = (float)a1[1][2] ;
	  V22[s][r][c] = (float)a1[2][2] ;
	  V23[s][r][c] = (float)a1[3][2] ;
	  V31[s][r][c] = (float)a1[1][3] ;
	  V32[s][r][c] = (float)a1[2][3] ;
	  V33[s][r][c] = (float)a1[3][3] ;
	} else {	      
	  eigval[0]=eigval[1]=eigval[2]=0.0;
	  l1=l2=l3= 0.0;
	  eigvect1[r][c][1]=eigvect1[r][c][2]=eigvect1[r][c][3]=0.0;
	      
	  V11[s][r][c]=V12[s][r][c]=V13[s][r][c]=0.0;
	  V21[s][r][c]=V22[s][r][c]=V23[s][r][c]=0.0;
	  V31[s][r][c]=V32[s][r][c]=V33[s][r][c]=0.0;
	}

	if (op==0 || op==1) {
	  if ((fwrite(&l1,sizeof(float),1,fp_l1))!=1) {
	    fprintf(stderr,"Error writing data to l1 at %d,%d,%d (1-offset)\n",c,r,s);exit(1);
	  }
	  if ((fwrite(&l2,sizeof(float),1,fp_l2))!=1) {
	    fprintf(stderr,"Error writing data to l2 at %d,%d,%d (1-offset)\n",c,r,s);exit(1);
	  }
	  if ((fwrite(&l3,sizeof(float),1,fp_l3))!=1) {
	    fprintf(stderr,"Error writing data to l3 at %d,%d,%d (1-offset)\n",c,r,s);exit(1);
	  }
	}
      }
  
  printf("Writing output files.\n");

  for(s=1;s<=zres;s++)
    for(r=1;r<=yres;r++)
      for(c=1;c<=xres;c++) {
	if (op==0 || op==2 || op==3) {
	  if ((fwrite(&V11[s][r][c],sizeof(float),1,fp_e1))!=1) {
	    fprintf(stderr,"Error writing component 1 to e1 at %d,%d,%d (1-offset)\n",c,r,s);exit(1);
	  }
	}
	if (op==0 || op==2) {
	  if ((fwrite(&V21[s][r][c],sizeof(float),1,fp_e2))!=1) {
	    fprintf(stderr,"Error writing component 1 to e2 at %d,%d,%d (1-offset)\n",c,r,s);exit(1);
	  }
	  if ((fwrite(&V31[s][r][c],sizeof(float),1,fp_e3))!=1) {
	    fprintf(stderr,"Error writing component 1 to e3 at %d,%d,%d (1-offset)\n",c,r,s);exit(1);
	  }
	}
      }


  for(s=1;s<=zres;s++)
    for(r=1;r<=yres;r++)
      for(c=1;c<=xres;c++) {
	if (op==0 || op==2 || op==3) {
	  if ((fwrite(&V12[s][r][c],sizeof(float),1,fp_e1))!=1) {
	    fprintf(stderr,"Error writing component 2 to e1 at %d,%d,%d (1-offset)\n",c,r,s);exit(1);
	  }
	}
	if (op==0 || op==2) {
	  if ((fwrite(&V22[s][r][c],sizeof(float),1,fp_e2))!=1) {
	    fprintf(stderr,"Error writing component 2 to e2 at %d,%d,%d (1-offset)\n",c,r,s);exit(1);
	  }
	  if ((fwrite(&V32[s][r][c],sizeof(float),1,fp_e3))!=1) {
	    fprintf(stderr,"Error writing component 2 to e3 at %d,%d,%d (1-offset)\n",c,r,s);exit(1);
	  }
	}
      }

  for(s=1;s<=zres;s++)
    for(r=1;r<=yres;r++)
      for(c=1;c<=xres;c++) {
	if (op==0 || op==2 || op==3) {
	  if ((fwrite(&V13[s][r][c],sizeof(float),1,fp_e1))!=1) {
	    fprintf(stderr,"Error writing component 3 to e1 at %d,%d,%d (1-offset)\n",c,r,s);exit(1);
	  }
	}
	if (op==0 || op==2) {
	  if ((fwrite(&V23[s][r][c],sizeof(float),1,fp_e2))!=1) {
	    fprintf(stderr,"Error writing component 3 to e2 at %d,%d,%d (1-offset)\n",c,r,s);exit(1);
	  }
	  if ((fwrite(&V33[s][r][c],sizeof(float),1,fp_e3))!=1) {
	    fprintf(stderr,"Error writing component 3 to e3 at %d,%d,%d (1-offset)\n",c,r,s);exit(1);
	  }
	}
      }

  free_vector(xc,1,6);
  free_vector(d1,1,3);
  free_vector(e1,1,3);
  free_matrix(a1,1,3,1,3);
  free_d3tensor(eigvect1,1,yres,1,xres,1,3);

  free_f3tensor(V11,1,zres,1,yres,1,xres);
  free_f3tensor(V12,1,zres,1,yres,1,xres);
  free_f3tensor(V13,1,zres,1,yres,1,xres);
  free_f3tensor(V21,1,zres,1,yres,1,xres);
  free_f3tensor(V22,1,zres,1,yres,1,xres);
  free_f3tensor(V23,1,zres,1,yres,1,xres);
  free_f3tensor(V31,1,zres,1,yres,1,xres);
  free_f3tensor(V32,1,zres,1,yres,1,xres);
  free_f3tensor(V33,1,zres,1,yres,1,xres);
  
  for (i=0;i<6;i++) free(data[i]);

  if (op==0 || op==2 || op==3) fclose(fp_e1);
  if (op==0 || op==2) {
    fclose(fp_e2);
    fclose(fp_e3);
  }
  if (op==0 || op==1) {
    fclose(fp_l1);
    fclose(fp_l2);
    fclose(fp_l3);
  }
  fclose(fp_tensor);

  return 0;
}



	
