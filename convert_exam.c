/**********************************************************/
/*                       9/18/2003                        */
/*                 Konstantinos Arfanakis                 */
/*                                                        */
/* This program will convert the I.xxx files that come    */
/* from the scanner to I.###.s$$ files, where ### is the  */
/* index number of the volume in which this file belongs, */
/* and $$ is the slice number for the file.The final files*/
/* include the same size of header as the original files  */
/* but the new files have the IMGF information at the     */
/* first bytes. This way the byte order of the files can  */
/* be detected automatically by all other programs.       */
/* This program is called with several arguments:         */
/* convert_exam 256 256 7904 13 31 2 1                    */
/* 256 256 are the x,y resolution respectively, 7904 is   */
/* the headersize, 13 is the number of volumes you will   */
/* end up with, 31 is the number of slices, 2 is the      */
/* number of repetitions of the first volume included in  */
/* I.xxx (the repetitions of the first volume will be     */
/* averaged by this program)(in case of DTI the first     */
/* volume is the T2 images), and 1 shows that the I.xxx   */
/* files were written in a different "endian" system than */
/* the current.                                           */ 
/*                                                        */
/**********************************************************/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stddef.h>
#include <math.h>

#define NR_END 1
#define FREE_ARG char* 

void writeheader(FILE *fpimage,int header_length,int xres,int yres,float sclfctr);
double ***d3tensor(long nrl, long nrh, long ncl, long nch, long ndl, long ndh);
void free_d3tensor(double ***t, long nrl, long nrh, long ncl, long nch, long ndl, long ndh);
void swap16(void* p);
void nrerror(char error_text[]);

int main(argc,argv)
int argc;
char **argv;
{
  FILE *fpin,*fpout;
  char filenamein[20],filenameout[20];
  int xres,yres,headersize,volumes,slices,t2vols,endian;
  int x,y,i,j,k,l;
  int count = 1;
  double ***image;
  char ctemp;
  short temp;

  if (argc!=8)
    {
      fprintf(stderr,"\nIncorrect Usage. Arguments are missing.");
      fprintf(stderr,"\nTo use this program correctly provide the following arguments");
      fprintf(stderr,"\nconvert_exam xresolution yresolution headersize numoffinalvolumes numofslices repsof1stvolume endianconflict\n");
      exit(1);
    }

  xres=atoi(argv[1]);
  yres=atoi(argv[2]);
  headersize=atoi(argv[3]);
  volumes=atoi(argv[4]);
  slices=atoi(argv[5]);
  t2vols=atoi(argv[6]);
  endian=atoi(argv[7]);

  image=d3tensor(1,xres,1,yres,1,slices);

  for(x=1;x<=xres;x++) /* initialize image matrix */
    for(y=1;y<=yres;y++)
      for(j=1;j<=slices;j++)
	image[x][y][j]=0.;

  /* calculate the average of the T2 volumes and produce 1 mean T2 volume */
  for(k=1;k<=t2vols;k++) 
    for(j=1;j<=slices;j++)
      {
	sprintf(filenamein,"I.%03d",count);
	fpin = fopen(filenamein,"rb") ; 
	count++ ;

	if (fpin == NULL)
	  {
	    fprintf(stderr,"Error opening file %s",filenamein);
	    exit(1);
	  }
		 
	for(l=0;l<headersize;l++)
	  fread(&ctemp,sizeof(char),1,fpin);

	for (x=1;x<=xres;x++)
	  for(y=1;y<=yres;y++)
	    {
	      fread(&temp,sizeof(short),1,fpin) ;

	      if(endian==1)
		swap16(&temp);
	      image[x][y][j]=image[x][y][j]+(double)temp;
	    }

	fclose(fpin);
      }

  for(x=1;x<=xres;x++)
    for(y=1;y<=yres;y++)
      for(j=1;j<=slices;j++)
	image[x][y][j]=image[x][y][j]/(double)t2vols;

  /* now, also read the rest of the files and write out the final images */
  
  for(i=1;i<=volumes;i++)
    for(j=1;j<=slices;j++)
      {
	if(i>1)
	  {
	    sprintf(filenamein,"I.%03d",count);
	    fpin = fopen(filenamein,"rb") ; 
	    count++ ;

	    if (fpin == NULL)
	      {
		fprintf(stderr,"Error opening file %s",filenamein);
		exit(1);
	      }
	    
	    for(l=0;l<headersize;l++)
	      fread(&ctemp,sizeof(char),1,fpin);
	  }

	sprintf(filenameout,"I.%03d.s%02d",i,j);
	fpout = fopen(filenameout,"wb") ; 
	
	if (fpout == NULL)
	  {
	    fprintf(stderr,"Error opening file %s",filenameout);
	    exit(1);
	  }
	
	writeheader(fpout,headersize,xres,yres,1.);
	
	if(i==1)
	  {
	    for(x=1;x<=xres;x++)
	      for(y=1;y<=yres;y++)
		{
		  temp=(short)image[x][y][j];
		  fwrite(&temp, sizeof(short),1,fpout);
		}
	  }
	else
	  {
	    for(x=1;x<=xres;x++)
	      for(y=1;y<=yres;y++)
		{
		  fread(&temp,sizeof(short),1,fpin) ;
		  if(endian==1)
		    swap16(&temp);
		  fwrite(&temp,sizeof(short),1,fpout) ;
		}
	  }

	if(i>1) 
	  fclose(fpin);

	fclose(fpout);
      }
	
	
  return 0;
}


void writeheader(FILE *fpimage,int header_length,int xres,int yres,float sclfctr)
{
	/*  char magicnum[] = "IMGF"; */
	  int bits_per_short=16;
	  short s1=044515,s2=043506;

      fseek(fpimage, (long) 0, (int) 0);

/*      fwrite(magicnum,sizeof(char),4,fpimage); */
      fwrite(&s1,sizeof(short),1,fpimage);
      fwrite(&s2,sizeof(short),1,fpimage);

      fwrite(&header_length,sizeof(int),1,fpimage);
      fwrite(&xres,sizeof(int),1,fpimage);
      fwrite(&yres,sizeof(int),1,fpimage);
      fwrite(&bits_per_short,sizeof(int),1,fpimage);
      fwrite(&sclfctr,sizeof(float),1,fpimage);
      fseek(fpimage, (long) header_length, (int) 0);
	  return;
}

double ***d3tensor(long nrl, long nrh, long ncl, long nch, long ndl, long ndh)
/* allocate a double 3tensor with range t[nrl..nrh][ncl..nch][ndl..ndh] */
{
	long i,j,nrow=nrh-nrl+1,ncol=nch-ncl+1,ndep=ndh-ndl+1;
	double ***t;

	/* allocate pointers to pointers to rows */
	t=(double ***) malloc((size_t)((nrow+NR_END)*sizeof(double**)));
	if (!t) nrerror("allocation failure 1 in s3tensor()");
	t += NR_END;
	t -= nrl;

	/* allocate pointers to rows and set pointers to them */
	t[nrl]=(double **) malloc((size_t)((nrow*ncol+NR_END)*sizeof(double*)));
	if (!t[nrl]) nrerror("allocation failure 2 in s3tensor()");
	t[nrl] += NR_END;
	t[nrl] -= ncl;

	/* allocate rows and set pointers to them */
	t[nrl][ncl]=(double *) malloc((size_t)((nrow*ncol*ndep+NR_END)*sizeof(double)));
	if (!t[nrl][ncl]) nrerror("allocation failure 3 in s3tensor()");
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
/* free a double d3tensor allocated by d3tensor() */
{
	free((FREE_ARG) (t[nrl][ncl]+ndl-NR_END));
	free((FREE_ARG) (t[nrl]+ncl-NR_END));
	free((FREE_ARG) (t+nrl-NR_END));
}

void swap16(void* p)
{
  ((unsigned char*)p)[0]^=((unsigned char*)p)[1];
  ((unsigned char*)p)[1]^=((unsigned char*)p)[0];
  ((unsigned char*)p)[0]^=((unsigned char*)p)[1];
}

void nrerror(char error_text[])
/* Numerical Recipes standard error handler */
{
        fprintf(stderr,"Numerical Recipes run-time error...\n");
        fprintf(stderr,"%s\n",error_text);
        fprintf(stderr,"...now exiting to system...\n");
        exit(1);
}
