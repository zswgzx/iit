/* SZ last updated: 150810
revised from subjtvdt-v1.c and should replace it
*/

#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include <math.h>

float varn(float *xc, int num)
{
  //find the standard deviation of a 1D array, ignore 0 in it!
  float avg,variance,sum=0.0;
  int i,count=0;

  //find average
  for (i=0;i<num;i++) {
    if (xc[i]!=0) { 
      count+=1;
      sum += xc[i];
    }
  }

  if (count!=0) {
    if (count==1) variance=0.0; //exclude the situation with only 1 nonzero value
    else {
      avg=sum/count;
      sum=0;
      //sum of square
      for (i=0;i<num;i++) {
	if (xc[i]!=0)      sum += pow(xc[i]-avg,2);
      }
      variance = sum/count;  //notice this is N variance, not N-1 !!
    }
  }
  else variance=0.0;

  return variance;
}

/*------------------------------------*/
int main(int argc, char **argv)
{ 
  char filename[100],subname[25];
  FILE *fp_tvdt,*fp_sub;
  int s,r,c,i,j,xres,yres,zres,subc=0,combi;
  float var,tvdt;
  const char *tsrFile = "reg%s.img";
 
  if (argc!=5) {
    fprintf(stderr,"%s - calculate total variance of diffusion tensor for raw tensor elements files (xx,yy,zz,xy,xz,yz order), unit is same as input\n",argv[0]);
    fprintf(stderr, "\nInput: subjects-??.txt,reg[subj].img\nOutput: tvdt.img (float)\n");
    fprintf(stderr,"\nUsage: %s <xdim> <ydim> <zdim> <# subj>\n", argv[0]);
    fprintf(stderr,"\ne.g. %s 181 217 181 18\n", argv[0]);
    exit(1);
  }
 
  xres = atoi(argv[1]);
  yres = atoi(argv[2]);
  zres = atoi(argv[3]);
  combi = atoi(argv[4]);

  char sub[combi][25];
  FILE *fp[combi];
  float tsr[combi];

//////////////////////////////////////////////////////////
  //open file that contains all subject names
  sprintf(filename,"subjects-%d.txt",combi);
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
    if (++subc == combi) break;
  }

  fclose(fp_sub);

  //open tensor files
  for (i=0;i<combi;i++) {
    sprintf(filename,tsrFile,sub[i],i);
    if((fp[i]=fopen(filename,"rb"))==NULL) {
      fprintf(stderr,"Error opening tensor file %s\n",filename);
      exit(1);
    }
  }

  //open file to write
  if((fp_tvdt=fopen("tvdt.img","wb"))==NULL) {
    fprintf(stderr,"Error opening file tvdt.img to write\n");
    exit(1);
  }

  //vox-wise loop
  for(s=0;s<zres;s++) {
    for(r=0;r<yres;r++) {
      for(c=0;c<xres;c++) {
	  tvdt=0.0;
	  //process 6 elements
	  for(i=0;i<6;i++) {
	    for (j=0;j<combi;j++) {
	      if((fread(&tsr[j],sizeof(float),1,fp[j]))!=1) {
		fprintf(stderr,"Error reading tensor file tensor.img\n");
		exit(1);
	      }
	    }

	    var=varn(tsr,combi);

	    if (i>2)	  tvdt+=2*var;   //off diag element add twice
	    else	  tvdt+=var;
	    /*
	    //debug
	    if (s==128 && r==128 && c==127) {
	      for (j=0;j<combi;j++) printf("tsr %d - %d =%f\n",j+1,i+1,tsr[j]);
	      printf("tvdt=%.9f\n",tvdt);
	    }
	    */
	  }

	  fwrite(&tvdt,sizeof(float),1,fp_tvdt);    
	}
      }
    }

  for (i=0;i<combi;i++) fclose(fp[i]);
  fclose(fp_tvdt);
  return 0;
}



	
