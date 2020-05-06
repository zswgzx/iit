/* SZ */
/* this version is much faster than BootstrapMeanTemplateByDTITK */

#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include <math.h>
#include <time.h>

#define HDR 352
typedef short sub_mask;

int main(int argc, char **argv)
{ 
  char filename[180],/*comm[40],*/subname[20],zerochar=0,*countmean;
  FILE *fp_txt,*fp_sub,*fp_meantemp,*fp_count,*fp_ssub,*fp_mask;
  unsigned int s,r,c,i,j,k,subc=0,N,XRES,YRES,ZRES,COMBI,NUMSUB,cnt,SSIZE;
  float zero=0.0,*meantsr,*subtsr;
  sub_mask *countsub;
  double duration,hr,mn,sec;
  time_t start, finish;

  /////////////////
  //check # of arguments
  if (argc!=7) {
      fprintf(stderr,"\n%s - make bootstrap mean templates\n",argv[0]);
      fprintf(stderr,"\nInput: subjects-[N].txt, [n]_[BS]rep.txt, subj_combined.nii, subj-mask.nii\nOutput: meantemp?.img (tsr elements volume by volume, order is the same as DTITK, no header)\n");
      fprintf(stderr,"\nUsage: %s <# subj chosen from N, i.e. n> <# bootstrap samples, i.e. BS> <xdim> <ydim> <zdim> <# subj, i.e. N>", argv[0]);
      fprintf(stderr,"\ne.g. %s 30 100 256 256 256 94.\n", argv[0]);
      exit(1);
    }

  N=atoi(argv[1]);
  COMBI=atoi(argv[2]);
  XRES=atoi(argv[3]);
  YRES=atoi(argv[4]);
  ZRES=atoi(argv[5]);
  NUMSUB=atoi(argv[6]);
  SSIZE=XRES*YRES*ZRES;

  char sub[NUMSUB][20];

  //////////////////////////////////////////////////////////
  //open subject namelist, store them to arrays
  sprintf(filename,"subjects-%d.txt",NUMSUB);
  if((fp_sub=fopen(filename,"r"))==NULL) {
    fprintf(stderr,"Error opening file %s\n",filename);
    exit(1);
  }

  while (fgets(subname, sizeof(subname), fp_sub)) {
    // cut newline from the end of subname
    char *tmp = subname;
    while (*tmp) {
      if (*tmp == '\n') *tmp = '\0';
      tmp++;
    }
    strcpy(sub[subc],subname);
    if (++subc == NUMSUB) break;
   }

  fclose(fp_sub);
  ////////////////////////////////////////////////
  //open list txt file
  int list[COMBI][N];

  sprintf(filename,"%d_%drep.txt",N,COMBI);
  if((fp_txt=fopen(filename,"r"))==NULL) {
    fprintf(stderr,"Error opening txt file %s\n",filename);
    exit(1);
  }

  for(i=0;i< COMBI;i++)
    for(j=0;j< N;j++) {
      fscanf(fp_txt,"%d",&list[i][j]); //scan # to matrix
    }

  fclose(fp_txt);
  ////////////////////////////////////////////////
  //main
  //allocate buffer to store 1 volume of data
  
  meantsr = malloc(sizeof(float)*SSIZE);
  if (meantsr == NULL) {
    printf("*** Failed to allocate buffer 1\n");
    exit(1);
  }

  subtsr = malloc(sizeof(float)*SSIZE);
  if (subtsr == NULL) {
    printf("*** Failed to allocate buffer 2\n");
    exit(1);
  }

  countmean = malloc(sizeof(char)*SSIZE);
  if (countmean == NULL) {
    printf("*** Failed to allocate buffer 3\n");
    exit(1);
  }

  countsub = malloc(sizeof(sub_mask)*SSIZE);
  if (countsub == NULL) {
    printf("*** Failed to allocate buffer 4\n");
    exit(1);
  }
  
  for(i=0;i< COMBI;i++) {
      printf("-------------------------------\nCreating mean temp %d...\n",i+1);
      ///////////////////////////////////////////
      //initialize count and template to all zero
      sprintf(filename,"meantemp%d.img",i+1);
      if((fp_meantemp=fopen(filename,"w+"))==NULL) {
	fprintf(stderr,"Error opening file %s\n",filename);
	exit(1);
      }

      for(s=1;s<=ZRES;s++)
	for(r=1;r<=YRES;r++)
	  for(c=1;c<=XRES;c++)
	    for(k=1;k<=6;k++) {
	      fwrite(&zero,sizeof(float),1,fp_meantemp);
	    }

      rewind(fp_meantemp);

      sprintf(filename,"count.img");
      if((fp_count=fopen(filename,"w+"))==NULL)	{
	fprintf(stderr,"Error opening file %s\n",filename);
	exit(1);
      }

      for(s=1;s<=ZRES;s++)
	for(r=1;r<=YRES;r++)
	  for(c=1;c<=XRES;c++) {
	    fwrite(&zerochar,sizeof(char),1,fp_count);
	  }

      rewind(fp_count);
      ////////////////////////////////
      // run thru all chosen subjects
      start = time(0);
      for(j=0;j<N;j++) {
	//determine a single # from row i col j
	k=list[i][j]-1;

	sprintf(filename,"../../tensors/registered/%s_combined.nii",sub[k]);
	if((fp_ssub=fopen(filename,"rb"))==NULL) {
	  fprintf(stderr,"Error opening file %s\n", filename);
	  exit(1);
	}
	if ((fseek(fp_ssub,HDR,SEEK_SET))!=0) {
	  fprintf(stderr,"Error fseeking file %s\n", filename);
	  exit(1);
	}

	sprintf(filename,"../../reg-masks/%s-mask.nii",sub[k]);
	if((fp_mask=fopen(filename,"rb"))==NULL) {
	  fprintf(stderr,"Error opening file %s\n", filename);
	  exit(1);
	}
	if ((fseek(fp_mask,HDR,SEEK_SET))!=0) {
	  fprintf(stderr,"Error fseeking file %s\n", filename);
	  exit(1);
	}

	//volume-vise loop
	if((fread(countsub,sizeof(sub_mask),SSIZE,fp_mask))!=SSIZE) {
	  fprintf(stderr,"Error reading mask file.\n");
	  exit(1);
	}

	if((fread(countmean,sizeof(char),SSIZE,fp_count))!=SSIZE) {
	  fprintf(stderr,"Error reading count file.\n");
	  exit(1);
	}
	rewind(fp_count); //important!

	//read from mask
	for(cnt=1;cnt<=SSIZE;cnt++) {
	  if (countsub[cnt-1] != 0) countmean[cnt-1]+=1;
	}
	if ((fwrite(countmean,sizeof(char),SSIZE,fp_count))!=SSIZE) {
	  fprintf(stderr,"Error writing temp mask file.\n");
	  exit(1);
	}

	for(s=1;s<=6;s++) {
	  //read from tensor
	  if((fread(subtsr,sizeof(float),SSIZE,fp_ssub))!=SSIZE) {
	    fprintf(stderr,"Error reading subject tensor file.\n");
	    exit(1);
	  }

	  if((fread(meantsr,sizeof(float),SSIZE,fp_meantemp))!=SSIZE) {
	    fprintf(stderr,"Error reading mean tensor file for adding in vol %d, subj %d.\n",s,j+1);
	    exit(1);
	  }
	  if ((fseek(fp_meantemp,(s-1)*sizeof(float)*SSIZE,SEEK_SET))!=0) {
	    fprintf(stderr,"fseek error in adding process vol %d\n",s);
	    exit(1);
	  } //important!

	  for(cnt=1;cnt<=SSIZE;cnt++) {
	    if (countsub[cnt-1] != 0)  meantsr[cnt-1]+=subtsr[cnt-1];
	  }

	  if ((fwrite(meantsr,sizeof(float),SSIZE,fp_meantemp))!=SSIZE) {
	    fprintf(stderr,"Error writing temp tensor file vol %d.\n",s);
	    exit(1);
	  }
	}

	rewind(fp_count);
	rewind(fp_meantemp);
	fclose(fp_ssub);
	fclose(fp_mask);
      }
    
    //calculate average and write data
      if ((fread(countmean,sizeof(char),SSIZE,fp_count))!=SSIZE) {
	fprintf(stderr,"Error reading temp count file.\n");
	exit(1);
      }

      for (s=1;s<=6;s++) {
	if((fread(meantsr,sizeof(float),SSIZE,fp_meantemp))!=SSIZE) {
	  fprintf(stderr,"Error reading mean tensor file.\n");
	  exit(1);
	}
	if ((fseek(fp_meantemp,(s-1)*sizeof(float)*SSIZE,SEEK_SET))!=0) {
	  fprintf(stderr,"fseek error in adding process vol %d\n",s);
	  exit(1);
	} //important!

	for (cnt=1;cnt<=SSIZE;cnt++) {
	  if (countmean[cnt-1] != 0) meantsr[cnt-1]/=countmean[cnt-1];
	  }
	if ((fwrite(meantsr,sizeof(float),SSIZE,fp_meantemp))!=SSIZE) {
	    fprintf(stderr,"Error writing average tensor file vol %d.\n",s);
	    exit(1);
	  }
      }
		      
      finish = time(0);
      fclose(fp_count);
      fclose(fp_meantemp);

      //show runtime
      duration = difftime(finish,start);
      hr=floor(duration/3600);
      mn=floor((duration-hr*3600)/60);
      sec=duration-hr*3600-mn*60;
      printf("Runtime: %.1lf hr %.1lf min and %.1lf sec.\n",hr,mn,sec);

      //sprintf(comm,"gzip meantemp%d.img",i+1);	  
      //system(comm);
    }
  
  free(meantsr);
  free(subtsr);
  free(countmean);
  free(countsub);

  return 0;
}



	
