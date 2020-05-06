/* SZ */
#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include <math.h>
#include <time.h>

#define COMMAND_LEN 256

int main(int argc, char **argv)
{ 
  char filename[180], comm[COMMAND_LEN],subname[20];
  FILE *fp_sub;
  unsigned int i,j,k,subc=0,nsubj,ncombi,nstart;

  double duration,hr,mn,sec;
  time_t start, finish;

  //check # of arguments
  if (argc!=4) {
    fprintf(stderr, "%s - make dti templates from nii tensors output by DTITK based on the combinations in rplist.txt.\n",argv[0]);
    fprintf(stderr, "\nInput: rplist.txt, subjects-??.txt, tensor and mask (nii.gz) in the same space, programs 'ExpLogTkNiiFloatTsr' and 'tkniitsr2dt'\nOutput: meantemplate???.img (depends on the last parameter in the Usage below and # rows in rplist.txt)\n");

    fprintf(stderr,"\nUsage: %s <# subjects> <# of combinations> <template # start from>", argv[0]);
    fprintf(stderr,"\ne.g. %s 72 100 1.\n", argv[0]);
    exit(1);
  }

  nsubj=atoi(argv[1]);
  ncombi=atoi(argv[2]);
  nstart=atoi(argv[3]);

  char sub[nsubj][30];

  //open file that contains all subject names
  sprintf(filename,"subjects-%d.txt",nsubj);
  if((fp_sub=fopen(filename,"r"))==NULL) {
    fprintf(stderr,"Error opening %s\n",filename);
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
    if (++subc == nsubj) break;
  }

  fclose(fp_sub);
  ////////////////////////////////////////////////
  //open list txt file
  unsigned int list[ncombi][nsubj];

  if((fp_sub=fopen("rplist.txt","r"))==NULL) {
    fprintf(stderr,"Error opening rplist.txt\n");
    exit(1);
  }

  for(i=0;i<ncombi;i++)
    for(j=0;j< nsubj;j++)
      fscanf(fp_sub,"%d",&list[i][j]); //scan # to matrix

  fclose(fp_sub);
  ////////////////////////////////////////////////
  //mainloop

  for (i=0;i<ncombi;i++) {
    start = time(0); // start timing
    printf("-------------------------------\nCreating mean temp %d...\n",nstart+i);

    // make empty mask 
    sprintf(comm,"fslmaths AM_6007_warped_mask -mul 0 sum-mask%d",nstart);
    system(comm);

    //make empty log tensor
    sprintf(comm,"TVtool -in log_AM_6007_final_masked.nii.gz -mask sum-mask%d.nii.gz -out sum-log%d.nii.gz",nstart,nstart);
    system(comm);

    // run thrught each subject in a selected row
    for (j=0;j<nsubj;j++) {
      // determine a single # from row i col j
      k=list[i][j]-1;

      // add log tensor
      sprintf(comm,"TVtool -in sum-log%d.nii.gz -add log_%s_final_masked.nii.gz -out sum-log%d.nii.gz",nstart,sub[k],nstart);
      system(comm);

      // add mask
      sprintf(comm,"fslmaths sum-mask%d -add %s_warped_mask sum-mask%d",nstart,sub[k],nstart);
      system(comm);
    }

    // average log tensor
    sprintf(comm,"TVDivideSV -tv sum-log%d.nii.gz -sv sum-mask%d.nii.gz -out meantemp%d.nii.gz",nstart,nstart,nstart+i);	  
    system(comm);

    sprintf(comm,"gunzip meantemp%d.nii.gz;",nstart+i);	  
    system(comm);

    // convert log tensor to normal tensor
    sprintf(comm,"./ExpLogTkNiiFloatTsr meantemp%d nii 256 256 256 1",nstart+i);	  
    system(comm);

    sprintf(comm,"mv exp_meantemp%d.nii meantemp%d.nii;",nstart+i,nstart+i);	  
    system(comm);
    
    // adjust diffusivity unit
    sprintf(comm,"TVtool -in meantemp%d.nii -scale 0.001 -out meantemp%d.nii",nstart+i,nstart+i);	  
    system(comm);

    //
    sprintf(comm,"./tkniitsr2dt meantemp%d 256 256 256",nstart+i);	  
    system(comm);

    sprintf(comm,"gzip meantemp%d.img;rm meantemp%d.nii",nstart+i,nstart+i);	  
    system(comm);
    
    finish = time(0);
    //show runtime
    duration = difftime(finish,start);
    hr=floor(duration/3600);
    mn=floor((duration-hr*3600)/60);
    sec=duration-hr*3600-mn*60;
    printf("Runtime: %.1lf hr %.1lf min and %.1lf sec.\n",hr,mn,sec);
  }

  // clean up unnecessary files
  sprintf(comm,"rm sum-mask%d.nii.gz sum-log%d.nii.gz",nstart,nstart);
  system(comm);
  //system("rm meanD*");

  return 0;
}



	
