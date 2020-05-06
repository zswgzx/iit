/* SZ 
program structure:
*/
#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include <math.h>
#include <time.h>

#define HDR 352
#define MAX_LEN 100
#define MAX_LABEL 30 // max number of labels in a voxel
#define TINY 1e-8

typedef int DATA_TYPE;

int main(int argc, char **argv)
{ 
  char filename[MAX_LEN],subname[MAX_LEN],mask;
  FILE *fpseg,*fpstat,*fpcount,*fpsublist;
  int s,r,c,i,j,subc=0,xres,yres,zres,numsub,max,max_count,max_loc,seg,stats,non_gm_cnt,total_gm;
  /*
  double duration,hr,mn,sec;
  time_t start, finish;
  */
  /////////////////
  //check # of arguments
  if (argc!=7) {
    fprintf(stderr, "%s- obtain the GM segmentation in final space from all subjects' GM segmentation by voting \n\nInput: subjects-?.txt,[subj]-iit3seg.nii (type int), GM mask (binary mask, char type)\nOutput: iit3seg.nii, seg-stats.nii(# subjects that contain that label in iit3seg.nii in each voxel),seg-stats-total.nii(# subjects that have been assigned a non-zero label in each voxel)\n\n",argv[0]);
    fprintf(stderr, "Usage: %s <xdim> <ydim> <zdim> <# subjects> <segmentation_suffix> <GM mask>\n",argv[0]);
    fprintf(stderr, "i.e. %s 256 256 256 71 -iit3seg.nii gm-mask.nii\n",argv[0]);
    exit(1);
  }

  xres=atoi(argv[1]);
  yres=atoi(argv[2]);
  zres=atoi(argv[3]);
  numsub=atoi(argv[4]);

  FILE *fpsub[numsub];
  char sub[numsub][MAX_LEN];
  DATA_TYPE label[numsub],zero_count=0,label_count[2][MAX_LABEL],nlabel_flag=0;

  //////////////////////////////////////////////////////////
  //open file that contains all subject names
  sprintf(filename,"subjects-%d.txt",numsub);
  if((fpsublist=fopen(filename,"r"))==NULL) {
    fprintf(stderr,"Error opening %s\n",filename);
    exit(1);
  }

  // store subject names to arrays
  while (fgets(subname, sizeof(subname), fpsublist)) {
    // cut newline from the end of subname
    char *tmp = subname;
    while (*tmp) {
      if (*tmp == '\n') *tmp = '\0';
      tmp++;
    }
    strcpy(sub[subc],subname);
    if (++subc == numsub) break;
  }
  fclose(fpsublist);

  //////////////////////////////////////////////////////////
  //open segmentation files of all subjects

  for (i=0;i!=numsub;++i) {
    sprintf(filename,"%s%s",sub[i],argv[5]);
    if((fpsub[i]=fopen(filename,"rb"))==NULL) {
      fprintf(stderr,"Error opening %s\n",filename);
      exit(1);
    }

    if((fseek(fpsub[i],HDR,SEEK_SET)) !=0) {
      fprintf(stderr,"fseek error in %s.\n",filename);
      exit(1);
    }
  }

  //open GM mask file
  sprintf(filename,"%s",argv[6]);
  if((fpsublist=fopen(filename,"rb"))==NULL) {
    fprintf(stderr,"Error opening %s\n",argv[6]);
    exit(1);
  }

  if((fseek(fpsublist,HDR,SEEK_SET)) !=0) {
    fprintf(stderr,"fseek error in %s.\n",argv[6]);
    exit(1);
  }

  ////////////////////////////////////////////////
  // open output files

  if((fpseg=fopen("iit3seg.nii","wb"))==NULL) {
    fprintf(stderr,"Error opening file iit3seg.nii to write\n");
    exit(1);
  }

  if((fpstat=fopen("seg-stats.nii","wb"))==NULL) {
    fprintf(stderr,"Error opening file seg-stats.nii to write\n");
    exit(1);
  }

  if((fpcount=fopen("seg-stats-total.nii","wb"))==NULL) {
    fprintf(stderr,"Error opening file seg-stats-total.nii to write\n");
    exit(1);
  }

  rewind(fpsub[0]);

  for (i=0;i!=HDR;++i) {
    subc=fgetc(fpsub[0]);
    fputc(subc,fpseg);
    fputc(subc,fpstat);
    fputc(subc,fpcount);
  }
  
  ////////////////////////////////////////////////
  //mainloop
  
  for (s=0;s!=zres;++s)
    for (r=0;r!=yres;++r)
      for (c=0;c!=xres;++c) {
	//read labels of all subjects
	for (i=0;i!=numsub;++i) {
	    if((fread(&label[i],sizeof(DATA_TYPE),1,fpsub[i]))!=1) {      
	      fprintf(stderr,"Error reading segmentation from %s at [%3d,%3d,%3d] (0-offset).\n",sub[i],c,r,s);
	      exit(1);
	    }
	}

	//read GM mask
	if((fread(&mask,sizeof(char),1,fpsublist))!=1) {      
	  fprintf(stderr,"Error reading GM mask from %s at [%3d,%3d,%3d] (0-offset).\n",argv[6],c,r,s);
	  exit(1);
	}

	//outside the mask
	if (!mask) { 
	  seg=stats=total_gm=0;
	  if ((fwrite(&seg,sizeof(int),1,fpseg))!=1) printf("fwrite error of label at [%3d,%3d,%3d] (0-offset)\n",c,r,s);
	  if ((fwrite(&stats,sizeof(int),1,fpstat))!=1) printf("fwrite error of stat at [%3d,%3d,%3d] (0-offset)\n",c,r,s);
	if ((fwrite(&total_gm,sizeof(int),1,fpcount))!=1) printf("fwrite error of stat at [%3d,%3d,%3d] (0-offset)\n",c,r,s);
	  continue;
	}

	// clear some counters
	for (i=0;i!=MAX_LABEL;++i) label_count[0][i]=label_count[1][i]=0;
	zero_count=nlabel_flag=non_gm_cnt=0;

	//find all possible non-zero labels
	for (i=0;i!=numsub;++i) {
	  if (label[i]==0) ++zero_count;
	  else {
	    label_count[0][nlabel_flag]=label[i];
	    if (nlabel_flag!=0) {
	      for (j=0;j!=nlabel_flag;++j) {
		if (label_count[0][nlabel_flag]!=label_count[0][j]) ;
		else break;
	      }
	      if (j==nlabel_flag) ++nlabel_flag;
	      if (nlabel_flag==MAX_LABEL) {
		printf("limit of max number of labels in a voxel is reached at [%3d,%3d,%3d] (0-offset)\n",c,r,s);
		break;
	      }
	    }
	    else ++nlabel_flag;
	  }
	}

	//find the number of subjects in each non-zero label
	for (i=0;i!=numsub;++i) {
	  if (label[i]) {
	    for (j=0;j!=nlabel_flag;++j) {
	      if (label[i]==label_count[0][j]) {
		++label_count[1][j];
		break;
	      }
	    }
	  }
	}

	//ignore non-GM labels by zeroing their occurances
	for (i=0;i!=nlabel_flag;++i) {
	  if (label_count[0][i]>=11100 || (label_count[0][i]>=47 && label_count[0][i]<=59) || (label_count[0][i]>=26 && label_count[0][i]<=27) || (label_count[0][i]>=17 && label_count[0][i]<=18) || (label_count[0][i]>=8 && label_count[0][i]<=12)) ;
	  else {
	    non_gm_cnt+=label_count[1][i];
	    label_count[1][i]=0;
	    //printf("non-GM label in voxel [%3d,%3d,%3d] (0-offset)\n",c,r,s);
	  }
	}

	//find max of non-zero label count and compare to zero count
	max=label_count[1][0];
	max_loc=0;
	max_count=1;
	if (nlabel_flag) {
	  for (i=1;i!=nlabel_flag;++i) {
	    if (label_count[1][i]>max) {
	      max=label_count[1][i];
	      max_loc=i;
	    }
	    else if (label_count[1][i]==max) ++max_count;
	  }
	  if (max_count>1) {
	    // random label assignment if there're more than 1 label 
	    srand(time(NULL));
	    subc=rand() % (max_count) +1;
	    for (j=0,i=0;i!=nlabel_flag;++i) {
	      if (label_count[1][i]==max) ++j;
	      if (j==subc) {
		max_loc=i;
		break;
	      }
	    }
	    //printf("%d labels have equal occurance in voxel [%3d,%3d,%3d] (from 0)\n",max_count,c,r,s);
	  }
	}

	total_gm=numsub-zero_count-non_gm_cnt;
	if (max==0) seg=stats=0;
	else {
	  //if (max==zero_count) printf("ambiguous label (equal vote for zero and non-zero label) in voxel [%3d,%3d,%3d] (0-offset)\n",c,r,s);
	  seg=label_count[0][max_loc];
	  stats=max;
	  /*
	  // account for cerebelum left/right issue on confidence index
	  subc=0;
	  for (i=0;i!=nlabel_flag;++i) {
	    if (label_count[0][i]==8 || label_count[0][i]==47) ++subc;
	  }
	  if (subc==2) total_gm=max; 
	  */
	}
	
	//write data
	if ((fwrite(&seg,sizeof(int),1,fpseg))!=1) printf("fwrite error of label at [%3d,%3d,%3d] (0-offset)\n",c,r,s);
	if ((fwrite(&stats,sizeof(int),1,fpstat))!=1) printf("fwrite error of stat at [%3d,%3d,%3d] (0-offset)\n",c,r,s);
	if ((fwrite(&total_gm,sizeof(int),1,fpcount))!=1) printf("fwrite error of stat-total at [%3d,%3d,%3d] (0-offset)\n",c,r,s);

	/*
	//debug 
	if (c==134 && r==69 && s==88) {
	  printf("labels at [%3d,%3d,%3d] (from 0) are:\n",c,r,s);
	  for (i=0;i!=numsub;++i) printf("%d ",label[i]);
	  printf("\n\n");

	  printf("zero_count=%d, non-gm_count=%d, nlabel_flag=%d\n\n",zero_count,non_gm_cnt,nlabel_flag);
	  for (i=0;i!=nlabel_flag;++i) printf("label %d count:%d\n",label_count[0][i],label_count[1][i]);
	  printf("max_loc=%d, max=%d, seg=%d\n",max_loc,max,seg);
	}	
	*/
      }
  
  fclose(fpseg);
  fclose(fpstat);
  fclose(fpcount);
  fclose(fpsublist);
  for (i=0;i!=numsub;++i) fclose(fpsub[i]);

  return 0;
}



	
