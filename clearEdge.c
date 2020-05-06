#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define ENDSIZE 100
#define MINGRPLNG 10
#define FAISTD 1
#define MAXGRPSLICE 20
#define EXTNEWHDR "_RimRemoved.hdr"
#define EXTNEWIMG "_RimRemoved.img"
#define TRHDR  "TR.hdr"
#define TRIMG  "TR.img"
#define FAHDR  "FA.hdr"
#define FAIMG  "FA.img"
#define EXTEDGEHDR "_EDGE.hdr"
#define EXTEDGEIMG "_EDGE.img"
#define EXTHDR ".hdr"
#define EXTIMG ".img"
#define HDRBYTES 348
#define MINTRPIX -3206.9250
typedef struct pix pixE;
typedef struct pixGroup pixG;

struct pixGroup {
   pixE *first;
   pixE *last;
   pixE *bLast;
   int length;
};

struct pix {
    int slice;
    int row;
    int column;
    pixE *next;
};

char *fileNames[8];
FILE *mainFAHdrFd = NULL;
FILE *mainFAImgFd = NULL;
FILE *newFAHdrFd = NULL;
FILE *newFAImgFd = NULL;
FILE *mainTRHdrFd = NULL;
FILE *mainTRImgFd = NULL;
FILE *newTRHdrFd = NULL;
FILE *newTRImgFd = NULL;
FILE *edgeLocations = NULL;
short int dim[5];
//vox_units[4]  
//cal_units[4]
short int datatype;
short int bitpix;
float vox_offset;
float calMax;
float calMin;
int glmax;
int glmin;
pixE *edgelist;  
pixE *end;
float maxpixfa;
int maxDepth;
int headersize;

void checkPix(float *pixels, float *current, int stdx, int cx, int cy, int depth,int z);
int nullpix(float pixel);
int evaluateFAH(void);
int evaluateFAI(void);
int evaluateTRH(void);
int evaluateTRI(void);
int checkFile(char *name, char *extension);
void nameFiles(char *argv[]);
void evaluate(void);
int interpHdr(char *bytes);
int checkHdrVals(void);
int interpTRI(float *pixels, int x, int y, int z);
void printEdges(void);
int interpFAI(float *pixels, int x, int y, int z, pixE **iterator);
int tooHigh(float pixel);
int samePix(pixE *pix1, pixE *pix2);

main(int argc, char *argv[]){
  if(argc!=7){
    printf("Please enter the following when executing:\nfilename_FA.hdr filename_FA.img filename_TR.hdr filename_TR.img maxIntensity maxDepthOfSearch\n");
    return;
  }
  else if(checkFile(argv[1],FAHDR)<0||checkFile(argv[2],FAIMG)<0||checkFile(argv[3],TRHDR)<0||checkFile(argv[4],TRIMG)<0)
   return;
    
  else if((maxpixfa = atof(argv[5]))==0){
    printf("Invalid Maximum Pixel Intensity For FA_NEW.img\n");
    printf("Recommended Intensity of 0.90 - 1.2\n");
    return;
  }
  maxDepth = atoi(argv[6]);
  if(maxDepth<0||maxDepth>1000){
    printf("Invalid Max Depth of Search\n");
    printf("Max depth of search must be greater than 0\n");
    return;
  }
  else{
    nameFiles(argv);
    evaluate();
  } 
 
}

int checkFile(char *name, char *extension){
  return 1;
  int location;
  char *base = basename(name);
  printf("%s\n",base);
  if((location = strcspn(base,"."))==strlen(base)){
    printf("Invalid %s File\n",extension);
    return -1;
  }
  if(strcmp((base+location-2),extension)!=0){
    printf("Invalid %s File\n",extension);
    return -1;
  }
  return 1;
}

void nameFiles(char *argv[]){
  char *temp = NULL;
  fileNames[0] = malloc(strlen(argv[1]));
  fileNames[1] = malloc(strlen(argv[2]));
  strcpy(fileNames[0],argv[1]);
  strcpy(fileNames[1],argv[2]);
  fileNames[2] = malloc(strlen(fileNames[0])+12);
  fileNames[3] = malloc(strlen(fileNames[1])+12);
  strcpy(fileNames[2],fileNames[0]);
  strcpy(fileNames[3],fileNames[1]);
  
  temp = fileNames[2];
  fileNames[2][strlen(fileNames[2])-sizeof(EXTIMG)+1] = NULL;  
  strcat(fileNames[2],EXTNEWHDR);
 

  temp = fileNames[3];
  fileNames[3][strlen(fileNames[3])-sizeof(EXTIMG)+1] = NULL;
  strcat(fileNames[3],EXTNEWIMG);
  
  fileNames[4] = malloc(strlen(argv[3]));
  fileNames[5] = malloc(strlen(argv[4]));
  strcpy(fileNames[4],argv[3]);
  strcpy(fileNames[5],argv[4]);
  fileNames[6] = malloc(strlen(argv[3])+5);
  fileNames[7] = malloc(strlen(argv[4])+5);
  strcpy(fileNames[6],fileNames[4]);
  strcpy(fileNames[7],fileNames[5]);
  
  temp = fileNames[6];
  fileNames[6][strlen(fileNames[6])-sizeof(EXTHDR)+1] = NULL;
  strcat(fileNames[6],EXTEDGEHDR);
  
  temp = fileNames[7];
  fileNames[7][strlen(fileNames[7])-sizeof(EXTIMG)+1] = NULL;
  strcat(fileNames[7],EXTEDGEIMG);

  //printf("0 %s\n1 %s\n2 %s\n3 %s\n4 %s\n5 %s\n6 %s\n7 %s\n",fileNames[0],fileNames[1],fileNames[2],fileNames[3],fileNames[4],fileNames[5],fileNames[6],fileNames[7]);

} 

void evaluate(void){
    if((mainFAHdrFd = fopen(fileNames[0],"r"))<0)
      printf("Error opening FA.hdr file\n");
    else if((mainFAImgFd = fopen(fileNames[1],"r"))<0)
      printf("Error opening FA.img file\n");
    else if((mainTRHdrFd = fopen(fileNames[4],"r"))<0)
      printf("Error opening TR.hdr file\n");
    else if((mainTRImgFd = fopen(fileNames[5],"r"))<0)
      printf("Error opening TR.img file\n");
    else{
      if(evaluateTRH()<0)
        return;
      else if(evaluateTRI()<0)
        return;
      else if(evaluateFAH()<0)
        return;
      else if(evaluateFAI()<0)
        return;
    }      
}

int evaluateFAH(void){
  char bytes[HDRBYTES];
  if(fread(bytes,1,HDRBYTES,mainFAHdrFd)<HDRBYTES){
    printf("Error reading FA.hdr file\n");
    return -1;
  }
  else{
    if(interpHdr(bytes)<0){
      return -1;
    }
    else{
      if((newFAHdrFd = fopen(fileNames[2],"wb"))<0){
        printf("Error generating new FA.hdr file\n");
        return -1;
      }
      else{
        if(fwrite(bytes,1,HDRBYTES,newFAHdrFd)<HDRBYTES){
          printf("Error writing new FA.hdr file\n");
          return -1;
        }
      }
    }
  }
  return 1;
}
int interpHdr(char *bytes){
  headersize = *((int *)&bytes[0]);
  dim[0] = *((short int *)&bytes[40]); 
  dim[1] = *((short int *)&bytes[42]);//dim 1 -3 correspond to x,y,z of image; may need to be changed
  dim[2] = *((short int *)&bytes[44]);
  dim[3] = *((short int *)&bytes[46]);
  //dim[2] = 160;
  //dim[3] = 133;
  dim[4] = *((short int *)&bytes[48]);
  //vox_units[4]  
  //cal_units[4]
  datatype = (short int)bytes[70];
  bitpix = (short int)bytes[72];
  vox_offset = bytes[108];
  calMax = (float)bytes[124];
  calMin = (float)bytes[128];
  glmax = (float)bytes[140];
  glmin = (float)bytes[144];
  if(checkHdrVals()<1){
    return -1;
  }
  return 1;
}

int checkHdrVals(void){
  //printf("HDR: %i\n%d\nx: %d\ny: %d\nz: %d\n%d\n%d\n%d\n%f\n%f\n%f\n%i\n%i\n\n",headersize,dim[0],dim[1],dim[2],dim[3],dim[4],datatype,bitpix,vox_offset,calMax,calMin,glmax,glmin);
return 1;
}

int tooHigh(float pixel){
  return pixel > maxpixfa;
}

void checkPix(float *pixels, float *current, int stdx, int cx, int cy, int depth,int z){
  if(depth==maxDepth||!tooHigh(*current))
   return;
  else{
    *current = 0;
    depth++;
    checkPix(pixels,&pixels[cy*stdx+cx+1],stdx,(cx+1),cy,depth,z);
    checkPix(pixels,&pixels[(cy-1)*stdx+cx],stdx,cx,(cy-1),depth,z);
    checkPix(pixels,&pixels[(cy+1)*stdx+cx],stdx,cx,(cy+1),depth,z);
    checkPix(pixels,&pixels[cy*stdx+cx-1],stdx,(cx-1),cy,depth,z);
    checkPix(pixels,&pixels[(cy-1)*stdx+cx+1],stdx,(cx+1),(cy-1),depth,z);
    checkPix(pixels,&pixels[(cy-1)*stdx+cx-1],stdx,(cx-1),(cy-1),depth,z);
    checkPix(pixels,&pixels[(cy+1)*stdx+cx+1],stdx,(cx+1),(cy+1),depth,z);
    checkPix(pixels,&pixels[(cy+1)*stdx+cx-1],stdx,(cx-1),(cy+1),depth,z);
  }
  return;
}  

int interpFAI(float *pixels,int x, int y, int z, pixE **iterator){
  float *current = NULL;
  while((*iterator)->slice==z&&(*iterator)!=end){
    current = &pixels[(*iterator)->row*x+(*iterator)->column];
    checkPix(pixels,current, x, (*iterator)->column, (*iterator)->row,0,(*iterator)->slice);
    (*iterator) = (*iterator)->next;
  }
  return 1;
}

int samePix(pixE *pix1, pixE *pix2){
  if(pix1->slice==pix2->slice&&pix1->row==pix2->row&&pix1->column==pix2->column)
    return 1;
return -1;
}

int evaluateFAI(){
  float blank [dim[1]*dim[2]];
  int read;
  float pixels[dim[1]*dim[2]];
  pixE *iterator = edgelist;
  int slice = 0,i;
  if((newFAImgFd = fopen(fileNames[3],"wb"))<0){
    printf("Error generating new FA.img file\n");
    return -1;
  }
  else{
    for(i = 0; i<dim[3]; i++){
      if((read = fread(pixels,sizeof(float),dim[1]*dim[2],mainFAImgFd))<dim[1]*dim[2]){
          printf("Error Reading Slice: %d of FA.img File\n",i);
          printf("%i\n",read);
        if(read == 0 && i>=ENDSIZE)
          fwrite(blank,sizeof(float),dim[1]*dim[2],newTRImgFd);
        else if(read != 0 && i>=ENDSIZE){
          fwrite(pixels,sizeof(float),read,newTRImgFd);
          fwrite(blank,sizeof(float),dim[1]*dim[2]-read,newTRImgFd);
        }
        else{
          return -1;
        }
      }
      else if(interpFAI(pixels,dim[1],dim[2],i,&iterator)<0){
	return -1;
      }
      else if(fwrite(pixels,sizeof(float),dim[1]*dim[2],newFAImgFd)<dim[1]*dim[2]){
        printf("Error Writing Slice: %d of FA_NEW.img File",i);
	return -1;
      }

    }	
  }
  return 1;
}

int evaluateTRI(){
  int read;
  float blank [dim[2]*dim[1]];
  edgelist = malloc(sizeof(pixE));
  edgelist->next = NULL;
  end = edgelist;
  float pixels[dim[2]*dim[1]];
  int rowOff = dim[1];
  int i = 0,b;
  if((newTRImgFd = fopen(fileNames[7],"wb"))<0){
    printf("Error generating new TR.img file\n");
    return -1;
  }
  else{
    for(i = 0; i<dim[3]; i++){
      if((read = fread(pixels,sizeof(float),dim[1]*dim[2],mainTRImgFd))<dim[1]*dim[2]){
	  printf("Error Reading Slice: %d of TR.img File\n",i);
	printf("%i\n",read);
	if(read == 0 && i>=ENDSIZE)
	  fwrite(blank,sizeof(float),dim[1]*dim[2],newTRImgFd);
	else if(read != 0 && i>=ENDSIZE){
	  fwrite(pixels,sizeof(float),read,newTRImgFd);
	  fwrite(blank,sizeof(float),dim[1]*dim[2]-read,newTRImgFd);
	}
	else{
          return -1; 
        }
      }
      else if(interpTRI(pixels,dim[1],dim[2],i)<0)
        return -1;
      else if(fwrite(pixels,sizeof(float),dim[1]*dim[2],newTRImgFd)<dim[2]*dim[1]){
        printf("Error Writing Slice: %d of TR_EDGE.img File",i);
        return -1;
      }
    }
  }
//printEdges();
return 1;
}

void printEdges(void){
  pixE *step = edgelist;
  while(step->next!=NULL){
  printf("Slice: %d, Row: %d, Column: %d\n",step->slice,step->row,step->column);
  step = step->next;
  }
}

int nullpix(float data){
   return (data<0.0001&&data>-.0001);
}

int interpTRI(float *pixels,int x,int y,int z){
  int i, b;
  pixE *current = end;
  for(i = 0; i<y; i++){
    for(b = 0; b<x-1; b++){
      if(!nullpix(pixels[i*x+b])){
	if(nullpix(pixels[i*x+b-1])||nullpix(pixels[i*x+b+1])){
            pixels[i*x+b] = MINTRPIX ;
            end->row = i;
	    end->column = b;
	    end->slice = z;
	    end->next = malloc(sizeof(pixE));
	    end = end->next;
	    end->next = NULL;
        }
	else if(nullpix(pixels[(i-1)*x+b])){
	    pixels[i*x+b] = MINTRPIX;
	    end->row = i;
            end->column = b;
            end->slice = z;
            end->next = malloc(sizeof(pixE));
            end = end->next;
            end->next = NULL;
        } 
        else if(nullpix(pixels[(i+1)*x+b])){
	    pixels[i*x+b] = MINTRPIX;
	    end->row = i;
            end->column = b;
            end->slice = z;
            end->next = malloc(sizeof(pixE));
            end = end->next;
            end->next = NULL;
        }
	else if(i<y-1&&i>0){
          if(nullpix(pixels[(i+1)*x+b+1])||nullpix(pixels[(i-1)*x+b+1])||nullpix(pixels[(i+1)*x+b-1])||nullpix(pixels[(i-1)*x+b-1])){
            pixels[i*x+b] = MINTRPIX;
	    end->row = i;
            end->column = b;
            end->slice = z;
            end->next = malloc(sizeof(pixE));
            end = end->next;
            end->next = NULL;
          }
        }
      }
    }
  }
  //checkHoles(current,pixels);
  return 1;
}

int evaluateTRH(){
  char bytes[HDRBYTES];
  if(fread(bytes,1,HDRBYTES,mainTRHdrFd)<HDRBYTES){
    printf("Error reading TR.hdr file\n");
    return -1;
  }
  else{
    if(interpHdr(bytes)<0){
      return -1;
    }
    else{
      if((newTRHdrFd = fopen(fileNames[6],"wb"))<0){
        printf("Error generating new TR.hdr file\n");
        return -1;
      }
      else{
        if(fwrite(bytes,1,HDRBYTES,newTRHdrFd)<HDRBYTES){
          printf("Error writing new TR.hdr file\n");
          return -1;
        }
      }
    }
  }
  return 1;
} 
