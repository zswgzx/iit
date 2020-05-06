#include <stdio.h>
#include <stdlib.h>

#define MYDATATYPE float
#define HEAD 352
#define MAX_FILENAME 60

int main(int argc,char *argv[])
{
  FILE *fpIn,*fpOut;
  MYDATATYPE *buf1,*buf2;
  char filename[MAX_FILENAME];
  int x,y,z,xres,yres,zres,nxres,nyres,nzres,cutL,cutR,cutT,cutB,cutD,cutU,res;

    //check # of arguments
    if (argc != 11) {
        printf("%s -- shrink 3D single volume image (default float) to desired dimension. apply to raw image regardless of orientation code in nii header, open the image in imageJ to check actual orientation issue\n", argv[0]);
        printf("Input: rawfile.img (float, no header)\nOutput:new_rawfile.img\n");
        printf("Usage: %s <rawfile> <XRES> <YRES> <ZRES> <XRES cut left> <XRES cut right> <YRES cut top> <YRES cut bottom> <ZRES cut down> <ZRES cut up>\n", argv[0]);
        printf("e.g.: %s rawfile.img 256 256 256 37 37 19 19 37 37\n", argv[0]);
        exit(1);
    }

    xres=atoi(argv[2]);
    yres=atoi(argv[3]);
    zres=atoi(argv[4]);
    cutL=atoi(argv[5]);cutR=atoi(argv[6]);
    cutT=atoi(argv[7]);cutB=atoi(argv[8]);
    cutD=atoi(argv[9]);cutU=atoi(argv[10]);

    nxres=xres-cutL-cutR;
    nyres=yres-cutT-cutB;
    nzres=zres-cutD-cutU;

    //didn't check if cutL etc. are no less than zero!

    //allocate buff for source and output data
    buf1 = (MYDATATYPE *) malloc(sizeof(MYDATATYPE)*xres*yres);
    if (buf1 == NULL) {
      printf("*** Failed to allocate buffer 1\n");
      exit(1);
    }

    buf2 = (MYDATATYPE *) malloc(sizeof(MYDATATYPE)*nxres*nyres);
    if (buf2 == NULL) {
      printf("*** Failed to allocate buffer 2\n");
      exit(1);
    }

    //open files to read/write
    fpIn=fopen(argv[1],"rb");
    if (fpIn == NULL) {
      printf("*** Failed to open %s\n",argv[1]);
      exit(1);
    }

    sprintf(filename,"new_%s",argv[1]);
    fpOut=fopen(filename,"wb");
    if (fpOut == NULL) {
      printf("*** Failed to open new_%s\n",argv[1]);
      exit(1);
    }

    //skip slices down
    res=fseek(fpIn,cutD*xres*yres*sizeof(MYDATATYPE),SEEK_SET);
    if (res) {
      printf("fseek error 1\n");
      exit(1);
    }

    //main part
    for (z=0;z<nzres;++z) {
      // read raw file
      if (fread(buf1,sizeof(MYDATATYPE),xres*yres,fpIn)!=xres*yres) {
	fprintf(stderr,"error reading slice %d(1 offset) of raw file\n",z+1+cutD);
	exit(1);
      }

      // cut edges
      for (y=0;y<nyres;y++)
	for (x=0;x<nxres;x++)
	  buf2[x+y*nxres]=buf1[(x+cutL)+(y+cutT)*xres];

      // write to output file
      if (fwrite(buf2,sizeof(MYDATATYPE),nxres*nyres,fpOut)!=nxres*nyres) {
	fprintf(stderr,"error writing slice %d(0 offset) to output\n",z);
	exit(1);
      }
    }

    free(buf1);
    free(buf2);
    fclose(fpIn);
    fclose(fpOut);

    return 0;
}



 

