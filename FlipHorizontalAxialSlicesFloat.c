/*SZ last updated: Oct 23,09 

Notice: . */ 

#include <stdio.h>
#include <stdlib.h>

#define HEAD 352
typedef float MYDATATYPE;

int main(int argc, char *argv[]) 
{
  FILE *fp1,*fp2;
  MYDATATYPE *buf1,*buf2;
  int x,y,j,XRES,YRES,ZRES,ret;

    //check # of arguments
    if (argc != 6) {
        printf("%s -- flip 3D char axial slice images (.img format) in left-right direction.\n", argv[0]);
        printf("Usage: %s <infile> <outfile> <xdim> <ydim> <zdim>\n", argv[0]);
        printf("e.g.: %s problem_vox.img new.img 256 256 256 \n", argv[0]);
        exit(1);
    }

    XRES=atoi(argv[3]);
    YRES=atoi(argv[4]);
    ZRES=atoi(argv[5]);

    //allocate buff from source and read data
    buf1 = (MYDATATYPE *) malloc(sizeof(MYDATATYPE)*XRES*YRES);
    if (buf1 == NULL) 
      {
	printf("*** Failed to allocate buffer1\n");
	exit(1);
      }

    buf2 = (MYDATATYPE *) malloc(sizeof(MYDATATYPE)*XRES*YRES);
    if (buf2 == NULL) 
      {
	printf("*** Failed to allocate buffer2\n");
	exit(1);
      }

    // open raw file
    fp1=fopen(argv[1],"rb");    
    if (fp1 == NULL) 
      {
	fprintf(stderr, "\nError opening file %s\n",argv[1]);
	exit(1);
      }
    // open output file
    fp2=fopen(argv[2],"wb");
    if (fp2 == NULL) 
      {
	fprintf(stderr, "\nError opening file %s\n",argv[2]);
	exit(1);
      }

    //main loop
    for (j=0;j<ZRES;j++)
      {
	ret=fread(buf1,sizeof(MYDATATYPE),XRES*YRES,fp1);
	if (ret != XRES*YRES) 
	  {
	    fprintf(stderr, "\nError reading slice %d from %s\n",j+1,argv[1]);
	    exit(1);
	  }

	// flip left-right of image
	for (y=0;y<YRES;y++)
	  for (x=0;x<XRES;x++)
	    buf2[x+YRES*y]=buf1[(XRES-1)-x+YRES*y];

	ret=fwrite(buf2,sizeof(MYDATATYPE),XRES*YRES,fp2);
	if (ret != XRES*YRES) 
	  {
	    fprintf(stderr, "\nError writing slice %d to %s\n",j+1,argv[2]);
	    exit(1);
	  }
      }
    free(buf1);
    free(buf2);
    //clean up
    fclose(fp1);
    fclose(fp2);
    return 0;
}
