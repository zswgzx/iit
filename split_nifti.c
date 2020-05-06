/*SZ last updated: 
 */ 

#include <stdio.h>
#include <stdlib.h>

#define HEAD 352

/*=================================================================*/
int main(int argc, char *argv[]) 
{
  FILE *fpraw,*fpsplit;
  char filename[60],ctmp,*buf;
  int i;
  long filesize;
  size_t res;

    //check # of arguments
    if (argc != 3) {
        printf("\n%s -- split nifti(.nii, hdr size: 352 byte) file to hdr and data pair.\nNOTICE:THE IMAGE PAIR WON'T BE DISPLAYED PROPERLY IN IMJ USING 'IMPORT' or 'OPEN' option!\n", argv[0]);
	printf("\nInput: filename.nii\n");
	printf("Output: user_specified.{hdr,img}, size of hdr=352 byte.\n\n");
        printf("Usage: %s <raw> <split>\n", argv[0]);
        printf("e.g.: %s filename.nii user_specified\n", argv[0]);
        exit(1);
    }
    
    //open source file to read
    sprintf(filename,"%s",argv[1]);
    fpraw = fopen(filename, "rb");
    if (fpraw == NULL) {
      printf("*** Failed to open %s\n", filename);
      exit(1);
    }

    fseek(fpraw,0,SEEK_END);
    filesize=ftell(fpraw);
    //printf("file size=%ld\n",filesize);
    rewind(fpraw);

    buf = (char*) malloc (sizeof(char)*(filesize-HEAD));
    if (buf == NULL) {fputs ("Memory error\n",stderr); exit (2);}
    
    //open header file to write
    sprintf(filename,"%s.hdr",argv[2]);
    fpsplit = fopen(filename, "wb");
    if (fpsplit == NULL) {
      printf("*** Failed to open %s\n", filename);
      exit(1);
    }

    for (i=1;i<=HEAD;i++) {
      if (fread(&ctmp,sizeof(char),1,fpraw) !=1) {
	fprintf(stderr,"error reading header from source img at %d\n",i);
	exit(1);
      }
      if (i==346) ctmp='i';
      fwrite(&ctmp,sizeof(char),1,fpsplit);
    }

    fclose(fpsplit);

    //open data file to write
    sprintf(filename,"%s.img",argv[2]);
    fpsplit = fopen(filename, "wb");
    if (fpsplit == NULL) {
      printf("*** Failed to open %s\n", filename);
      exit(1);
    }

    res = fread (buf,1,filesize-HEAD,fpraw);
    if (res != (filesize-HEAD)) {fputs ("Reading error",stderr); exit (3);}
    fwrite(buf,1,filesize-HEAD,fpsplit);
 
    fclose(fpsplit);
    free(buf);
    
  return 0;
}
