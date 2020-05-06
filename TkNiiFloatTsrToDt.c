/* SZ */

#include <stdio.h>
#include <stdlib.h>

#define HEAD 352
typedef float DATA_TYPE;

int main(int argc, char *argv[])
 {
   FILE *fpraw,*fpout;
   DATA_TYPE *buf,temp;
   char filename[100];
   unsigned int x,y,z,res,xres,yres,zres,volsize;

   if (argc != 5) {
     printf("%s -- Convert symmetric tensor volumes from DTI-TK to DTI-GUI format, default datatype: float. Adjust diffusivity back before running this program. Background tsr may still be isotropic without correction.\n", argv[0]);
     printf("\nInput: tensor.nii (lower triangle order, a.k.a xx yx yy zx zy zz)\n");
     printf("Output: tensor.img (xx yy zz xy xz yz order)\n");
     printf("\nUsage: %s <vol> <xres> <yres> <zres>\n", argv[0]);
     printf("e.g.: %s tensor 182 218 182.\n", argv[0]);
     exit(1);
   }

   xres=atoi(argv[2]);
   yres=atoi(argv[3]);
   zres=atoi(argv[4]);
   volsize=xres*yres*zres;

   buf = malloc(6*sizeof(DATA_TYPE)*volsize);
   if (buf == NULL) {
     printf("*** Failed to allocate buffer\n");
     exit(1);
   }
    
   sprintf(filename,"%s.nii",argv[1]);
   fpraw = fopen(filename, "rb");
   if (fpraw == NULL) {
     printf("*** Failed to open %s", filename);
     exit(1);
   }

   res = fseek(fpraw, HEAD, SEEK_SET);
   if (res) {
     printf("*** Failed to seek1 in %s", filename);
     exit(1);
   }

   res = fread(buf, sizeof(DATA_TYPE), 6*volsize, fpraw);
   if (res != 6*volsize) {
     printf("*** Failed to read in %s.nii\n", argv[1]);
     printf("Got %d elements\n", res);
     exit(1);
   }

   sprintf(filename,"%s.img",argv[1]);
   fpout = fopen(filename, "wb");
   if (fpout == NULL) {
     printf("*** Failed to open %s", filename);
     exit(1);
   }
   
   for (z=0;z<zres;z++)
     for (y=0;y<yres;y++)
       for(x=0;x<xres;x++) {
	 temp=buf[z*xres*yres+y*xres+x];
	 fwrite(&temp, sizeof(DATA_TYPE), 1, fpout);
	 
	 temp=buf[z*xres*yres+y*xres+x+2*volsize];
	 fwrite(&temp, sizeof(DATA_TYPE), 1, fpout);

	 temp=buf[z*xres*yres+y*xres+x+5*volsize];
	 fwrite(&temp, sizeof(DATA_TYPE), 1, fpout);

	 temp=buf[z*xres*yres+y*xres+x+volsize];
	 fwrite(&temp, sizeof(DATA_TYPE), 1, fpout);

	 temp=buf[z*xres*yres+y*xres+x+3*volsize];
	 fwrite(&temp, sizeof(DATA_TYPE), 1, fpout);

	 temp=buf[z*xres*yres+y*xres+x+4*volsize];
	 fwrite(&temp, sizeof(DATA_TYPE), 1, fpout);
	 
       }
   
   fclose(fpraw);
   fclose(fpout);
   free(buf);
   return 0;
}
