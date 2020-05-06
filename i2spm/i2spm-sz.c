/*

i2spm.c

082001 BPR:  Added TIGNORE option.
030701 BPR

Most recent change: fixed to accommodate such things as 5.0S-125.0S in
the XRANGE, YRANGE, and ZRANGE.

This program will convert single-slice GE image files (assuming a filename of
prefix.nnn.snn) to a series of spm-compatible 3D Analyze-type .img/.hdr files,
one per time point.

This code assumes that the order of coordinates, from fastest-changing
to slowest-changing, is X Y Z in the input data files. Also, it writes
this order (LPI) for SPM:

X increases from left to right
Y increases from posterior to anterior
Z increases from inferior to superior

The origin written is always (0,0,0) -- this means that, in general,
anatomical images like the usual spgr and functional data from epi
images will NOT be aligned in spm unless the Talairach normalization
is made. AFNI is the way to go for viewing unnormalized data.

Data is written out in the byte order native to the machine running the
program.

The input file is i2spm.in, should be located in the current directory,
and should contain this stuff:

# Comments

# i2spm.in
#
# This file contains information vital for the i2spm program,
# which converts GE I.* files and epirecon I.* files to SPM-
# compatible Analyze format.
#
# Note: Data files are written in the native byte order for SPM
# compatibility. The Analyze file format contains no byte order
# info in the header, and SPM expects native byte order when it
# opens files.

# Image slice file header size
HDRSIZE 80

# Slice files; include path if necessary (usually not). Typically,
#   I.%03i.s%02i  for EPI  (TIMEPTS must be > 1 for 2 %i specifiers;
#                           the first is time, the second slice)
#   I.%03i        for anat (TIMEPTS must be 1 for 1 %i specifier,
#                           which is slice number)
# Examples:
# INFILEMASK /atlas2/rogers/spm/jdesousa/S7/I.%03i.s%02i
# INFILEMASK I.%03i.s%02i
INFILEMASK I.%03i.s%02i

# Output file(s); as above. Typically,
#   S5-motor.%03i  for EPI
#   S2-spgr       for anat
# The extensions .img and .hdr are always automatically added.
OUTFILEMASK motor-S7.%03i

# Number of slices (also change the appropriate one of PIXEL_X,

# PIXEL_Y, or PIXEL_Z if you change this!)
SLICES 20

# Total number of time points, including the ones you'll ignore
TIMEPTS 134

# Ignore this many time points at the beginning; 0 to use them all
TIGNORE 5

# Slice axis (COR, AXI, SAG)
AXIS COR

# Number of pixels in each direction:
#  X is right to left
#  Y is anterior to posterior
#  Z is inferior to superior 
# One of these should be the same as SLICES.
PIXEL_X 64
PIXEL_Y 20
PIXEL_Z 64

# Orientation and origin info; extent is from voxel edge
# to voxel edge in the slice plane (i.e., FOV), and from
# center to center along the slice axis.
# Numbers are mm.
#
# The conversion to spm format uses this info to find
# the voxel size and the data file byte sequence; the
# conversion to afni format uses it to determine the origin.
#
# Order matters! For example, putting 120R-120L when it
# should be 120L-120R flips the image left-right. Use R-L,
# A-P, and S-I for the two axes you don't know from the
# scan info sheet.
#
# X is the sagittal axis, Y the coronal, and Z the axial.
XRANGE 120R-120L
YRANGE 69.8P-82.2A
ZRANGE 120S-120I

# SPM file scale factor. Generally, use 1.
SCALE 1

# Set OVERWRITE to YES to overwrite existing output files
# if you're running i2spm a second time. Set to NO to avoid
# overwriting. Do not comment this line out, though.
OVERWRITE YES

# Data type (only SHORT).
DATUM SHORT

*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
/* #include "nrc_bpr/nrutil.h" */

#define MAXLEN 255
#define INPUTS 16
#define NR_END 1
#define FREE_ARG char* 

short ***s3tensor(long nrl, long nrh, long ncl, long nch, long ndl, long ndh);
void free_s3tensor(short ***t, long nrl, long nrh, long ncl, long nch, long ndl, long ndh);

void
swap16(void* p) {
  /* Swap two bytes */
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

int
main(int argc, char **argv) {

  FILE *infp, *datafp, *outfp, *logfp;
  char *infname = "i2spm.in";
  char *logfname = "i2spm.log";
  char datafname[MAXLEN];
  char outfname[MAXLEN];
  char hdrfname[MAXLEN];
  char hdrfprefix[MAXLEN];
  char line[MAXLEN];
  char item[15], value[MAXLEN];
  int inflags[INPUTS];
  int incount=0;
  int i;

  char axis[4];
  int hdrsize;
  char inmask[MAXLEN], outmask[MAXLEN];
  int slices;
  int timepts;
  short pixel_x, pixel_y, pixel_z;
  float voxel_x, voxel_y, voxel_z;
  char xrange[20], yrange[20], zrange[20];
  float ex1, ex2;
  char dir1, dir2;
  float exr, exl, exa, exp, exi, exs;
  char orient[4] = "   ";
  float scale;
  char datum[6];
  char overwrite;

  short ***data;
  short twobytes;
  signed int swapflag;
  int x,y,z,t,s;
  int refx, refy, refz;
  int hdrfilesize;
  int len, pos;
  int tignore;

  char zeros[84], auxfile[24], short_outfname[18];
  char descrip[80];
  char datatype[10];
  short origin[5];
  const short four=4, one=1, sixteen=16;
  const int glmax=32767, glmin=0;
  const float fzero = 0.0;

  for( i=0 ; i<INPUTS ; i++ )
    inflags[i] = 0;

  // SZ added these below, 10292014
  if (argc!=2) {
      fprintf(stderr,"\nUsage: %s <byte order swap? Y:-1, N:otherwise >", argv[0]);
      fprintf(stderr,"\ne.g. %s -1.\n", argv[0]);
      exit(1);
    }

  swapflag=atoi(argv[1]);
  if (swapflag == -1) printf("Swapping byte order.\n");
  else printf("Keeping original byte order.\n");
  // SZ additions end, 10292014

  /* Open input file */
  infp = fopen(infname, "r");
  if( infp == NULL ) {
    fprintf(stderr, "Could not open input file i2spm.in.\n");
    fprintf(stderr, "Is it in the current directory?\n");
    exit(1);
  }

  /* Open log file */
  logfp = fopen(logfname, "w");
  if( logfp == NULL ) {
    fprintf(stderr, "Could not open log file i2spm.log.\n");
    exit(1);
  }

  /* Parse input file, a line at a time, flagging each item as
     we read it in. Loop until we've got them all */
  while( incount < INPUTS ) {

    if( fgets(line, MAXLEN, infp) == NULL) {
      fprintf(stderr,"File i2spm.in ended before all input was found.\n");
      exit(1);
    }

    if( strncmp(line, "#", 1) == 0 ) continue;
    if( strlen(line) < 4 ) continue;

    if( sscanf(line, "%s %s", item, value) == 2 ) {
      
      if( strcmp(item, "OVERWRITE") == 0 ) {
	overwrite = value[0];
	inflags[14] = 1;
      }
      if( strcmp(item, "AXIS") == 0 ) {
	inflags[0] = 1;
	strcpy(axis, value);
      }
      else if( strcmp(item, "HDRSIZE") == 0 ) {
	inflags[1] = 1;
	hdrsize = atoi(value);
      }
      else if( strcmp(item, "INFILEMASK") == 0 ) {
	inflags[2] = 1;
	strcpy(inmask, value);
      }
      else if( strcmp(item, "OUTFILEMASK") == 0 ) {
	inflags[3] = 1;
	strcpy(outmask, value);
      }
      else if( strcmp(item, "SLICES") == 0 ) {
	inflags[4] = 1;
	slices = atoi(value);
      }
      else if( strcmp(item, "TIMEPTS") == 0 ) {
	inflags[5] = 1;
	timepts = atoi(value);
      }
      else if( strcmp(item, "TIGNORE") == 0 ) {
	inflags[15] = 1;
	tignore = atoi(value);
      }
      else if( strcmp(item, "PIXEL_X") == 0 ) {
	inflags[6] = 1;
	pixel_x = atoi(value);
      }
      else if( strcmp(item, "PIXEL_Y") == 0 ) {
	inflags[7] = 1;
	pixel_y = atoi(value);
      }
      else if( strcmp(item, "PIXEL_Z") == 0 ) {
	inflags[8] = 1;
	pixel_z = atoi(value);
      }
      else if( strcmp(item, "XRANGE") == 0 ) {
	inflags[9] = 1;
	strcpy(xrange,value);
      }
      else if( strcmp(item, "YRANGE") == 0 ) {
	inflags[10] = 1;
	strcpy(yrange,value);
      }
      else if( strcmp(item, "ZRANGE") == 0 ) {
	inflags[11] = 1;
	strcpy(zrange,value);
      }
      else if( strcmp(item, "SCALE") == 0 ) {
	inflags[12] = 1;
	scale = atof(value);
      }
      else if( strcmp(item, "DATUM") == 0 ) {
	inflags[13] = 1;
	if( strcmp(value,"SHORT") && strcmp(value,"FLOAT") ) {
	  fprintf(stderr, "Invalid data type. Must be SHORT or FLOAT.\n");
	  exit(1);
	}
	else
	  strcpy(datum, value);

	if( strcmp(datum, "FLOAT") == 0 ) {
	  fprintf(stderr,
		  "i2spm can't cope with float data right now.\n");
	  exit(1);
	}

      }
    }
    incount = 0;
    for( i=0; i<INPUTS ; i++ )
      incount += inflags[i];

  }

  /* Calculate voxel size */
  /* exr is right extent in mm, exl is left extent in mm,
     ex1 and ex2 are just temporary vars for parsing. orient[0]
     contains R if the first data point in the file is right of
     the last, but L if the first data point is left of the last.
     Similar stuff is true for the other two axes. Recall that 
     xrange is of the form '120R-120L'. */

  sscanf(xrange, "%f%c-%f%c", &ex1, &dir1, &ex2, &dir2);

  if( dir1=='R' || dir1=='r' ) {   /* xrange is R-? */
    if( dir2=='L' || dir2=='l' ) {  /* range is R-L */
      exr = ex1;  /* first number given is R extent */
      exl = ex2;  /* second number is L extent */
      orient[0] = 'R';
    }
    if( dir2=='R' || dir2=='r' ) {  /* xrange is R-R */
      if( ex1<ex2 ) {  /* first number is LEFTmost value */
	exl = -ex1;
	exr = ex2;
	orient[0] = 'L';
      }
      else if( ex1>=ex2 ) {  /* first num is RIGHTmost value; */
	exl = -ex2;          /* lump in the ex1==ex2 case here */
	exr = ex1;
	orient[0] = 'R';
      }
    }
  }
  else if( dir1=='L' || dir1=='l' ) {  /* xrange is L-? */
    if( dir2=='R' || dir2=='r' ) {  /* range is L-R */
      exl = ex1;  /* first number given is L extent */
      exr = ex2;  /* second number is R extent */
      orient[0] = 'L';
    }
    if( dir2=='L' || dir2=='l' ) {  /* xrange is L-L */
      if( ex1<ex2 ) {  /* first number is RIGHTmost value */
	exr = -ex1;
	exl = ex2;
	orient[0] = 'R';
      }
      else if( ex1>=ex2 ) {  /* first num is LEFTmost value; */
	exr = -ex2;           /* lump in the ex1==ex2 case here */
	exl = ex1;
	orient[0] = 'L';
      }
    }
  }
  else {
    fprintf(stderr, "Could not parse XRANGE.");
    exit(1);
  }
  if( axis[0]=='S' )
    voxel_x = (exr+exl) / (slices-1);
  else
    voxel_x = (exr+exl) / pixel_x;


  /* Now for YRANGE */

  sscanf(yrange, "%f%c-%f%c", &ex1, &dir1, &ex2, &dir2);

  if( dir1=='A' || dir1=='a' ) {   /* range is A-? */
    if( dir2=='P' || dir2=='p' ) {  /* range is A-P */
      exa = ex1;  /* first number given is A extent */
      exp = ex2;  /* second number is P extent */
      orient[1] = 'A';
    }
    if( dir2=='A' || dir2=='a' ) {  /* range is A-A */
      if( ex1<ex2 ) {  /* first number is Pmost value */
	exp = -ex1;
	exa = ex2;
	orient[1] = 'P';
      }
      else if( ex1>=ex2 ) {  /* first num is Amost value; */
	exp = -ex2;          /* lump in the ex1==ex2 case here */
	exa = ex1;
	orient[1] = 'A';
      }
    }
  }
  else if( dir1=='P' || dir1=='p' ) {  /* xrange is P-? */
    if( dir2=='A' || dir2=='a' ) {  /* range is P-A */
      exp = ex1;  /* first number given is P extent */
      exa = ex2;  /* second number is A extent */
      orient[1] = 'P';
    }
    if( dir2=='P' || dir2=='p' ) {  /* xrange is P-P */
      if( ex1<ex2 ) {  /* first number is Amost value */
	exa = -ex1;
	exp = ex2;
	orient[1] = 'A';
      }
      else if( ex1>=ex2 ) {  /* first num is Pmost value; */
	exa = -ex2;           /* lump in the ex1==ex2 case here */
	exp = ex1;
	orient[1] = 'P';
      }
    }
  }
  else {
    fprintf(stderr, "Could not parse YRANGE.");
    exit(1);
  }
  if( axis[0]=='C' )
    voxel_y = (exa+exp) / (slices-1);
  else
    voxel_y = (exa+exp) / pixel_y;


  /* Finally, ZRANGE */

  sscanf(zrange, "%f%c-%f%c", &ex1, &dir1, &ex2, &dir2);

  if( dir1=='S' || dir1=='s' ) {   /* range is S-? */
    if( dir2=='I' || dir2=='i' ) {  /* range is S-I */
      exs = ex1;  /* first number given is S extent */
      exi = ex2;  /* second number is I extent */
      orient[2] = 'S';
    }
    if( dir2=='S' || dir2=='s' ) {  /* range is S-S */
      if( ex1<ex2 ) {  /* first number is Imost value */
	exi = -ex1;
	exs = ex2;
	orient[2] = 'I';
      }
      else if( ex1>=ex2 ) {  /* first num is Smost value; */
	exi = -ex2;          /* lump in the ex1==ex2 case here */
	exs = ex1;
	orient[2] = 'S';
      }
    }
  }
  else if( dir1=='I' || dir1=='i' ) {  /* xrange is I-? */
    if( dir2=='S' || dir2=='s' ) {  /* range is I-S */
      exi = ex1;  /* first number given is I extent */
      exs = ex2;  /* second number is S extent */
      orient[2] = 'I';
    }
    if( dir2=='I' || dir2=='i' ) {  /* xrange is I-I */
      if( ex1<ex2 ) {  /* first number is Smost value */
	exs = -ex1;
	exi = ex2;
	orient[2] = 'S';
      }
      else if( ex1>=ex2 ) {  /* first num is Imost value; */
	exs = -ex2;           /* lump in the ex1==ex2 case here */
	exi = ex1;
	orient[2] = 'I';
      }
    }
  }
  else {
    fprintf(stderr, "Could not parse ZRANGE.");
    exit(1);
  }
  if( axis[0]=='A' )
    voxel_z = (exs+exi) / (slices-1);
  else
    voxel_z = (exs+exi) / pixel_z;


  /* Log the values we just read in, so they can be checked */
  fprintf(logfp, "\nFinished reading input file.\n");
  fprintf(logfp, "   Slice axis: %s\n", axis);
  fprintf(logfp, "   Header size: %i\n", hdrsize);
  fprintf(logfp, "   Input file mask: %s\n", inmask);
  fprintf(logfp, "   Output file mask: %s\n", outmask);
  fprintf(logfp, "   Slices: %i\n", slices);
  fprintf(logfp, "   Time points: %i\n", timepts);
  fprintf(logfp, "   Ignoring %i time points at the beginning.\n", tignore);
  fprintf(logfp, "   Number of pixels (x,y,z): %i, %i, %i\n",
	 pixel_x, pixel_y, pixel_z);
  fprintf(logfp, "   L-R extent: %fL, %fR\n", exl, exr);
  fprintf(logfp, "   P-A extent: %fP, %fA\n", exp, exa);
  fprintf(logfp, "   S-I extent: %fS, %fI\n", exs, exi);
  fprintf(logfp, "   Calculated voxel size in mm (x,y,z): %f, %f, %f\n",
	  voxel_x, voxel_y, voxel_z);
  fprintf(logfp, "   Input file data orientation: %s\n", orient);
  fprintf(logfp, "   Scale factor: %f\n", scale);
  fprintf(logfp, "   Data type: %s\n\n", datum);


  /* Make one 3D file for each time point */
  for( t=1+tignore ; t<=timepts ; t++ ) {

    /* Allocate memory for this time point's chunk o' data --
       at the moment, it has to be short data */
    data = s3tensor(1,pixel_x, 1,pixel_y, 1,pixel_z);
    
    /* Read in each slice at this time point */
    for( s=1; s<=slices ; s++ ) {

      if( timepts == 1)
	sprintf(datafname, inmask, s);
      else
	sprintf(datafname, inmask, t, s);

      datafp = fopen(datafname, "rb");
      if( datafp == NULL ) {
	fprintf(stderr, "Could not open file %s.\n",datafname);
	exit(1);
      }

      /* See if we need to swap bytes as data is read */
      /* SZ commented these below, 10292014
      fread(&twobytes,sizeof(short),1,datafp);
      if(twobytes==044515) {
	swapflag=1; 
	fprintf(logfp, "%s is apparently in native byte order. Not swapping.\n",
		datafname);
      }
      else if(twobytes==046511) {
	swapflag=-1;
	fprintf(logfp, "%s is apparently in foreign byte order. Swapping.\n",
		datafname);
      }
      else {
	swapflag=0;
	fprintf(logfp, "Byte order undetermined for %s. Not swapping.\n",
		datafname);
      }
          SZ comments end, 10292014 */
      /* Jump to the beginning of the data */
      fseek(datafp, hdrsize, SEEK_SET);
      
      /* What to do now depends on the slice axis */
      switch( axis[0] ) {
      case 'A':
	/* axial slices; read X,Y */
	for( y=1 ; y<=pixel_y ; y++ ) {
	  for( x=1 ; x<=pixel_x ; x++ ) {
	    if( fread(&twobytes,sizeof(short),1,datafp) != 1 ) {
	      fprintf(stderr, "Error reading %s.\n", datafname);
	      exit(1);
	    }
	    if(swapflag==-1) swap16(&twobytes);
	    refx = x; refy = y; refz = s;
	    if(orient[0]=='R') refx = pixel_x-x+1;
	    if(orient[1]=='A') refy = pixel_y-y+1;
	    if(orient[2]=='S') refz = slices-s+1;
	    data[refx][refy][refz] = twobytes;
	  }
	}
	break;

      case 'C':
	/* coronal slices; read X,Z */
	for( z=1 ; z<=pixel_z ; z++ ) {
	  for( x=1 ; x<=pixel_x ; x++ ) {
	    if( fread(&twobytes,sizeof(short),1,datafp) != 1 ) {
	      fprintf(stderr, "Error reading %s.\n", datafname);
	      exit(1);
	    }
	    if(swapflag==-1) swap16(&twobytes);
	    refx = x; refy = s; refz = z;
	    if(orient[0]=='R') refx = pixel_x-x+1;
	    if(orient[1]=='A') refy = slices-s+1;
	    if(orient[2]=='S') refz = pixel_z-z+1;
	    data[refx][refy][refz] = twobytes;
	  }
	}
	break;

      case 'S':
	/* sagittal slices; read Y,Z */
	for( z=1 ; z<=pixel_z ; z++ ) {
	  for( y=1 ; y<=pixel_y ; y++ ) {
	    if( fread(&twobytes,sizeof(short),1,datafp) != 1 ) {
	      fprintf(stderr, "Error reading %s.\n", datafname);
	      exit(1);
	    }
	    if(swapflag==-1) swap16(&twobytes);
	    refx = s; refy = y; refz = z;
	    if(orient[0]=='R') refx = slices-s+1;
	    if(orient[1]=='A') refy = pixel_y-y+1;
	    if(orient[2]=='S') refz = pixel_z-z+1;
	    data[refx][refy][refz] = twobytes;
	  }
	}
	break;

      default:
	fprintf(stderr, "Incorrect axis specified in i2spm.in\n");
	exit(1);
      }

      fclose(datafp);

    }


    /* Now that we've read all slices at one time point into data,
       write data out to the appropriate .img file */
    
    if( timepts == 1) {
      sprintf(outfname, outmask);
      sprintf(hdrfprefix, outmask);
    }
    else {
      sprintf(outfname, outmask, t);
      sprintf(hdrfprefix, outmask, t);
    }
    strcat(outfname, ".img");
    strcpy(hdrfname, hdrfprefix);
    strcat(hdrfname, ".hdr");
    
    /* Check if the file exists -- this doesn't cut it if the
       file has write permission but not read permission */
    if( ((outfp = fopen(outfname, "rb")) != NULL) && (overwrite=='N') ) {
      fprintf(stderr, "File exists: %s\n", outfname);
      exit(1);
    }
    outfp = fopen(outfname, "wb");
    if( outfp == NULL ) {
      fprintf(stderr, "Error opening %s for write.\n", outfname);
      exit(1);
    }

    /* Write the data */
    for( z=1 ; z<=pixel_z ; z++)
      for( y=1 ; y<=pixel_y ; y++ )
	for( x=1 ; x<=pixel_x ; x++ )
	  fwrite(&data[x][y][z], sizeof(short), 1, outfp);
    fclose(outfp);

    /* Write out a header for this file, in imitation
       of the spm_hwrite matlab script. Don't check whether
       the file exists, on the theory that a .hdr without an
       .img is useless anyway. Types: int must be 4 bytes,
       short must be 2 bytes */
    outfp = fopen(hdrfname, "wb");
    if( outfp == NULL ) {
      fprintf(stderr, "Error opening %s for write.\n", hdrfname);
      exit(1);
    }

    /* Set the origin */
    for( i=0 ; i<=4 ; i++ )
      origin[i] = 0;
    
    /* Get the bare filename for the db_name field */
    len = strlen(hdrfprefix);
    for( pos=len ; pos>=0 ; pos-- )
      if( hdrfprefix[pos] == '/' )
	break;
    pos++;
    if( len-pos >= 16 ) {
      for( i=0 ; i<=16 ; i++)
	short_outfname[i] = hdrfprefix[pos+i];
      short_outfname[17] = 0;
    }
    else {
      for( i=0 ; i<=len-pos ; i++ )
	short_outfname[i] = hdrfprefix[pos+i];
      for( i=len-pos+1 ; i<=17 ; i++)
	short_outfname[i] = 0;
    }

    /* More stuff to put in the header */
    for( i=0 ; i<85 ; i++ )
      zeros[i] = 0;
    for( i=0 ; i<80 ; i++ )
      descrip[i] = 0;
    strcpy(auxfile, "none                   ");
    hdrfilesize = 348;
    strcpy(datatype, "dsr      ");
    strcpy(descrip, "spm compatible");

    /* header_key */
    fseek(outfp, 0, SEEK_SET);

    fwrite(&hdrfilesize, sizeof(int), 1, outfp);
    fwrite(datatype, sizeof(char), 10, outfp);
    fwrite(short_outfname, sizeof(char), 18, outfp);
    fwrite(zeros, sizeof(int), 1, outfp);
    fwrite(zeros, sizeof(short), 1, outfp);
    fwrite("r", sizeof(char), 1, outfp);
    fwrite("0", sizeof(char), 1, outfp);

    /* image_dimension */
    fseek(outfp, 40, SEEK_SET);

    fwrite(&four, sizeof(short), 1, outfp);        /* start dim */
    fwrite(&pixel_x, sizeof(short), 1, outfp);
    fwrite(&pixel_y, sizeof(short), 1, outfp);
    fwrite(&pixel_z, sizeof(short), 1, outfp);
    fwrite(&one, sizeof(short), 1, outfp);
    fwrite(zeros, sizeof(short), 1, outfp); /* writing shorts from char */
    fwrite(zeros, sizeof(short), 1, outfp); /*   array messy but OK */
    fwrite(zeros, sizeof(short), 1, outfp);        /* end dim */
    fwrite("mm", sizeof(char), 2, outfp);
    fwrite(zeros, sizeof(char), 1, outfp);
    fwrite(zeros, sizeof(char), 1, outfp);
    
    fwrite(zeros, sizeof(char), 8, outfp);
    fwrite(zeros, sizeof(short), 1, outfp);
    fwrite(&four, sizeof(short), 1, outfp);        /* data type: 4 = short  */
    fwrite(&sixteen, sizeof(short), 1, outfp);       /* bits per datum */
    fwrite(zeros, sizeof(short), 1, outfp);
    fwrite(&fzero, sizeof(float), 1, outfp);        /* start pixdim */
    fwrite(&voxel_x, sizeof(float), 1, outfp);
    fwrite(&voxel_y, sizeof(float), 1, outfp);
    fwrite(&voxel_z, sizeof(float), 1, outfp);
    fwrite(&fzero, sizeof(float), 1, outfp);
    fwrite(&fzero, sizeof(float), 1, outfp);
    fwrite(&fzero, sizeof(float), 1, outfp);
    fwrite(&fzero, sizeof(float), 1, outfp);        /* end pixdim */
    fwrite(&fzero, sizeof(float), 1, outfp);        /* OFFSET */
    fwrite(&scale, sizeof(float), 1, outfp);    /* scale - funused1 */
    fwrite(&fzero, sizeof(float), 1, outfp);
    fwrite(&fzero, sizeof(float), 1, outfp);
    fwrite(&fzero, sizeof(float), 1, outfp);
    fwrite(&fzero, sizeof(float), 1, outfp);
    fwrite(zeros, sizeof(int), 1, outfp);
    fwrite(zeros, sizeof(int), 1, outfp);
    fwrite(&glmax, sizeof(int), 1, outfp);     /* global max */
    fwrite(&glmin, sizeof(int), 1, outfp);     /* global min */

    /* image_dimension */
    fwrite(descrip, sizeof(char), 80, outfp);   /* description */
    fwrite(auxfile, sizeof(char), 24, outfp);
    fwrite(zeros, sizeof(char), 1, outfp);
    fwrite(origin, sizeof(short), 5, outfp);
    fwrite(zeros, sizeof(char), 85, outfp);

    if( ftell(outfp) != hdrfilesize ) {
      fprintf(stderr, "Error writing %s: wrong length.\n", hdrfname);
      exit(1);
    }

    fclose(outfp);

  }

  /* Clean up */
  fclose(logfp);
  free_s3tensor(data, 1,pixel_x, 1,pixel_y, 1,pixel_z);
  return(0);
}

short ***s3tensor(long nrl, long nrh, long ncl, long nch, long ndl, long ndh)
/* allocate a short 3tensor with range t[nrl..nrh][ncl..nch][ndl..ndh] */
{
	long i,j,nrow=nrh-nrl+1,ncol=nch-ncl+1,ndep=ndh-ndl+1;
	short ***t;

	/* allocate pointers to pointers to rows */
	t=(short ***) malloc((size_t)((nrow+NR_END)*sizeof(short**)));
	if (!t) nrerror("allocation failure 1 in s3tensor()");
	t += NR_END;
	t -= nrl;

	/* allocate pointers to rows and set pointers to them */
	t[nrl]=(short **) malloc((size_t)((nrow*ncol+NR_END)*sizeof(short*)));
	if (!t[nrl]) nrerror("allocation failure 2 in s3tensor()");
	t[nrl] += NR_END;
	t[nrl] -= ncl;

	/* allocate rows and set pointers to them */
	t[nrl][ncl]=(short *) malloc((size_t)((nrow*ncol*ndep+NR_END)*sizeof(short)));
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

void free_s3tensor(short ***t, long nrl, long nrh, long ncl, long nch, long ndl, long ndh)
/* free a short s3tensor allocated by s3tensor() */
{
	free((FREE_ARG) (t[nrl][ncl]+ndl-NR_END));
	free((FREE_ARG) (t[nrl]+ncl-NR_END));
	free((FREE_ARG) (t+nrl-NR_END));
}
