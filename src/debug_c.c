/*                             debug_c.c
                              -----------
*/

/* Date:   15-Feb-2002
   Author: wd (Wolfgang.Dobler@kis.uni-freiburg.de)
   Description:
 Define function used for debugging and diagnostics of the F90 code.
 Written in C because some things (in particular IO) are difficult to code
 in Fortran.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>


/* Beware: You need to choose single or double precision here */
#define REAL float
#define FINT int
#define NBYTES 4
/*  #define REAL double */
/*  #define NBYTES 8 */

/* Pick correct number of underscores here (2 for g77 without
   `-fno-second-underscore', 1 for most other compilers).
   Use the `-DFUNDERSC=1' option in the makefile to set this.
*/

#if (FUNDERSC == 2)
#  define FTNIZE(name) name##__
#elif (FUNDERSC == 0)
#  define FTNIZE(name) name
#else
#  define FTNIZE(name) name##_
#endif

/*
  do imn=1,ny*nz
    n=nn(imn)
    m=mm(imn)
    if (necessary(imn)) call finalise_isendrcv_bdry(f)
*/

/* ---------------------------------------------------------------------- */

int output_stenciled_c_(char *filename, REAL *stenc, FINT *ndim,
			FINT *i, FINT *iy, FINT *iz, REAL *t,
			FINT *nx, FINT *ny, FINT *nz, FINT *nghost,
			FINT *fnlen)
/* Writes a scalar field to a file mimicking the Fortran record structure
   of the file. This subroutine is called once for each stencil.
   fnlen   -- length of filename (needed due to unclear Fortran--C mapping of
              strings)
   i,ix,iy -- position of stencil. i is the loop index.
*/
{
  REAL zero=0;
  static FILE *file;
  static char *fname;
  int nstenc,ilast,bcount,j;
  long int mx,my,mz;
  long int datasize,pos;

  mx = *nx + 2*(*nghost);
  my = *ny + 2*(*nghost);
  mz = *nz + 2*(*nghost);

  nstenc = (*ny)*(*nz);		/* Number of stencils (sans ghosts) */
  ilast = nstenc;

  if (*i == 1) {
    /* Called for the first time:
       - open file
       - calculate and write byte count
    */
    /* Extract filename as C string */
    fname = (char *)calloc(*fnlen+1, sizeof(char));
    strncpy(fname,filename,*fnlen);
    fname[*fnlen] = 0;
    /* Open file */
    file = fopen(fname, "w");
    if (file == NULL) {
      fprintf(stderr, "debug_c.c: Can't open file %s\n", fname);
      abort;
    }
    datasize = mx*my*mz*sizeof(REAL);
    bcount = datasize;
    fwrite(&bcount, sizeof(bcount), 1, file);
  }

  /* Neither first nor last call:
     - position appropriately
     - write nghosts zeros, one stencil, nghost zeros
  */
  pos = sizeof(FINT) + mx*(*iy-1 + my*(*iz-1))*sizeof(REAL);
if (*i == 1){
  fprintf(stderr,"i=%d, iy=%d, iz=%d, pos=%d\n", *i, *iy,*iz,pos);
}
  fseek(file, pos, SEEK_SET);
  for (j=0;j<*nghost;j++) {
    fwrite(&zero, sizeof(REAL), 1, file);
  }
  fwrite(stenc, sizeof(REAL), *nx, file);
  for (j=0;j<*nghost;j++) {
    fwrite(&zero, sizeof(REAL), 1, file);
  }



  if (*i == ilast) {
    /* Last call:
       - position after data block
       - write byte count
       - write time as short record
       - close file
    */
    datasize = mx*my*mz*sizeof(REAL);
    pos = (long int)(datasize+sizeof(FINT));
    fseek(file, pos, SEEK_SET);
    bcount = datasize;
    fwrite(&bcount, sizeof(bcount), 1, file);
    /* Write time record */
    bcount = sizeof(REAL);
    fwrite(&bcount, sizeof(bcount), 1, file);
    fwrite(t, sizeof(REAL), 1, file);
    fwrite(&bcount, sizeof(bcount), 1, file);
    fclose(file);
    free(fname);
  }

  return 1;
}

/* ---------------------------------------------------------------------- */


/* End of file debug_c.c */
