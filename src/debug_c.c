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

int FTNIZE(output_penciled_c) (char *filename, REAL *pencil, FINT *ndim,
			FINT *i, FINT *iy, FINT *iz, REAL *t,
			FINT *nx, FINT *ny, FINT *nz, FINT *nghost,
			FINT *fnlen)
/* Writes a scalar field to a file mimicking the Fortran record structure
   of the file. This subroutine is called once for each pencil.
   ndim           -- 1 for scalar, 3 for vector field
   i,ix,iy        -- position of pencil. i is the loop index.
   n[x-z],nghosts -- data layout
   fnlen          -- length of filename (needed due to obscure Fortran--C
                     mapping of strings)

   Using linked list data structure to allow for calling this routine with
   different variables in parallel. This assumes that the calls occur in
   the same order throughought the whole tencil loop -- the opposite is
   very implausible.
*/
{
  typedef struct node{
    int depth;
    struct node *next;
    struct node *prev;
    FILE *file;
  } filenode;
  static filenode root_node, *current_node;
  filenode *new_node;

  REAL zero=0.;
  static char *fname;
  int npencil,ilast,bcount;
  int j,k,l,m;
  long int mx,my,mz;
  long int datasize,pos;
  int first_pencil,last_pencil;
  static int first_var=1,max_depth,toggle;	/* Processing first variable */

  mx = *nx + 2*(*nghost);
  my = *ny + 2*(*nghost);
  mz = *nz + 2*(*nghost);

  npencil = (*ny)*(*nz);		/* Number of pencils (excluding ghosts) */
  ilast = npencil;

  first_pencil = (*i == 1);
  last_pencil = (*i == ilast);


    /* Called for the first time:
       - open file
       - calculate and write byte count
    */
  if (first_pencil) {
    if (first_var) {		/* First pencil, first variable: need to
				   initialise list structure */
      root_node.depth=0;
      root_node.prev=0;		/* just to be safe */
      current_node = &root_node;
      max_depth = 0;
      toggle = 1;
    } else {			/* Append new node */
      new_node = (filenode *) malloc(sizeof(filenode));
      new_node->prev = current_node;
      new_node->depth = current_node->depth + 1;
      max_depth++;
      current_node->next = new_node;
      current_node = new_node;	/* move on in the linked list */
    }
    /* Extract filename as C string */
    fname = (char *)calloc(*fnlen+1, sizeof(char));
    strncpy(fname,filename,*fnlen);
    fname[*fnlen] = 0;
    /* Open file */
    current_node->file = fopen(fname, "w");
    if (current_node->file == NULL) {
      fprintf(stderr, "debug_c.c: Can't open file %s\n", fname);
      abort();
    }
    free(fname);
    datasize = mx*my*mz*(*ndim)*sizeof(REAL);
    bcount = datasize;
    fwrite(&bcount, sizeof(bcount), 1, current_node->file);
  }


  /* Second pencil called for first variable:
     - reset first_var and toggle
     (needs to be done just once; for the following calls, depth
     vs. maxdepth is the correct criterion)
  */
  if (!first_pencil && toggle) {
    toggle = 0;
    first_var = 1;
  }


  /* Any call:
     - position appropriately
     - write nghosts zeros, one pencil, nghost zeros
  */
  if (!first_pencil) {
    if (first_var) {		/* reset current_node to start from top */
      current_node = &root_node; 
    } else {			/* move on in linked list */
      current_node = current_node->next;
    }
  }
  for (m=0;m<*ndim;m++) { 
    pos = sizeof(FINT) + mx*(*iy-1 + my*(*iz-1 + mz*m))*sizeof(REAL);
    fseek(current_node->file, pos, SEEK_SET);
    for (j=0;j<*nghost;j++) {
      fwrite(&zero, sizeof(REAL), 1, current_node->file);
    }
    fwrite(pencil+m*(*nx), sizeof(REAL), *nx, current_node->file);
    for (j=0;j<*nghost;j++) {
      fwrite(&zero, sizeof(REAL), 1, current_node->file);
    }
    if (first_pencil || (current_node->depth < max_depth)) {
      first_var = 0;
    } else {			/* processing last variable */
      first_var = 1;		/* reset for next call */
    }
  }
  

  /* Last call:
     - fill remaining ghost zones with zeros
     - position after data block
     - write byte count
     - write time as short record
     - close file
  */
  if (last_pencil) {
    /* Zero out remaining ghost zones */
    for (m=0;m<*ndim;m++) { 
      for (k=0;k<*nghost;k++) {
	for (l=0;l<mz;l++) {
	  pos = sizeof(FINT) + mx*(k + my*(l + mz*m))*sizeof(REAL);
	  fseek(current_node->file, pos, SEEK_SET);
	  for (j=0;j<mx;j++) {
	    fwrite(&zero, sizeof(REAL), 1, current_node->file);
	  }
	  pos = sizeof(FINT) + mx*((my-k-1) + my*(l + mz*m))*sizeof(REAL);
	  fseek(current_node->file, pos, SEEK_SET);
	  for (j=0;j<mx;j++) {
	    fwrite(&zero, sizeof(REAL), 1, current_node->file);
	  }
	}
      }
      for (k=(*nghost);k<my-(*nghost);k++) {
	for (l=0;l<*nghost;l++) {
	  pos = sizeof(FINT) + mx*(k + my*(l + mz*m))*sizeof(REAL);
	  fseek(current_node->file, pos, SEEK_SET);
	  for (j=0;j<mx;j++) {
	    fwrite(&zero, sizeof(REAL), 1, current_node->file);
	  }
	  pos = sizeof(FINT) + mx*(k + my*(mz-l-1 + mz*m))*sizeof(REAL);
	  fseek(current_node->file, pos, SEEK_SET);
	  for (j=0;j<mx;j++) {
	    fwrite(&zero, sizeof(REAL), 1, current_node->file);
	  }
	}
      }
    }
    /* Write byte count */
    datasize = mx*my*mz*(*ndim)*sizeof(REAL);
    pos = (long int)(datasize+sizeof(FINT));
    fseek(current_node->file, pos, SEEK_SET);
    bcount = datasize;
    fwrite(&bcount, sizeof(bcount), 1, current_node->file);
    /* Write time record */
    bcount = sizeof(REAL);
    fwrite(&bcount, sizeof(bcount), 1, current_node->file);
    fwrite(t, sizeof(REAL), 1, current_node->file);
    fwrite(&bcount, sizeof(bcount), 1, current_node->file);
    fclose(current_node->file);
    /* Set first_var to zero -- unless this was the very last cycle, in
       which case we set it to 1 for the next round */
    if (current_node->depth < max_depth) { /* Is this the last variable? */
      first_var = 0;
    } else {
      first_var = 1;
      /* Walk through list and free memory */
      while (current_node->depth >0) {
	current_node = current_node->prev;
	free(current_node->next);
      }
    }
  }

  return 1;
}

/* ---------------------------------------------------------------------- */


/* End of file debug_c.c */
