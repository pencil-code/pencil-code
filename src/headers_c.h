/*                             headers_c.h
                              ------------
*/

/* $Id$
   Description:
     Common headers for all of our C files. Mostly required to get the
     number of underscores and single vs. double precision right.
*/

/* Choose single or double precision here (typically done from the Makefile) */
#ifdef DOUBLE_PRECISION
#  define REAL double
#  define FINT int
#  define NBYTES 8
#  define GSL_PREC GSL_PREC_DOUBLE
#else
#  define REAL float
#  define FINT int
#  define NBYTES 4
#  define GSL_PREC  GSL_PREC_SINGLE
#endif

/* Pick correct number of underscores here (2 for g95 without
   `-fno-second-underscore', 1 for most other compilers).
   Use the `-DFUNDERSC=1' option in the Makefile to set this.
*/
#if (FUNDERSC == 0)
#  define FTNIZE(name) name
#elif (FUNDERSC == 1)
#  define FTNIZE(name) name##_
#else
#  define FTNIZE(name) name##__
#endif

/* End of file */
