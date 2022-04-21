/* $Id: vis5d.h,v 1.8 1997/01/02 17:25:29 billh Exp $ */
/* $Id: v5d.h,v 1.16 1996/08/23 13:03:12 billh Exp $ */

/* Vis5D version 5.0 */
/* Vis5D version 4.3 */
/* Vis5D version 4.2 */

#ifndef V5D_H
#define V5D_H

/*
 * A numeric version number which we can test for in utility programs which
 * use the v5d functions.  For example, we can do tests like this:
 * #if V5D_VERSION > 42
 *      do something
 * #else
 *      do something else
 * #endif
 *
 * If V5D_VERSION is not defined, then its value is considered to be zero.
 */

#define V5D_VERSION 42

/*
 * Define our own 1 and 2-byte data types.  We use these names to avoid
 * collisions with types defined by the OS include files.
 */
typedef unsigned char V5Dubyte;     /* Must be 1 byte, except for cray */
typedef unsigned short V5Dushort;   /* Must be 2 byte, except for cray */

#define MISSING 1.0e35
#define IS_MISSING(X)  ( (X) >= 1.0e30 )

/* Limits on 5-D grid size:  (must match those in v5df.h!!!) */
#define MAXVARS      30
#define MAXTIMES    200
#ifdef MCIDAS_SIDECAR
#define MAXROWS     500
#define MAXCOLUMNS  500
#else
#define MAXROWS     500
#define MAXCOLUMNS  500
#endif
#define MAXLEVELS   200

/************************************************************************/
/***                                                                  ***/
/*** Functions for writing v5d files.  See README file for details.   ***/
/*** These are the functions user's will want for writing file        ***/
/*** converters, etc.                                                 ***/
/***                                                                  ***/
/************************************************************************/

extern int v5dCreate( const char *name,
                      int numtimes, int numvars,
                      int nr, int nc, const int nl[],
                      const char varname[MAXVARS][10],
                      const int timestamp[],
                      const int datestamp[],
                      int compressmode,
                      int projection,
                      const float proj_args[],
                      int vertical,
                      const float vert_args[] );

extern int v5dWrite( int time, int var, const float data[] );

extern int v5dClose( void );

extern int v5dSetLowLev( int lowlev[] );

extern int v5dSetUnits( int var, const char *units );

/************************************************************************/
/***                                                                  ***/
/*** Definition of v5d struct and function prototypes.                ***/
/*** These functions are used by vis5d and advanced v5d utilities.    ***/
/***                                                                  ***/
/************************************************************************/

#define MAXPROJARGS 100
#define MAXVERTARGS (MAXLEVELS+1)

/*
 * This struct describes the structure of a .v5d file.
 */
typedef struct {
    /* PUBLIC (user can freely read, sometimes write, these fields) */
        int NumTimes;                   /* Number of time steps */
        int NumVars;                    /* Number of variables */
        int Nr;                         /* Number of rows */
        int Nc;                         /* Number of columns */
        int Nl[MAXVARS];                /* Number of levels per variable */
        int LowLev[MAXVARS];            /* Lowest level per variable */
        char VarName[MAXVARS][10];      /* 9-character variable names */
        char Units[MAXVARS][20];        /* 19-character units for variables */
        int TimeStamp[MAXTIMES];        /* Time in HHMMSS format */
        int DateStamp[MAXTIMES];        /* Date in YYDDD format */
        float MinVal[MAXVARS];          /* Minimum variable data values */
        float MaxVal[MAXVARS];          /* Maximum variable data values */

        /* This info is used for external function computation */
        short McFile[MAXTIMES][MAXVARS];/* McIDAS file number in 1..9999 */
        short McGrid[MAXTIMES][MAXVARS];/* McIDAS grid number in 1..? */

        int VerticalSystem;             /* Which vertical coordinate system */
        float VertArgs[MAXVERTARGS];    /* Vert. Coord. Sys. arguments... */

        /*
        IF VerticalSystem==0 THEN
                -- Linear scale, equally-spaced levels in generic units
                VertArgs[0] = Height of bottom-most grid level in generic units
                VertArgs[1] = Increment between levels in generic units
        ELSE IF VerticalSystem==1 THEN
                -- Linear scale, equally-spaced levels in km
                VertArgs[0] = Height of bottom grid level in km
                VertArgs[1] = Increment between levels in km
        ELSE IF VerticalSystem==2 THEN
                -- Linear scale, Unequally spaced levels in km
                VertArgs[0] = Height of grid level 0 (bottom) in km
                ...                ...
                VertArgs[n] = Height of grid level n in km
        ELSE IF VerticalSystem==3 THEN
                -- Linear scale, Unequally spaced levels in mb
                VertArgs[0] = Pressure of grid level 0 (bottom) in mb
                ...             ...
                VertArgs[n] = Pressure of grid level n in mb
        ENDIF
        */

        int Projection;                     /* Which map projection */
        float ProjArgs[MAXPROJARGS];        /* Map projection arguments... */

        /*
        IF Projection==0 THEN
                -- Rectilinear grid, generic units
                ProjArgs[0] = North bound, Y coordinate of grid row 0
                ProjArgs[1] = West bound, X coordiante of grid column 0
                ProjArgs[2] = Increment between rows
                ProjArgs[3] = Increment between colums
                NOTES: X coordinates increase to the right, Y increase upward.
                NOTES: Coordinate system is right-handed.
        ELSE IF Projection==1 THEN
                -- Cylindrical equidistant (Old VIS-5D)
                -- Rectilinear grid in lat/lon
                ProjArgs[0] = Latitude of grid row 0, north bound, in degrees
                ProjArgs[1] = Longitude of grid column 0, west bound, in deg.
                ProjArgs[2] = Increment between rows in degrees
                ProjArgs[3] = Increment between rows in degrees
                NOTES: Coordinates (degrees) increase to the left and upward.
        ELSE IF Projection==2 THEN
                -- Lambert conformal
                ProjArgs[0] = Standared Latitude 1 of conic projection
                ProjArgs[1] = Standared Latitude 2 of conic projection
                ProjArgs[2] = Row of North/South pole
                ProjArgs[3] = Column of North/South pole
                ProjArgs[4] = Longitude which is parallel to columns
                ProjArgs[5] = Increment between grid columns in km
        ELSE IF Projection==3 THEN
                -- Polar Stereographic
                ProjArgs[0] = Latitude of center of projection
                ProjArgs[1] = Longitude of center of projection
                ProjArgs[2] = Grid row of center of projection
                ProjArgs[3] = Grid column of center of projection
                ProjArgs[4] = Increment between grid columns at center in km
        ELSE IF Projection==4 THEN
                -- Rotated
                ProjArgs[0] = Latitude on rotated globe of grid row 0
                ProjArgs[1] = Longitude on rotated globe of grid column 0
                ProjArgs[2] = Degrees of latitude on rotated globe between
                                grid rows
                ProjArgs[3] = Degrees of longitude on rotated globe between
                                grid columns
                ProjArgs[4] = Earth latitude of (0, 0) on rotated globe
                ProjArgs[5] = Earth longitude of (0, 0) on rotated globe
                ProjArgs[6] = Clockwise rotation of rotated globe in degrees
        ENDIF
        */

        int CompressMode;        /* 1, 2 or 4 = # bytes per grid point */
        char FileVersion[10];    /* 9-character version number */

    /* PRIVATE (not to be touched by user code) */
        unsigned int FileFormat; /* COMP5D file version or 0 if .v5d */
        int FileDesc;            /* Unix file descriptor */
        char Mode;               /* 'r' = read, 'w' = write */
        int CurPos;              /* current position of file pointer */
        int FirstGridPos;        /* position of first grid in file */
        int GridSize[MAXVARS];   /* size of each grid */
        int SumGridSizes;        /* sum of GridSize[0..NumVars-1] */
} v5dstruct;

extern float pressure_to_height( float pressure);

extern float height_to_pressure( float height );

extern int v5dYYDDDtoDays( int yyddd );

extern int v5dHHMMSStoSeconds( int hhmmss );

extern int v5dDaysToYYDDD( int days );

extern int v5dSecondsToHHMMSS( int seconds );

extern void v5dPrintStruct( const v5dstruct *v );

extern v5dstruct *v5dNewStruct( void );

extern void v5dFreeStruct( v5dstruct* v );

extern void v5dInitStruct( v5dstruct *v );

extern int v5dVerifyStruct( const v5dstruct *v );

extern void v5dCompressGrid( int nr, int nc, int nl, int compressmode,
                             const float data[], void *compdata,
                             float ga[], float gb[],
                             float *minval, float *maxval );

extern int v5dSizeofGrid( const v5dstruct *v, int time, int var );

extern v5dstruct *v5dOpenFile( const char *filename, v5dstruct *v );

extern int v5dCreateFile( const char *filename, v5dstruct *v );

extern v5dstruct *v5dUpdateFile( const char *filename, v5dstruct *v );

extern int v5dCloseFile( v5dstruct *v );

extern int v5dWriteCompressedGrid( const v5dstruct *v,
                                   int time, int var,
                                   const float *ga, const float *gb,
                                   const void *compdata );

extern int v5dWriteGrid( v5dstruct *v, int time, int var, const float data[] );

#endif

/*
 * Functions to do binary I/O of floats, ints, etc. with byte swapping
 * as needed.
 */

#ifndef BINIO_H
#define BINIO_H

/* Include files which define SEEK_SET, O_RD_ONLY, etc. */
/* and prototype open(), close(), lseek(), etc. */
#include <unistd.h>
#include <fcntl.h>

extern void flip4( const unsigned int *src, unsigned int *dest, int n );

extern void flip2( const unsigned short *src, unsigned short *dest, int n );

#ifdef _CRAY
  extern void cray_to_ieee_array( long *dest, const float *source, int n );
  extern void ieee_to_cray_array( float *dest, const long *source, int n );
#endif

/**********************************************************************/
/*****                         Write Functions                    *****/
/**********************************************************************/

extern int write_bytes( int f, const void *b, int n );

extern int write_int2_array( int f, const short *iarray, int n );

extern int write_uint2_array( int f, const unsigned short *iarray, int n );

extern int write_int4( int f, int i );

extern int write_int4_array( int f, const int *iarray, int n );

extern int write_float4( int f, float x );

extern int write_float4_array( int f, const float *x, int n );

extern int write_block( int f, const void *data, int elements, int elsize );

#endif

/* Vis5D version 4.3
 * This configuration file contains options which can be safely
 * changed by the user.
 */

#ifndef VIS5D_H
#define VIS5D_H

/*
 * Amount of physical RAM in megabytes:
 * vis5d normally uses a bounded amount of memory to avoid swapping.
 * When the limit is reached, the least-recently-viewed graphics will
 * be deallocated.  If MBS is set to 0, however, vis5d will use ordinary
 * malloc/free and not deallocate graphics (ok for systems with a lot
 * of memory (>=128MB)).
 */
#define MBS 32

/* Default topography file: */
#define TOPOFILE "EARTH.TOPO"

/* Default map lines files: */
#define WORLDFILE "OUTLSUPW"
#define USAFILE "OUTLUSAM"

/* Default filename of Tcl startup commands: */
#define TCL_STARTUP_FILE "vis5d.tcl"

/* Default directory to search for user functions: */
#define FUNCTION_PATH "userfuncs"

/* Default animation rate in milliseconds: */
#define ANIMRATE 100

/* Default scale and exponent values for logrithmic vertical coordinate system: */
#define DEFAULT_LOG_SCALE  1012.5
#define DEFAULT_LOG_EXP  -7.2

/***  USERS:  DON'T CHANGE ANYTHING BEYOND THIS POINT  ***/
/*
 * Define BIG_GFX to allow larger isosurfaces, contour slices, etc. if
 * there's enough memory.
#if MBS==0 || MBS>=128
#  define BIG_GFX
#endif
 */

#define BIG_GFX

/*
 * Shared by code above and below API:
 */
#define MAX_LABEL   1000
#define MAX_FUNCS   100

#endif
