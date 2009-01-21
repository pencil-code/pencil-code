### $Id$
### Machine Specific Makefile for modular pencil code
### 
### The appropriate Makefile will be chosen from the config tree and included during at compile time 
### 

## Uncomment this line to inherit the Makefile in the parent directory
include ../Makefile


FC=mpif90
F77=$(FC)
FFLAGS=-lowercase -I/usr/mpi/include -O4 -pvec -parallel=0
FFLAGS=-lowercase -Os -noparallel -O4
CFLAGS=-O4 -DFUNDERSC=0
LDMPI=-L/usr/mpi/lib/lib32 -lfmpi -lmpi

