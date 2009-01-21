### $Id$
### Machine Specific Makefile for modular pencil code
### 
### The appropriate Makefile will be chosen from the config tree and included during at compile time 
### 

## Uncomment this line to inherit the Makefile in the parent directory
include ../Makefile

FC=ifc_6.0 -I/usr/local/lib/LAM/include 
FC=pgf90 
LDMPI=-static -L/usr/local/lib/LAM/lib -llamf77mpi -lmpi -llam -lnsl -lutil -L/usr/local/lib/PGI/linux86/lib -lpgc 
FFLAGS=-fast -lm 
F90FLAGS= 
F77=$(FC) -lm  
