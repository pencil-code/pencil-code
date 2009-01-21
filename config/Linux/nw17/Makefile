### $Id$
### Machine Specific Makefile for modular pencil code
### 
### The appropriate Makefile will be chosen from the config tree and included during at compile time 
### 

## Uncomment this line to inherit the Makefile in the parent directory
include ../Makefile

FFLAGS= -O3 -w90 
FC=ifc -I/usr/lib/lam/include 
LDMPI=-L/usr/lib/lam/lib -llammpio -lpmpi -llamf77mpi -lmpi -llam -lnsl -lutil 
F77=$(FC) 
