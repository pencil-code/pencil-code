### $Id$
### Machine Specific Makefile for modular pencil code
### 
### The appropriate Makefile will be chosen from the config tree and included during at compile time 
### 

## Uncomment this line to inherit the Makefile in the parent directory
include ../Makefile

FC=f95 -I/usr/lib/lam/include -mismatch 
FC=ifc -I/usr/lib/lam/include 
LDMPI=-L/usr/lib/lam/lib -llam -lmpi -lpmpi -llamf77mpi -llammpi++ -llammpio 
F77=$(FC) -dusty 
