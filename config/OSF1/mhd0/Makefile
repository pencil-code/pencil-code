### $Id: Makefile,v 1.1 2003-08-25 13:05:41 mee Exp $
### Machine Specific Makefile for modular pencil code
### 
### The appropriate Makefile will be chosen from the config tree and included during at compile time 
### 

## Uncomment this line to inherit the Makefile in the parent directory
include ../Makefile

FC=f95 
FC=f90 
OPTFLAGS=-fast -O5 -notransform_loops 
FFLAGS=$(OPTFLAGS) -tune ev6 -arch ev6 
