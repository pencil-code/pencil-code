### $Id$
### Machine Specific Makefile for modular pencil code
### 
### The appropriate Makefile will be chosen from the config tree and included during at compile time 
### 

## Uncomment this line to inherit the Makefile in the parent directory
include ../Makefile

FFLAGS= -64 -O3 -mips4 -r12000 -LNO:cs1=32k:cs2=8M -I/share/apps/unsupported/pacx/standard/include -L/share/apps/unsupported/pacx/standard/lib 
LDMPI= -L/share/apps/unsupported/pacx/standard/lib -lpacxf -lpacx -lmpi 
