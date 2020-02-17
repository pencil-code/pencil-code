### $Id: Makefile,v 1.20 2008-07-02 16:48:47 brandenb Exp $
### Makefile for modular pencil code

## Default settings: they can be overwritten by src/Makefile.local
## which should be present in each local to run directory.

LDMPI=-lmpi

MPICOMM   = nompicomm
GLOBAL    = noglobal
IO        =   io_dist
FFT       = nofft
DEBUG     =   debug_c

HYDRO     =   hydro
DENSITY   =   density
ENTROPY   = noentropy
MAGNETIC  = nomagnetic
RADIATION = noradiation
PSCALAR   = nopscalar

GRAVITY   = nograv
FORCING   = noforcing
SHEAR     = noshear
TIMEAVG   = notimeavg

FFLAGS_GENERAL=
# 
# Default double flags (can be overwritten in machine-specific section).
#
FFLAGS_DOUBLE = -r8
CFLAGS_DOUBLE = -DDOUBLE_PRECISION

## Overwrite settings that are local to each run directory
include Makefile.local

MPICOMM   = mpicomm_
LDMPI=

include Makefile.inc  ## needed ? 

## Default settings
FC=f90
FFLAGS=-O1
CFLAGS=-O3 -DFUNDERSC=1

## The following block (`Begin/End machine dependent') specifies compilers
## and switches for given architectures and machine names. See
## `perldoc adapt-mkfile' for full documentation of the syntax.

### Begin machine dependent

## IRIX64:
#FC=f90
#FFLAGS= -64 -O3 -mips4
#FFLAGS= -auto_use mpi_interface -64 -O3 -C -macro_expand  #(gridur|embla)
#FFLAGS= -auto_use mpi_interface -64 -O3 -C -macro_expand  #(Antares)
#FFLAGS= -pfalist -64 -O3 -mips4 -r10000 -C -macro_expand  #(Grand)
#FFLAGS= -64 -O3 -mips4 -r12000 #(Ukaff)
#FFLAGS= -64 -O3 -auto_use mpi_interface -mips4 #(ukaff|origin)
#OPTFLAGS= -64
#F77=$(FC)
#CFLAGS=-O3 -64 -DFUNDERSC=1

## OSF1:
#FC=f90 #(Mhd.)
#OPTFLAGS=
#OPTFLAGS=-fast -O5
#FFLAGS=$(OPTFLAGS) -tune ev6 -arch ev6 #(Mhd.)
#F77=$(FC)
#LIBSIG=-lsig

## Linux:
## (AB/22-06-03: Here I've cut+pasted the current block from $PENCIL_HOME/src)
## The -tpp7 flag leads to a 10% speedup on nl6
#FC=f95
#FC=mpif90
#FFLAGS=-O3
#F90FLAGS=-u -tpp7
#F90FLAGS= # (fire)
#F90FLAGS=
#LDMPI=
##LDMPI=-lfmpich -lmpich -lmpe -lpmpich -llmpe # (copson) 
##LDMPI=-L/usr/X11R6/lib/ -lfmpich -lampe -lmpe -lpmpich -lX11 # (copson)
#FC=mpif90 # (Samos|fire|copson)
#FC=/usr/lib/lam/bin/mpif95 #(Cincinnatus|Owen|Master)
##FC=/home/brandenb/pacx/bin/pacxfc #PACX #(fe[0-3])
##FC=mpif77 #(fe[0-3])
#FC=mpif77 #(ns0)
#FC=ifc_6.0 -I/usr/local/lib/LAM/include #(fe[0-9])
#FC=/people/disk2/brandenb/bin/ifort -I/usr/local/lib/LAM/include #(fe[0-9])
#FC=~dobler/opt/intel/compiler70/ia32/bin/ifc -I/usr/local/lib/LAM/include #(fe[0-9])
#FC=f95 -I/usr/lib/lam/include -mismatch #(morvern)      
#FC=pgf90 #(sintef*|vsl176)
#FFLAGS_DOUBLE= -r8 #(vsl176)
#FFLAGS= -O3 -I/cvos/shared/apps/mvapich/pgi/1.2rc2/include/ #(vsl176)
#FFLAGS= -O4 -Munroll -Mvect=prefetch -I/cvos/shared/apps/mvapich/pgi/1.2rc2/include/ #(vsl176)
##FFLAGS= -O0 -C -Mbounds -Mdwarf2 -Melf  -I/cvos/shared/apps/mvapich/pgi/1.2rc2/include/ #(vsl176)
##FFLAGS= -O3 -I/cvos/shared/apps/mvapich/pgi/1.2rc2/include/ -Mrecursive #(vsl176)
#LD_MPI= -L/cvos/shared/apps/mvapich/pgi/1.2rc2/lib -lmpich -L/cvos/shared/apps/ofed/1.3.1/lib64/ -libumad -libverbs -lrdmacm #(vsl176) 
#FC=ifort #(vsl21*|vsl22*|vsl230)
#FFLAGS= -I./ -D_REENTRANT -I/opt/scali/include #(vsl21*|vsl22*|vsl230)
#LD_MPI= -L/opt/scali/lib -lfmpi -lmpi -lpthread #(vsl21*|vsl22*|vsl230)
#F77=$(FC) #(sintef*|vmware|erland|bohm|igor|onsager|vsl21*|vsl22*|vsl230)
#FC=ifc -I/usr/lib/lam/include #(morvern)  
##FC=/home/brandenb/pacx/bin/pacxfc #(ns0)
#FC=ifc #(nl6|mmtf|nl19|hamlet)
#FC=f95 #(karbala)
#FC=pgf90 #(erland|igor|tycho|vmware|linux29|onsager|fe4)
#FFLAGS=-O # (Samos)
#FFLAGS= # (fire)
#LDMPI=-static -L/usr/local/lib/LAM/lib -llamf77mpi -lmpi -llam -lnsl -lutil -L/usr/local/lib/PGI/linux86/lib -lpgc #(fe[0-4])
##AB: next line uncommented on Tobi's advice
#LDMPI=-static -L/usr/local/lib/LAM/lib -llamf77mpi -lmpi -llam -lnsl -lutil -pthread #(fe[0-4])
#LDMPI=-L/usr/lib/lam/lib -llam -lmpi -lpmpi -llamf77mpi -llammpi++ -llammpio #(morvern)
#LDMPI=-static -L/opt/lam/lib -lmpi -llamf77mpi -llam  #(nw50)
##FFLAGS=-O3 -ipo #PACX #(fe[0-3])
#FFLAGS=-fast -lm #(erland|igor|tycho|vmware|linux29|onsager|fe4)
#F90FLAGS= #(erland|igor|tycho|vmware|linux29|onsager|fe4)
##FFLAGS=-O3 -ipo #(fe[0-3])
#FFLAGS= #(fe[0-3])  ##(compiles faster/axel)
#FFLAGS=-O4 -C -gline -Wc,-malign-double #(Mephisto)
#FFLAGS=-O4 #(nl3)
#FFLAGS= -O3 -w90 #(cosmo)
#FC=ifc -I/usr/include #(cosmo)  
#FC=ifc -I/usr/include #(nl6)  
#LDMPI=-L/usr/lib -lmpi -llamf77mpi -llam -lpmpi -llammpi++ -llammpio #(cosmo)
#LDMPI= # (Samos|fire)
#CFLAGS=-O3 -DFUNDERSC=1 -Wall
#F77=$(FC)
#F77=$(FC) -dusty #(nl3|Mephisto|mmtf|karbala|morvern)
#F77=$(FC) -lm  #(tycho|vmware|erland|igor|onsager|fe4)
#F77=$(FC) #(Cincinnatus|Owen|Master|cosmo|nl6)
##
## DCSC cluster - (Copenhagen, not Odense!)
## www.dcsc.ku.dk
#FC=pgf90 #(fend0[1-3])
#CC=pgcc #(fend0[1-3])
#FFLAGS=-fastsse -tp k8-64e -I/usr/local/topspin/mpi/mpich/include #(fend0[1-3])
#CFLAGS=-fastsse -tp k8-64e -DFUNDERSC=1 #(fend0[1-3])
#LD_MPI=-L/usr/local/topspin/mpi/mpich/lib64 -lmpich -lpthread -lmpich_p #(fend0[1-3])
#FC=gfortran -I. #(copernicus)
#FFLAGS= -O3 -I/usr/local/include #(copernicus)
#FFLAGS_DOUBLE=-fdefault-real-8 #(copernicus)
##LD_MPI=-L/usr/local/lib -lfmpich -lmpich -pthread #(copernicus)
##
## Axel's laptop (nl6, 2.0 GHz)
##
#FC=ifort #(nl6)
#FFLAGS=-O3 -I/usr/include/lam -CB #(nl6)
#FFLAGS=-O3 -I/usr/include/lam #(nl6)
#LD_MPI=-I/usr/include/lam -L/usr/lib/lam/lib -llammpio -llamf77mpi -lmpi -llam -lutil #(nl6)
#LIBMPISUFFIX=_#(nl6)
#CFLAGS=-O3 -DFUNDERSC=1 -Wall #(nl6)
#F77FLAGS=-w90 #(nl6)
##
## Axel on Neolith (ifc 7.1)
##
#FC=ifort #(neolith)
#FFLAGS=-O3 -Nmpi#(neolith)
#CC=icc #(neolith)
#CFLAGS=-O3 -DFUNDERSC=1 #(neolith)
#LD_MPI=-Nscampi #(neolith)
#LD_MPI=-L/opt/scali/lib -lfmpi -lmpi -lpthread #(neolith)
##LIBMPISUFFIX=_#(neolith)
##
## Tegner (Stockholm)                
## 1516 24-nodes per cabinet 2.1GHz Opterons                      
##                                                                        
#FC=ftn #(tegner*)
#F90=$(FC) #(tegner*)
#CC=cc #(tegner*)
#FFLAGS=-O3 #(tegner*) 
#CFLAGS=-DFUNDERSC=1 -O3 #(tegner*)
##
## Lindgren - CrayXE6 - (Stockholm)                
## 1516 24-nodes per cabinet 2.1GHz Opterons                      
##                                                                        
#FC=ftn #(emil*)
#F90=$(FC) #(emil*)
#CC=cc #(emil*)
#FFLAGS=-O3 #(emil*) 
#CFLAGS=-DFUNDERSC=1 -O3 #(emil*)
## SunOS:
#FC=pghpf #(Lomond)
#FFLAGS= -Msmp -O4 #(Lomond)
#FC=f95 #(Lomond)
#FFLAGS=-O3 -I/opt/local/mpich/include #(Lomond)
#LDMPI=-L/opt/local/mpich/lib -lmpich -lsocket -lnsl -laio #(Lomond)
#F77=$(FC) #(Lomond)
##
## Sisu - CrayXC30 - (Kajaani)
##
## do first "module load intel" 
##
#FC=ifort #(sisu*)
#FFLAGS=-O1 #(sisu*)
#CC=icc #(sisu*)
#CFLAGS=  #(sisu*)
#LD_MPI=  #(sisu*)
##
### old
##FC=ftn #(sisu*)
##F90=$(FC) #(sisu*)
##CC=cc #(sisu*)
##FFLAGS=-O3 -ffree -e m -J experimental -J magnetic #(sisu*)
##FFLAGS=-O3 -ffree -e m  #(sisu*)
##FFLAGS_DOUBLE=-s default64 #(sisu*)
##FFLAGS_DOUBLE=-r8 #(sisu*)
##CFLAGS=-DFUNDERSC=1 -O3 #(sisu*)
##
## Hermit - CrayXE6 - (Stuttgart)
##
## do first "module load intel"
#FC=ifort #(eslogin0*)
#FFLAGS=-O3 #(eslogin0*)
#CC=icc #(eslogin0*)
#CFLAGS=  #(eslogin0*)
#LD_MPI=  #(eslogin0*)
##                                      
## GWDG Cluster in Goettingen
##
#FC = mpif90 -std95 # (gwdu101)
#F77 = mpif90 # (gwdu101)
#FFLAGS = -implicitnone # (gwdu101)
#FFLAGS_DOUBLE = -r8 # (gwdu101)
#CC = mpicc # (gwdu101)
#CFLAGS = -DFUNDERSC=1 # (gwdu101)
#CFLAGS_DOUBLE = -DDOUBLE_PRECISION # (gwdu101)
#LDMPI=  # (gwdu101)
##
#FC = mpif90 -std95 # (gwdu102)
#F77 = mpif90 # (gwdu102)
#FFLAGS = -implicitnone # (gwdu102)
#FFLAGS_DOUBLE = -r8 # (gwdu102)
#CC = mpicc # (gwdu102)
#CFLAGS = -DFUNDERSC=1 # (gwdu102)
#CFLAGS_DOUBLE = -DDOUBLE_PRECISION # (gwdu102)
#LDMPI=  # (gwdu102)
##

## UNICOS/mk:
#FC=f90
#FFLAGS=-O3 -em

## HI-UX/MPP:
#FC=f90
#FFLAGS=-I/usr/mpi/include -O4 -pvec -parallel=0
#LDMPI=-L/usr/mpi/lib/lib32 -lfmpi -lmpi

## AIX:
#FC=xlf90
#FC=mpxlf90 #(Sp[12])
#FFLAGS=-O5 -qsuffix=f=f90

### End machine dependent


ifeq ($(REAL_PRECISION),double)
  FFLAGS := $(FFLAGS) $(FFLAGS_DOUBLE) 
endif

LD=$(FC)
LDFLAGS=$(FFLAGS) $(FFLAGS_GENERAL)

.SUFFIXES: .f90 .mod

.f.o:
	$(F77) $(OPTFLAGS) $(FFLAGS_GENERAL) -c $*.f
.f90.o:
	$(FC) $(FFLAGS) $(FFLAGS_GENERAL) -c $*.f90
.mod.o:
	$(FC) $(FFLAGS) $(FFLAGS_GENERAL) -c $*.f90
.c.o:
	$(CC) $(CFLAGS) -c $*.c

#####

default: makefile cparam.inc code

makefile: Makefile
	$(PENCIL_HOME)/src/scripts/adapt-mkfile Makefile makefile
code:
	make remesh.x
#	make remesh.x remesh_ext.x

technical=cparam.o cdata.o common.o syscalls.o syscalls_ansi.o general.o
remesh=remesh.o $(technical)
#remesh_ext=remesh_ext.o $(technical)

remesh.x: $(remesh) Makefile.local
	$(LD) $(LDFLAGS) $(remesh) $(LDMPI) -o remesh.x 

#remesh_ext.x: $(remesh_ext) Makefile.local
#	$(LD) $(LDFLAGS) $(remesh_ext) $(LDMPI) -o remesh_ext.x 

remesh.o:    remesh.f90 common.local $(technical)
#remesh_ext.o:    remesh_ext.f90 common.local $(technical)
#
cparam.o:    cparam.f90 cparam.local cparam.inc
cdata.o:     cdata.f90 cparam.o
common.o:    common.f90 common.local cdata.o cparam.o
syscalls.o:  syscalls.f90 syscalls_ansi.o
syscalls_ansi.o:  syscalls_ansi.c
general.o:   general.f90 cdata.o cparam.o syscalls.o
#deriv.o:     deriv.f90 cdata.o

# Clean up, but leave makefile and TAGS
clean:
	rm -f *.x *.o *.mod *.d *.inc
	rm -f work.pc work.pcl

# Clean up all
cleann: clean
	rm -f makefile TAGS
