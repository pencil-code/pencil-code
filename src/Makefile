### Makefile for modular pencil code

## Default settings: they can be overwritten by src/Makefile.local
## which should be present in each local to run directory.

LIBSIG=-L$(HOME)/f90/lib  ##(Wolfgang, should we take this out?)
LDMPI=-lmpi
MPICOMM=nompicomm
GLOBAL=noglobal
GRAVITY=nograv
FORCING=noforcing
ENTROPY=noentropy
MAGNETIC=nomagnetic
HYDRO=hydro
DENSITY=density
IO=io_dist
FFT=nofft
DEBUG=debug_c

## Overwrite settings that are local to each run directory
include Makefile.local

## Default settings
FC=f90
FFLAGS=-O3
CFLAGS=-O3 -DFUNDERSC=1


### Begin machine dependent

## IRIX64:
#FC=f90
#FFLAGS= -64 -O3 -mips4
#FFLAGS= -auto_use mpi_interface -64 -O3 -C -macro_expand  #(gridur|embla)
#FFLAGS= -auto_use mpi_interface -64 -O3 -C -macro_expand  #(Antares)
#FFLAGS= -pfalist -64 -O3 -mips4 -r10000 -C -macro_expand  #(Grand)
#FFLAGS= -64 -O3 -mips4 -r12000 #(Ukaff)
#FFLAGS= -64 -O3 -mips4 #(ukaff|origin)
#F77=$(FC)
CFLAGS=-O3 -64 -DFUNDERSC=1

## OSF1:
#FC=f90 #(Mhd.)
#OPTFLAGS=
#OPTFLAGS=-fast -O5
#FFLAGS=$(OPTFLAGS) -tune ev6 -arch ev6 #(Mhd.)
#F77=$(FC)
#LIBSIG=-lsig

## Linux:
#FC=f95
#FFLAGS=-O3 -u
#LDMPI=
#FC=mpif90 # (Samos)
#FC=/usr/lib/lam/bin/mpif95 #(Cincinnatus|Owen|Master)
#FC=ifc	#(Nq[0-9]+)
#FC=/home/brandenb/bin/f95 #(Nq[0-9]+)
#FC=pgf90 #(Nq[0-9]+)
#FC=pgf90 #(erland|igor|tycho)
#FFLAGS=-O # (Samos)
##FFLAGS=-O3 -I/usr/lib/lam/include -static #(Nq[0-9]+)
#FFLAGS=-O3 -static #(Nq[0-9]+)
#FFLAGS=-fast -lm #(erland|igor|tycho)
#FFLAGS=-O4 -u -Wc,-malign-double	#(Nl3|Mephisto)
#FFLAGS= -fast -DPGF90 #(Nq[0-9]+)
#LDMPI=-llammpi++ -llammpio -lpmpi -llamf77mpi -lmpi -llam -lnsl -lutil #(Nq[0-9]+)
#LDMPI= # (Samos)
#CFLAGS=-O3 -DFUNDERSC=1 -Wall
#F77=$(FC) -O4 -dusty
#F77=$(FC) -O3 #(Cincinnatus|Owen|Master)

## SunOS:
#FC=pghpf #(Lomond)
#FFLAGS= -Msmp -O4  #(Lomond)
#FC=f95
#FFLAGS=-O3 -I/opt/local/mpich/include
#LDMPI=-L/opt/local/mpich/lib -lmpich -lsocket -lnsl -laio
#F77=$(FC)

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


LD=$(FC)
LDFLAGS=$(FFLAGS)

.SUFFIXES: .f90 .mod

.f.o:
	$(F77) $(OPTFLAGS) -c $*.f
.f90.o:
	$(FC) $(FFLAGS) -c $*.f90
.mod.o:
	$(FC) $(FFLAGS) -c $*.f90
.c.o:
	$(CC) $(CFLAGS) -c $*.c

#####

default: makefile code


makefile: Makefile
	adapt-mkfile Makefile makefile
code:
	make start.x run.x

technical=cparam.o cdata.o general.o slices.o deriv.o boundcond.o register.o prints.o
physics=$(GLOBAL).o $(HYDRO).o $(DENSITY).o $(FORCING).o $(GRAVITY).o $(ENTROPY).o $(MAGNETIC).o
comm_deps=$(technical) $(IO).o sub.o param_io.o $(MPICOMM).o $(physics)
start_deps_f90=$(comm_deps)
start_deps=$(start_deps_f90) $(DEBUG).o $(FFT).o
run_deps_f90=$(comm_deps) equ.o timestep.o wsnaps.o
run_deps=$(run_deps_f90) $(DEBUG).o $(FFT).o
start=start.o $(start_deps)
run=run.o $(run_deps)

start.x: $(start) Makefile.local
	$(LD) $(LDFLAGS) $(start) $(LDMPI) -o start.x 

run.x: $(run) Makefile.local
	$(LD) $(LDFLAGS) $(run) $(LDMPI) $(LIBSIG) -o run.x 

start.o:     start.f90 $(start_deps)
run.o:       run.f90 $(run_deps)
#
cparam.o:    cparam.f90 cparam.local
general.o:   general.f90
cdata.o:     cdata.f90 cparam.o
slices.o:    slices.f90 cdata.o
deriv.o:     deriv.f90 cdata.o
sub.o:       sub.f90 cdata.o deriv.o $(MPICOMM).o general.o $(IO).o
param_io.o:  param_io.f90 timestep.o $(physics)
register.o:  register.f90 cdata.o $(MPICOMM).o sub.o $(physics)
prints.o:    prints.f90 cdata.o sub.o $(physics)
equ.o:       equ.f90 cdata.o $(MPICOMM).o sub.o $(physics) slices.o
boundcond.o: boundcond.f90 cdata.o $(ENTROPY).o
timestep.o:  timestep.f90 cdata.o $(MPICOMM).o equ.o $(FORCING.o)
wsnaps.o:    wsnaps.f90 cdata.o $(MPICOMM).o boundcond.o
#
$(MPICOMM).o:   $(MPICOMM).f90 cdata.o general.o
$(GLOBAL).o:    $(GLOBAL).f90 cdata.o $(IO).o $(MPICOMM).o
$(GRAVITY).o:   $(GRAVITY).f90 cdata.o $(MPICOMM).o sub.o slices.o $(GLOBAL).o
$(FORCING).o:   $(FORCING).f90 cdata.o $(MPICOMM).o sub.o
$(ENTROPY).o:   $(ENTROPY).f90 cdata.o $(MPICOMM).o sub.o $(GLOBAL).o $(HYDRO).o $(DENSITY).o $(IO).o slices.o
$(MAGNETIC).o:  $(MAGNETIC).f90 cdata.o $(MPICOMM).o sub.o slices.o $(FFT).o
$(HYDRO).o:     $(HYDRO).f90 cdata.o $(MPICOMM).o sub.o $(GRAVITY).o $(DENSITY).o $(MAGNETIC).o slices.o
$(DENSITY).o:   $(DENSITY).f90 cdata.o $(MPICOMM).o sub.o $(GRAVITY).o slices.o
$(IO).o:        $(IO).f90 cdata.o $(MPICOMM).o $(DEBUG).o
# $(FFT).o:       $(FFT).f
nodebug.o:   nodebug.f90 cdata.o $(MPICOMM).o


TAG_SOURCE=$(run_deps_f90:.o=.f90) start.f90 run.f90 debug_c.c
TAGS: $(TAG_SOURCE)
	etags $(TAG_SOURCE)

# Clean up, but leave makefile and TAGS
clean:
	rm -f *.x *.o *.mod *.d
	rm -f work.pc work.pcl

# Clean up all
cleann: clean
	rm -f makefile TAGS
