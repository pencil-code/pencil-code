### Makefile for modular pencil code

LIBSIG=-L$(HOME)/f90/lib
LDMPI=-lmpi
MPICOMM=nompicomm
FORCING=forcing
ENTROPY=entropy
MAGNETIC=magnetic
IO=io_dist

### Begin machine dependent

## IRIX64:
#FC=f90
#FFLAGS= -64 -O3 -C -macro_expand  #(Antares)
#FFLAGS= -pfalist -64 -O3 -mips4 -r10000 -C -macro_expand  #(Grand)
#FFLAGS= -64 -O3 -mips4 -r12000 #(Ukaff)
#FFLAGS= -64 -O3 -mips4 #(ukaff|origin)
#F77=$(FC)

## OSF1:
#FC=f90 #(Mhd.)
#OPTFLAGS=
#OPTFLAGS=-fast -O5
#FFLAGS=$(OPTFLAGS) -tune ev6 -arch ev6 #(Mhd.)
#F77=$(FC)
#LIBSIG=

## Linux:
#FC=f95		#(norlap2)
#FC=pgf90	#(norhpc)
#FFLAGS=-O2 -Minfo 	#(norhpc)
#FFLAGS= -fast -DPGF90 #(Nq[0-9]+)
#FFLAGS=-O4 -Wc,-malign-double	#(norlap2)
#LDMPI=-lmpich #(Nq[0-9]+)
#LDMPI= #(norlap2)
#MPICOMM=nompicomm#(norlap2)
#FC=mpif90 # (Samos)
#FFLAGS=-O # (Samos)
#LDMPI= # (Samos)

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

# Switch off MPI for all machines (unless commented out)
#LDMPI=
#MPICOMM=nompi

LD=$(FC)
LDFLAGS=$(FFLAGS)

.SUFFIXES: .f90 .mod

.f90.o:
	$(FC) $(FFLAGS) -c $*.f90
.mod.o:
	$(FC) $(FFLAGS) -c $*.f90

#####

default: makefile code


makefile: Makefile
	adapt-mkfile Makefile makefile
code:
	make start.x run.x

comm_deps=cparam.o cdata.o general.o slices.o deriv.o $(IO).o sub.o $(MPICOMM).o register.o $(FORCING).o $(ENTROPY).o $(MAGNETIC).o
start_deps=$(comm_deps)
run_deps=$(comm_deps) equ.o timestep.o

start=start.o $(start_deps)
start.x: $(start)
	$(LD) $(LDFLAGS) $(start) $(LDMPI) -o start.x 

run=run.o $(run_deps)
run.x: $(run)
	$(LD) $(LDFLAGS) $(run) $(LDMPI) $(LIBSIG) -lsig -o run.x 

start.o:     start.f90 $(start_deps)
run.o:       run.f90 $(run_deps)
#
cparam.o:    cparam.f90
general.o:   general.f90
cdata.o:     cdata.f90 cparam.o
slices.o:    slices.f90 cdata.o
deriv.o:     deriv.f90 cdata.o
$(IO).o:     $(IO).f90 cdata.o
sub.o:       sub.f90 cdata.o deriv.o $(IO).o $(MPICOMM).o
register.o:  register.f90 cdata.o sub.o $(MPICOMM).o $(ENTROPY).o $(MAGNETIC).o
equ.o:       equ.f90 cdata.o general.o sub.o deriv.o $(ENTROPY).o
timestep.o:  timestep.f90 cdata.o $(FORCING.o) $(MPICOMM).o
#
$(MPICOMM).o:  $(MPICOMM).f90 cdata.o general.o
$(FORCING).o:  $(FORCING).f90 cdata.o $(MPICOMM).o sub.o
$(ENTROPY).o:  $(ENTROPY).f90 cdata.o deriv.o $(MPICOMM).o sub.o slices.o
$(MAGNETIC).o: $(MAGNETIC).f90 cdata.o deriv.o $(MPICOMM).o sub.o slices.o

TAG_SOURCE=$(run_deps:.o=.f90) start.f90 run.f90
TAGS: $(TAG_SOURCE)
	etags $(TAG_SOURCE)

clean:
	rm -f *.x *.o *.mod
	rm -f makefile
	rm -f TAGS
