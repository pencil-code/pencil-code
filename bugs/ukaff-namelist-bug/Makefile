### Makefile for modular pencil code



## Default settings
FC=f90
CFLAGS=-O3 -DFUNDERSC=1

### Begin machine dependent

## IRIX64:
#FC=f90
#FFLAGS= -auto_use mpi_interface -64 -O3 -C -macro_expand  #(Antares)
#FFLAGS= -pfalist -64 -O3 -mips4 -r10000 -C -macro_expand  #(Grand)
#FFLAGS= -64 -O3 -mips4 -r12000 #(Ukaff)
#FFLAGS= -auto_use mpi_interface -64 -O0 -fullwarn -C -mips4 #(ukaff|origin)
#FFLAGS= -64 -O0 -DEBUG:verbose_runtime=ON -C -g -mips4 #(ukaff|origin)
#FFLAGS= -n32
#FFLAGS= -64
#F77=$(FC)
CFLAGS=-O3 -64 -DFUNDERSC=1

## OSF1:
#FC=f90 #(Mhd.)
#OPTFLAGS=
#OPTFLAGS=-fast -O5
#FFLAGS=$(OPTFLAGS) -tune ev6 -arch ev6 #(Mhd.)
#F77=$(FC)
#LIBSIG=

## Linux:
#FC=f95
#FFLAGS=-O3 -u
#LDMPI=
#FC=mpif90 # (Samos)
#FC=/usr/lib/lam/bin/mpif95 #(Cincinnatus|Owen)
#FC=pgf90 #(Nq[0-9]+)
#FC=ifc	#(Nq[0-9]+)
#FC=/home/brandenb/bin/f95 #(Nq[0-9]+)
#FFLAGS=-O # (Samos)
#FFLAGS= -fast -DPGF90 #(Nq[0-9]+) ##(old)
#FFLAGS=-O3 -I/usr/lib/lam/include -static #(Nq[0-9]+)
#FFLAGS=-O4 -u -Wc,-malign-double	#(Nl3|Mephisto)
#LDMPI=-llammpi++ -llammpio -lpmpi -llamf77mpi -lmpi -llam -lnsl -lutil #(Nq[0-9]+)
#LDMPI= # (Samos)
#CFLAGS=-O3 -DFUNDERSC=1 -Wall

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
	make run.x


run_deps=cdata.o magnetic.o hydro.o param_io.o
run=run.o $(run_deps)

run.x: $(run)
	$(LD) $(LDFLAGS) $(run) $(LDMPI) $(LIBSIG) -o run.x 
run.o:       run.f90 $(run_deps)
cdata.o:     cdata.f90
param_io.o:  param_io.f90 magnetic.o hydro.o cdata.o
magnetic.o:  magnetic.f90 cdata.o
hydro.o:     hydro.f90 cdata.o magnetic.o
#


# Clean up, but leave makefile and TAGS
clean:
	rm -f *.x *.o *.mod *.d
	rm -f work.pc work.pcl

# Clean up all
cleann: clean
	rm -f makefile TAGS
