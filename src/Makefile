### $Id: Makefile,v 1.186 2003-06-16 16:29:56 dobler Exp $
### Makefile for modular pencil code

## Default settings: they can be overwritten by src/Makefile.local
## which should be present in each local to run directory.

LIBSIG=
LDMPI=-lmpi



MPICOMM      = nompicomm
GLOBAL       = noglobal
IO           =   io_dist
FFT          = fft
FFTPACK      = fftpack
POWER        = nopower_spectrum
DEBUG        =   debug_c

HYDRO        =   hydro
VISCOSITY    =   visc_const
DENSITY      =   density
ENTROPY      = noentropy
MAGNETIC     = nomagnetic
RADIATION    = noradiation
PSCALAR      = nopscalar
DUSTVELOCITY = nodustvelocity
DUSTDENSITY  = nodustdensity
INTERSTELLAR = nointerstellar
IONIZATION   = noionization
STRUCT_FUNC  = nostruct_func
INITCOND_SPEC= noinitcond_spec

GRAVITY      = nograv
FORCING      = noforcing
SHEAR        = noshear
TIMEAVG      = notimeavg

## This allows us to add -r8 for double precision in Makefile.local
FFLAGS_GENERAL=

## Overwrite settings that are local to each run directory
include Makefile.local

## Default settings
FC=f90
FFLAGS=-O3
FFLAGS=
F90FLAGS=
F77FLAGS=
CFLAGS=-O3 -DFUNDERSC=1

## The following block (`Begin/End machine dependent') specifies compilers
## and switches for given architectures and machine names. See
## `perldoc adapt-mkfile' for full documentation of the syntax.
## Compiler flags:
##   FFLAGS are used by both the F90 compiler $(FC) and $(F77);
##   F90FLAGS are only used by $(FC)
##   F77FLAGS are only used by $(F77)

### Begin machine dependent

## IRIX64:
#FC=f90
#FFLAGS= -64 -O3 -mips4
#FFLAGS= -auto_use mpi_interface -64 -O3 -C -macro_expand  #(gridur|embla|balder)
#FFLAGS= -auto_use mpi_interface -64 -O3 -C -macro_expand  #(Antares)
#FFLAGS= -pfalist -64 -O3 -mips4 -r10000 -C -macro_expand  #(Grand)
#FFLAGS= -64 -O3 -mips4 -r12000 #(Ukaff)
#FFLAGS= -64 -O3 -auto_use mpi_interface -mips4 #(ukaff|origin)
#FFLAGS= -64 -O3 -mips4 -r12000 -flist #(green)
#FFLAGS= -64 -O3 -mips4 -r12000 #(wrenstandard)
#FFLAGS= -64 -O3 -mips4 -r12000 -LNO:cs1=32k:cs2=8M #(wrencache)
#FFLAGS= -64 -O3 -mips4 -r12000 -LNO:cs1=32k:cs2=8M -I/share/apps/unsupported/pacx/standard/include -L/share/apps/unsupported/pacx/standard/lib #(wrenpacx)
#LDMPI= -L/share/apps/unsupported/pacx/standard/lib -lpacxf -lpacx -lmpi #(wrenpacx)
#FFLAGS= -64 -O3 -mips4 -r12000 -LNO:cs1=32k:cs2=8M -I/share/apps/unsupported/pacx/standard/include -L/share/apps/unsupported/pacx/standard/lib #(wren)
#LDMPI= -L/share/apps/unsupported/pacx/standard/lib -lpacxf -lpacx -lmpi #(wren)
##FC=/c/brandenb/pacx/bin/pacxfc #(tacsg2)
#FFLAGS=-64 -O3 -mips4 #(tacsg2)
#OPTFLAGS= -64
#F77=$(FC)
#CFLAGS=-O3 -64 -DFUNDERSC=1

## OSF1:
#FC=f95 #(Mhd.)
#FC=f90 #(Mhd.)
#OPTFLAGS=-fast -O5 
#OPTFLAGS=-fast -O5 -notransform_loops #(Mhd.)
#FFLAGS=$(OPTFLAGS) -tune ev6 -arch ev6 #(Mhd.)
#F77=$(FC)
#LIBSIG=-lsig

## Linux:
## The -tpp7 flag leads to a 10% speedup on nl6
#FC=f95
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
#FC=f95 -I/usr/lib/lam/include -mismatch #(morvern)      
#FC=ifc -I/usr/lib/lam/include #(morvern)  
##FC=/home/brandenb/pacx/bin/pacxfc #(ns0)
#FC=ifc #(nl6|mmtf|nl19)
#FC=f95 #(karbala)
#FC=pgf90 #(erland|igor|tycho|linux29|onsager|fe4)
#FFLAGS=-O # (Samos)
#FFLAGS= # (fire)
#LDMPI=-static -L/usr/local/lib/LAM/lib -llamf77mpi -lmpi -llam -lnsl -lutil -L/usr/local/lib/PGI/linux86/lib -lpgc #(fe[0-4])
#LDMPI=-L/usr/lib/lam/lib -llam -lmpi -lpmpi -llamf77mpi -llammpi++ -llammpio #(morvern)
#LDMPI=-static -L/opt/lam/lib -lmpi -llamf77mpi -llam  #(nw50)
##FFLAGS=-O3 -ipo #PACX #(fe[0-3])
#FFLAGS=-fast -lm #(erland|igor|tycho|linux29|onsager|fe4)
#F90FLAGS= #(erland|igor|tycho|linux29|onsager|fe4)
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
#F77=$(FC) -lm  #(tycho|erland|igor|onsager|fe4)
#F77=$(FC) #(Cincinnatus|Owen|Master|cosmo|nl6)


## SunOS:
#FC=f90
#FC=pghpf #(Lomond)
#FFLAGS=-O4
#FFLAGS= -Msmp -O4  #(Lomond)
#CFLAGS=-O4 -DFUNDERSC=1
#CFLAGS=-O4 -DFUNDERSC=0 #(ursa)
#LDMPI=-L/opt/local/mpich/lib -lmpich -lsocket -lnsl -laio
#F77=$(FC)

## UNICOS/mk:
#FC=f90
#FFLAGS=-O3 -em

## HI-UX/MPP:
#FC=mpif90
#F77=$(FC)
#FFLAGS=-lowercase -I/usr/mpi/include -O4 -pvec -parallel=0
#FFLAGS=-lowercase -OSS -noparallel -O4
#CFLAGS=-O4 -DFUNDERSC=0
#LDMPI=-L/usr/mpi/lib/lib32 -lfmpi -lmpi

## SUPER-UX:
#FC=f90 -I/usr/include
#F77=$(FC)
#FFLAGS=-C hopt

## AIX:
#FC=xlf90
#FC=mpxlf90 #(Sp[12])
#FFLAGS=-O5 -qsuffix=f=f90

### End machine dependent
LD=$(FC)
LDFLAGS=$(FFLAGS) $(FFLAGS_GENERAL)
# -Mpencil.map

## Generic rules:
.SUFFIXES: .f90 .mod

.f.o:
	$(F77) $(FFLAGS) $(FFLAGS_GENERAL) $(F77FLAGS) -c $*.f
.f90.o:
	$(FC) $(FFLAGS) $(FFLAGS_GENERAL) $(F90FLAGS) -c $*.f90
.mod.o:
	$(FC) $(FFLAGS) $(FFLAGS_GENERAL) $(F90FLAGS) -c $*.f90
.c.o:
	$(CC) $(CFLAGS) -c $*.c

#####

default: makefile cparam.inc code

pacx: pacx-makefile cparam.inc code

makefile: Makefile
	adapt-mkfile Makefile makefile
pacx-makefile: Makefile
	pacx-adapt-mkfile Makefile pacx-makefile
	cp pacx-makefile makefile
cparam.inc: Makefile Makefile.local
	mkcparam Makefile Makefile.local -o cparam.inc
mpicomm_double.f90: mpicomm.f90
	sed -e 's/REAL/DOUBLE_PRECISION/g' mpicomm.f90 > mpicomm_double.f90
mpicomm_.f90: mpicomm.f90
	sed -e 's/\(MPI_[A-Z_]*\)\([ ]*(\)/\1_\2/' mpicomm.f90 | sed -e 's/\(MPI_[A-Z_]*\)\([ ]*!\)/\1_\2/' > mpicomm_.f90
mpicomm_double_.f90: mpicomm_double.f90
	sed -e 's/\(MPI_[A-Z_]*\)\([ ]*(\)/\1_\2/' mpicomm_double.f90 | sed -e 's/\(MPI_[A-Z_]*\)\([ ]*!\)/\1_\2/' > mpicomm_double_.f90

code:
	make start.x run.x read_videofiles.x

technical=cparam.o cdata.o general.o slices.o deriv.o boundcond.o $(INITCOND_SPEC).o initcond.o register.o prints.o $(TIMEAVG).o
physics=$(GLOBAL).o $(VISCOSITY).o $(HYDRO).o $(DENSITY).o $(FORCING).o $(GRAVITY).o $(INTERSTELLAR).o $(ENTROPY).o $(IONIZATION).o $(MAGNETIC).o $(SHEAR).o $(PSCALAR).o $(RADIATION).o $(DUSTVELOCITY).o $(DUSTDENSITY).o
comm_deps=$(technical) sub.o $(IO).o wsnaps.o param_io.o $(MPICOMM).o $(physics) $(FFT).o $(FFTPACK).o $(POWER).o $(STRUCT_FUNC).o
start_deps_f90=$(comm_deps) $(DEBUG).o
start_deps=$(start_deps_f90)
run_deps_f90=$(comm_deps) equ.o timestep.o
run_deps=$(run_deps_f90) $(DEBUG).o
start=start.o $(start_deps)
run=run.o $(run_deps)
read_videofiles_deps=cparam.o general.o
read_videofiles=read_videofiles.o $(read_videofiles_deps)

start.x: $(start) Makefile.local
	$(LD) $(LDFLAGS) $(start) $(LDMPI) -o start.x 

run.x: $(run) Makefile.local
	$(LD) $(LDFLAGS) $(run) $(LDMPI) $(LIBSIG) -o run.x 

read_videofiles.x: $(read_videofiles)
	$(LD) $(LDFLAGS) $(read_videofiles) -o read_videofiles.x 

start.o:     start.f90 $(start_deps)
run.o:       run.f90 $(run_deps)
read_videofiles.o: read_videofiles.f90 $(read_videofiles_deps)
#
cparam.o:    cparam.f90 cparam.local cparam.inc
general.o:   general.f90 cparam.o
cdata.o:     cdata.f90 cparam.o
slices.o:    slices.f90 cdata.o sub.o
deriv.o:     deriv.f90 cdata.o
fft.o:	     fft.f
$(FFTPACK).o:   $(FFTPACK).f
#fft_nr.o:    fft_nr.f90
sub.o:       sub.f90 cdata.o deriv.o $(MPICOMM).o general.o
param_io.o:  param_io.f90 timestep.o $(TIMEAVG).o $(physics)
register.o:  register.f90 cdata.o $(MPICOMM).o sub.o $(IO).o param_io.o $(physics) $(TIMEAVG).o prints.o
prints.o:    prints.f90 cdata.o sub.o $(physics)
equ.o:       equ.f90 cdata.o $(MPICOMM).o sub.o $(IO).o $(physics) boundcond.o \
             $(GLOBAL).o slices.o
boundcond.o: boundcond.f90 cdata.o $(ENTROPY).o $(MAGNETIC).o $(DENSITY).o $(RADIATION).o $(HYDRO).o
initcond.o:  initcond.f90 cdata.o general.o
kepvor.o : kepvor.f90 cdata.o general.o
noinitcond_spec.o: noinitcond_spec.f90 cdata.o general.o
timestep.o:  timestep.f90 cdata.o $(MPICOMM).o equ.o $(FORCING).o
wsnaps.o:    wsnaps.f90 cdata.o $(MPICOMM).o $(IO).o boundcond.o $(POWER).o $(STRUCT_FUNC).o $(RADIATION).o $(IONIZATION).o
#
$(MPICOMM).o:      $(MPICOMM).f90 cdata.o general.o $(FFT).o 
$(GLOBAL).o:       $(GLOBAL).f90 cdata.o $(IO).o $(MPICOMM).o
$(VISCOSITY).o:    $(VISCOSITY).f90 cdata.o cparam.o $(DENSITY).o sub.o
$(SHEAR).o:	   $(SHEAR).f90 $(MPICOMM).o cdata.o sub.o cparam.o $(HYDRO).o
$(GRAVITY).o:      $(GRAVITY).f90 cdata.o $(MPICOMM).o sub.o slices.o $(GLOBAL).o
$(FORCING).o:      $(FORCING).f90 cdata.o $(MPICOMM).o sub.o $(HYDRO).o
$(ENTROPY).o:      $(ENTROPY).f90 cdata.o $(MPICOMM).o sub.o initcond.o $(INITCOND_SPEC).o $(GLOBAL).o $(VISCOSITY).o $(HYDRO).o $(IO).o slices.o $(INTERSTELLAR).o $(IONIZATION).o
$(MAGNETIC).o:     $(MAGNETIC).f90 cdata.o $(MPICOMM).o sub.o $(IO).o initcond.o $(INITCOND_SPEC).o $(DENSITY).o slices.o $(FFT).o $(GRAVITY).o
$(RADIATION).o:	   $(RADIATION).f90 cdata.o $(MPICOMM).o sub.o initcond.o $(INITCOND_SPEC.o) $(DENSITY).o $(IONIZATION).o
$(HYDRO).o:        $(HYDRO).f90 cdata.o $(MPICOMM).o sub.o initcond.o $(INITCOND_SPEC).o $(VISCOSITY).o $(GRAVITY).o $(DENSITY).o $(MAGNETIC).o slices.o
$(DENSITY).o:      $(DENSITY).f90 cdata.o $(MPICOMM).o sub.o $(IO).o initcond.o $(INITCOND_SPEC).o $(GRAVITY).o slices.o
$(PSCALAR).o:      $(PSCALAR).f90 cdata.o $(MPICOMM).o sub.o initcond.o $(INITCOND_SPEC).o
$(DUSTVELOCITY).o: $(DUSTVELOCITY).f90 cdata.o $(MPICOMM).o sub.o initcond.o $(INITCOND_SPEC).o $(HYDRO).o $(GRAVITY).o
$(DUSTDENSITY).o:  $(DUSTDENSITY).f90 cdata.o $(MPICOMM).o sub.o initcond.o $(DUSTVELOCITY).o
$(INTERSTELLAR).o: $(INTERSTELLAR).f90 cdata.o $(MPICOMM).o sub.o initcond.o $(INITCOND_SPEC).o $(IONIZATION).o #Expected deps!
$(IONIZATION).o:   $(IONIZATION).f90 cdata.o $(MPICOMM).o sub.o initcond.o $(INITCOND_SPEC).o $(DENSITY).o
$(IO).o:           $(IO).f90 cdata.o sub.o $(MPICOMM).o $(DEBUG).o
#$(TIMEAVG).o:      $(TIMEAVG).f90 cdata.o sub.o $(IO).o $(MPICOMM).o ctimeavg.local
$(TIMEAVG).o:      $(TIMEAVG).f90 cdata.o sub.o $(IO).o $(MPICOMM).o
$(FFT).o:	   cdata.o
nodebug.o:         nodebug.f90 cdata.o $(MPICOMM).o
$(POWER).o:        $(POWER).f90 cdata.o sub.o $(MPICOMM).o general.o
$(STRUCT_FUNC).o:  $(STRUCT_FUNC).f90 cdata.o sub.o $(MPICOMM).o general.o
#axeltest $(INITCOND_SPEC).o:$(INITCOND_SPEC).f90

TAG_SOURCE=$(run_deps_f90:.o=.f90) start.f90 run.f90 debug_c.c
TAGS: $(TAG_SOURCE)
	etags $(TAG_SOURCE)

# Clean up, but leave makefile and TAGS
clean:
	rm -f *.x *.o *.mod *.d *.inc *.il
	rm -f work.pc work.pcl

# Clean up all
cleann: clean
	rm -f makefile TAGS

# Remove all build files, but keep executables
run:
	rm -f *.o *.mod *.d *.inc
	rm -f work.pc work.pcl


# Alternative to make run. To copy executables to parent directory.
# Make clean(n) could then be run on the source.
install:
	cp *.x ../
