# Makefile

.SUFFIXES: .o .f90 .x

G95=g95
FFLAGS=-r8 -freal=nan


# Generic rule:
.f90.o:
	$(G95) $(FFLAGS) -c $*.f90

# Specific rules:
default : test

test: start.x
	env G95_FPU_INVALID=Yes ./start.x

start.x: general.o start.o
	$(G95) $(FFLAGS) general.o start.o -o start.x

start.o: start.f90 general.o

tarball:
	(cd ..; tar czf g95-fp-invalid.tgz g95-fp-invalid/*)

clean:
	rm -f *.o *.x *.mod

# End of file
