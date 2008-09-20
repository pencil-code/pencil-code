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
	./start.x

start.x: general.o start.o
	$(G95) $(FFLAGS) general.o start.o -o start.x

start.o: start.f90 general.o

clean:
	rm -f *.o *.x *.mod

# End of file
