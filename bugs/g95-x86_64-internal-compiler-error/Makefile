# Makefile

FC=g95
FFLAGS=-O4

## Generic rules:
.SUFFIXES:                      # get rid of that annoying Modula rule
.SUFFIXES: .f90 .o
.f90.o:
	$(FC) $(FFLAGS) -c $*.f90 



default: general.o

clean:
	rm -f *.o *.mod


# End of file
