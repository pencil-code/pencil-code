# Makefile
#
# Usage:
#   make FC=ifort8.0 test


.SUFFIXES: .f90 .x

.f90.x:
	$(FC) $(FFLAGS) $*.f90 -o $*.x

PROG=minmask

code: $(PROG).x

test: $(PROG).x
	$(PROG).x

clean:
	rm -rf *.x work.pc work.pcl

cleann: clean

# End of file
