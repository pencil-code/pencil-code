# Makefile
#
# Usage:
#   make f90
#   make FC=ifort test


.SUFFIXES: .f90 .x

.f90.x:
	$(FC) $(FFLAGS) $*.f90 -o $*.x

code: seed.x

test: seed.x
	seed.x

clean:
	rm -rf *.x work.pc work.pcl

cleann: clean

# End of file
