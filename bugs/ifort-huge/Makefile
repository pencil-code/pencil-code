# Makefile
#
# Usage:
#   make f90
#   make FC=ifort test


.SUFFIXES: .f90 .x

.f90.x:
	$(FC) $(FFLAGS) $*.f90 -o $*.x

code: huge.x

test: huge.x
	huge.x
	@echo "Should give:"
	@echo "   1.000000      3.4028235E+38"

clean:
	rm -rf *.x work.pc work.pcl

cleann: clean

# End of file
