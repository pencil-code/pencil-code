# Makefile
#
# Usage:
#   make f90
#   make FC=ifort test


.SUFFIXES: .f90 .x

.f90.x:
	$(FC) $(FFLAGS) $*.f90 -o $*.x

XFILES=spread-tmp.x spread_255.x spread_256.x spread_1791.x spread_1792.x

code: $(XFILES)

test: $(XFILES)
	for x in $(XFILES); do ./$$x; done

f90: spread.f90
	perl \
	  -e '@files=qw($(XFILES));' \
	  -e '@nums=map /^spread_([0-9]+).x$$/, @files;' \
	  -e 'foreach $$n (@nums) {' \
	  -e '  system("sed s/\\<N\\>/$$n/ < spread.f90 > spread_$$n.f90");' \
	  -e '}'

clean:
	rm -rf *.x work.pc work.pcl

cleann:
	make clean
	rm -rf spread_[0-9]*.f90
	

# End of file
