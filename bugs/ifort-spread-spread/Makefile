# Makefile

.SUFFIXES: .f90 .x

.f90.x:
	$(FC) $(FFLAGS) $*.f90 -o $*.x

XFILES=spread-tmp.x spread_255.x spread_256.x

code: $(XFILES)

test: $(XFILES)
	spread-tmp.x
	spread_255.x
	spread_256.x

clean:
	rm -rf *.x


# End of file
