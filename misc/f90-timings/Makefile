# Makefile

# To do:
# - Generate different param.inc dynamically, as advertised
# - Use Pencil's adapt-mkfile framework (also as advertised)

.SUFFIXES:			# clear builtin suffixes
.SUFFIXES: .f90 .o .x


FC=g95
FFLAGS=-O3

MOD_O = timings_module.o array_valued_module.o
EXE = array-val-funcs.x

## Generic rules:
##
.f90.x:
	$(FC) $(FFLAGS) $*.f90 $(MOD_O) -o $*.x

.f90.o:
	$(FC) $(FFLAGS) $*.f90 -c -o $*.o


## Explicit rules:
##

default: test

code: $(EXE)

test: code
	for x in $(EXE); do \
	    echo "$${x}:"; \
	    ./$${x}; \
	done

# Dependencies of executables [most or all will depend on $(MOD_O)]:
array-val-funcs.x: array-val-funcs.f90 $(MOD_O)


# Other stuff
$(MOD_O): param.inc

param.inc: Makefile
	@touch param.inc
	@printf "%s\n" \
	    '! param.inc' \
	    '!' \
	    '! Written by Makefile; if you edit this file, you get what you deserve.' \
	    '' \
	    'integer, parameter :: nx=1024, niter=1000000' \
	    '' \
	    '! End of file' \
	    > param.inc



clean:
	rm -rf *.o *.x *.mod
	rm -rf param.inc

cleanall: clean
cleann: cleanall

# End of file
