# Makefile
#
# Usage:
#   make f90
#   make FC=ifort test

.SUFFIXES: .f90 .x

.f90.x:
        ftn allocatable.f90 -o allocatable.x

test: allocatable.x
        ./allocatable.x
        @echo "   "
        @echo "Should give on both: 2"

clean:
        rm -rf *.x

# End of file
