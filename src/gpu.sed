/^ *$/ b end
/^ *#/ b end
/^ *[A-Z0-9_]* *= *-/ b end
/io_/ b end
/DEBUG/ b end
/DERIV/ b end
/GHOSTFOLD/ b end
/GSL/ b end
/INITIAL_CONDITION/ b end
/FIXED_POINT/ b end
/LSODE/ b end
#/MPICOMM/ b end
/NSCBC/ b end
/POWER/ b end
/SIGNAL/ b end
/SLICES/ b end
/SOLID_CELLS/ b end
/STRUCT_FUNC/ b end
/SYSCALLS/ b end
/TESTPERTURB/ b end
/TIMEAVG/ b end
/VISCOSITY/ b end
/WENO_TRANSPORT/ b end
/STREAMLINES/ b end
/YINYANG/ b end
/PARTICLES/ b end
/MPICOMM/ b end
/MULTITHREADING/ b end
/GPU_VENDOR/ b end
/RUNTIME_COMPILATION/ b end
/TRANSPILATION/ b end
/GENERATE_DSL_CODE/ b end
/LIBRARIES/ b end
/FARRAY/ b end
/^ *[A-Z0-9_]* *= *no/ b end
s/^ *REAL_PRECISION *= *double *$/PRECISION=DOUBLE/
s/^ *REAL_PRECISION *= *8 *$/PRECISION=DOUBLE/
s/^ *REAL_PRECISION *= *single *$/PRECISION=SINGLE/
s/^ *REAL_PRECISION *= *4 *$/PRECISION=SINGLE/
t store
b cont
: store
H
b end 
: cont
s/^ *GPU *= *\([A-Za-z0-9_]*\) *$/\U\1=\U\1/
t store1
b cont1
: store1
#H
b end
: cont1
/IO[_ ]/ b end
#s/^.*= *\([A-Za-z0-9_][A-Za-z0-9_]*\)/#define \U\1/ 
s/\([a-zA-Z_0-9]\)  *\([a-zA-Z_0-9]\)/\1.f90 ..\/\2/g
s/^ *\([A-Z_][A-Z0-9_]*\) *= *\([A-Za-z0-9_][A-Za-z0-9_/. ]*\) *$/#define L\1 1 \n#define L\2_MODULE 1 \/\/ ..\/\2.f90/ 
p
s/.*\/\/ *\([a-zA-Z_0-9\.].*[a-zA-Z_0-9]\) *$/\1 \\/
H
: end
$! d
${
g
#if there is a PRECISION=... entry, put it as the first line, add \n MODULESOURCES= \ and append the parts
#before and after the PRECISION entry with an additional \ in between to fill the empty line
/DOUBLE/ {
          s/^\(.*\)PRECISION=DOUBLE\(.*\)$/PRECISION=DOUBLE\nMODULESOURCES= \\\1\\\2/
          t out
         }
/SINGLE/ {
          s/^\(.*\)PRECISION=SINGLE\(.*\)$/PRECISION=SINGLE\nMODULESOURCES= \\\1\\\2/
          t out
         }
#if there is no PRECISION=... entry, put as first line MODULESOURCES= \ and append the remaining part
s/^\(.*\)$/MODULESOURCES= \\\1/
: out
#write pattern space to PC_modulesources.h
w CUDA_MAKEDIR/PC_modulesources.h
#w PC_modulesources.h
#delete pattern space
d
}
