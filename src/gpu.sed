/^ *$/ b end
/^ *#/ b end
/^ *[A-Z0-9_]* *= *-/ b end
/io_/ b end
s/^ *[A-Z0-9_]* *= *noviscosity *$/#undef VISCOSITY/ 
t prin
b cont0
: prin
p
t end
: cont0
/^ *[A-Z0-9_]* *= *no/ b end
s/^ *REAL_PRECISION *= *double *$/DOUBLE_PRECISION=DOUBLE_PRECISION/
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
s/^ *\([A-Z0-9_][A-Z0-9_]*\) *= *\([a-z0-9_][a-z0-9_]*\) *$/#define \1 \/\/ ..\/..\/\2.f90/ 
p
s/.*\/\/ *\(\.\.\/\.\.\/[a-z0-9_][a-z0-9_]*\.f90\) *$/\1 \\/
H
: end
$! d
${
g
/DOUBLE_PRECISION/ {
                    s/^\(.*\)DOUBLE_PRECISION=DOUBLE_PRECISION\(.*\)$/DOUBLE_PRECISION=DOUBLE_PRECISION\nMODULESOURCES= \\\1\2\\/
                    t out
                   }
s/^\(.*\)$/MODULESOURCES= \\\1/
: out
w cuda/src_new/PC_modulesources
d
}
