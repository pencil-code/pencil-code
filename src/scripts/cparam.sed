/^ *type *[a-zA-Z0-9_]* *$/,/^ *end *type *[a-zA-Z0-9_]* *$/ d
/^ *!/ d
s/.*/\L&/g
#s/huge *(0)/std::numeric_limits<FINT>::max()/g 
s/huge *(0)/INT_MAX/g 
s/huge *(int.*)/INT_MAX/g 
s/huge *(0\.0*)/REAL_MAX/g 
s/huge *(0\.*0*d0)/DBL_MAX/g 
s/^\(.*rkind8.*\)huge *( *[a-zA-Z0-9_]* *)/\1HUGE_VAL/g
s/\([0-9.]\) *[dD] *\([-0-9]\)/\1E\2/g
/tiny *(/ d
/selected[a-z_]* *(/ d
/implicit  *none/ d
/cparam_pencils\.inc/ d
/dimension *(/ d
/integer, *parameter *:: *[ir]kind[0-9]* *=/ d
s/, *bind *( *c *)//
s/\([a-zA-Z0-9_]*\) *\*\* *\([^,]*\)\([,$]\)/cpu_pow(\1,\2)\3/g
s/\([a-zA-Z0-9_]*\) *\*\* *\([^,]*\) *$/cpu_pow(\1,\2)/g
#s/\([a-zA-Z0-9_]*\) *\*\* *\([(]?[-]?[a-zA-Z0-9_]*[)]?\)/cpu_pow(\1,\2)/
s/^ *include *.\([a-z]*\.inc\). *$/#if IN_DSL\n  #include "..\/..\/..\/\1_c.h"\n#else\n  #include "\1_c.h"\n#endif/
s/^ *include *.\([a-z]*\.local\). *$/#if IN_DSL\n  #include "..\/..\/..\/\1_c.h"\n#else\n  #include "\1_c.h"\n#endif/
s/\([^ ]\) *!.*$/\1/
s/^ *module .*$/\n#if IN_DSL\n #include "..\/..\/..\/headers_c.h"\n#else\n  \n#include <float.h>\n#include <limits.h>\n #include "headers_c.h"\n#endif\n#define y0 y0_\n/ 
/end *module / d
s/integer *( *kind *= *ikind8 *) *, *parameter *::/const long long /
s/integer *( *kind *= *ikind4 *) *, *parameter *::/const long /
s/real *( *kind *= *rkind8 *) *, *parameter *::/const double /
s/int *(\([a-zA-Z0-9_]*\) *, *kind *= *ikind8 *)/(long long)\1/g
s/integer *, *parameter *::/const int /
s/logical *, *parameter *::/const bool /
s/(\/ '/\ \/\*\* /g
s/' \/)/\*\*\//g
s|^    /)|**/|g
s/real *, *parameter *::/const real /
s/\([^0-9a-zA-Z_]\)min(/\1MIN(/g
s/\([^0-9a-zA-Z_]\)max(/\1MAX(/g
s/double  *precision *, *parameter *::/const double /
s/\.true\./true/g
s/\.false\./false/g
s/\(^[^#].*[^&">]\) *$/\1;/
s/& *$//
s/dbl_max/DBL_MAX/g 
s/\.or\./ | /g
s/\.and\./ & /g
/const *int *nghost *=/ {
h
s/const *int *nghost *= *\([0-9]*\) *;/#undef NGHOST\n#define NGHOST \1\n#undef STENCIL_ORDER\nhostdefine STENCIL_ORDER (2*NGHOST)\n/
t cont
b end
: cont
w astaroth/DSL/local/PC_nghost.h
: end
g
}
