/^ *type *[a-zA-Z0-9_]* *$/,/^ *end *type *[a-zA-Z0-9_]* *$/ d
/^ *!/ d
s/.*/\L&/g
#s/huge *(0)/std::numeric_limits<FINT>::max()/g 
s/huge *(0)/INT_MAX/g 
s/huge *(0\.0*)/AC_REAL_MAX/g 
s/huge *(0\.*0*d0)/DBL_MAX/g 
s/\([0-9.]\) *[dD] *\([-0-9]\)/\1E\2/g
/tiny *(/ d
/selected[a-z_]* *(/ d
/implicit  *none/ d
/cparam_pencils\.inc/ d
/dimension *(/ d
s/\([a-zA-Z0-9_]*\) *\*\* *\([^,]*\)\([,$]\)/pow(\1,\2)\3/g
s/\([a-zA-Z0-9_]*\) *\*\* *\([^,]*\) *$/pow(\1,\2)/g
#s/\([a-zA-Z0-9_]*\) *\*\* *\([(]?[-]?[a-zA-Z0-9_]*[)]?\)/pow(\1,\2)/
s/include *.\([a-z]*\.inc\). *$/# include "\1_c.h"/
s/include *.\([a-z]*\.local\). *$/# include "\1_c.h"/
s/\([^ ]\) *!.*$/\1/
s/^ *module .*$/# pragma once \n# include <limits.h> \n# include "headers_c.h"/ 
/end *module / d
s/integer *( *kind *= *ikind8 *) *, *parameter *::/const long long /
s/integer *, *parameter *::/const FINT /
s/logical *, *parameter *::/const int /
s/real *, *parameter *::/const REAL /
s/\([^0-9a-zA-Z_]\)min(/\1MIN(/g
s/\([^0-9a-zA-Z_]\)max(/\1MAX(/g
s/double  *precision *, *parameter *::/const double /
/,kind=ikind8/ d
s/\.true\./true/g
s/\.false\./false/g
s/\(^[^#].*[^&">]\) *$/\1;/
s/const  *FINT  *nghost *= *\([0-9]*\) *;/const FINT nghost=\1;\n#ifdef NGHOST \n#undef NGHOST \n#endif \n#define NGHOST \1\n#ifdef STENCIL_ORDER\n#undef STENCIL_ORDER\n#endif \n#define STENCIL_ORDER (2\*\1)/
s/& *$//

