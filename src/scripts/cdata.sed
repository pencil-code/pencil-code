#remove user defined type definitions
/^ *type *[a-zA-Z0-9_]* *$/,/^ *end *type *[a-zA-Z0-9_]* *$/ d
#remove comment lines
/^ *!/ d
#remove lines containing use, public, equivalence, character or target
/^ *use/ d
/^ *public/ d
/^ *equivalence/ d
/^ *character/ d
/, *target/ d
#remove quantities used in an equivalence
s/, *z_allprocs\([^0-9a-zA-Z_]\)/\1/
s/\([^0-9a-zA-Z_]\)z_allprocs[ ,]/\1/
s/\([^0-9a-zA-Z_]\)z_allprocs$/\1/
s/, *zgrid\([^0-9a-zA-Z_]\)/\1/
s/\([^0-9a-zA-Z_]\)zgrid[ ,]/\1/
s/\([^0-9a-zA-Z_]\)zgrid$/\1/
#remove allocatable quantities
/( *:/ d
#make everything lowercase
#!!!s/.*/\L&/g
#remove array initializations
s/(\/.*\/)//g
#replace double precision exponent symbol by single precision one
s/\([0-9.]\) *[dD] *\([-0-9]\)/\1E\2/g
#remove lines containing implicit none
/implicit  *none/ d
s/integer *( *kind *= *ikind8 *) *::/long long/
s/real *( *kind *= *rkind8 *) *::/ double/ 
s/integer *( *KIND*= *ikind8 *) *::/long long/
s/real *( *KIND*= *rkind8 *) *::/ double/ 
#remove volatile
s/, *volatile//
#remove comment at line end
s/\([^ ]\) *!.*$/\1/
#insert pragma and includes instead of module
#s/^ *module .*$/  # pragma once \n  namespace PC\n{\n  # include "headers_c.h" \n  # include "defines_cdata.h"/ 
s/^ *module .*$/  #pragma once \n  #include "headers_c.h" \n  #include "defines_cdata.h" \n  #include "cdata_refs.h" \n/ 
#replace dimension...:: by *::
s/, *dimension.*::/*::/
#transform parameter attribute to const
s/integer *, *parameter *\([\*]*\):: */const FINT \1/
#transform integer, real, logical, double precision with and without * into extern <type>, replace :: by ,. * remains
s/integer *\([*]*\):: */extern FINT \1,/
s/real *\([*]*\):: */extern REAL \1,/
s/logical *\([*]*\):: */extern int \1,/
s/double  *precision *\([*]*\):: */extern double \1,/
s/integer *( *kind *= *ikind8 *) *\([*]*\):: */extern long long \1,/
s/real *( *kind *= *rkind8 *) *\([*]*\):: */extern double \1,/ 
#remove initializations for non-constant items
/const/ !s/ *= *[a-zA-Z0-9_]* *([^)][^)]*)[^,]*,/,/g
/const/ !s/ *= *[a-zA-Z0-9_]* *([^)][^)]*)[^,]*$//
/const/ !s/ *= *[^,]*,/,/g        
/const/ !s/ *= *[^,]* *$//g
#insert pointer assignment for array quantities
/^ *extern  *[a-zA-Z][a-zA-Z]* *\*,/ {
                s/, *\([a-zA-Z0-9_][a-zA-Z0-9_]*\)/,*\1=\&\1_/g
               }
/extern/ {
      h
      s/^ *extern  *[a-zA-Z][a-zA-Z]* *[\*]* *,/,/
      # remove terminating , and & in continuation lines
      s/, *& *$//
      s/, *\*\([a-zA-Z0-9_][a-zA-Z0-9_]*\) *= *\&[^,]*/#define \1_ MODULE_PREFIXcdataMODULE_INFIX\L\1\UMODULE_SUFFIX \n/g
      s/, *\([a-zA-Z0-9_][a-zA-Z0-9_]*\)/#define \1 MODULE_PREFIXcdataMODULE_INFIX\L\1\UMODULE_SUFFIX \n/g
      /^ *, *$/ d
      w tmp
      g
      /^ *extern  *[a-zA-Z][a-zA-Z]* *\* *,/ {
        s/, *\*\([a-zA-Z0-9_][a-zA-Z0-9_]*\) *= *\&[^,]*/,\1_/g
        s/\* *,//
        s/\& *$//
        s/[,]* *$/;/
        w cdata_refs.h
      }
      g
     }
#remove first comma if not in a continuation line 
/^ *extern[^,]*,/ s/,/ /
#set static attribute for array pointers
/^ *extern[^\*]*\*/ {
         s/extern/static/
         s/\*//
        }
#set ; when there is no continuation symbol " or > (for include statements)
s/\([^&">]\) *$/\1;/
#remove continuation symbol
s/& *$//
#remove endmodule
/^ *endmodule.*$/ d
