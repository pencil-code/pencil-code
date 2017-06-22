/^ *type *[a-zA-Z0-9_]* *$/,/^ *end *type *[a-zA-Z0-9_]* *$/ d
/^ *!/ d
/^ *use/ d
/^ *public/ d
/^ *equivalence/ d
/^ *character/ d
/, *target/ d
/( *: *)/ d
s/.*/\L&/g
s/(\/.*\/)//g
s/\([0-9.]\) *[dD] *\([-0-9]\)/\1E\2/g
/implicit  *none/ d
/kind=ikind8/ d
s/\([^ ]\) *!.*$/\1/
//s/^ *module .*$/  # pragma once \n  namespace PC\n{\n  # include "headers_c.h" \n  # include "defines_cdata.h"/ 
s/^ *module .*$/  # pragma once \n  # include "headers_c.h" \n  # include "defines_cdata.h"/ 
s/, *dimension.*::/*::,/
s/integer *, *parameter *::/const FINT /
s/integer\([*]*\) *::/extern FINT \1/
s/real\([*]*\) *::/extern REAL \1/
s/logical\([*]*\) *::/extern bool \1/
s/double  *precision\([*]*\) *::/extern double \1/
/const/ !s/ *= *[^,]*,/,/g
/const/ !s/ *= *[^,]* *$//g
/extern/ {
      h
      s/^ *extern  *[a-zA-Z][a-zA-Z]* *\* *,/,/
      s/^ *extern  *[a-zA-Z][a-zA-Z]* */,/
      s/, *& *$//
      s/, *\([a-zA-Z0-9_][a-zA-Z0-9_]*\)/#define \1 MPREcdataMINF\1MSUF \n/g
      /^ *, *$/ d
      w tmp
      g
      s/\* *,/**,/
     }
/extern.*\*\*/ {
                s/, *\([a-zA-Z][a-zA-Z]*\)/,*\1/g
                s/\*\*,/ /
               }
s/\([^&">]\) *$/\1;/
s/& *$//
/^ *extern.*\*.*, *$/ {
                      : loop
                      n
                      s/ *= *[^,]*\([,;]\)/\1/g
                      s/ *= *[^,]* *$//g
                      s/\([^&">,]\) *$/\1;/
                      s/&.*$//
                      s/\([a-zA-Z0-9_][a-zA-Z0-9_]*\)/*\1/g
                      t cont
                      : cont
                      s/; *$/; /
                      t
                      b loop        
	             }
s/\*\([a-zA-Z0-9_][a-zA-Z0-9_]*\)/\1[1]/g
/^ *endmodule.*$/ d
