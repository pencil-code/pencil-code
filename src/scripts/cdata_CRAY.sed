#restrict to necessary variables
1,/BEGIN C BINDING/ d
/END C BINDING/, $ d
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
#use 3-vector types
s/integer *, *dimension *( *3 *) *::/extern int3arr,/
s/logical *, *dimension *( *3 *) *::/extern bool3arr,/
s/real *, *dimension *( *3 *) *::/extern real3arr,/
#remove array initializations
s/(\/.*\/)//g
#replace double precision exponent symbol by single precision one
s/\([0-9.]\) *[dD] *\([-0-9]\)/\1E\2/g
#remove lines containing implicit none
/implicit  *none/ d
#make KIND lowercase
s/KIND *=/kind=/
s/integer *( *kind *= *ikind8 *) *::/long long/
s/real *( *kind *= *rkind8 *) *::/ double/
#remove volatile
s/, *volatile//
#remove comment at line end
s/\([^ ]\) *!.*$/\1/
#insert pragma and includes instead of module
#s/^ *module .*$/  # pragma once \n  namespace PC\n{\n  # include "headers_c.h" \n  # include "defines_cdata.h"/ 
s/^ *module .*$/\n/
#replace dimension...:: by *::
s/, *dimension(\([^,:]*\)) *:: *\([a-zA-Z0-9_]*\)/:: \2[\1]/
s/, *dimension.*::/*::/
#transform parameter attribute to const
s/integer *, *parameter *\([\*]*\):: *//
s/real *, *parameter *\([\*]*\):: *//
#transform integer, real, logical, double precision with and without * into extern <type>, replace :: by ,. * remains
#s/integer *\([*]*\):: */extern int \1,/
#s/real *\([*]*\):: */extern real \1,/
#s/logical *\([*]*\):: */extern bool \1,/
#s/double  *precision *\([*]*\):: */extern double \1,/
#s/integer *( *kind *= *ikind8 *) *\([*]*\)::/extern long long \1,/
#s/real *( *kind *= *rkind8 *) *\([*]*\)::/extern double \1,/
##remove initializations for non-constant items
#/const/ !s/ *= *[a-zA-Z0-9_]* *([^)][^)]*)[^,]*,/,/g
#/const/ !s/ *= *[a-zA-Z0-9_]* *([^)][^)]*)[^,]*$//
#/const/ !s/ *= *[^,]*,/,/g        
#/const/ !s/ *= *[^,]* *$//g
#insert pointer assignment for array quantities
#/^ *extern  *[a-zA-Z][a-zA-Z]* *\*,/ {
#                s/, *\([a-zA-Z0-9_][a-zA-Z0-9_]*\)/,*\1/g
#               }
#s/\(^ *extern  *[a-zA-Z][a-zA-Z]* *\)[\*], *[\*]/\1,*/
#/extern/ {
      # remove terminating , and & in continuation lines
      ##s/, *& *$//
      ##/^ *, *$/ d
#     }
#remove first comma if not in a continuation line 
/^ *extern[^,]*,/ s/,/ /
#set ; when there is no continuation symbol " or > (for include statements)
s/\([^&">]\) *$/\1;/
#remove continuation symbol
s/& *$//
#remove endmodule
/^ *endmodule.*$/ d
#change all variable names to lowercase and create defines for those, which contain uppercase letters
#/[A-Z]/ {:loop
#         s/\([A-Za-z_0-9]*[A-Z][A-Za-z_0-9]*\)\([,;]\)/#define \1 \L\1\2\1\2/
#         t cont
#         b
#         :cont
#         h
#         s/^.*\(#[^,;]*\)[;,].*$/\1/
#         w tmp
#         x
#         s/\(.*\)#define[^;,]*[;,]\(.*\)/\1\2/
#         b loop
#        }
