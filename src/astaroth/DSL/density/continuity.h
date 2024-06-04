rhs=0.
glnrho = gradient(LNRHO)
#include "../density/diffusivity.h"
return rhs - dot(vecvalue(UU), glnrho) - divergence(UU)

