// checked 17.6
rhs=real3(0.,0.,0.)
#include "../magnetic/magnetic_diffusivity.h"
#if LHYDRO
  rhs += cross(vecvalue(UU), curl(AA))
#endif
return rhs

