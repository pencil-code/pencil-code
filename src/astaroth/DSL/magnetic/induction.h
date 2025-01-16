// checked 17.6
rhs=real3(0.,0.,0.)
#include "../magnetic/magnetic_diffusivity.h"
#if LHYDRO
  if (linduction){
    if (lupw_aa){
      rhs += UU*gradient_tensor(AA) - ugrad_upw(AA,UU)
    } else {
      rhs += cross(UU, curl(AA))
    }
  }
#endif
return rhs
