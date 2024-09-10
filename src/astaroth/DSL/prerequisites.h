// Provides all declarations and functions needed for the formulation of the PDEs' rhss by DSL code
// and finally for the definition of the solve kernel.
#define IN_DSL 1

#include "PC_nghost.h"
#include "../../PC_moduleflags.h"
#include "../stdlib/math"
#include "../stdlib/utils/intrinsics.h"

#include "../../../cparam.inc_c.h"
const real impossible=3.9085e37;
const int mcom = mvar+maux_com
#define mcom 9  //!!!
//#include "../../../cparam_c.h"

#include "fieldecs.h"
#include "../stdlib/operators.h"
#include "../stdlib/integrators.h"
#include "../stdlib/units.h"
#include "../stdlib/utils/kernels.h"
#include "../phys_consts.h"
#include "PC_modulepardecs.h"

// declare here reduction results needed for the timestep

#ifdef LHYDRO
  output real AC_maxadvec    // for all velocities - fluid and wave
#endif
#ifdef LENTROPY
  output real AC_maxchi
#endif

#ifdef LDENSITY
  int AC_ldensity_nolog
  #define LNRHO RHO
  #define ldensity_nolog AC_ldensity_nolog
#endif
#include "../bcs/funcs.h"
#include "../bcs/funcs_overload.h"

#ifdef LFORCING
  #include "../forcing/pcstyleforcing.h"
  #define ADDFORCE + force
#else
  #define ADDFORCE
#endif

gmem real AC_fbcx[mcom][2], AC_fbcx_2[mcom][2]
gmem real AC_fbcy[mcom][2], AC_fbcy_1[mcom][2], AC_fbcy_2[mcom][2]
gmem real AC_fbcz[mcom][2], AC_fbcz_1[mcom][2], AC_fbcz_2[mcom][2]

#include "equations.h"
