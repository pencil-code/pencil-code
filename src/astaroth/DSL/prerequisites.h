// Provides all declarations and functions needed for the formulation of the PDEs' rhss by DSL code
// and finally for the definition of the solve kernel.

real AC_xorig
real AC_yorig
real AC_zorig

#include "PC_nghost.h"
#include "../../PC_moduleflags.h"
#include "../stdlib/math"
#include "../stdlib/utils/intrinsics.h"

// declare here reduction results needed for the timestep

#ifdef LHYDRO
output real AC_maxadvec    // for all velocities - fluid and wave
#endif
#ifdef LENTROPY
output real AC_maxchi
#endif

#include "../../../cparam.inc_c.h"

#include "fieldecs.h"
#include "../stdlib/operators.h"
#include "../stdlib/integrators.h"
#include "../stdlib/units.h"
#include "../stdlib/utils/kernels.h"
#include "../phys_consts.h"
#include "PC_modulepardecs.h"
#ifdef LDENSITY
  int AC_ldensity_nolog
  #define LNRHO RHO
  #define ldensity_nolog AC_ldensity_nolog
#endif
#ifdef LFORCING
  #include "../forcing/pcstyleforcing.h"
  #define ADDFORCE if (step_num==2) {uu += forcing()}
#else
  #define ADDFORCE 
#endif
//int mcom = mvar+maux_com
/*real fbcx[mcom][2], fbcx_2[mcom][2]
real fbcy[mcom][2], fbcy_1[mcom][2], fbcy_2[mcom][2]
real fbcz[mcom][2], fbcz_1[mcom][2], fbcz_2[mcom][2]
*/  
#include "equations.h"
