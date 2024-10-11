// Provides all declarations and functions needed for the formulation of the PDEs' rhss by DSL code
// and finally for the definition of the solve kernel.
#define IN_DSL 1

#define double real
#define cpu_pow pow
#define REAL_MAX AC_REAL_MAX

#include "PC_nghost.h"
#include "../../PC_moduleflags.h"
#include "../stdlib/math"
#include "../stdlib/utils/intrinsics.h"

#include "../../../cparam_c.h"

#define AC_n1 n1
#define AC_m1 m1
#define AC_l1 l1

#define AC_mxgrid mxgrid
#define AC_mygrid mygrid
#define AC_mzgrid mzgrid

#include "fieldecs.h"
#include "../stdlib/grid.h"
#include "../stdlib/derivs.h"
#include "../stdlib/operators.h"
#include "../stdlib/integrators.h"
#include "../stdlib/units.h"
#include "../stdlib/utils/kernels.h"
#include "PC_modulepardecs.h"
//TP: x,y and z are too generic macros
#undef x
#undef y
#undef z
#define AC_NGHOST NGHOST

// declare here reduction results needed for the timestep

#ifdef LHYDRO
  output real AC_maxadvec    // for all velocities - fluid and wave
#endif
#ifdef LENTROPY
  output real AC_maxchi
  #define LENERGY 1       // a hack for the moment
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
