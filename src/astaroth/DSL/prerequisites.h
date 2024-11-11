// Provides all declarations and functions needed for the formulation of the PDEs' rhss by DSL code
// and finally for the definition of the solve kernel.
#define IN_DSL 1

#define double real
#define cpu_pow pow
#define REAL_MAX AC_REAL_MAX

#include "../../PC_moduleflags.h"
#include "../stdlib/math"
#include "../stdlib/utils/intrinsics.h"

#include "../../../cparam_c.h"
#include "PC_nghost.h"

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
#include "../stdlib/map.h"
#include "PC_modulepardecs.h"
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
  #define LNRHO RHO
#endif
#include "../bcs/funcs.h"
#include "../bcs/funcs_overload.h"

#ifdef LFORCING
  #include "../forcing/pcstyleforcing.h"
  #define ADDFORCE + force
#else
  #define ADDFORCE
#endif

#include "equations.h"
