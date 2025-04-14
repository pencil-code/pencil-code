// Provides all declarations and functions needed for the formulation of the PDEs' rhss by DSL code
// and finally for the definition of the solve kernel.
#define IN_DSL 1
#define double real
#define cpu_pow pow
#define REAL_MAX AC_REAL_MAX

const int rkind8 = 0
const int prof_nz = 150
struct PC_rhs_update
{
	real3 dt
	real max_advec
}

#include "PC_nghost.h"

#define AC_mx AC_mlocal.x
#define AC_my AC_mlocal.y
#define AC_mz AC_mlocal.z

#define AC_nx AC_nlocal.x
#define AC_ny AC_nlocal.y
#define AC_nz AC_nlocal.z

#define AC_nxgrid AC_ngrid.x
#define AC_nygrid AC_ngrid.y
#define AC_nzgrid AC_ngrid.z

#define AC_dsx AC_ds.x
#define AC_dsy AC_ds.y
#define AC_dsz AC_ds.z

#include "../../PC_moduleflags.h"
#include "../stdlib/math"
#include "../stdlib/utils/intrinsics.h"

#include "../../../cparam_c.h"
#include "../../../cparam_pencils.inc_c.h"

#define AC_n1 n1
#define AC_m1 m1
#define AC_l1 l1

#define AC_mxgrid mxgrid
#define AC_mygrid mygrid
#define AC_mzgrid mzgrid

#include "fieldecs.h"
#include "../stdlib/grid.h"
#include "../stdlib/optimized_integrators.h"
//#include "../stdlib/units.h"
#include "../stdlib/utils/kernels.h"
//#include "../stdlib/map.h"
#include "PC_modulepardecs.h"
#include "../stdlib/derivs.h"
#include "../stdlib/operators.h"
#define AC_NGHOST NGHOST

// declare here reduction results needed for the timestep

#ifdef LHYDRO
  output real AC_maxadvec    // for all velocities - fluid and wave
#endif
#ifdef LENTROPY
  output real AC_maxchi
  #define LENERGY 1       // a hack for the moment
#endif
output real AC_maximum_error
output real AC_dt1_max

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

enum PC_SUB_STEP_NUMBER
{
	PC_FIRST_SUB_STEP,
	PC_SECOND_SUB_STEP,
	PC_THIRD_SUB_STEP,
	PC_FOURTH_SUB_STEP,
	PC_FIFTH_SUB_STEP,
}
