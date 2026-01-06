// Provides all declarations and functions needed for the formulation of the PDEs' rhss by DSL code
// and finally for the definition of the solve kernel.
#include "../typedefs.h"
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
struct PC_gridpars
{
 	real3 dline_1
       	real dxyz_2
        real dxyz_4
        real dxyz_6
}
enum PC_SUB_STEP_NUMBER
{
	PC_FIRST_SUB_STEP,
	PC_SECOND_SUB_STEP,
	PC_THIRD_SUB_STEP,
	PC_FOURTH_SUB_STEP,
	PC_FIFTH_SUB_STEP,
}

#include "PC_nghost.h"

#define n__mod__cdata (vertexIdx.z+1)
#define m__mod__cdata (vertexIdx.y+1)
#define AC_n__mod__cdata (vertexIdx.z+1)
#define AC_m__mod__cdata (vertexIdx.y+1)

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
#include "../stdlib/grid/funcs.h"
#include "../stdlib/utils/intrinsics.h"
#ifndef LTRAINING
#include "../stdlib/slope_limited_diffusion.h"
#endif

#include "../../../cparam_c.h"
#include "../../../cparam_pencils.inc_c.h"

#define AC_n1 n1
#define AC_m1 m1
#define AC_l1 l1

#define AC_mxgrid mxgrid
#define AC_mygrid mygrid
#define AC_mzgrid mzgrid

#include "fieldecs.h"
#include "../stdlib/optimized_integrators.h"
//#include "../stdlib/units.h"
#include "../stdlib/utils/kernels.h"
//#include "../stdlib/map.h"
#include "PC_modulepardecs.h"

#include "../stdlib/general_derivs.h"
#include "../stdlib/general_operators.h"
#define AC_NGHOST NGHOST

// declare here reduction results needed for the timestep

#ifdef LHYDRO
  output real AC_maxadvec    // for all velocities - fluid and wave
#endif
#ifdef LENTROPY
  output real AC_maxdiffchi
  #define LENERGY 1          // a hack for the moment
#endif
global output real AC_maximum_error
#ifdef LVISCOSITY
  output real AC_maxdiffnu
#endif
#ifdef LMAGNETIC
  output real AC_maxdiffeta
#endif

output real AC_dt1_max
global output real AC_Arms
const int AC_xbot__mod__equationofstate=1
const int AC_xtop__mod__equationofstate=nx

#ifdef LDENSITY
  #define LNRHO RHO
#endif
#include "../bcs/funcs.h"
#include "../bcs/funcs_overload.h"

#ifdef LFORCING
  #include "../forcing/pcstyleforcing.h"
#endif

#include "../get_grid_mn.h"
#include "../steps_two_full.h"
#include "equations.h"
