// Provides all declarations and functions needed for the formulation of the PDEs' rhss by DSL code
// and finally for the definition of the solve kernel.

#include "PC_nghost.h"
#define nghost NGHOST

// declare here reduction results needed for the timestep
output real AC_maxadvec    // for all velocities - fluid and wave

#include "fieldecs.h"
#include "../stdlib/operators.h"
#include "../stdlib/integrators.h"
#include "../stdlib/units.h"
#include "../stdlib/utils.h"
#include "../../PC_moduleflags.h"
#include "../phys_consts.h"
#include "PC_modulepardecs.h"
#ifdef LFORCING
  #include "../forcing/pcstyleforcing.h"
#endif
#ifdef LDENSITY
  #define LNRHO RHO
#endif
#include "equations.h"