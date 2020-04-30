
#include <stdio.h>
#include <math.h>
#include <algorithm>

#include "common/qualify.h"
#include "utils/utils.h"
#include "gpu/cuda/cuda_generic.cuh"

#include "../../cparam_c.h"
#include "../../cdata_c.h"
#include "../../diagnostics_c.h"
#define save_name qualified(diagnostics,save_name,MODPRE,MODIN,MODSUF)

// Contains declarations if idiag_* variables and definition of function init_diagnostics tb called by gpu_astaroth.
#define EXTERN
#include "diagnostics/PC_module_diagfuncs.h"

static real diag; 

void init_diagnostics()
{
  #include "diagnostics/PC_modulediags_init.h"
  //printf("idiag_urms %d\n", idiag_urms);
}

void timeseries_diagnostics(const Grid & h_grid)
{
  // Calculate and save all of the diagnostic variables calculated within the CUDA devices. 

  // Contains automatically generated calls to reduce_cuda_PC according to required diagnostics of the different modules.
  #include "diagnostics/PC_modulediags.h"

  if (idiag_mass){
    diag=reduce_cuda_PC(ldensity_nolog ? SUM_SCAL : SUM_EXP, h_grid.LNRHO)*box_volume/nw;
    save_name(diag,idiag_mass);
  }
/* not yet automatically generated density diagnostics calls */
                if (idiag_rhomax>0) {
                        diag=reduce_cuda_PC(MAX_SCAL, h_grid.LNRHO);
                        if (!ldensity_nolog) diag = exp(diag);       //Change away from the logarithmic form
                        save_name(diag,idiag_rhomax);
                }
                if (idiag_rhomin>0) {
                        diag=reduce_cuda_PC(MIN_SCAL, h_grid.LNRHO);
                        if (!ldensity_nolog) diag = exp(diag);       //Change away from the logarithmic form
                        diag=-diag;
                        save_name(diag,idiag_rhomin);
                }
                if (idiag_rhom){
                        diag=reduce_cuda_PC(ldensity_nolog ? SUM_SCAL : SUM_EXP, h_grid.LNRHO);
                        save_name(diag,idiag_rhom);
                }
 		if (idiag_rhorms){
                        diag=reduce_cuda_PC(ldensity_nolog ? RMS_SCAL : RMS_EXP, h_grid.LNRHO);
    	        	save_name(diag,idiag_rhorms);
		}
}
