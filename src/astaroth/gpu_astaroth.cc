/*                             gpu_astaroth_ansi.cu
                               --------------------

   Date:   8-Feb-2017
   Author: M. Rheinhardt & J. Pekkilae
   Description:
           ANSI C and standard library callable function wrappers for ASTAROTH-nucleus functions to be called from Fortran.
  Comments: 
*/

//General
#include <cmath>
#include <algorithm>
#include <limits>
using namespace std;
#include "submodule/include/astaroth.h"

//PC interface
#define real AcReal
#define CUDA_ERRCHK(ans) {}             // suppress calls not needed at the moment
#include "../sub_c.h"                   // provides set_dt
#include "../boundcond_c.h"             // provides boundconds[xyz] etc.
//#include "diagnostics/diagnostics.h"
#include "loadStore.h"
#define EXTERN 
#include "PC_module_parfuncs.h"

static AcMesh mesh;
static Node *node_handle;
//static ForcingParams forcing_params;

/***********************************************************************************************/
AcReal max_advec()
{
        AcReal uxmax, uymax, uzmax, maxadvec_=0.;
        return maxadvec_;
}
/***********************************************************************************************/
AcReal max_diffus()
{
        AcReal maxdiffus_=0.;
#if LVISCOSITY
        maxdiffus_=nu*dxyz_2[nghost];
#endif
#if LMAGNETIC
        maxdiffus_=max(maxdiffus_,eta*dxyz_2[nghost]);
#endif
        return maxdiffus_;
}
/***********************************************************************************************/
//Do the 'isubstep'th integration step on all GPUs on the node and handle boundaries
//TODO: note, isubstep starts from 0 on the GPU side (i.e. rk3 does substeps 0, 1, 2)
//
extern "C" void substepGPU(int isubstep, bool full=false, bool early_finalize=false)
{
#if LFORCING
    //Update forcing params
    if (isubstep == itorder) { 
        // forcing_params.Update();    	                // calculate on CPU // %JP: TODO not available in astaroth.h yet
        // GPUUpdateForcingCoefs(&forcing_params);		// load into GPU    // %JP: TODO not available in astaroth.h yet
    }
#endif  
    if (lfirst && ldt) {
         AcReal dt1_advec  = max_advec()/cdt;
         AcReal dt1_diffus = max_diffus()/cdtv;
         AcReal dt1_=sqrt(pow(dt1_advec,2) + pow(dt1_diffus,2));
         set_dt(dt1_);
    }
    //Transfer the updated ghost zone to the device(s) in the node 

    if (full)
    {
        acLoad(mesh);
    }

    //if (ldiagnos) timeseries_diagnostics(h_grid);

    //Integrate on the GPUs in this node
    //NOTE: In Astaroth, isubstep is {0,1,2}, in PC it is {1,2,3}

    if (early_finalize) {
    // MPI communication has already finished, hence the full domain can be advanced.
      if (!full)
      {
          //loadOuterHalos(mesh);
          acLoad(mesh);
      }
      acIntegrateStep(isubstep-1, dt);

    } else {
    // MPI communication has not yet finished, hence only the inner domain can be advanced.

      int3 start=(int3){l1i,m1i,n1i}, end=(int3){l2i,m2i,n2i};
      acIntegrateStepWithOffset(isubstep-1, dt,start,end);
  
      acSynchronize();

      finalize_isendrcv_bdry((AcReal*)mesh.vertex_buffer[0]);

  // TODO: call for update of ghosts!
      loadOuterFront(mesh);
      start=(int3){l1,m1,n1}; end=(int3){l2,m2,n1i-1};     // front plate
      acIntegrateStepWithOffset(isubstep-1, dt,start,end);
  
      loadOuterBack(mesh);
      start=(int3){l1,m1,n2i+1}; end=(int3){l2,m2,n2};     // back plate
      acIntegrateStepWithOffset(isubstep-1, dt,start,end);
  
      loadOuterBottom(mesh);
      start=(int3){l1,m1,n1i}; end=(int3){l2,m1i-1,n2i};   // bottom plate
      acIntegrateStepWithOffset(isubstep-1, dt,start,end);
  
      loadOuterTop(mesh);
      start=(int3){l1,m2i+1,n1i}; end=(int3){l2,m2,n2i};   // top plate
      acIntegrateStepWithOffset(isubstep-1, dt,start,end);
  
      loadOuterLeft(mesh);
      start=(int3){l1,m1i,n1i}; end=(int3){l1i-1,m2i,n2i};   // left plate
      acIntegrateStepWithOffset(isubstep-1, dt,start,end);
  
      loadOuterRight(mesh);
      start=(int3){l2i+1,m1i,n1i}; end=(int3){l2,m2i,n2i};   // right plate
      acIntegrateStepWithOffset(isubstep-1, dt,start,end);
  
      acSynchronize();
    }

    //storeInnerHalos(mesh);
    acStore(&mesh);
}
/* ---------------------------------------------------------------------- */
extern "C" void registerGPU(AcMesh& mesh, AcReal* farray)
{
    long long offset=0;
//printf("farray-address = %d \n", (long)(farray));
    for (int i = 0; i < NUM_VTXBUF_HANDLES; ++i) {
//printf("farray-address offset,i= %d %d \n", offset, (long)((AcReal*)farray+offset));
        mesh.vertex_buffer[VertexBufferHandle(i)] = (AcReal*)farray+offset;
        offset+=mw;
    }
}
/* ---------------------------------------------------------------------- */
//Setup the GPUs in the node to be ready for computation
extern "C" void initGPU()
{
    //Initialize GPUs in the node
    AcResult res=acCheckDeviceAvailability();
}
void setupConfig(AcMeshInfo & config){
     config.int_params[AC_nx]=nx;
     config.int_params[AC_ny]=ny;
     config.int_params[AC_nz]=nz;

     config.real_params[AC_dsx]=dx;
     config.real_params[AC_dsy]=dy;
     config.real_params[AC_dsz]=dz;
//printf("nx etc. %d %d %d %f %f %f \n",nx,ny,nz,dx,dy,dz);
     config.real_params[AC_xlen]=lxyz[0];
     config.real_params[AC_ylen]=lxyz[1];
     config.real_params[AC_zlen]=lxyz[2];
     config.real_params[AC_xorig]=xyz0[0];
     config.real_params[AC_yorig]=xyz0[1];
     config.real_params[AC_zorig]=xyz0[2];
//printf("lxyz etc. %f %f %f %f %f %f \n",lxyz[0],lxyz[1],lxyz[2],xyz0[0],xyz0[1],xyz0[2]);
     config.real_params[AC_unit_density]=unit_density;
     config.real_params[AC_unit_velocity]=unit_velocity;
     config.real_params[AC_unit_length]=unit_length;
//printf("units etc. %f %f %f \n", unit_density, unit_velocity, unit_length);
#if LVISCOSITY
#include "../viscosity_pars_c.h"
     config.real_params[AC_nu_visc]=nu;
     config.real_params[AC_zeta]=zeta;
//printf("nu etc. %f %f \n", nu, zeta);
#endif
#if LMAGNETIC
     config.real_params[AC_eta]=eta;
#endif
     config.real_params[AC_mu0]=mu0;
#include "../equationofstate_pars_c.h"
     config.real_params[AC_cs_sound]=sqrt(cs20);
     config.real_params[AC_gamma]=gamma;
     config.real_params[AC_cv_sound]=cv;
     config.real_params[AC_cp_sound]=cp;
     config.real_params[AC_lnT0]=lnt0;
     config.real_params[AC_lnrho0]=lnrho0;
//printf("eos etc. %f %f %f %f %f %f \n",sqrt(cs20),gamma,cv,cp,lnt0,lnrho0);
}
extern "C" void initializeGPU()
{
    //Setup configurations used for initializing and running the GPU code

    	setupConfig(mesh.info);
        AcResult res=acInit(mesh.info);         //Allocs memory on the GPU and loads device constants
        res=acNodeCreate(iproc, mesh.info, node_handle);
        initLoadStore();
    // initialize diagnostics
       //init_diagnostics();
}
/***********************************************************************************************/
extern "C" void copyFarray() 
{
       AcResult res=acStore(&mesh);
}
/***********************************************************************************************/
//Destroy the GPUs in the node (not literally hehe)
extern "C" void finalizeGPU()
{
       finalLoadStore();
    //Deallocate everything on the GPUs and reset
       AcResult res=acQuit();
}
/***********************************************************************************************/
