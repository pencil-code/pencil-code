/*                             gpu_astaroth_ansi.cu
                               --------------------

   Date:   8-Feb-2017
   Author: M. Rheinhardt & J. Pekkilae
   Description:
           ANSI C and standard library callable function wrappers for ASTAROTH-nucleus functions to be called from Fortran.
  Comments: 
*/

//General
//#include "common/errorhandler.h"
//#include "common/forcing.h"

//PC interface
#include "submodule/include/PC_moduleflags.h"
#include "../../cparam_c.h"
#include "../../cdata_c.h"
#include "../../sub_c.h"
#include "diagnostics/diagnostics.h"
#define EXTERN extern
#include "submodule/include/PC_module_parfuncs.h"
#include "submodule/include/astaroth.h"

static real* halo_buffer;
static AcMesh mesh;
static ForcingParams forcing_params;

/***********************************************************************************************/
real max_advec()
{
        real uxmax, uymax, uzmax, maxadvec_=0.;
        return maxadvec_;
}
/***********************************************************************************************/
real max_diffus()
{
        real maxdiffus_=0.;
#ifdef VISCOSITY
        maxdiffus_=nu*dxyz_2[nghost];
#endif
#ifdef MAGNETIC
        maxdiffus_=max(maxdiffus_,eta*dxyz_2[nghost]);}
#endif
        return maxdiffus_;
}
/***********************************************************************************************/
//Do the 'isubstep'th integration step on all GPUs on the node and handle boundaries
//TODO: note, isubstep starts from 0 on the GPU side (i.e. rk3 does substeps 0, 1, 2)
//
extern "C" void substepGPU(int isubstep, bool full=false)
{
#ifdef FORCING
    //Update forcing params
    if (isubstep == itorder) { 
        forcing_params.Update();    	                // calculate on CPU
        GPUUpdateForcingCoefs(&forcing_params);		// load into GPU
    }
#endif  
    if (lfirst && ldt) {
         real dt1_advec  = max_advec()/cdt;
         real dt1_diffus = max_diffus()/cdtv;
         real dt1_=sqrt(pow(dt1_advec,2) + pow(dt1_diffus,2));
         set_dt(dt1_);
    }
    //Transfer the updated ghost zone to the device(s) in the node 

    if (full)
    {
        acLoad(mesh);
    }

    if (ldiagnos) timeseries_diagnostics(h_grid);

    //Integrate on the GPUs in this node
    //NOTE: In Astaroth, isubstep is {0,1,2}, in PC it is {1,2,3}

    int3 start=(int3){l1i,m1i,n1i}, end=(int3){l2i,m2i,n2i};
    acIntegrateStepBlocked(isubstep-1, dt,start,end);

    acSynchronize();
// call for update_ghosts!
    if (!full)
    {
        //!!!GPULoadOuterHalos(&h_grid,halo_buffer);
        acLoad(mesh);
    }
    start=(int3){l1,m1,n1}; end=(int3){l2,m2,n1i-1};     // front plate
    acIntegrateStepBlocked(isubstep-1, dt,start,end);

    start=(int3){l1,m1,n2i+1}; end=(int3){l2,m2,n2};     // back plate
    acIntegrateStepBlocked(isubstep-1, dt,start,end);

    start=(int3){l1,m1,n1i}; end=(int3){l2,m1i-1,n2i};   // bottom plate
    acIntegrateStepBlocked(isubstep-1, dt,start,end);

    start=(int3){l1,m2i+1,n1i}; end=(int3){l2,m2,n2i};   // top plate
    acIntegrateStepBlocked(isubstep-1, dt,start,end);

    start=(int3){l1,m1i,n1i}; end=(int3){l1i-1,m2i,n2i};   // left plate
    acIntegrateStepBlocked(isubstep-1, dt,start,end);

    start=(int3){l2i+1,m1i,n1i}; end=(int3){l2,m2i,n2i};   // right plate
    acIntegrateStepBlocked(isubstep-1, dt,start,end);

    acSynchronize();

    //!!!GPUStoreInternalHalos(&h_grid,halo_buffer);
    acStore(mesh);
}
/* ---------------------------------------------------------------------- */
extern "C" void registerGPU(AcMesh& mesh, real* farray)
{
    long long offset=0;

    for (int i = 0; i < NUM_VTXBUF_HANDLES; ++i) {
        mesh.vertex_buffer[VertexBufferHandle(i)] = farray+offset;
        offset+=mw;
    }
}
/* ---------------------------------------------------------------------- */
//Setup the GPUs in the node to be ready for computation
extern "C" void initGPU()
{
    //Initialize GPUs in the node
    AcResult res=acGetDevice();
}
extern "C" void initializeGPU()
{
    //Setup configurations used for initializing and running the GPU code
#ifdef DOUBLE_PRECISION
printf("DOUBLE PRECISION!\n");
#endif
    	SetupConfig(mesh.config);
        AcResult res=acInitialize(mesh.config);         //Allocs memory on the GPU and loads device constants

    // initialize diagnostics
       init_diagnostics();

       //cparams.print();
}
void SetupConfig(acMeshInfo & config){

     config.int_params[AC_nx]=nx;
     config.int_params[AC_ny]=ny;
     config.int_params[AC_nz]=nz;

     config.real_params[AC_dsx]=dx;
     config.real_params[AC_dsy]=dy;
     config.real_params[AC_dsz]=dz;

     config.real_params[AC_xlen]=Lxyz[0];
     config.real_params[AC_ylen]=Lxyz[1];
     config.real_params[AC_zlen]=Lxyz[2];
     config.real_params[AC_xorig]=xyz0[0];
     config.real_params[AC_yorig]=xyz0[1];
     config.real_params[AC_zorig]=xyz0[2];
     config.real_params[AC_unit_density]=unit_density;
     config.real_params[AC_unit_velocity]=unit_velocity;
     config.real_params[AC_unit_length]=unit_length;
#ifdef LVISCOSITY
     config.real_params[AC_nu_visc]=nu;
     config.real_params[AC_zeta]=zeta;
#endif
     config.real_params[AC_cs_sound]=sqrt(cs20);
#ifdef LMAGNETIC
     config.real_params[AC_eta]=eta;
#endif;
     config.real_params[AC_mu0]=mu0;
     config.real_params[AC_gamma]=gamma;
     config.real_params[AC_cv_sound]=cv;
     config.real_params[AC_cp_sound]=cp;
     //config.real_params[AC_lnT0]=;
     config.real_params[AC_lnrho0]=lnrho0;
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
    //Deallocate everything on the GPUs and reset
       AcResult res=acQuit();
}
/***********************************************************************************************/
