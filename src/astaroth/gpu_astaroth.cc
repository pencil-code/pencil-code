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
#include <stdio.h>
//using namespace std;

#define CUDA_ERRCHK(X)
#include "submodule/src/core/math_utils.h"
#include "submodule/include/astaroth.h"

//PC interface
#define real AcReal
#define EXTERN 
#define FINT int

#include "/homeappl/home/mreinhar/git/pencil-code/samples/gputest/src/cparam_c.h"
#include "/homeappl/home/mreinhar/git/pencil-code/samples/gputest/src/cdata_c.h"

#include "../sub_c.h"                   // provides set_dt
//#include "../boundcond_c.h"             // provides boundconds[xyz] etc.
#include "../mpicomm_c.h"               // provides finalize_sendrecv_bdry
//#include "diagnostics/diagnostics.h"
#include "loadStore.h"
#if LFORCING
  #include "forcing.h"
#endif

#include "PC_module_parfuncs.h"

static AcMesh mesh;
Node node;
DeviceConfiguration devConfig;
int halo_yz_size=0;
#if LFORCING
static ForcingParams forcing_params;
#endif
/***********************************************************************************************/
AcReal max_advec()
{
        AcReal maxadvec_=0.;
#if LHYDRO
        AcReal umax=acReduceVec(RTYPE_MAX,VTXBUF_UUX,VTXBUF_UUY,VTXBUF_UUZ);
#endif
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
        maxdiffus_=std::max(maxdiffus_,eta*dxyz_2[nghost]);
#endif
        return maxdiffus_;
}
/***********************************************************************************************/
//Do the 'isubstep'th integration step on all GPUs on the node and handle boundaries
//
extern "C" void substepGPU(int isubstep, bool full=false, bool early_finalize=false)
{
#if LFORCING
    //Update forcing params
    if (isubstep == itorder) { 
         forcing_params.Update();  // calculate on CPU and load into GPU
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
    //if (1) {
    // MPI communication has already finished, hence the full domain can be advanced.
      if (!full)
      {
          loadOuterHalos(mesh);
          //acLoad(mesh);
      }
      acSynchronizeMesh();
      acIntegrateStep(isubstep-1, dt);

    } else {
    // MPI communication has not yet finished, hence only the inner domain can be advanced.

      int3 start=(int3){l1i+2,m1i+2,n1i+2}-1, end=(int3){l2i-2,m2i-2,n2i-2}-1;   // -1 shift because of C indexing convention
      acIntegrateStepWithOffset(isubstep-1,dt,start,end);
      acSynchronize();
 
      int iarg1=1, iarg2=NUM_VTXBUF_HANDLES; 
      printf("mesh.vertex_buffer, iargs= %p %d %d \n",mesh.vertex_buffer[0], iarg1, iarg2);
      finalize_isendrcv_bdry((AcReal*) mesh.vertex_buffer[0], &iarg1, &iarg2);

      loadOuterFront(mesh,STREAM_0);
      start=(int3){l1i+2,m1i+2,n1}; end=(int3){l2i-2,m2i-2,n1i+1};     // integrate inner front plate
      acDeviceIntegrateSubstep(devConfig.devices[0], STREAM_0, isubstep-1, start, end, dt);
 
      loadOuterBack(mesh,STREAM_1);
      start=(int3){l1i+2,m1i+2,n2i-1}; end=(int3){l2i-2,m2i-2,n2};     // integrate inner back plate
      acDeviceIntegrateSubstep(devConfig.devices[0], STREAM_1, isubstep-1, start, end, dt);
  
      loadOuterBot(mesh,STREAM_2);
      start=(int3){l1i+2,m1,n1}; end=(int3){l2i-2,m1i+1,n2};   // integrate inner bottom plate
      acDeviceIntegrateSubstep(devConfig.devices[0], STREAM_2, isubstep-1, start, end, dt);
  
      loadOuterTop(mesh,STREAM_3);
      start=(int3){l1i+2,m2i-1,n1}; end=(int3){l2i-2,m2,n2};   // integrate inner top plate
      acDeviceIntegrateSubstep(devConfig.devices[0], STREAM_3, isubstep-1, start, end, dt);
  
      loadOuterLeft(mesh,STREAM_4);
      start=(int3){l1,m1,n1}; end=(int3){l1i+1,m2,n2};   // integrate inner left plate
      acDeviceIntegrateSubstep(devConfig.devices[0], STREAM_4, isubstep-1, start, end, dt);
  
      loadOuterRight(mesh,STREAM_5);
      start=(int3){l2i-1,m1,n1}; end=(int3){l2,m2,n2};   // integrate inner right plate
      acDeviceIntegrateSubstep(devConfig.devices[0], STREAM_5, isubstep-1, start, end, dt);
  
      acSynchronize();
    }

    //storeInnerHalos(mesh);
    acStore(&mesh);
}
/***********************************************************************************************/
extern "C" void registerGPU(AcReal* farray)
{
    size_t offset=0;

    for (int i = 0; i < NUM_VTXBUF_HANDLES; ++i) {
        //mesh.vertex_buffer[VertexBufferHandle(i)] = (AcReal*)farray+offset;
        mesh.vertex_buffer[VertexBufferHandle(i)] = &farray[offset];
//printf("&farray[offset]= %d \n", (size_t)&farray[offset]);
        offset+=mw;
    }
}
/***********************************************************************************************/
extern "C" void initGPU()
{
    //Initialize GPUs in the node
    AcResult res=acCheckDeviceAvailability();
}
/***********************************************************************************************/
void setupConfig(AcMeshInfo & config){

     config.int_params[AC_nx]=nx;
     config.int_params[AC_ny]=ny;
     config.int_params[AC_nz]=nz;

     config.int_params[AC_mx] = mx;
     config.int_params[AC_my] = my;
     config.int_params[AC_mz] = mz;
     config.int_params[AC_nx_min] = l1;
     config.int_params[AC_nx_max] = l2;
     config.int_params[AC_ny_min] = m1;
     config.int_params[AC_ny_max] = m2;
     config.int_params[AC_nz_min] = n1;
     config.int_params[AC_nz_max] = n2;
     config.int_params[AC_mxy]  = mx*my;
     config.int_params[AC_nxy]  = nx*ny;
     config.int_params[AC_nxyz] = nw;
     config.int_params[AC_yz_plate_bufsize] = halo_yz_size;

     config.real_params[AC_dsx]=dx;
     config.real_params[AC_dsy]=dy;
     config.real_params[AC_dsz]=dz;
printf("nx etc. %d %d %d %f %f %f \n",nx,ny,nz,dx,dy,dz);
printf("dxmin, dxmax= %f %f \n", dxmin,dxmax); //it, isubstep);
     config.real_params[AC_inv_dsx] = 1./dx;
     config.real_params[AC_inv_dsy] = 1./dy;
     config.real_params[AC_inv_dsz] = 1./dz;
     config.real_params[AC_dsmin]   = std::min(dx,std::min(dy,dz));
     config.real_params[AC_xlen]=lxyz[0];
     config.real_params[AC_ylen]=lxyz[1];
     config.real_params[AC_zlen]=lxyz[2];
     config.real_params[AC_xorig]=xyz0[0];
     config.real_params[AC_yorig]=xyz0[1];
     config.real_params[AC_zorig]=xyz0[2];
printf("lxyz etc. %f %f %f %f %f %f \n",lxyz[0],lxyz[1],lxyz[2],xyz0[0],xyz0[1],xyz0[2]);
     config.real_params[AC_unit_density]=unit_density;
     config.real_params[AC_unit_velocity]=unit_velocity;
     config.real_params[AC_unit_length]=unit_length;
printf("units etc. %lf %lf %lf \n", unit_density, unit_velocity, unit_length);

#include "PC_modulepars.h"
#if LVISCOSITY
#include "../viscosity_pars_c.h"
     config.real_params[AC_nu_visc]=nu;
     config.real_params[AC_zeta]=zeta;
printf("nu etc. %f %f \n", nu, zeta);
#endif
#if LMAGNETIC
#include "../magnetic_pars_c.h"
     config.real_params[AC_eta]=eta;
printf("eta etc. %f \n", eta);
#endif
     config.real_params[AC_mu0]=mu0;
#include "../equationofstate_pars_c.h"
     config.real_params[AC_cs_sound]=sqrt(cs20);
     config.real_params[AC_gamma]=gamma;
     config.real_params[AC_cv_sound]=cv;
     config.real_params[AC_cp_sound]=cp;
     config.real_params[AC_lnT0]=lntt0;
     config.real_params[AC_lnrho0]=lnrho0;
printf("eos etc. %f %f %f %f %f %f \n",sqrt(cs20),gamma,cv,cp,lntt0,lnrho0);
#if LENTROPY
#include "../energy_pars_c.h"
     config.real_params[AC_chi]=chi;
printf("chi %f \n", chi);
#endif
}

#define PUT(ptr,n_elem) \
  acDeviceLoadScalarArray(devConfig.devices[0], STREAM_DEFAULT, AC_##ptr, 0, ptr, n_elem);

void loadProfiles(AcMeshInfo & config){
#if LFORCING
#include "../forcing_pars_c.h"
     config.int_params[AC_iforcing_zsym]=iforcing_zsym;
     config.real_params[AC_k1_ff]=k1_ff;
     printf("profx_ampl, val= %d %lf %lf\n", profx_ampl, profx_ampl[0], profx_ampl[nx-1]);
     PUT(profx_ampl,nx)
     PUT(profy_ampl,my)
     PUT(profz_ampl,mz)
     PUT(profx_hel,nx)
     PUT(profy_hel,my)
     PUT(profz_hel,mz)
#endif
}
extern "C" void initializeGPU()
{
    //Setup configurations used for initializing and running the GPU code
        initLoadStore();
    	setupConfig(mesh.info);
        AcResult res=acInit(mesh.info);
        res=acGetNode(&node);
        acNodeQueryDeviceConfiguration(node, &devConfig);
        loadProfiles(mesh.info);

    // initialize diagnostics
       //init_diagnostics();
}
/***********************************************************************************************/
extern "C" void copyFarray() 
{
       AcResult res=acStore(&mesh);
}
/***********************************************************************************************/
extern "C" void finalizeGPU()
{
       finalLoadStore();

    //Deallocate everything on the GPUs and reset
       AcResult res=acQuit();
}
/***********************************************************************************************/
