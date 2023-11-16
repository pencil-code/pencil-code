/*
                               --------------------

   Date:   8-Feb-2017
   Author: M. Rheinhardt & J. Pekkilae
   Description:
           ANSI C and standard library callable function wrappers for ASTAROTH-nucleus functions to be called from Fortran.
  Comments: 
*/

//General
//#include <cmath>
#include <algorithm>
#include <chrono>
#include <iostream>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <mpi.h>

#define CUDA_ERRCHK(X)

#include "submodule/acc-runtime/api/errchk.h"
#include "submodule/acc-runtime/api/math_utils.h"
#include "submodule/include/astaroth.h"
#include "submodule/build/acc-runtime/api/user_defines.h"
#include "submodule/src/core/kernels/kernels.h"
#include "submodule/src/core/task.h"
#include "submodule/acc-runtime/api/math_utils.h"
#include "submodule/include/astaroth_utils.h"
#define real AcReal
#define EXTERN 
#define FINT int
__thread int tp_int;

//PC interface
#include "PC_moduleflags.h"
#include "../cparam_c.h"
#include "../cdata_c.h"
#include "../sub_c.h"                   // provides set_dt
#include "../boundcond_c.h"             // provides boundconds[xyz] etc.
#include "../mpicomm_c.h"               // provides finalize_sendrcv_bdry
//#include "diagnostics/diagnostics.h"
#if PACKED_DATA_TRANSFERS
  //#include "loadStore.h"
#endif
#if LFORCING
//   #include "forcing.h"
#endif

#include "PC_module_parfuncs.h"

static bool with_boundconds = true;
static AcMesh mesh;
static AcMesh test_mesh;
static AcTaskGraph* graph_1;
static AcTaskGraph* graph_2;
static AcTaskGraph* graph_3;
static AcTaskGraph* randomize_graph;
static int pid;
static int rank;
DeviceConfiguration devConfig;
int halo_xz_size[2]={0,0}, halo_yz_size[2]={0,0};
static int l0=1;

static AcReal* xtop_buffer;
static AcReal* xbot_buffer;

static AcReal* ytop_buffer;
static AcReal* ybot_buffer;
#if LFORCING
// static ForcingParams forcing_params;
// #include "../forcing_c.h"
#endif
int
DCONST(const AcIntParam param)
{
  return mesh.info.int_params[param];
}
int3
DCONST(const AcInt3Param param)
{
  return mesh.info.int3_params[param];
}
AcReal
DCONST(const AcRealParam param)
{
  return mesh.info.real_params[param];
}
AcReal3
DCONST(const AcReal3Param param)
{
  return mesh.info.real3_params[param];
}
int DEVICE_VTXBUF_IDX(const int x, const int y, const int z){
    return x+ mx*(y+my*(z));
}
/***********************************************************************************************/
void boundcond_y_interfaced(AcReal *f, int *ivar1, int *ivar2, bool *finished){
    printf("Hi :). Successful call from taskgraph to boundcond_interfaced \n");
    boundconds_y_c(f,ivar1,ivar2);
    *finished = true;

} 
/***********************************************************************************************/
void boundcond_z_interfaced(AcReal *f, int *ivar1, int *ivar2, bool *finished){
    printf("Hi :). Successful call from taskgraph to boundcond_z_interfaced \n");
    boundconds_z_c(f,ivar1,ivar2);
    *finished = true;

} 
/***********************************************************************************************/
void boundcond_x_interfaced(AcReal *f, int *ivar1, int *ivar2, bool *finished){
    printf("Hi :). Successful call from taskgraph to boundcond_x_interfaced \n");
    boundconds_x_c(f,ivar1,ivar2);
    *finished = true;

} 
/***********************************************************************************************/
void
print_diagnostics(const int pid, const int step, const AcReal dt, const AcReal simulation_time,
                  FILE* diag_file, const AcReal sink_mass, const AcReal accreted_mass,
                  int* found_nan)
{

    AcReal buf_rms, buf_max, buf_min;
    const int max_name_width = 16;

    // Calculate rms, min and max from the velocity vector field
    acGridReduceVec(STREAM_DEFAULT, RTYPE_MAX, VTXBUF_UUX, VTXBUF_UUY, VTXBUF_UUZ, &buf_max);
    acGridReduceVec(STREAM_DEFAULT, RTYPE_MIN, VTXBUF_UUX, VTXBUF_UUY, VTXBUF_UUZ, &buf_min);
    acGridReduceVec(STREAM_DEFAULT, RTYPE_RMS, VTXBUF_UUX, VTXBUF_UUY, VTXBUF_UUZ, &buf_rms);

    if(pid==0)(pid, "Step %d, t_step %.3e, dt %e s\n", step, double(simulation_time),
                      double(dt));
    if(pid==0)(pid, "  %*s: min %.3e,\trms %.3e,\tmax %.3e\n", max_name_width, "uu total",
                      double(buf_min), double(buf_rms), double(buf_max));
    if (pid == 0) {
        fprintf(diag_file, "%d %e %e %e %e %e ", step, double(simulation_time), double(dt),
                double(buf_min), double(buf_rms), double(buf_max));
    }

#if LBFIELD
    acGridReduceVec(STREAM_DEFAULT, RTYPE_MAX, BFIELDX, BFIELDY, BFIELDZ, &buf_max);
    acGridReduceVec(STREAM_DEFAULT, RTYPE_MIN, BFIELDX, BFIELDY, BFIELDZ, &buf_min);
    acGridReduceVec(STREAM_DEFAULT, RTYPE_RMS, BFIELDX, BFIELDY, BFIELDZ, &buf_rms);

    acLogFromRootProc(pid, "  %*s: min %.3e,\trms %.3e,\tmax %.3e\n", max_name_width, "bb total",
                      double(buf_min), double(buf_rms), double(buf_max));
    if (pid == 0) {
        fprintf(diag_file, "%e %e %e ", double(buf_min), double(buf_rms), double(buf_max));
    }

    acGridReduceVecScal(STREAM_DEFAULT, RTYPE_ALFVEN_MAX, BFIELDX, BFIELDY, BFIELDZ, VTXBUF_LNRHO,
                        &buf_max);
    acGridReduceVecScal(STREAM_DEFAULT, RTYPE_ALFVEN_MIN, BFIELDX, BFIELDY, BFIELDZ, VTXBUF_LNRHO,
                        &buf_min);
    acGridReduceVecScal(STREAM_DEFAULT, RTYPE_ALFVEN_RMS, BFIELDX, BFIELDY, BFIELDZ, VTXBUF_LNRHO,
                        &buf_rms);

    acLogFromRootProc(pid, "  %*s: min %.3e,\trms %.3e,\tmax %.3e\n", max_name_width, "vA total",
                      double(buf_min), double(buf_rms), double(buf_max));
    if (pid == 0) {
        fprintf(diag_file, "%e %e %e ", double(buf_min), double(buf_rms), double(buf_max));
    }
#endif

    // Calculate rms, min and max from the variables as scalars
    for (int i = 0; i < NUM_VTXBUF_HANDLES; ++i) {
        acGridReduceScal(STREAM_DEFAULT, RTYPE_MAX, VertexBufferHandle(i), &buf_max);
        acGridReduceScal(STREAM_DEFAULT, RTYPE_MIN, VertexBufferHandle(i), &buf_min);
        acGridReduceScal(STREAM_DEFAULT, RTYPE_RMS, VertexBufferHandle(i), &buf_rms);

        acLogFromRootProc(pid, "  %*s: min %.3e,\trms %.3e,\tmax %.3e\n", max_name_width,
                          vtxbuf_names[i], double(buf_min), double(buf_rms), double(buf_max));
        if (pid == 0) {
            fprintf(diag_file, "%e %e %e ", double(buf_min), double(buf_rms), double(buf_max));
        }

        if (isnan(buf_max) || isnan(buf_min) || isnan(buf_rms)) {
            *found_nan = 1;
        }
    }

    if ((sink_mass >= AcReal(0.0)) || (accreted_mass >= AcReal(0.0))) {
        if (pid == 0) {
            fprintf(diag_file, "%e %e ", double(sink_mass), double(accreted_mass));
        }
    }

    if (pid == 0) {
        fprintf(diag_file, "\n");
    }

#if LSINK
    acLogFromRootProc(pid, "sink mass is: %.15e \n", double(sink_mass));
    acLogFromRootProc(pid, "accreted mass is: %.15e \n", double(accreted_mass));
#endif

    fflush(diag_file);
    fflush(stdout);
}
/***********************************************************************************************/
AcReal max_advec()
{
        AcReal maxadvec_=0.;
#if LHYDRO
	AcReal umax;
        acGridReduceVec(STREAM_DEFAULT,RTYPE_MAX,VTXBUF_UUX,VTXBUF_UUY,VTXBUF_UUZ,&umax);
#endif
        return maxadvec_;
}
/***********************************************************************************************/
AcReal max_diffus()
{
        AcReal maxdiffus_=0.;
#if LVISCOSITY
        maxdiffus_=nu*dxyz_2[nghost-1];
#endif
#if LMAGNETIC
        maxdiffus_=std::max(maxdiffus_,eta*dxyz_2[nghost-1]);
#endif
        return maxdiffus_;
}
int
id_to_tag(int3 id)
{
    return ((3 + id.x) % 3) * 9 + ((3 + id.y) % 3) * 3 + (3 + id.z) % 3;
}
int3
tag_to_id(int _tag)
{
    int3 _id = (int3){(_tag) / 9, ((_tag) % 9) / 3, (_tag) % 3};
    _id.x    = _id.x == 2 ? -1 : _id.x;
    _id.y    = _id.y == 2 ? -1 : _id.y;
    _id.z    = _id.z == 2 ? -1 : _id.z;
    ERRCHK_ALWAYS(id_to_tag(_id) == _tag);
    return _id;
}
void loadBoundary()
{
    Device device = acGridGetDevice();
    int3 nn = (int3){
        device->local_config.int_params[AC_nx],
        device->local_config.int_params[AC_ny],
        device->local_config.int_params[AC_nz],
    };
    int min_halo_tag = 1;
    int max_halo_tag = 27;
    uint3_64 decomp = {1,2,1};
    for (int tag = min_halo_tag; tag < max_halo_tag; tag++) {
        Region region = Region(RegionFamily::Exchange_output, tag, {mx,my,mz},{},0);
        if (Region::is_on_boundary(decomp, rank, tag, BOUNDARY_Z_BOT)) {
            for(int vtxbuf=0;vtxbuf<NUM_VTXBUF_HANDLES;vtxbuf++) acDeviceLoadVertexBufferWithOffset(device,STREAM_2,mesh,(VertexBufferHandle) vtxbuf,{0,0,0},{0,0,0},mx*my*NGHOST);
            // acDeviceLoadMeshWithOffset(device,STREAM_2,mesh, region.position,region.position,region.dims.x*region.dims.y*region.dims.z);
        }
        else if(Region::is_on_boundary(decomp, rank, tag, BOUNDARY_Z_TOP)){
            for(int vtxbuf=0;vtxbuf<NUM_VTXBUF_HANDLES;vtxbuf++) acDeviceLoadVertexBufferWithOffset(device,STREAM_2,mesh,(VertexBufferHandle) vtxbuf,{0,0,0},{0,0,0},mx*my*NGHOST);
        }
        else if(Region::is_on_boundary(decomp, rank, tag, BOUNDARY_Y_BOT)){
            for(int vtxbuf=0;vtxbuf<NUM_VTXBUF_HANDLES;vtxbuf++) acDeviceLoadVertexBufferWithOffset(device,STREAM_2,mesh,(VertexBufferHandle) vtxbuf,{0,0,0},{0,0,0},mx*my*NGHOST);
        }
        else if(Region::is_on_boundary(decomp, rank, tag, BOUNDARY_Y_TOP)){
            for(int vtxbuf=0;vtxbuf<NUM_VTXBUF_HANDLES;vtxbuf++) acDeviceLoadVertexBufferWithOffset(device,STREAM_2,mesh,(VertexBufferHandle) vtxbuf,{0,0,0},{0,0,0},mx*my*NGHOST);
        }
        else if(Region::is_on_boundary(decomp, rank, tag, BOUNDARY_X_TOP)){
            for(int vtxbuf=0;vtxbuf<NUM_VTXBUF_HANDLES;vtxbuf++) acDeviceLoadVertexBufferWithOffset(device,STREAM_2,mesh,(VertexBufferHandle) vtxbuf,{0,0,0},{0,0,0},mx*my*NGHOST);
        }
        else if(Region::is_on_boundary(decomp, rank, tag, BOUNDARY_X_BOT)){
            for(int vtxbuf=0;vtxbuf<NUM_VTXBUF_HANDLES;vtxbuf++) acDeviceLoadVertexBufferWithOffset(device,STREAM_2,mesh,(VertexBufferHandle) vtxbuf,{0,0,0},{0,0,0},mx*my*NGHOST);
        }
    }
    acGridSynchronizeStream(STREAM_ALL);
}
/***********************************************************************************************/
//Do the 'isubstep'th integration step on all GPUs on the node and handle boundaries.
//
extern "C" void substepGPU(int isubstep, bool full=false, bool early_finalize=false)
{
#if LFORCING
    // //Update forcing params

    // if (isubstep == itorder) 
    //      forcing_params.Update();  // calculate on CPU and load into GPU
#endif
    if (lfirst && ldt) {
         AcReal dt1_advec  = max_advec()/cdt;
         AcReal dt1_diffus = max_diffus()/cdtv;
         AcReal dt1_=sqrt(pow(dt1_advec,2) + pow(dt1_diffus,2));
         set_dt(dt1_);
    }
    acGridLoadScalarUniform(STREAM_DEFAULT, AC_dt, dt);
    acGridSynchronizeStream(STREAM_DEFAULT);
    //Transfer the updated ghost zone to the device(s) in the node 

    if (full){
        //normalize
        for(int i=0;i<mx*my*mz;i++){
          for(int vtxbufidx = 0; vtxbufidx< NUM_VTXBUF_HANDLES; vtxbufidx++)
          {
            mesh.vertex_buffer[vtxbufidx][i] = mesh.vertex_buffer[vtxbufidx][i]*0.1;
          }
        }
        acDeviceLoadMesh(acGridGetDevice(),STREAM_DEFAULT,mesh);
        // acGridLoadMesh(STREAM_DEFAULT,mesh);
        acGridSynchronizeStream(STREAM_ALL);
        // acGridExecuteTaskGraph(randomize_graph,1);
        acGridSynchronizeStream(STREAM_ALL);
	      //example of loading profile
        //acDeviceLoadProfile(acGridGetDevice(), STREAM_DEFAULT, mesh,PROFILE_X);
    }
    if(!early_finalize){

    }
        acGridSynchronizeStream(STREAM_ALL);
    // printf("TP: done randomize\n");
    // fflush(stdout);
    if (full){
        // acGridLoadMesh(STREAM_DEFAULT,mesh);
        // acGridLoadProfile(STREAM_DEFAULT, PROFILE_X, mesh);
    }
    acGridSynchronizeStream(STREAM_ALL);
    if(!with_boundconds){
        loadBoundary();
    }
    acGridLoadScalarUniform(STREAM_DEFAULT, AC_dt, dt);
    if(isubstep == 1){
        acGridExecuteTaskGraph(graph_1,1);
    }
    if(isubstep == 2){
        acGridExecuteTaskGraph(graph_2,1);
    }
    if(isubstep == 3){
        acGridExecuteTaskGraph(graph_3,1);
    }
    acGridSynchronizeStream(STREAM_ALL);
    if(!with_boundconds){
        loadBoundary();
    }
    acGridSynchronizeStream(STREAM_ALL);
    // printf("Done substep: %d\n",isubstep);
    // fflush(stdout);
    // int found_nan;
    // FILE* diag_file = fopen("astaroth_timeseries.ts", "a");
    // int rank;
    // print_diagnostics(rank, 1, dt, dt,diag_file, 0.0001, 0.0001,
    //               &found_nan);
    return;
}
extern "C" void testBcKernel(AcReal* farray_in, AcReal* farray_truth){
    AcMesh mesh_true;
    AcMesh mesh_test;
    AcReal epsilon = 0.00001;

    size_t offset=0;
    for (int i = 0; i < NUM_VTXBUF_HANDLES; ++i) {
        mesh_test.vertex_buffer[VertexBufferHandle(i)] = &farray_in[offset];
        offset+=mw;
    }
    offset = 0;
    for (int i = 0; i < NUM_VTXBUF_HANDLES; ++i) {
        mesh_true.vertex_buffer[VertexBufferHandle(i)] = &farray_truth[offset];
        offset+=mw;
    }

    // AcReal cv1 = 1.66666663;
    // lnrho0 = 0;
    // cp = 1.0;
    // AcReal gamma_m1 = 0.666666627;
    //Run the gpu code serially
    int3 dims = {mx,my,1};
    for(int i=0;i<dims.x;i++){
        for(int j=0;j<dims.y;j++){
            for(int k=0;k<dims.z;k++){
                //Remember to convert to zero based index :))))))
                // AcReal ftopktop = 1.0;
                // AcReal rho_xy=exp(mesh.vertex_buffer[ilnrho-1][DEVICE_VTXBUF_IDX(i,j,n2-1)]);
                // AcReal cs2_xy = mesh.vertex_buffer[iss-1][DEVICE_VTXBUF_IDX(i,j,n2-1)];
                // cs2_xy=cs20*exp(gamma_m1*(mesh.vertex_buffer[ilnrho-1][DEVICE_VTXBUF_IDX(i,j,n2-1)]-lnrho0)+cv1*cs2_xy);
                // AcReal tmp_xy=ftopktop/cs2_xy;
                // for(int z=1;z<=3;z++){
                //     rho_xy = mesh.vertex_buffer[ilnrho-1][DEVICE_VTXBUF_IDX(i,j,n2+z-1)]-mesh.vertex_buffer[ilnrho-1][DEVICE_VTXBUF_IDX(i,j,n2-z-1)];
                //     mesh.vertex_buffer[iss-1][DEVICE_VTXBUF_IDX(i,j,n2+z-1)] = mesh.vertex_buffer[iss-1][DEVICE_VTXBUF_IDX(i,j,n2-z-1)] + cp*(cp-cv)*(-rho_xy-1.0*tmp_xy);
                // }
                //#include "res.cuh"

            }
        }
    }

    bool passed = true;
    for(int i=0;i<mx;i++){
        for(int j=0;j<my;j++){
            for(int k=0;k<mz;k++){
                for(int ivar=0;ivar<mfarray;ivar++){
                    AcReal out_val = mesh_test.vertex_buffer[ivar][DEVICE_VTXBUF_IDX(i,j,k)];
                    AcReal true_val= mesh_true.vertex_buffer[ivar][DEVICE_VTXBUF_IDX(i,j,k)];
                    if(fabs(out_val - true_val) > epsilon){
                        passed = false;
                        printf("C val wrong at %d,%d,%d\n",i,j,k);
                        printf("field = %d",ivar);
                        printf("C val: %f\tF val: %f\n",out_val,true_val);
                    }
                }
            }
        }
    }
    if(passed){
        printf("Passed C test :)\n");
    } else{
        printf("Did not pass C test :(\n");
    }
    if(!passed) return;
    printf("Starting GPUtest\n");
    fflush(stdout);
    acGridSynchronizeStream(STREAM_ALL);
    // acGridLoadMesh(STREAM_DEFAULT,mesh);
    // printf("loaded mesh\n");
    // acGridTestBCKernel({mx,my,1});


    acGridSynchronizeStream(STREAM_ALL);
    printf("after bc kernel\n");
    fflush(stdout);
    acGridStoreMesh(STREAM_DEFAULT,&mesh);
    acGridSynchronizeStream(STREAM_ALL);
    printf("after store\n");
    fflush(stdout);
    for(int i=0;i<mx;i++){
        for(int j=0;j<my;j++){
            for(int k=0;k<mz;k++){
                for(int ivar=0;ivar<mfarray;ivar++){
                    AcReal out_val = mesh.vertex_buffer[ivar][DEVICE_VTXBUF_IDX(i,j,k)];
                    AcReal true_val= mesh_true.vertex_buffer[ivar][DEVICE_VTXBUF_IDX(i,j,k)];
                    if(fabs(out_val - true_val) > epsilon){
                        passed = false;
                        printf("GPU val wrong at %d,%d,%d\n",i,j,k);
                        printf("field = %d",ivar);
                        printf("GPU val: %f\tTRUE val: %f\n",out_val,true_val);
                    }
                }
            }
        }
    }
    if(passed){
        printf("Passed GPU test :)\n");
    }
    else{
        printf("Did not pass GPU test :(\n");
    }
    return;

}
/***********************************************************************************************/
extern "C" void registerGPU(AcReal* farray)
{

    // AcReal* profile_x_host = (AcReal*)malloc(sizeof(AcReal)*mx);
    // for(int i=0;i<mx;i++){
    //     profile_x_host[i] = (AcReal)i;
    //     printf("profile_x_host[%d]=%f\n",i,profile_x_host[i]);
    // }
    // mesh.profiles[PROFILE_X] = profile_x_host; 


    size_t offset=0;
    for (int i = 0; i < NUM_VTXBUF_HANDLES; ++i) {
        mesh.vertex_buffer[VertexBufferHandle(i)] = &farray[offset];
        test_mesh.vertex_buffer[VertexBufferHandle(i)] = (AcReal*)malloc(sizeof(AcReal)*mw);
        offset+=mw;
    }
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
}
/***********************************************************************************************/
extern "C" void initGPU()
{
    // Initialize GPUs in the node
    AcResult res=acCheckDeviceAvailability();
}
/***********************************************************************************************/
void setupConfig(AcMeshInfo & config){

printf("nx etc. %d %d %d %.14f %.14f %.14f \n",nxgrid,nygrid,nzgrid,dx,dy,dz);
     config.int_params[AC_nx]=nxgrid;
     config.int_params[AC_ny]=nygrid;
     config.int_params[AC_nz]=nzgrid;
     config.int_params[AC_mx] = mxgrid;
     config.int_params[AC_my] = mygrid;
     config.int_params[AC_mz] = mzgrid;
      config.int_params[AC_nx_min] = l1;
      config.int_params[AC_nx_max] = l2;
      config.int_params[AC_ny_min] = m1;
      config.int_params[AC_ny_max] = m2;
      config.int_params[AC_nz_min] = n1;
      config.int_params[AC_nz_max] = n2;
     config.real_params[AC_dsx]=dx;
     config.real_params[AC_dsy]=dy;
     config.real_params[AC_dsz]=dz;
printf("%d: l1i etc. %d %d %d %d %d %d \n", pid, l1i,l2i,n1i,n2i,m1i,m2i);
printf("%d: l1 etc. %d %d %d %d %d %d \n", pid, l1,l2,n1,n2,m1,m2);
     config.real_params[AC_dsmin]   = std::min(dx,std::min(dy,dz));
     config.real_params[AC_xlen]=6.28318548;
     config.real_params[AC_ylen]=6.28318548;
     config.real_params[AC_zlen]=6.28318548;
     config.real_params[AC_xorig]=-3.14159274;
     config.real_params[AC_yorig]=-3.14159274;
  config.real_params[AC_zorig]=-3.14159274;
     config.real_params[AC_mu0]=mu0;
    printf("Done setUpConfig\n");

#include "PC_modulepars.h"

}
/***********************************************************************************************/
void checkConfig(AcMeshInfo & config){
printf("check that config is correct\n");
//printf("setupConfig:mesh.info.real_params[AC_k1_ff]= %f \n",config.real_params[AC_k1_ff]);
#if LENTROPY
     printf("lpressuregradientgas= %d %d \n", lpressuregradient_gas, config.int_params[AC_lpressuregradient_gas]);
#endif
#if LENTROPY
     printf("chi= %f %f \n", chi, config.real_params[AC_chi]);
#endif
#if LVISCOSITY
     printf("nu= %f %f \n", nu, config.real_params[AC_nu_visc]);
     printf("zeta= %f %f \n", zeta, config.real_params[AC_zeta]);
#endif
#if LMAGNETIC
     printf("eta= %f %f \n", eta, config.real_params[AC_eta]);
#endif
#if LEOS
     printf("cs20= %f %f \n", cs20, config.real_params[AC_cs20]);
     printf("gamma= %f %f \n", gamma, config.real_params[AC_gamma]);
     printf("cv= %f %f \n", cv, config.real_params[AC_cv_sound]);
     printf("cp= %f %f \n", cp, config.real_params[AC_cp_sound]);
     printf("lnT0= %f %f \n", lnTT0, config.real_params[AC_lnTT0]);
     printf("lnrho0= %f %f \n", lnrho0, config.real_params[AC_lnrho0]);
#endif
#if LFORCING
    //  printf("iforcing_zsym= %f %f \n", iforcing_zsym, config.int_params[AC_iforcing_zsym]);
    //  printf("k1_ff= %f %f \n", k1_ff, config.real_params[AC_k1_ff]);
    //  printf("tforce_stop= %f %f \n", tforce_stop, config.real_params[AC_tforce_stop]);
     //printf("k1_ff,profx_ampl, val= %f %d %lf %lf\n", k1_ff, profx_ampl, profx_ampl[0], profx_ampl[nx-1]);
#endif
}
/***********************************************************************************************/

void loadProfiles(AcMeshInfo & config){
// #if LFORCING
//      PUT(profx_ampl,nx,0,0)
//      PUT(profy_ampl,0,my,0)
//      PUT(profz_ampl,0,0,mz)
//      PUT(profx_hel,nx,0,0)
//      PUT(profy_hel,0,my,0)
//      PUT(profz_hel,0,0,mz)
// #endif
    // acGridLoadProfile(STREAM_DEFAULT, PROFILE_X, config);
}
/***********************************************************************************************/
extern "C" void copyVBApointers(AcReal **in, AcReal **out){
    Device device = acGridGetDevice();
    *in = device->vba.in[0];
    *out = device->vba.out[0];
}

extern "C" void getFArrayIn(AcReal **p_f_in){
    Device device = acGridGetDevice();
    *p_f_in = device->vba.in[0];
}
extern "C" void initializeGPU(AcReal **farr_GPU_in, AcReal **farr_GPU_out)
{
    //Setup configurations used for initializing and running the GPU code
#if PACKED_DATA_TRANSFERS
        // initLoadStore();
#endif
        setupConfig(mesh.info);



printf("Before set domain decomp\n");
printf("nprocy_c: 2:%d\n",nprocy);
printf("nprocx_c: 1:%d\n",nprocx);
printf("nprocz_c: 1:%d\n",nprocz);
printf("hi :(\n");
fflush(stdout);
printf("Setting domain int3\n");
int3 decomp = {nprocx,nprocy,nprocx};
fflush(stdout);
acCheckDeviceAvailability();
printf("Setted domain int3\n");
fflush(stdout);

acGridSetDomainDecomposition(decomp);
checkConfig(mesh.info);
fflush(stdout);
    acGridInit(mesh.info);

// printf("Done grid init\n");
//         AcReal *p[2];
// printf("Before acGridGetVBApointers\n");
//         if (acGridGetVBApointers(p)==AC_SUCCESS) {
// printf("After acGridGetVBApointers\n");
//           *farr_GPU_in=p[0];
//           *farr_GPU_out=p[1];
// printf("From grid layer: vbapointer= %p %p \n", *farr_GPU_in, *farr_GPU_out);
//         } else {
//           *farr_GPU_in=NULL;
//           *farr_GPU_out=NULL;
//         }
    VertexBufferHandle all_fields[NUM_VTXBUF_HANDLES];
    for(int i=0;i<NUM_VTXBUF_HANDLES;i++){
        all_fields[i] = (VertexBufferHandle)i;
    }

    randomize_graph = acGridBuildTaskGraph
        ({

        acHaloExchange(all_fields),
        acBoundaryCondition(BOUNDARY_XYZ, BOUNDCOND_PERIODIC, all_fields),
        acCompute(KERNEL_randomize, all_fields),
        });
    if(with_boundconds){
        graph_1 = acGridBuildTaskGraph({
            acHaloExchange(all_fields),
            acBoundaryCondition(BOUNDARY_XYZ, BOUNDCOND_PERIODIC, all_fields),
            acCompute(KERNEL_twopass_solve_intermediate, all_fields),
            acCompute(KERNEL_twopass_solve_final, all_fields),
        });
    }else{
        graph_1 = acGridBuildTaskGraph({
            acHaloExchange(all_fields),
            acCompute(KERNEL_twopass_solve_intermediate, all_fields),
            acCompute(KERNEL_twopass_solve_final, all_fields),
        });

    }
    if(with_boundconds){
        graph_2 = acGridBuildTaskGraph(
        {
            acHaloExchange(all_fields),
            acBoundaryCondition(BOUNDARY_XYZ, BOUNDCOND_PERIODIC, all_fields),
            acCompute(KERNEL_twopass_solve_intermediate, all_fields),
            acCompute(KERNEL_twopass_solve_final, all_fields),
        });
    }else{
        graph_2 = acGridBuildTaskGraph(
        {
            acHaloExchange(all_fields),
            acCompute(KERNEL_twopass_solve_intermediate, all_fields),
            acCompute(KERNEL_twopass_solve_final, all_fields),
        });
    }
    if(with_boundconds){
        graph_3 = acGridBuildTaskGraph(
        {
            acHaloExchange(all_fields),
            acBoundaryCondition(BOUNDARY_XYZ, BOUNDCOND_PERIODIC, all_fields),
            acCompute(KERNEL_twopass_solve_intermediate, all_fields),
            acCompute(KERNEL_twopass_solve_final, all_fields),
        });
    }else{
        graph_3 = acGridBuildTaskGraph(
        {
            acHaloExchange(all_fields),
            acCompute(KERNEL_twopass_solve_intermediate, all_fields),
            acCompute(KERNEL_twopass_solve_final, all_fields),
        });

    }
    printf("BUILD graphs\n");
    fflush(stdout);
    // initialize diagnostics
       //init_diagnostics();
}
/***********************************************************************************************/
extern "C" void copyFarray() 
{
    //    AcResult res=acGridStoreMesh(STREAM_DEFAULT, &mesh);
    // acGridGetDevice();
    // AcReal* vba_in = acGridGetVBApointerIn();
    // acGridSynchronizeStream(STREAM_ALL);
    // // hipMemcpy(&mesh.vertex_buffer[0],vba_in,mw*NUM_VTXBUF_HANDLES*sizeof(AcReal),hipMemcpyDeviceToHost);
    acGridSynchronizeStream(STREAM_ALL);
    acDeviceStoreMesh(acGridGetDevice(), STREAM_DEFAULT, &mesh);
    // acDeviceStoreMesh(acGridGetDevice(), STREAM_DEFAULT, &mesh);
    acGridSynchronizeStream(STREAM_ALL);
// printf("store all %d \n",res); fflush(stdout);
}
/***********************************************************************************************/
extern "C" void finalizeGPU()
{
#if PACKED_DATA_TRANSFERS
   acDeviceStoreMesh(acGridGetDevice(), STREAM_DEFAULT, &mesh);
#endif
    // Deallocate everything on the GPUs and reset
       AcResult res=acGridQuit();
}
/***********************************************************************************************/
