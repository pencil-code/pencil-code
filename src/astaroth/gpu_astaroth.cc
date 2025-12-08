/*
   Date:   8-Feb-2017
   Author: M. Rheinhardt & j. Pekkilae
   Description:
           ANSI C and standard library callable function wrappers for ASTAROTH-nucleus functions to be called from Fortran.
  Comments: 
*/
// General headers.
#include <sstream>
#include <math.h>
#include <algorithm>
#include <chrono>
#include <iostream>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <mpi.h>
#include <sys/resource.h>
#include <fstream>
#include <dlfcn.h>

//TP: defined here since mpi.h can have its own definition of DOUBLE_PRECISION
//    and we don't want to conflict it with it. This is at least true on my laptop
#if AC_DOUBLE_PRECISION
#define DOUBLE_PRECISION 1
#endif

#define CUDA_ERRCHK(X)
int nt = 0;
int counter = 0;
int train_counter = 0;
int done = 0;
int my_rank;
bool calculated_coeff_scales = false;

// Astaroth headers
#include "astaroth.h"
#include "submodule/stdlib/fft.h"

#define real AcReal
#include "math_utils.h"
#define EXTERN
#define FINT int

#if DOUBLE_PRECISION
  #define REAL_MAX DBL_MAX
#else
  #define REAL_MAX FLT_MAX
#endif
#include "../cparam_c.h"

//TP: these are ugly but for the moment we live with these
#if TRANSPILATION
  #define limplicit_diffusion_with_fft limplicit_diffusion_with_fft__mod__implicitdiffusion
  #define limplicit_diffusion_with_cg  limplicit_diffusion_with_cg__mod__implicitdiffusion
  #define limplicit_resistivity limplicit_resistivity__mod__magnetic
  #define limplicit_viscosity   limplicit_viscosity__mod__viscosity
  #define lcumulative_df_on_gpu lcumulative_df_on_gpu__mod__cdata
  #define lbidiagonal_derij lbidiagonal_derij__mod__cdata
  #define nu nu__mod__viscosity
  #define nu_hyper2 nu_hyper2__mod__viscosity
  #define nu_hyper3 nu_hyper3__mod__viscosity
  #define n_odevars n_odevars__mod__cdata 
  #define AC_f_ode AC_f_ode__mod__cdata 
  #define gamma gamma__mod__equationofstate

  #define eta eta__mod__magnetic
  #define eta_hyper2 eta_hyper2__mod__magnetic
  #define eta_hyper3 eta_hyper3__mod__magnetic

  #define chi chi__mod__energy
  #define chi_hyper2 chi_hyper2__mod__energy
  #define chi_hyper3 chi_hyper3__mod__energy

  #define lmorton_curve lmorton_curve__mod__cdata
  #define ltest_bcs     ltest_bcs__mod__cdata
  #define num_substeps  num_substeps__mod__cdata
  #define maux_vtxbuf_index maux_vtxbuf_index__mod__cdata
  #define read_vtxbuf_from_gpu read_vtxbuf_from_gpu__mod__cdata 
  #define ldt ldt__mod__cdata
  #define dt dt__mod__cdata
  #define it it__mod__cdata
  #define dx dx__mod__cdata
  #define dy dy__mod__cdata
  #define dz dz__mod__cdata
  #define mu0 mu0__mod__cdata
  #define ldensity_nolog ldensity_nolog__mod__cdata 
  #define cdt            cdt__mod__cdata 
  #define itorder        itorder__mod__cdata
  #define dtinc          dtinc__mod__cdata
  #define dtdec          dtdec__mod__cdata

  #define AC_ldensity_nolog AC_ldensity_nolog__mod__cdata
  #define AC_mu0  AC_mu0__mod__cdata

  #define cdtv  cdtv__mod__cdata
  #define cdtv2 cdtv2__mod__cdata
  #define cdtv3 cdtv3__mod__cdata

  #define lcylindrical_coords lcylindrical_coords__mod__cdata
  #define lspherical_coords   lspherical_coords__mod__cdata
  #define lcartesian_coords   lcartesian_coords__mod__cdata

  #define lequidist lequidist__mod__cdata

  #define dx_1 dx_1__mod__cdata
  #define dy_1 dy_1__mod__cdata
  #define dz_1 dz_1__mod__cdata

  #define dx_tilde dx_tilde__mod__cdata
  #define dy_tilde dy_tilde__mod__cdata
  #define dz_tilde dz_tilde__mod__cdata
  
  #define rcyl_mn1 rcyl_mn1__mod__cdata
  #define r1_mn    r1_mn__mod__cdata
  #define sin1th   sin1th__mod__cdata
  #define cotth    cotth__mod__cdata

  #define luse_trained_tau luse_trained_tau__mod__training
  #define lcpu_timestep_on_gpu   lcpu_timestep_on_gpu__mod__cdata
  #define lperi                  lperi__mod__cdata
  #define lxyz                   lxyz__mod__cdata
  #define lac_sparse_autotuning  lac_sparse_autotuning__mod__cdata
  #define ldebug ldebug__mod__cdata
  
  #define deltay  deltay__mod__cdata
  #define eps_rkf eps_rkf__mod__cdata
  #define itauxx itauxx__mod__training
  #define itauxy itauxy__mod__training
  #define itauxz itauxz__mod__training
  #define itauyy itauyy__mod__training
  #define itauyz itauyz__mod__training
  #define itauzz itauzz__mod__training

  #define lread_all_vars_from_device lread_all_vars_from_device__mod__cdata
  #define lcuda_aware_mpi            lcuda_aware_mpi__mod__cdata
  #define lsecond_force lsecond_force__mod__forcing
  #define lforce_helical lforce_helical__mod__forcing

  #define lrmv lrmv__mod__cdata
  #define lconserve_total_mass lconserve_total_mass__mod__density
  #define tstart_selfgrav tstart_selfgrav__mod__selfgravity
#endif

AcReal dt1_interface;
static int rank;
static MPI_Comm comm_pencil = MPI_COMM_NULL;

// Astaroth objects instantiation.
static AcMesh mesh = acInitMesh();
//static AcMesh test_mesh;

void torch_trainCAPI(int sub_dims[3], AcReal* input, AcReal* label, AcReal* loss_val);
void torch_inferCAPI(int sub_dims[3], AcReal* input, AcReal* label);
void scaling();
void print_debug();
float MSE();
//void torch_createmodel(const char* name, const char* config_fname, MPI_Comm mpi_comm, int device);

extern "C" void copyFarray(AcReal* f);    // forward declaration
extern "C" void loadFarray(); // forward declaration
void denormalize(std::string filename, AcRealSymmetricTensor &tau_means, AcRealSymmetricTensor &tau_stds);

/***********************************************************************************************/
AcReal cpu_pow(AcReal const val, AcReal exponent)
{
// Masks hip GPU power function.
    return std::pow(val, exponent);
}
/***********************************************************************************************/
AcReal sign(const AcReal a, const AcReal b)
{
	return (b < (AcReal)0.0) ? -abs(a) : abs(a);
}
/***********************************************************************************************/
bool has_nans(AcMesh mesh_in)
{
  bool res = false;
  AcMeshDims dims = acGetMeshDims(acGridGetLocalMeshInfo());
  const auto DEVICE_VTXBUF_IDX = [&](const int x, const int y, const int z)
  				{
					return acVertexBufferIdx(x,y,z,mesh.info);
  				};
  for (size_t i = dims.n0.x; i < dims.n1.x; i++)
  {
    for (size_t j = dims.n0.y; j < dims.n1.y; j++)
    {
      for (size_t k = dims.n0.z; k < dims.n1.z; k++)
      {
        for (size_t ivar = 0; ivar < acGetNumFields(); ivar++)
        {
          if (mesh_in.vertex_buffer[ivar] == NULL) continue;
          if (isnan(mesh_in.vertex_buffer[ivar][DEVICE_VTXBUF_IDX(i, j, k)]))
          {
            res = true;
            acLogFromRootProc(rank,"nan at %d,%d,%d\n", i, j, k);
            acLogFromRootProc(rank,"field = %d", ivar);
          }
        }
      }
    }
  }
  return res;
}
/***********************************************************************************************/
extern "C" void testRHS(AcReal *farray_in, AcReal *dfarray_truth)
{
  const auto DEVICE_VTXBUF_IDX = [&](const int x, const int y, const int z)

  				{
					return acVertexBufferIdx(x,y,z,mesh.info);
  				};
  // __energy_MOD_pushpars2c(p_par_energy);
  // mesh.info.profiles[AC_dlnhcond_prof]=dlnhcond_prof; // [2-1] dlnhcond_prof real(nz)
  // for (int i=20;i<30;++i)
  //   acLogFromRootProc(rank,"C: dlnhcond_prof %d=%f\n",i,mesh.info.profiles[AC_dlnhcond_prof][i]);
  // fflush(stdout);
  // return;
  // make_tasks(diagnostics_func, reduction_func, finalize_func, write_func);
  // thread_pool.WaitAll();
  constexpr int num_of_steps = 100;

  AcMesh mesh_true;
  AcMesh mesh_test;
  const AcReal epsilon = (AcReal)pow(0.1,12);
  // constexpr AcReal epsilon = 0.0;
  constexpr AcReal local_dt = (AcReal)0.001;
  Device dev = acGridGetDevice();
  acDeviceSetInput(acGridGetDevice(),AC_dt,local_dt);
  acGridSynchronizeStream(STREAM_ALL);

  size_t offset = 0;
  for (int i = 0; i < acGetNumFields(); ++i)
  {
    mesh_test.vertex_buffer[VertexBufferHandle(i)] = &farray_in[offset];
    offset += mw;
  }
  offset = 0;
	for (int i = 0; i < acGetNumFields(); ++i)
  {
    mesh_true.vertex_buffer[VertexBufferHandle(i)] = &dfarray_truth[offset];
    offset += mw;
  }
  AcMeshDims dims = acGetMeshDims(acGridGetLocalMeshInfo());
  acLogFromRootProc(rank,"n0: %d,%d,%d\trank: %d\n", dims.n0.x, dims.n0.y, dims.n0.z, rank);
  acLogFromRootProc(rank,"n1: %d,%d,%d\trank: %d\n", dims.n1.x, dims.n1.y, dims.n1.z, rank);

  //dryrun
  // acGridLaunchKernel(STREAM_DEFAULT, twopass_solve_intermediate_step0, dims.n0,dims.n1);
  // acGridSynchronizeStream(STREAM_ALL);
  // acGridSwapBuffers();
  // acGridSynchronizeStream(STREAM_ALL);
  // acGridLaunchKernel(STREAM_DEFAULT, twopass_solve_final_step0, dims.n0,dims.n1);
  // acGridSynchronizeStream(STREAM_ALL);
  // acGridSwapBuffers();
  // acGridSynchronizeStream(STREAM_ALL);
  // acGridLaunchKernel(STREAM_DEFAULT, twopass_solve_intermediate_step1, dims.n0,dims.n1);
  // acGridSynchronizeStream(STREAM_ALL);
  // acGridSwapBuffers();
  // acGridSynchronizeStream(STREAM_ALL);

  // acGridExecuteTaskGraph(rhs_test_graph,1);
  AcTaskGraph *rhs =  acGetOptimizedDSLTaskGraph(AC_rhs);
  for (int i = 0; i < 3; ++i)
  {
  	acDeviceSetInput(acGridGetDevice(), AC_step_num, (PC_SUB_STEP_NUMBER) i);
  	acGridExecuteTaskGraph(rhs,1);
  	acGridSynchronizeStream(STREAM_ALL);
  }

  acGridSynchronizeStream(STREAM_ALL);
  acDeviceLoadMesh(acGridGetDevice(), STREAM_DEFAULT, mesh_test);
  acGridSynchronizeStream(STREAM_ALL);
  acDeviceStoreMesh(acGridGetDevice(), STREAM_DEFAULT, &mesh);
  acGridSynchronizeStream(STREAM_ALL);

  bool loaded_correct = true;
  for (size_t i = 0; i < mx; i++)
  {
    for (size_t j = 0; j < my; j++)
    {
      for (size_t k = 0; k < mz; k++)
      {
        for (size_t ivar = 0; ivar < mfarray; ivar++)
        {
          AcReal out_val = mesh.vertex_buffer[ivar][DEVICE_VTXBUF_IDX(i, j, k)];
          AcReal true_val = mesh_test.vertex_buffer[ivar][DEVICE_VTXBUF_IDX(i, j, k)];
          if (out_val != true_val)
          {
            loaded_correct = false;
            acLogFromRootProc(rank,"Loaded val wrong at %d,%d,%d\n", i, j, k);
            acLogFromRootProc(rank,"field = %d", ivar);
            acLogFromRootProc(rank,"Loaded val: %f\tTRUE val: %f\n", (double)out_val, (double)true_val);
          }
        }
      }
    }
  }
  //AcReal derx_ux;
  //AcReal derx_normal;
  //AcReal dery_uy;
  //AcReal derz_uz;
  //test calculating divu different ways
  //const int y = 42;
  //const int x = 42;
  //const int z = 42;
  //const AcReal AC_inv_dsx = 2.0;

    // [0][0][-3] = -AC_inv_dsx * DER1_3,
    // [0][0][-2] = -AC_inv_dsx * DER1_2,
    // [0][0][-1] = -AC_inv_dsx * DER1_1,
    // [0][0][1]  = AC_inv_dsx * DER1_1,
    // [0][0][2]  = AC_inv_dsx * DER1_2,
    // [0][0][3]  = AC_inv_dsx * DER1_3

  // derx_normal = -AC_inv_dsx*DER1_3*(mesh.vertex_buffer[0][DEVICE_VTXBUF_IDX(x-3,y,z)]);
  // derx_normal += -AC_inv_dsx*DER1_2*(mesh.vertex_buffer[0][DEVICE_VTXBUF_IDX(x-2,y,z)]);
  // derx_normal += -AC_inv_dsx*DER1_1*(mesh.vertex_buffer[0][DEVICE_VTXBUF_IDX(x-1,y,z)]);
  // derx_normal += AC_inv_dsx*DER1_1*(mesh.vertex_buffer[0][DEVICE_VTXBUF_IDX(x+1,y,z)]);
  // derx_normal += AC_inv_dsx*DER1_2*(mesh.vertex_buffer[0][DEVICE_VTXBUF_IDX(x+2,y,z)]);
  // derx_normal += AC_inv_dsx*DER1_3*(mesh.vertex_buffer[0][DEVICE_VTXBUF_IDX(x+3,y,z)]);
  // acLogFromRootProc(rank,"test derx_ux: %.7e\n",derx_ux);
  // acLogFromRootProc(rank,"normal derx_ux. %.7e\n",derx_normal);
  fflush(stdout);
  // return;
  if (loaded_correct)
  {
    acLogFromRootProc(rank,"loaded correct data\n");
  }
  else
  {
    acLogFromRootProc(rank,"loaded incorrect data :(\n");
  }

  //set output buffer to 0 since if we are reading from it we don't want NaNs
  // acGridExecuteTaskGraph(rhs_test_graph, 1);
  // acGridSynchronizeStream(STREAM_ALL);
  // acGridExecuteTaskGraph(rhs_test_rhs_1, 1);
  // acGridSynchronizeStream(STREAM_ALL);
  acGridLaunchKernel(STREAM_DEFAULT, AC_BUILTIN_RESET, dims.n0, dims.n1);
  acGridSynchronizeStream(STREAM_ALL);

  //actual run
  //TP: works for only rk3 but is anyone anyways using this? Probably not
  for (int i=0;i<num_of_steps;i++){
    for (int substep = 0; substep <3; ++substep)
    {
    	acDeviceSetInput(acGridGetDevice(), AC_step_num, (PC_SUB_STEP_NUMBER)substep);
    	acGridExecuteTaskGraph(acGetOptimizedDSLTaskGraph(AC_rhs),1);
    	acGridSynchronizeStream(STREAM_ALL);
    }
    acGridSynchronizeStream(STREAM_ALL);
    acDeviceStoreMesh(acGridGetDevice(), STREAM_DEFAULT, &mesh);
  //  acGridSynchronizeStream(STREAM_ALL);
  //  AcReal max_uux2 = 0.0;
  //  AcReal max_uuy2 = 0.0;
  //  AcReal max_uuz2 = 0.0;
  //for (int i = dims.n0.x; i < dims.n1.x; i++)
  //{
  //  for (int j = dims.n0.y; j < dims.n1.y; j++)
  //  {
  //    for (int k = dims.n0.z; k < dims.n1.z; k++)
  //    {
  //        max_uux2 = max(max_uux2,pow(mesh.vertex_buffer[0][DEVICE_VTXBUF_IDX(i,j,k)],2));
  //        max_uuy2 = max(max_uuy2,pow(mesh.vertex_buffer[1][DEVICE_VTXBUF_IDX(i,j,k)],2));
  //        max_uuz2 = max(max_uuz2,pow(mesh.vertex_buffer[2][DEVICE_VTXBUF_IDX(i,j,k)],2));
  //    }
  //  }
  //}
  //acLogFromRootProc(rank,"GPU: uumax: %.7e\n", pow(max_uux2+max_uuy2+max_uuz2,0.5));
  //  acGridSynchronizeStream(STREAM_ALL);
  }
    acGridSynchronizeStream(STREAM_ALL);
  // acGridLaunchKernel(STREAM_DEFAULT, twopass_solve_intermediate_step0, dims.n0,dims.n1);
  // acGridSynchronizeStream(STREAM_ALL);
  // acGridSwapBuffers();
  // acGridSynchronizeStream(STREAM_ALL);
  // acGridLaunchKernel(STREAM_DEFAULT, twopass_solve_final_step0, dims.n0,dims.n1);
  // acGridSynchronizeStream(STREAM_ALL);
  // acGridSwapBuffers();
  // acGridSynchronizeStream(STREAM_ALL);
  // acGridLaunchKernel(STREAM_DEFAULT, twopass_solve_intermediate_step1, dims.n0,dims.n1);
  // acGridSynchronizeStream(STREAM_ALL);
  // acGridSwapBuffers();
  // acGridSynchronizeStream(STREAM_ALL);

  // acGridLaunchKernel(STREAM_DEFAULT, twopass_solve_final_step1, dims.n0,dims.n1);
  // acGridSynchronizeStream(STREAM_ALL);
  // acGridSwapBuffers();
  // acGridSynchronizeStream(STREAM_ALL);

  // acGridLaunchKernel(STREAM_DEFAULT, twopass_solve_intermediate_step2, dims.n0,dims.n1);
  // acGridSynchronizeStream(STREAM_ALL);
  // acGridSwapBuffers();
  // acGridSynchronizeStream(STREAM_ALL);

  // acGridLaunchKernel(STREAM_DEFAULT, twopass_solve_final_step2, dims.n0,dims.n1);
  // acGridSynchronizeStream(STREAM_ALL);
  // acGridSwapBuffers();
  // acGridSynchronizeStream(STREAM_ALL);

  acGridSynchronizeStream(STREAM_ALL);
  acDeviceStoreMesh(acGridGetDevice(), STREAM_DEFAULT, &mesh);

  bool passed = true;
  AcReal max_abs_not_passed_val=-1.0;
  AcReal true_pair{};
  AcReal max_abs_relative_difference =-1.0;
  AcReal max_abs_value = -1.0;
  AcReal min_abs_value = 1.0;
  AcReal gpu_val_for_largest_diff{};
  AcReal true_val_for_largest_diff{};
  int num_of_points_where_different[NUM_VTXBUF_HANDLES] = {0};

  for (size_t i = dims.n0.x; i < dims.n1.x; i++)
  {
    for (size_t j = dims.n0.y; j < dims.n1.y; j++)
    {
      for (size_t k = dims.n0.z; k < dims.n1.z; k++)
      {
        for (size_t ivar = 0; ivar < acGetNumFields(); ivar++)
        {
          AcReal out_val = mesh.vertex_buffer[ivar][DEVICE_VTXBUF_IDX(i, j, k)];
          AcReal true_val = mesh_true.vertex_buffer[ivar][DEVICE_VTXBUF_IDX(i, j, k)];
          AcReal abs_diff = fabs(out_val - true_val);
          if (fabs(true_val) > max_abs_value) max_abs_value = fabs(true_val);
          if (fabs(true_val) < min_abs_value) min_abs_value = fabs(true_val);
          if ((AcReal)(abs_diff/true_val) > epsilon || (true_val == (AcReal)0.0 && (AcReal)fabs(out_val) > (AcReal)pow(0.1,13)) || (epsilon == (AcReal)0.0 && true_val != out_val))
          {
            passed = false;
            num_of_points_where_different[ivar]++;
            // acLogFromRootProc(rank,"rhs val wrong at %d,%d,%d\n", i, j, k);
            // acLogFromRootProc(rank,"field = %d", ivar);
            // acLogFromRootProc(rank,"GPU val: %.7e\tTRUE val: %.7e\n", out_val, true_val);
            // acLogFromRootProc(rank,"PID: %d\n", rank);
            if (max_abs_not_passed_val<abs(out_val)){
              max_abs_not_passed_val = abs(out_val);
              true_pair = true_val;
            }
            if (true_val != (AcReal)0.0){
              if (max_abs_relative_difference<(abs_diff/true_val)){
                max_abs_relative_difference=(abs_diff/true_val);
                gpu_val_for_largest_diff = out_val;
                true_val_for_largest_diff = true_val;
              }
            }  
          }
          if (isnan(out_val))
          {
            acLogFromRootProc(rank,"nan before at %d,%d,%d,%d!\n!",i,j,k,ivar);
            acLogFromRootProc(rank,"%.7e\n",(double)out_val);
          }
        }
      }
    }
  }
  auto volume_size = [](const auto& a)
  {
	  return a.x*a.y*a.z;
  };
  for (int ivar=0;ivar<acGetNumFields();ivar++)
    acLogFromRootProc(rank,"ratio of values wrong for field: %d\t %f\n",ivar,(double)num_of_points_where_different[ivar]/volume_size(dims.n1-dims.n0));
  passed &= !has_nans(mesh);
  if (passed)
  {
    acLogFromRootProc(rank,"Passed GPU test :)\n");
  }
  else
  {
    acLogFromRootProc(rank,"Did not pass GPU test :(\n");
  }
  acLogFromRootProc(rank,"max abs not passed val: %.7e\t%.7e\n",(double)max_abs_not_passed_val, (double)fabs(true_pair));
  acLogFromRootProc(rank,"max abs relative difference val: %.7e\n",(double)max_abs_relative_difference);
  acLogFromRootProc(rank,"largest difference: %.7e\t%.7e\n",(double)gpu_val_for_largest_diff, (double)true_val_for_largest_diff);
  acLogFromRootProc(rank,"abs range: %.7e-%7e\n",(double)min_abs_value,(double)max_abs_value);
  fflush(stdout);
}
/***********************************************************************************************/
AcReal to_real(void* param, const char* name)
{
	if (param == NULL)
	{
		fprintf(stderr,"Passed NULL to pushparsc, expected real: %s!!\n",name);
		abort();
	}
	return *((AcReal*)param);
}
/***********************************************************************************************/
int to_int(void* param, const char* name)
{
	if (param == NULL)
	{
		fprintf(stderr,"Passed NULL to pushparsc, expected int: %s!!\n",name);
		abort();
	}
	return *((int*)param);
}
/***********************************************************************************************/
bool to_bool(void* param, const char* name)
{
	if (param == NULL)
	{
		fprintf(stderr,"Passed NULL to pushparsc, expected bool: %s!!\n",name);
		abort();
	}
	return *((bool*)param);
}
/***********************************************************************************************/
int3 to_int3(void* param, const char* name)
{
	if (param == NULL)
	{
		fprintf(stderr,"Passed NULL to pushparsc, expected int3: %s!!\n",name);
		abort();
	}
        int* arr = (int*)param;
        return (int3){arr[0],arr[1],arr[2]};
}
/***********************************************************************************************/
AcReal3 to_real3(void* param, const char* name)
{
	if (param == NULL)
	{
		fprintf(stderr,"Passed NULL to pushparsc, expected real3: %s!!\n",name);
		abort();
	}
        AcReal* arr = (AcReal*)param;
        return (AcReal3){arr[0],arr[1],arr[2]};
}
/***********************************************************************************************/
AcBool3 to_bool3(void* param, const char* name)
{
	if (param == NULL)
	{
		fprintf(stderr,"Passed NULL to pushparsc, expected bool3: %s!!\n",name);
		abort();
	}
        bool* arr = (bool*)param;
        return (AcBool3){arr[0],arr[1],arr[2]};
}
/***********************************************************************************************/
//typedef void (*rangefunc)(const int a, const int b);
/***********************************************************************************************/
// PC interface headers.
#include "PC_moduleflags.h"

#if TRANSPILATION
#if LFORCING
torus_rect to_torus_rect(void* param, const char* name)
{
       if (param == NULL)
       {
               fprintf(stderr,"Passed NULL to pushparsc: %s!!\n",name);
               abort();
       }
       //TP: placeholder for now before testing is torus_react being POD sufficient for save access here
       return (torus_rect){};
}
#endif
#endif
/***********************************************************************************************/
//#include "../cdata_c.h"
#include "../sub_c.h"           // provides set_dt
#include "../boundcond_c.h"     // provides boundconds[xyz] etc.
#include "PC_module_parfuncs.h" // provides stuff from physics modules

#if LFORCING
  #include "../forcing_c.h"     // provides forcing_pars_hel
  #include "forcing.h"
#endif

/***********************************************************************************************/
int DCONST(const AcIntParam param)
{
  return mesh.info[param];
}
/***********************************************************************************************/
int3 DCONST(const AcInt3Param param)
{
  return mesh.info[param];
}
/***********************************************************************************************/
AcReal DCONST(const AcRealParam param)
{
  return mesh.info[param];
}
/***********************************************************************************************/
AcReal3 DCONST(const AcReal3Param param)
{
  return mesh.info[param];
}
/***********************************************************************************************/
template <typename T>
void
pushpars_special(const char* func_name, T p_pars)
{
 	void* handle = dlopen("src/special.so",RTLD_NOW | RTLD_LOCAL);
	auto sub = reinterpret_cast<void (*)(T)>(dlsym(handle,func_name));
	sub(p_pars);
	dlclose(handle);
}
/***********************************************************************************************/
#define PCLoad acPushToConfig
/***********************************************************************************************/
void setupConfig(AcMeshInfo& config)
{
  config = acInitInfo();
 
  #include "PC_modulepars.h"

  PCLoad(config, AC_use_cuda_aware_mpi,lcuda_aware_mpi);
  PCLoad(config, AC_bidiagonal_derij,lbidiagonal_derij);
  //TP: loads for non-Cartesian derivatives

#if TRANSPILATION
  PCLoad(config, AC_x,x__mod__cdata);
  PCLoad(config, AC_y,y__mod__cdata);
  PCLoad(config, AC_z,z__mod__cdata);

  PCLoad(config, AC_inv_cyl_r,rcyl_mn1);
  PCLoad(config, AC_inv_r,r1_mn);
  PCLoad(config, AC_inv_sin_theta,sin1th);
  PCLoad(config, AC_cot_theta,cotth);

//  loads for non-equidistant grids
  PCLoad(config,AC_nonequidistant_grid, (AcBool3){!lequidist.x,!lequidist.y,!lequidist.z});
  PCLoad(config,AC_inv_mapping_func_derivative_x,dx_1);
  PCLoad(config,AC_inv_mapping_func_derivative_y,dy_1);
  PCLoad(config,AC_inv_mapping_func_derivative_z,dz_1);

  PCLoad(config,AC_mapping_func_tilde_x,dx_tilde);
  PCLoad(config,AC_mapping_func_tilde_y,dy_tilde);
  PCLoad(config,AC_mapping_func_tilde_z,dz_tilde);

//  loads for slope-limited-diffusion related arrays
  PCLoad(config, AC_x12,x12__mod__cdata);
  PCLoad(config, AC_y12,y12__mod__cdata);
  PCLoad(config, AC_sinth12,sinth12__mod__cdata);
  PCLoad(config, AC_z12,z12__mod__cdata);
#endif

  PCLoad(config, AC_rk_order, itorder);
  PCLoad(config, AC_shear,lshear);
#if TRANSPILATION
  PCLoad(config, AC_rk_cumulative_df,lcumulative_df_on_gpu);
#endif

  if (lcartesian_coords)
  {
          PCLoad(config, AC_coordinate_system, AC_CARTESIAN_COORDINATES);
  }
  if (lspherical_coords)
  {
          PCLoad(config, AC_coordinate_system, AC_SPHERICAL_COORDINATES);
  }
  if (lcylindrical_coords)
  {
          PCLoad(config, AC_coordinate_system, AC_CYLINDRICAL_COORDINATES);
  }
  //TP: this is needed to not run out of memory for spherical-gdisk-planet-thermo on norlx51
#if LNEWTON_COOLING
  PCLoad(config,AC_only_default_stream_for_taskgraphs, lcpu_timestep_on_gpu);
#endif
  PCLoad(config, AC_domain_decomposition, (int3) {nprocx,nprocy,nprocz});
  PCLoad(config, AC_ngrid, (int3){nxgrid,nygrid,nzgrid});
  PCLoad(config, AC_skip_single_gpu_optim, true);

  PCLoad(config,AC_decompose_strategy,AC_DECOMPOSE_STRATEGY_EXTERNAL);
  PCLoad(config,AC_proc_mapping_strategy,lmorton_curve ? AC_PROC_MAPPING_STRATEGY_MORTON : AC_PROC_MAPPING_STRATEGY_LINEAR);

  PCLoad(config, AC_include_3d_halo_corners, ltraining);
  PCLoad(config,AC_MPI_comm_strategy,AC_MPI_COMM_STRATEGY_DUP_USER);
  config.comm->handle = comm_pencil;

// Grid and geometry related parameters

  PCLoad(config,AC_ds,(AcReal3){dx,dy,dz});
  PCLoad(config,AC_periodic_grid,lperi);
//  Overwrites Astaroth's default formula for AC_len in case it does not cover everything like non-equidistant grids
  PCLoad(config,AC_len,lxyz);

  #if LTRAINING
  AcRealSymmetricTensor tau_means{};
  AcRealSymmetricTensor tau_stds{};
	
  //The statistics for doing the inverse are only needed when using the trained tau
  if(luse_trained_tau) denormalize("normalizer.bin", tau_means, tau_stds);

	

  PCLoad(config,AC_tau_means,tau_means);
  PCLoad(config,AC_tau_stds,tau_stds);
  #endif
  PCLoad(config,AC_sparse_autotuning,lac_sparse_autotuning);

//  Enter physics related parameters in config.

  #if LDENSITY
    PCLoad(config,AC_ldensity_nolog,ldensity_nolog);
    //printf("ldensity_nolog is %d \n",config[AC_ldensity_nolog]);//ldensity_nolog);
  #endif
  //Device dev = acGridGetDevice();
#if LHYDRO
  //dev->output.real_outputs[AC_maxadvec]=0.;
#endif
#if LENTROPY
  //dev->output.real_outputs[AC_maxchi]=0.;
#endif
  acHostUpdateParams(&config); 
  if (!ltraining) config.runtime_compilation_log_dst = "ac_compilation.log";
  char cwd[9000];
  cwd[0] = '\0';
  const char* err = getcwd(cwd, sizeof(cwd));
  if (err == NULL) 
  {
	  fprintf(stderr,"Was not able to get cwd!\n");
	  exit(EXIT_FAILURE);
  }
  char build_path[18000];
  sprintf(build_path,"%s/src/astaroth/submodule/build",cwd);
  config.runtime_compilation_build_path = strdup(build_path);
}
/***********************************************************************************************/
//TP: x,y and z macros are too general
#undef x
#undef y
#undef z
/***********************************************************************************************/
std::array<AcReal,3> visc_get_max_diffus()
{
#if LVISCOSITY
	 return {nu,nu_hyper2,nu_hyper3};
#else
	 return {0.0,0.0,0.0};
#endif
}
/***********************************************************************************************/
std::array<AcReal,3> magnetic_get_max_diffus()
{
#if Lmagnetic_MODULE
	return {eta,eta_hyper2,eta_hyper3};
#else
	return {0.0,0.0,0.0};
#endif
}
/***********************************************************************************************/
std::array<AcReal,3> energy_get_max_diffus()
{
#if TRANSPILATION
	return {0.0,0.0,0.0};
#else
 	constexpr AcReal chi_hyper2=0.;
#if LENTROPY
	return {gamma*chi,gamma*chi_hyper2,gamma*chi_hyper3};
#else
	return {0.0,0.0,0.0};
#endif
#endif
}
/***********************************************************************************************/
std::array<AcReal,3> elem_wise_max(const std::array<AcReal,3>& a,const std::array<AcReal,3>& b,const std::array<AcReal,3>& c)
{
 return {
	  std::max(std::max(a[0],b[0]),c[0]),
	  std::max(std::max(a[1],b[1]),c[1]),
	  std::max(std::max(a[2],b[2]),c[2])
	};
}
/***********************************************************************************************/
AcReal max_diffus(AcReal maxnu_dyn, AcReal maxchi_dyn)
{
  AcReal3 dxyz_vals = get_dxyzs();

  auto max_diffusions = elem_wise_max(visc_get_max_diffus(), magnetic_get_max_diffus(), energy_get_max_diffus());
#if LENTROPY
  max_diffusions[0] = std::max(max_diffusions[0],maxchi_dyn);
#endif
#if LVISCOSITY
  max_diffusions[0] = std::max(max_diffusions[0],maxnu_dyn);
#endif

  return max_diffusions[0]*dxyz_vals.x/cdtv + max_diffusions[1]*dxyz_vals.y/cdtv2 + max_diffusions[2]*dxyz_vals.z/cdtv3;
}
/***********************************************************************************************/
AcReal calc_dt1_courant(const AcReal t)
{
#if TRANSPILATION
      const AcReal gpu_max_dt1 = acDeviceGetOutput(acGridGetDevice(),AC_dt1_max);
      if (lfractional_tstep_advance__mod__cdata)
      {
	      return std::max((double)gpu_max_dt1,1./(dt_incr__mod__cdata*t));
      }
      return gpu_max_dt1;
#else
      AcReal maxadvec = 0.;
#if LHYDRO
      maxadvec = acDeviceGetOutput(acGridGetDevice(), AC_maxadvec)/cdt;
      //if (rank==0) printf("rank, maxadvec= %d %e \n", rank, maxadvec);
#endif
      AcReal maxchi_dyn = 0.;
#if LENTROPY
      maxchi_dyn = acDeviceGetOutput(acGridGetDevice(), AC_maxchi);
#endif
      AcReal maxnu_dyn = 0.;
#if LVISCOSITY
      maxnu_dyn = acDeviceGetOutput(acGridGetDevice(), AC_maxnu);
//if (rank==0) printf("maxnu_dyn= %e \n", maxnu_dyn);
#endif
      fprintf(stderr,"Maxadvec: %14e\n",maxadvec);
      return (AcReal)sqrt(pow(maxadvec, 2) + pow(max_diffus(maxnu_dyn,maxchi_dyn), 2));
#endif
}
/***********************************************************************************************/
AcReal GpuCalcDt(const AcReal t)
{
	acGridSynchronizeStream(STREAM_ALL);
  	acDeviceSetInput(acGridGetDevice(), AC_step_num, (PC_SUB_STEP_NUMBER) 0);
	const auto graph = acGetOptimizedDSLTaskGraph(AC_calculate_timestep);
	acGridExecuteTaskGraph(graph,1);
	acGridSynchronizeStream(STREAM_ALL);
	acDeviceSwapBuffers(acGridGetDevice());
	return calc_dt1_courant(t);
}
/***********************************************************************************************/
extern "C" void sourceFunctionAndOpacity(int inu)
{
#if Lradiation_ray_MODULE
  	acDeviceSetInput(acGridGetDevice(), AC_frequency_bin,inu);

//TP: need some fields (like TT) on the boundaries to calculate Srad and kappa_rho
	acGridHaloExchange();
  	auto bcs = acGetOptimizedDSLTaskGraph(boundconds);	
	acGridExecuteTaskGraph(bcs,1);

	const auto info = acDeviceGetLocalConfig(acGridGetDevice());
	const Volume offsets = (Volume)
	{
		(size_t)nxgrid == 1 ? 0 : 1,
		(size_t)nygrid == 1 ? 0 : 1,
		(size_t)nzgrid == 1 ? 0 : 1
	};
	const Volume start = acGetMinNN(info)-offsets;
	const Volume end   = acGetMaxNN(info)+offsets;
	acGridExecuteTaskGraph(acGetOptimizedDSLTaskGraph(get_source_function_and_opacity,start,end),1);
//TP: should not do it here but for now for testing do everything together (works only if there is only one frequency)
//    but this explicit call to calculate sources for radiation is anyways too granular: have just one call to interface to calculate Qrad
        const auto Qintrinsic_graph = acGetOptimizedDSLTaskGraph(Qintrinsic_steps,true,Qintrinsic_bcs);
	acGridExecuteTaskGraph(Qintrinsic_graph,1);
        const auto Qextrinsic_graph = acGetOptimizedDSLTaskGraph(Qextrinsic_steps);
	acGridExecuteTaskGraph(bcs,1);
	acGridExecuteTaskGraph(Qextrinsic_graph,1);
#endif
}
/***********************************************************************************************/
extern "C" void splitUpdate(const AcReal error_tolerance, const int max_steps)
{
#if LIMPLICIT_DIFFUSION
	if(limplicit_diffusion_with_fft)
	{
#if LMAGNETIC
		if(limplicit_resistivity && lbfield)
		{
			ac_fft_split_diffusion_update(acGetF_BX(),dt,eta,acGetSPLIT_DIFFUSION_UPDATE_BUFFER_REAL(),acGetSPLIT_DIFFUSION_UPDATE_BUFFER_IMAG(),acGetSPLIT_DIFFUSION_UPDATE_BUFFER_REAL(),acGetSPLIT_DIFFUSION_UPDATE_BUFFER_IMAG());
			ac_fft_split_diffusion_update(acGetF_BY(),dt,eta,acGetSPLIT_DIFFUSION_UPDATE_BUFFER_REAL(),acGetSPLIT_DIFFUSION_UPDATE_BUFFER_IMAG(),acGetSPLIT_DIFFUSION_UPDATE_BUFFER_REAL(),acGetSPLIT_DIFFUSION_UPDATE_BUFFER_IMAG());
			ac_fft_split_diffusion_update(acGetF_BZ(),dt,eta,acGetSPLIT_DIFFUSION_UPDATE_BUFFER_REAL(),acGetSPLIT_DIFFUSION_UPDATE_BUFFER_IMAG(),acGetSPLIT_DIFFUSION_UPDATE_BUFFER_REAL(),acGetSPLIT_DIFFUSION_UPDATE_BUFFER_IMAG());
		}
#endif
#if LVISCOSITY
		if(limplicit_viscosity)
		{
			ac_fft_split_diffusion_update(acGetUUX(),dt,nu,acGetSPLIT_DIFFUSION_UPDATE_BUFFER_REAL(),acGetSPLIT_DIFFUSION_UPDATE_BUFFER_IMAG(),acGetSPLIT_DIFFUSION_UPDATE_BUFFER_REAL(),acGetSPLIT_DIFFUSION_UPDATE_BUFFER_IMAG());
			ac_fft_split_diffusion_update(acGetUUY(),dt,nu,acGetSPLIT_DIFFUSION_UPDATE_BUFFER_REAL(),acGetSPLIT_DIFFUSION_UPDATE_BUFFER_IMAG(),acGetSPLIT_DIFFUSION_UPDATE_BUFFER_REAL(),acGetSPLIT_DIFFUSION_UPDATE_BUFFER_IMAG());
			ac_fft_split_diffusion_update(acGetUUZ(),dt,nu,acGetSPLIT_DIFFUSION_UPDATE_BUFFER_REAL(),acGetSPLIT_DIFFUSION_UPDATE_BUFFER_IMAG(),acGetSPLIT_DIFFUSION_UPDATE_BUFFER_REAL(),acGetSPLIT_DIFFUSION_UPDATE_BUFFER_IMAG());
		}
#endif
	}
	else if(limplicit_diffusion_with_cg)
	{
		acGridExecuteTaskGraph(acGetOptimizedDSLTaskGraph(implicit_diffusion_init_cg),1);
		acDeviceSetInput(acGridGetDevice(),AC_implicit_diffusion_coefficient,eta);
		for(int field = 0; field < 3; ++field)
		{
			int step_num = 0;
			acDeviceSetInput(acGridGetDevice(),AC_CG_FIELD,CG_FIELD(field));
			AcReal residual = 1e100;
			bool over_max_steps = false;
			acGridExecuteTaskGraph(acGetOptimizedDSLTaskGraph(implicit_diffusion_get_residual),1);
			residual = sqrt(acDeviceGetOutput(acGridGetDevice(),AC_implicit_diffusion_residual));
			while(residual > error_tolerance && !over_max_steps)
			{
				acGridExecuteTaskGraph(acGetOptimizedDSLTaskGraph(implicit_diffusion_cg_step),1);
				acGridExecuteTaskGraph(acGetOptimizedDSLTaskGraph(implicit_diffusion_get_residual),1);
				residual = sqrt(acDeviceGetOutput(acGridGetDevice(),AC_implicit_diffusion_residual));
				++step_num;
				over_max_steps |= (max_steps > 0 && step_num >= max_steps);
			}
		}
		acDeviceSetInput(acGridGetDevice(),AC_implicit_diffusion_coefficient,nu);
		for(int field = 3; field < 6; ++field)
		{
			acDeviceSetInput(acGridGetDevice(),AC_CG_FIELD,CG_FIELD(field));
			AcReal residual = 1e100;
			acGridExecuteTaskGraph(acGetOptimizedDSLTaskGraph(implicit_diffusion_get_residual),1);
			residual = sqrt(acDeviceGetOutput(acGridGetDevice(),AC_implicit_diffusion_residual));
			bool over_max_steps = false;
			int step_num = 0;
			while(residual > error_tolerance && !over_max_steps)
			{
				acGridExecuteTaskGraph(acGetOptimizedDSLTaskGraph(implicit_diffusion_cg_step),1);
				acGridExecuteTaskGraph(acGetOptimizedDSLTaskGraph(implicit_diffusion_get_residual),1);
				residual = sqrt(acDeviceGetOutput(acGridGetDevice(),AC_implicit_diffusion_residual));
				++step_num;
				over_max_steps |= (max_steps > 0 && step_num >= max_steps);
			}
		}
	}
	else
	{
		fprintf(stderr,"Unknown method for implicit diffusion update on GPU!\n");
		fflush(stderr);
		exit(EXIT_FAILURE);
	}

#endif
}
/***********************************************************************************************/
extern "C" void registerGPU()
{
  // AcReal* profile_x_host = (AcReal*)malloc(sizeof(AcReal)*mx);
  // for (int i=0;i<mx;i++){
  //     profile_x_host[i] = (AcReal)i;
  //     acLogFromRootProc(rank,"profile_x_host[%d]=%f\n",i,profile_x_host[i]);
  // }
  // mesh.profiles[PROFILE_X] = profile_x_host;

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  mesh.info = acInitInfo();
#if AC_RUNTIME_COMPILATION
#else
  AcResult res = acCheckDeviceAvailability();

  if (res == AC_FAILURE) 
  {
	  fprintf(stderr,"No devices!\n");
	  exit(EXIT_FAILURE);
  }
#endif
}
/***********************************************************************************************/
// used as a flag to check if we need to calcualte the taus & uumean again
bool called_training = false;

// flags for printing if we are training or infering
bool calling_train = false;
bool calling_infer = false;

int randomNumber;
std::vector<float>val_loss;
std::vector<float>train_loss;
std::vector<double>train_time;
std::vector<double>val_time;

bool loaded_stats = false;
/***********************************************************************************************/
void denormalize(std::string filename, AcRealSymmetricTensor &tau_means, AcRealSymmetricTensor &tau_stds){
	if(!loaded_stats){
		loaded_stats = true;
		std::ifstream f(filename, std::ios::binary);
		if (!f.is_open()){
			fprintf(stderr, "Could not open stats file");
			fflush(stderr);
			exit(EXIT_FAILURE);
		}

		float acc_count;
		float num_acc;
		std::vector<float> acc_sum;
		std::vector<float> acc_sum_squared;
		
		while(f.peek() != EOF){
			uint32_t name_len;
			f.read(reinterpret_cast<char*>(&name_len), sizeof(name_len));
			if (f.eof()) break;
			std::string name;
			name.resize(name_len);
			f.read(&name[0], name_len);

			std::vector<int> shape;
			std::vector<float> data;

			uint32_t ndim = 1;
			f.read(reinterpret_cast<char*>(&ndim), sizeof(ndim));
			shape.resize(ndim);
			f.read(reinterpret_cast<char*>(shape.data()), ndim * sizeof(uint32_t));

			size_t num_elem = 1;
			for (auto s : shape) num_elem *=s;
			data.resize(num_elem);
			f.read(reinterpret_cast<char*>(data.data()), num_elem * sizeof(float));

			if (name == "acc_count"){
				acc_count = data[0];
			}
			
			if (name == "num_acc"){
				num_acc = data[0];
			}

			if(name == "acc_sum"){
				acc_sum = data;
			}

			if (name == "acc_sum_squared"){
				acc_sum_squared = data;
			}
		}
	
	float safe_count = std::max(acc_count, 1.0f);
	std::vector<float> means(acc_sum.size());
	for (int index = 0; index < acc_sum.size(); index++){
		means[index] = acc_sum[index] / safe_count;
	}
	std::vector<float> stds(acc_sum_squared.size());
	for (int index = 0; index < acc_sum.size(); index++){
		float var = (acc_sum_squared[index]/ safe_count) - (means[index] * means[index]);
		stds[index] = std::max(std::sqrt(var), 1e-8f);
	}

	tau_means.xx = means[0];
	tau_means.yy = means[1];
	tau_means.zz = means[2];
	tau_means.xy = means[3];
	tau_means.yz = means[4];
	tau_means.xz = means[5];


	tau_stds.xx = stds[0];
	tau_stds.yy = stds[1];
	tau_stds.zz = stds[2];
	tau_stds.xy = stds[3];
	tau_stds.yz = stds[4];
	tau_stds.xz = stds[5];

	/*
	std::cout << "printing acc_count: " << acc_count <<  "\n" << std::flush;

	std::cout << "printing num_acc: " << num_acc <<  "\n" << std::flush;

	std::cout << "printing acc_sum:" << std::flush;
	for (auto &i : acc_sum){
		std::cout << i << "";
	}
	std::cout << "\n" << std::flush;


	std::cout << "printing acc_sum_squared:" << std::flush;
	for (auto &i : acc_sum_squared){
		std::cout << i << "";
	}
	std::cout << "\n" << std::flush;
	*/
	


	}
}



/***********************************************************************************************/
extern "C" void torch_infer_c_api(int itstub){	
#if TRAINING
	#include "user_constants.h"
	if(!luse_trained_tau && itstub!=1) return;
	if(!calling_infer){
		AcRealSymmetricTensor tau_means = mesh.info[AC_tau_means];
		AcRealSymmetricTensor tau_stds  = mesh.info[AC_tau_stds];
		fprintf(stderr,"Calling infer\n");
		fprintf(stderr,"means xx: %f, yy: %f, zz: %f, xy: %f, yz: %f, xz: %f\n", tau_means.xx, tau_means.yy, tau_means.zz, tau_means.xy, tau_means.yz, tau_means.xz);
		fprintf(stderr,"stds xx: %f, yy: %f, zz: %f, xy: %f, yz: %f, xz: %f\n", tau_stds.xx, tau_stds.yy, tau_stds.zz, tau_stds.xy, tau_stds.yz, tau_stds.xz);
		fflush(stderr);
	}
	calling_infer = true;

	if (!called_training){
		randomNumber = 5;
		acDeviceSetInput(acGridGetDevice(), AC_ranNum, randomNumber);

		auto calc_uumean_tau = acGetOptimizedDSLTaskGraph(initialize_uumean_tau);
		acGridExecuteTaskGraph(calc_uumean_tau, 1);

  	auto bcs = acGetOptimizedDSLTaskGraph(boundconds);	
		acGridExecuteTaskGraph(bcs,1);

	}
		

	AcReal* out = NULL;

	AcReal* uumean_ptr = NULL;
	AcReal* tau_infer_ptr = NULL;
	
	/*
	if (randomNumber == 0) acDeviceGetVertexBufferPtrs(acGridGetDevice(), UUMEANBatch[0].x, &uumean_ptr, &out);
	if (randomNumber == 1) acDeviceGetVertexBufferPtrs(acGridGetDevice(), UUMEANBatch[1].x, &uumean_ptr, &out);
	if (randomNumber == 2) acDeviceGetVertexBufferPtrs(acGridGetDevice(), UUMEANBatch[2].x, &uumean_ptr, &out);
	if (randomNumber == 3) acDeviceGetVertexBufferPtrs(acGridGetDevice(), UUMEANBatch[3].x, &uumean_ptr, &out);
	if (randomNumber == 4) acDeviceGetVertexBufferPtrs(acGridGetDevice(), UUMEANBatch[4].x, &uumean_ptr, &out);
	if (randomNumber == 5) acDeviceGetVertexBufferPtrs(acGridGetDevice(), UUMEANBatch[5].x, &uumean_ptr, &out);
	*/

	acDeviceGetVertexBufferPtrs(acGridGetDevice(), uumean.x, &uumean_ptr, &out);

	acDeviceGetVertexBufferPtrs(acGridGetDevice(), TAU_INFERRED.xx, &tau_infer_ptr, &out);

	acGridHaloExchange();

	double start;
	double end;

	start = MPI_Wtime();

	torch_inferCAPI((int[]){mx,my,mz}, uumean_ptr, tau_infer_ptr);
	
	end = MPI_Wtime();
	auto descale_inf = acGetOptimizedDSLTaskGraph(descale_inferred_taus);
	acGridExecuteTaskGraph(descale_inf, 1);

  auto bcs = acGetOptimizedDSLTaskGraph(boundconds);	
	acGridExecuteTaskGraph(bcs,1);

	float vloss = MSE();
 	
	if(itstub == 1){
		//fprintf(stderr, "Validation error is: %.50f\n", vloss);
		//fprintf(stderr, "Validation took %f seconnds\n", (end-start));
		//fflush(stderr);
  	val_time.push_back((end-start));
  	val_loss.push_back(vloss);
		print_debug();
	}
#endif
}
/***********************************************************************************************/
extern "C" void torch_train_c_api(AcReal *loss_val, int itstub) {
#if TRAINING
	#include "user_constants.h"
	#include <stdlib.h>
	if(itstub != 1) return;

	if(!calling_train){
		fprintf(stderr,"Calling training\n");
		fflush(stderr);
	}

	//fprintf(stderr,"The value of it: %d", it);
	//fflush(stderr);

	called_training = true;

	calling_train = true;

	/*
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

	if (my_rank == 0){
		//fprintf(stderr, "The iteration number on c++ is: %d, ranNum is %d\n", it, randomNumber);	
		//fflush(stderr);
		MPI_Bcast(&randomNumber, 1, MPI_INT, 0, MPI_COMM_WORLD);
		MPI_Barrier(MPI_COMM_WORLD);
		acDeviceSetInput(acGridGetDevice(), AC_ranNum, randomNumber);
	}
	*/

	auto calc_uumean_tau = acGetOptimizedDSLTaskGraph(initialize_uumean_tau);
	acGridExecuteTaskGraph(calc_uumean_tau, 1);

  auto bcs = acGetOptimizedDSLTaskGraph(boundconds);	
  acGridExecuteTaskGraph(bcs,1);
  	

	AcReal* out = NULL;
  
  AcReal* uumean_ptr = NULL;
  AcReal* TAU_ptr = NULL;
  *loss_val = 0.1;
  
  acDeviceGetVertexBufferPtrs(acGridGetDevice(), tau.xx, &TAU_ptr, &out);
  acDeviceGetVertexBufferPtrs(acGridGetDevice(), uumean.x, &uumean_ptr, &out);
  
  acGridHaloExchange();
  
  double start, end;
  
  start = MPI_Wtime();
  
  torch_trainCAPI((int[]){mx,my,mz}, uumean_ptr, TAU_ptr, loss_val);

  end = MPI_Wtime();

  //fprintf(stderr,"Time for one time step: %f\n", (end-start));
  //fprintf(stderr,"Loss after training: %.7f\n", *loss_val);
  //fflush(stderr);
  train_time.push_back((end-start));
  train_loss.push_back(*loss_val);
  train_counter++;

	if (it==nt){
		double total_train_time = 0.0;
		double total_val_time = 0.0;

		double train_variance = 0.0;
		double val_variance = 0.0;

		for (int i=0;i<train_time.size();i++){
			total_train_time += train_time[i];
		}

		for (int i=0;i<val_time.size();i++){
			total_val_time += val_time[i];
		}
		
		
		float avg_train_loss = 0.0;
		float avg_val_loss = 0.0;

		for(int i=0;i<train_loss.size();i++){
			avg_train_loss += train_loss[i];
		}

		for(int i=0;i<val_loss.size();i++){
			avg_val_loss += val_loss[i];
		}

		avg_train_loss = avg_train_loss / train_loss.size();

		avg_val_loss = avg_val_loss / val_loss.size();

	
		for(int i=0;i<train_loss.size();i++){
			train_variance += (train_loss[i] - avg_train_loss) * (train_loss[i] - avg_train_loss);
		}

		for(int i=0;i<val_loss.size();i++){
			val_variance += (val_loss[i] - avg_val_loss) * (val_loss[i] - avg_val_loss);
		}
			
		train_variance = train_variance / train_loss.size();

		val_variance = val_variance / val_loss.size();
		
	
		fprintf(stderr,"The total time taken for training: %f.7\n", total_train_time);
		fprintf(stderr,"The average time taken for training: %f.7\n", total_train_time/train_time.size());

		fprintf(stderr,"The total time taken for validation: %f.7\n", total_val_time);
		fprintf(stderr,"The average time taken for validation: %f.7\n", total_val_time/val_time.size());

		
		fprintf(stderr,"Training loss variance: %f.7\n", train_variance);
		fprintf(stderr,"Validation loss variance: %f.7\n", val_variance);

		fflush(stderr);

		std::ofstream myFile;
		std::string fileString = "train_loss_" + std::to_string(my_rank)  + ".csv";	

		myFile.open(fileString);

    myFile << "epoch,train_loss\n";

		for (int i=0;i<train_loss.size();i++){
			myFile << i << "," << train_loss[i] << "\n";
		}
		
		myFile.close();


		fileString = "val_loss_" + std::to_string(my_rank)  + ".csv";	

		myFile.open(fileString);


    myFile << "epoch,val_loss\n";

		for (int i=0;i<val_loss.size();i++){
			myFile << i << "," << val_loss[i] << "\n";
		}

		myFile.close();
	}

#endif
}
/***********************************************************************************************/
float MSE(){
#if TRAINING
	#include "user_constants.h"

	auto calc_validation = acGetOptimizedDSLTaskGraph(calc_validation_loss);
	acGridExecuteTaskGraph(calc_validation, 1);

 	auto bcs = acGetOptimizedDSLTaskGraph(boundconds);	
	acGridExecuteTaskGraph(bcs,1);

	calculated_coeff_scales = true;
	return (acDeviceGetOutput(acGridGetDevice(), AC_l2_sum))/(6*nxgrid*nygrid*nzgrid);
#else
        return 0;
#endif
}
/***********************************************************************************************/
void load_f_ode()
{
        //TP: this is simply the initial implementation
        //TODO: benchmark what is the most efficient way of getting ode array to the GPU each substep
#if TRANSPILATION
        if (n_odevars > 0)
        {
                acDeviceSynchronizeStream(acGridGetDevice(),STREAM_DEFAULT);
                acDeviceLoad(acGridGetDevice(), STREAM_DEFAULT, mesh.info, AC_f_ode);
                acDeviceSynchronizeStream(acGridGetDevice(),STREAM_DEFAULT);
        }
#endif
}
/***********************************************************************************************/
extern "C" void beforeBoundaryGPU(bool lrmv, int isubstep, double t)
{
// Load values of ODE variables to GPU since before boundary may use them
	load_f_ode();

// Load those dynamical parameters which depend on time to GPU

 	acDeviceSetInput(acGridGetDevice(), AC_lrmv, lrmv);
 	acDeviceSetInput(acGridGetDevice(), AC_t, AcReal(t));

// Execute all "before-boundary-actions", which do not update the halos, by separate task graph

	acGridExecuteTaskGraph(acGetOptimizedDSLTaskGraph(AC_before_boundary_steps),1);

//Some Fields are directly calculated on the halos like yH in ioncalc.
//Could reformulate the kernels in a way that the bc is simply the same kernel as the normal calculation
//But don't want to repeat calc too often so this is an somewhat easy to way to do it
#if TRANSPILATION
	const auto steps_updating_halos = acGetOptimizedDSLTaskGraph(AC_before_boundary_steps_including_halos,
					          (Volume){0,0,0},acGetLocalMM(acGridGetLocalMeshInfo()));
	if (!acGridTaskGraphIsEmpty(steps_updating_halos))
	{
		AcTaskGraph* bcs = acGetOptimizedDSLTaskGraph(boundconds);
		acGridExecuteTaskGraph(bcs,1);                            // apply boundconds
		acGridHaloExchange();                                     // halo communication
		acGridExecuteTaskGraph(steps_updating_halos,1);           // f-array update
	}
#endif
#if LSELFGRAVITY
	if (t>=tstart_selfgrav)
	{
		acGridExecuteTaskGraph(acGetOptimizedDSLTaskGraph(AC_calc_selfgravity_rhs),1);
		acDeviceFFTR2C(acGridGetDevice(),acGetRHS_POISSON(),RHS_POISSON_COMPLEX);                     // FFT of rhs of Poisson eq.
    		AcMeshDims dims = acGetMeshDims(acGridGetLocalMeshInfo());
  		acGridLaunchKernel(STREAM_DEFAULT, selfgravity_poisson_solve, dims.n0, dims.n1);
		acGridSynchronizeStream(STREAM_ALL);
		acDeviceFFTC2R(acGridGetDevice(),SELFGRAVITY_POTENTIAL_COMPLEX,acGetSELFGRAVITY_POTENTIAL()); // FFT backtransform. of solution
		//TP: A placeholder for iterative solvers, in the future choose the number of solving steps based on the norm of the residual
		/**
		for (int i = 0; i < 100; ++i)
		{
			acGridExecuteTaskGraph(acGetOptimizedDSLTaskGraph(AC_sor_step),1);
		}
		**/
		acGridExecuteTaskGraph(acGetOptimizedDSLTaskGraph(AC_calc_final_potential),1);
	}

#endif
#if LNEWTON_COOLING
	acGridExecuteTaskGraph(acGetOptimizedDSLTaskGraph(AC_integrate_tau),1);
#endif
}
/***********************************************************************************************/
bool idx_init = false;
std::vector<size_t> idx_cache;
std::vector<double> buffer;


void print_debug() {
if (it % 5 !=0) return;
#if TRAINING
    #include "user_constants.h"
		
		std::string fname = "snapshots/snapshot_rank_" + std::to_string(my_rank) + "_it_" + std::to_string(it) + ".bin";
		std::ifstream infile(fname, std::ios::binary);
		if (infile.good()) return;
		
		counter = it;
		
		/*
		MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
		acDeviceSetInput(acGridGetDevice(), AC_ranNum, randomNumber);
		*/

	  acGridHaloExchange();
    copyFarray(NULL);

    AcMeshDims dims = acGetMeshDims(acGridGetLocalMeshInfo());

		/*
  	const auto DEVICE_VTXBUF_IDX = [&](const int x, const int y, const int z)
  				{
					return acVertexBufferIdx(x,y,z,mesh.info);
  				};
		*/


		int x_size = (dims.m1.x - dims.m0.x);
		int y_size = (dims.m1.y - dims.m0.y);
		int z_size = (dims.m1.z - dims.m0.z);

		const size_t n_points =  x_size * y_size * z_size;
		const int n_fields = 22;
		
		if(!idx_init){
			idx_cache.reserve(n_points);
			buffer.reserve(n_points * n_fields);

    	for (size_t i = dims.m0.x; i < dims.m1.x; i++) {
      	for (size_t j = dims.m0.y; j < dims.m1.y; j++) {
        	for (size_t k = dims.m0.z; k < dims.m1.z; k++) {
						idx_cache.push_back(acVertexBufferIdx(i,j,k,mesh.info));
					}
				}
			}
			idx_init = true;
		}
		
		buffer.clear();

		const AcReal* tau_xx_buf = mesh.vertex_buffer[tau.xx];
    const AcReal* tau_inferred_xx_buf = mesh.vertex_buffer[TAU_INFERRED.xx];
    const AcReal* tau_yy_buf = mesh.vertex_buffer[tau.yy];
    const AcReal* tau_inferred_yy_buf = mesh.vertex_buffer[TAU_INFERRED.yy];
    const AcReal* tau_zz_buf = mesh.vertex_buffer[tau.zz];
    const AcReal* tau_inferred_zz_buf = mesh.vertex_buffer[TAU_INFERRED.zz];
    const AcReal* tau_xy_buf = mesh.vertex_buffer[tau.xy];
    const AcReal* tau_inferred_xy_buf = mesh.vertex_buffer[TAU_INFERRED.xy];
    const AcReal* tau_yz_buf = mesh.vertex_buffer[tau.yz];
    const AcReal* tau_inferred_yz_buf = mesh.vertex_buffer[TAU_INFERRED.yz];
    const AcReal* tau_xz_buf = mesh.vertex_buffer[tau.xz];
    const AcReal* tau_inferred_xz_buf = mesh.vertex_buffer[TAU_INFERRED.xz];
    const AcReal* uumean_x_buf = mesh.vertex_buffer[uumean.x];
    const AcReal* uumean_y_buf = mesh.vertex_buffer[uumean.y];
    const AcReal* uumean_z_buf = mesh.vertex_buffer[uumean.z];
    const AcReal* uux_buf = mesh.vertex_buffer[UUX];
    const AcReal* uuy_buf = mesh.vertex_buffer[UUY];
    const AcReal* uuz_buf = mesh.vertex_buffer[UUZ];

    const double it_double = static_cast<double>(it);
    const int yz_size = y_size * z_size;

		for(size_t idx_i=0; idx_i < idx_cache.size(); ++idx_i){
			
			const size_t idx = idx_cache[idx_i];

			const int i = dims.m0.x + (idx_i / (yz_size));
			const int j = dims.m0.y + ((idx_i / z_size) % y_size);
			const int k = dims.m0.z + (idx_i % z_size);
			
			 buffer.push_back(it_double);

       buffer.push_back(tau_xx_buf[idx]);
       buffer.push_back(tau_inferred_xx_buf[idx]);

       buffer.push_back(tau_yy_buf[idx]);
       buffer.push_back(tau_inferred_yy_buf[idx]);

       buffer.push_back(tau_zz_buf[idx]);
       buffer.push_back(tau_inferred_zz_buf[idx]);

       buffer.push_back(tau_xy_buf[idx]);
       buffer.push_back(tau_inferred_xy_buf[idx]);

       buffer.push_back(tau_yz_buf[idx]);
       buffer.push_back(tau_inferred_yz_buf[idx]);

       buffer.push_back(tau_xz_buf[idx]);
       buffer.push_back(tau_inferred_xz_buf[idx]);


       buffer.push_back(uumean_x_buf[idx]);
       buffer.push_back(uumean_y_buf[idx]);
       buffer.push_back(uumean_z_buf[idx]);

       buffer.push_back(uux_buf[idx]);
       buffer.push_back(uuy_buf[idx]);
       buffer.push_back(uuz_buf[idx]);

       buffer.push_back(static_cast<double>(i));
       buffer.push_back(static_cast<double>(j));
       buffer.push_back(static_cast<double>(k));
		}
	

		//std::ostringstream fname;
    //fname << "snapshots/snapshot_rank_" + std::to_string(my_rank) + "_it_" + std::to_string(it) + ".bin";
    std::ofstream out(fname, std::ios::binary);
    out.write(reinterpret_cast<const char*>(buffer.data()), buffer.size() * sizeof(double));
    out.close();

	#endif
	if(called_training) called_training = false;
}
/***********************************************************************************************/
extern "C" void afterSubStepGPU()
{
	if (acDeviceGetInput(acGridGetDevice(), AC_step_num) == PC_FIRST_SUB_STEP)
	{
#if LGRAVITATIONAL_WAVES_HTXK
	    acDeviceFFTR2PlanarBatched(acGridGetDevice(), acGetF_STRESS_0(),acGetAC_tpq_re__mod__gravitational_waves_htxk_0(),acGetAC_tpq_im__mod__gravitational_waves_htxk_0(),6);
	    acGridExecuteTaskGraph(acGetOptimizedDSLTaskGraph(AC_gravitational_waves_solve_and_stress),1);
#endif
	}
	acGridExecuteTaskGraph(acGetOptimizedDSLTaskGraph(AC_after_timestep),1);
}
/***********************************************************************************************/
extern "C" void substepGPU(int isubstep, double t)
//
//  Do the 'isubstep'th integration step on all GPUs on the node and handle boundaries.
//
{
   //TP: logs performance metrics of Astaroth
   const bool log = false;
#if LFORCING
  //Update forcing params
   if (lsecond_force) 
   {
	   fprintf(stderr,"Second forcing force not yet implemented on GPU!\n");
	   exit(EXIT_FAILURE);
   }
   if (isubstep == num_substeps) forcing_params.Update();  // calculate on CPU and load into GPU
#endif
  acDeviceSetInput(acGridGetDevice(), AC_step_num,(PC_SUB_STEP_NUMBER) (isubstep-1));
  if (lshear) 
  {
	  acDeviceSetInput(acGridGetDevice(), AC_shear_delta_y, deltay);
  }
  Device dev = acGridGetDevice();
  //TP: done in this more complex manner to ensure the actually integrated time and the time reported by Pencil agree
  //if we call set_dt after the first timestep there would be slight shift in dt what Pencil sees and what is actually used for time integration
  
  if (isubstep == 1) 
  {
          static bool lfirst_timestep_calculated = false;
	  //TP: lcpu_timestep_on_gpu enables the same timestep as PC when testing
	  if (ldt && lcourant_dt && (!lfirst_timestep_calculated || lcpu_timestep_on_gpu)) 
	  {
		dt1_interface = GpuCalcDt(AcReal(t));
	  	lfirst_timestep_calculated = true;
	  }
	  if (ldt) set_dt(dt1_interface);
	  acDeviceSetInput(acGridGetDevice(), AC_dt,dt);
  }
  acDeviceSetInput(acGridGetDevice(), AC_t,(AcReal)t);
  //fprintf(stderr,"before acGridExecuteTaskGraph");
  AcTaskGraph *rhs =  acGetOptimizedDSLTaskGraph(AC_rhs);
  auto start = MPI_Wtime();
  acGridExecuteTaskGraph(rhs, 1);
  auto end = MPI_Wtime();
  if (log && !rank) fprintf(stderr,"RHS TOOK: %14e\n",end-start);
  if (ldt && (   (isubstep == 5 && !lcourant_dt) 
              || (isubstep == 1 &&  lcourant_dt)
             )
     )
  {
    constexpr AcReal unit = 1.0;
    AcReal dt1_;
    if (!lcourant_dt)
    {
      const AcReal maximum_error = acDeviceGetOutput(acGridGetDevice(), AC_maximum_error)/eps_rkf;
      AcReal dt_;
      const AcReal dt_increase=-unit/(itorder+dtinc);
      const AcReal dt_decrease=-unit/(itorder-dtdec);
      constexpr AcReal safety=(AcReal)0.95;
      if (maximum_error > 1)
      {
      	// Step above error threshold so decrease the next time step
      	const AcReal dt_temp = safety*dt*pow(maximum_error,dt_decrease);
      	// Don't decrease the time step by more than a factor of ten
        constexpr AcReal decrease_factor = (AcReal)0.1;
      	dt_ = sign(max(abs(dt_temp), decrease_factor*abs(dt)), dt);
      } 
      else
      {
      	dt_ = dt*pow(maximum_error,dt_increase);
      }
      set_next_dt(dt_);
      dt1_ = unit/dt_;
    }
    else 
    {
      dt1_ = calc_dt1_courant(AcReal(t));
    }
    dt1_interface = dt1_;
  }
  return;
}
/***********************************************************************************************/
void copyFarray(AcReal* f)
{
  #include "user_constants.h"

  acGridSynchronizeStream(STREAM_ALL);
  //TP: for now only copy the advanced fields back
  //TODO: should auxiliaries be needed on the GPU like e.g. Shock be copied? They can always be recomputed on the host if needed
  //If doing training we read all since we might want TAU components to calculate e.g. validation error
  const int all_vtxbufs = acGetNumFields();
  const int end = ltraining ? acGetNumFields(): 
	  	  lread_all_vars_from_device ? mfarray : mvar;

  AcMesh* dst = &mesh;
  AcMesh tmp;
  if (dimensionality <= 1 || (dimensionality == 2 && nygrid == 1))
  {
  	acHostMeshCopy(mesh, &tmp);
	dst = &tmp;
  }
  for (int i = 0; i < all_vtxbufs; ++i)
  {
	  //Have to specialize for training for now since we are reading fields that do not exist in Fortran: TODO: allocate fields on the CPU in Fortran
	  const int index = ltraining ? i : 
		  	    i < mvar ? i : maux_vtxbuf_index[i];
	  bool read_var = (i < end);
	  read_var &= (index != -1);
	  read_var |= (read_vtxbuf_from_gpu[i]);
	  read_var |= ltraining;
          if (!read_var) continue;
	  acDeviceStoreVertexBuffer(acGridGetDevice(),STREAM_DEFAULT,VertexBufferHandle(index),dst);
  }
  acGridSynchronizeStream(STREAM_ALL);
  //TP: Astaroth does not allocate ghost zones for inactive dimensions for 1d and 2d simulations, and unlike for xy  we cannot simply offset into the farray so have to manually copy the values
  //    This is fine since 1d simulations are anyways mainly for testing
  if (dimensionality == 1)
  {
  	for (int i = 0; i < end; ++i)
  	{
		if (i >= mvar && maux_vtxbuf_index[i] == -1) continue;
 		if (nxgrid != 1)
		{
			for (int x = 0; x < nx; ++x)
			{
				const size_t f_index = (NGHOST+x)+ mx*(NGHOST + my*(NGHOST));
				mesh.vertex_buffer[i][f_index] = dst->vertex_buffer[i][NGHOST+x];
			}
		}
 		if (nygrid != 1)
		{
			for (int y = 0; y < ny; ++y)
			{
				const size_t f_index = (NGHOST)+ mx*((y+NGHOST) + my*(NGHOST));
				mesh.vertex_buffer[i][f_index] = dst->vertex_buffer[i][NGHOST+y];
			}
		}
 		if (nzgrid != 1)
		{
			for (int z = 0; z < nz; ++z)
			{
				const size_t f_index = (NGHOST)+ mx*(NGHOST + my*(z+NGHOST));
				mesh.vertex_buffer[i][f_index] = dst->vertex_buffer[i][NGHOST+z];
			}
		}
	}
  	acHostMeshDestroy(&tmp);
  }
  if (dimensionality == 0)
  {
    	for (int i = 0; i < mfarray; ++i)
  	{
		const int index = (i < mvar) ? i : maux_vtxbuf_index[i];
		if (index == -1) continue;
		const size_t f_index = (NGHOST) + mx*(NGHOST + my*(NGHOST));
		const size_t ac_index = 0;
 		mesh.vertex_buffer[i][f_index] = dst->vertex_buffer[i][ac_index];
	}
  }
  if (dimensionality == 2 && nygrid == 1)
  {
    	for (int i = 0; i < mfarray; ++i)
  	{
		const int index = (i < mvar) ? i : maux_vtxbuf_index[i];
		if (index == -1) continue;
		for (int x = 0; x < nx; ++x)
		{
			for (int z = 0; z < nz; ++z)
			{
				const size_t f_index = (x+NGHOST) + mx*(NGHOST + my*(z+NGHOST));
				const size_t ac_index = (x+NGHOST) + mx*(z+NGHOST);
 				mesh.vertex_buffer[index][f_index] = dst->vertex_buffer[index][ac_index];
			}
		}
	}
  	acHostMeshDestroy(&tmp);
  }
}
/***********************************************************************************************/
void checkConfig(AcMeshInfo &config)
{
 acLogFromRootProc(rank,"Check that config is correct\n");
 acLogFromRootProc(rank,"d[xyz]: %.14f %.14f %.14f \n", dx, dy, dz);
// acLogFromRootProc(rank,"rank= %d: l1, l2, n1, n2, m1, m2= %d %d %d %d %d %d \n", rank, l1, l2, n1, n2, m1, m2);
// acLogFromRootProc(rank,"zlen= %.14f %.14f \n", config[AC_len].z, lxyz[2]);
 /*
#if LHYDRO
 acLogFromRootProc(rank,"lpressuregradientgas= %d %d \n", lpressuregradient_gas, config[AC_lpressuregradient_gas]);
#endif
#if LENTROPY
 acLogFromRootProc(rank,"chi= %f %f \n", chi, config[AC_chi]);
 acLogFromRootProc(rank,"nkramers= %f %f \n", nkramers, config[AC_nkramers]);
 acLogFromRootProc(rank,"hcond0_kramers= %f %f \n", hcond0_kramers, config[AC_hcond0_kramers]);
 acLogFromRootProc(rank,"hcond_Kconst= %f %f \n", hcond_Kconst, config[AC_hcond_Kconst]);
 //acLogFromRootProc(rank,"Fbot= %f %f \n", Fbot, config[AC_Fbot]);
 acLogFromRootProc(rank,"chi_t= %f %f \n", chi_t, config[AC_chi_t]);
#endif
#if LVISCOSITY
 acLogFromRootProc(rank,"nu= %f %f \n", nu, config[AC_nu]);
 acLogFromRootProc(rank,"zeta= %f %f \n", zeta, config[AC_zeta]);
#endif
#if LMAGNETIC
  acLogFromRootProc(rank,"eta= %f %f \n", eta, config[AC_eta]);
#endif
#if LEOS
  acLogFromRootProc(rank,"cs20= %f %f \n", cs20, config[AC_cs20]);
  //  acLogFromRootProc(rank,"gamma= %f %f \n", gamma, get_real_param(config,comp_infoAC_gamma));
  acLogFromRootProc(rank,"gamma_m1= %f %f \n", gamma_m1, config[AC_gamma_m1]);
  acLogFromRootProc(rank,"gamma1= %f %f \n", gamma1, config[AC_gamma1]);
  acLogFromRootProc(rank,"cv= %f %f \n", cv, config[AC_cv]);
  acLogFromRootProc(rank,"cp= %f %f \n", cp, config[AC_cp]);
  acLogFromRootProc(rank,"lnT0= %f %f \n", lnTT0, config[AC_lnTT0]);
  acLogFromRootProc(rank,"lnrho0= %f %f \n", lnrho0, config[AC_lnrho0]);
#endif
#if LFORCING
  acLogFromRootProc(rank,"iforcing_zsym= %f %f \n", iforcing_zsym, config[AC_iforcing_zsym]);
  acLogFromRootProc(rank,"k1_ff= %f %f \n", k1_ff, config[AC_k1_ff]);
  acLogFromRootProc(rank,"tforce_stop= %f %f \n", tforce_stop, config[AC_tforce_stop]);
  acLogFromRootProc(rank,"k1_ff,profx_ampl, val= %f %d %lf %lf\n", k1_ff, profx_ampl, profx_ampl[0], profx_ampl[nx-1]);
#endif
  acLogFromRootProc(rank,"mu0= %f %f \n", (double)mu0, (double)config[AC_mu0]);
*/
}
/***********************************************************************************************/
extern "C" void getFArrayIn(AcReal **p_f_in)
{
  return;   //!!!
  AcReal* out = NULL;

  AcReal* uux_ptr = NULL;
  AcReal* uuy_ptr = NULL;
  AcReal* uuz_ptr = NULL;

  acDeviceGetVertexBufferPtrs(acGridGetDevice(),UUX,&uux_ptr,&out);
  acDeviceGetVertexBufferPtrs(acGridGetDevice(),UUY,&uuy_ptr,&out);
  acDeviceGetVertexBufferPtrs(acGridGetDevice(),UUZ,&uuz_ptr,&out);
  if (uux_ptr + mw != uuy_ptr) fprintf(stderr, "UU not contiguous\n");
  if (uuy_ptr + mw != uuz_ptr) fprintf(stderr, "UU not contiguous\n");
  acDeviceGetVertexBufferPtrs(acGridGetDevice(),VertexBufferHandle(0),p_f_in,&out);
}
/***********************************************************************************************/
extern "C" void copyVBApointers(AcReal **in, AcReal **out)
{
  AcReal* uux_ptr = NULL;
  AcReal* uuy_ptr = NULL;
  AcReal* uuz_ptr = NULL;

  acDeviceGetVertexBufferPtrs(acGridGetDevice(),UUX,&uux_ptr,out);
  acDeviceGetVertexBufferPtrs(acGridGetDevice(),UUY,&uuy_ptr,out);
  acDeviceGetVertexBufferPtrs(acGridGetDevice(),UUZ,&uuz_ptr,out);
  if (uux_ptr + mw != uuy_ptr) fprintf(stderr, "UU not contiguous\n");
  if (uuy_ptr + mw != uuz_ptr) fprintf(stderr, "UU not contiguous\n");

  acDeviceGetVertexBufferPtrs(acGridGetDevice(),VertexBufferHandle(0),in,out);
}
/***********************************************************************************************/
void autotune_all_integration_substeps()
{
  for (int i = 0; i < num_substeps; ++i)
  {
  	acDeviceSetInput(acGridGetDevice(), AC_step_num,(PC_SUB_STEP_NUMBER)i);
        if (rank==0 && ldebug) printf("memusage before GetOptimizedDSLTaskGraph= %f MBytes\n", acMemUsage()/1024.);
  	acDeviceSetInput(acGridGetDevice(), AC_lrmv,false);
	acGetOptimizedDSLTaskGraph(AC_rhs);
  	acDeviceSetInput(acGridGetDevice(), AC_lrmv,true);
	acGetOptimizedDSLTaskGraph(AC_rhs);
        if (rank==0 && ldebug) printf("memusage after GetOptimizedDSLTaskGraph= %f MBytes\n", acMemUsage()/1024.);
 	beforeBoundaryGPU(true,i,0.0);
        beforeBoundaryGPU(false,i,0.0);
  }
  sourceFunctionAndOpacity(0);
  splitUpdate(1e11,1);
}
/***********************************************************************************************/
extern "C" void loadFarray()
{
  AcMesh src = mesh;
  AcMesh tmp;
  if (dimensionality == 1)
  {
  	acHostMeshCopy(mesh, &tmp);
	src = tmp;
    	for (int i = 0; i < mfarray; ++i)
  	{
		const int index = (i < mvar) ? i : maux_vtxbuf_index[i];
		if (index == -1) continue;
 		if (nxgrid != 1)
		{
			for (int x = 0; x < nx; ++x)
			{
				const size_t f_index = (NGHOST+x)+ mx*(NGHOST + my*(NGHOST));
				src.vertex_buffer[index][NGHOST+x]  = mesh.vertex_buffer[index][f_index];
			}
		}
		else if (nygrid != 1)
		{
			for (int y = 0; y < ny; ++y)
			{
				const size_t f_index = (NGHOST)+ mx*((NGHOST+y) + my*(NGHOST));
				src.vertex_buffer[index][NGHOST+y]  = mesh.vertex_buffer[index][f_index];
			}
		}
		else if (nzgrid != 1)
		{
			for (int z = 0; z < nz; ++z)
			{
				const size_t f_index = (NGHOST)+ mx*(NGHOST + my*(z+NGHOST));
				src.vertex_buffer[index][NGHOST+z]  = mesh.vertex_buffer[index][f_index];
			}
		}
	}
  }
  if (dimensionality == 2 && nygrid == 1)
  {
  	acHostMeshCopy(mesh, &tmp);
	src = tmp;
    	for (int i = 0; i < mfarray; ++i)
  	{
		const int index = (i < mvar) ? i : maux_vtxbuf_index[i];
		if (index == -1) continue;
		for (int x = 0; x < nx; ++x)
		{
			for (int z = 0; z < nz; ++z)
			{
				const size_t f_index = (x+NGHOST) + mx*(NGHOST + my*(z+NGHOST));
				const size_t ac_index = (x+NGHOST) + mx*(z+NGHOST);
				src.vertex_buffer[index][ac_index] = mesh.vertex_buffer[index][f_index];
			}
		}
	}
  }
  if (dimensionality == 0)
  {
  	acHostMeshCopy(mesh, &tmp);
	src = tmp;
    	for (int i = 0; i < mfarray; ++i)
  	{
		const int index = (i < mvar) ? i : maux_vtxbuf_index[i];
		if (index == -1) continue;
		const size_t f_index = (NGHOST) + mx*(NGHOST + my*(NGHOST));
		const size_t ac_index = 0;
		src.vertex_buffer[index][ac_index] = mesh.vertex_buffer[index][f_index];
	}
  }
  acGridSynchronizeStream(STREAM_ALL);
  {
    for (int i = 0; i < mvar; ++i)
  	acDeviceLoadVertexBuffer(acGridGetDevice(), STREAM_DEFAULT, src, VertexBufferHandle(i));

    int n_aux_on_gpu = 0;
    for (int i = 0; i < mfarray; ++i)
      if (maux_vtxbuf_index[i] != -1)
      {
	      n_aux_on_gpu++;
  		acDeviceLoadVertexBuffer(acGridGetDevice(), STREAM_DEFAULT, src, VertexBufferHandle(maux_vtxbuf_index[i]));
      }
    for (int i = 0; i < mfarray-mvar-maux; ++i)
    {
  	acDeviceLoadVertexBuffer(acGridGetDevice(), STREAM_DEFAULT, src, VertexBufferHandle(mvar+n_aux_on_gpu+i));
    }
  }
  acGridSynchronizeStream(STREAM_ALL);
  if (dimensionality == 1) acHostMeshDestroy(&tmp);
}
/***********************************************************************************************/
void generate_bcs()
{
	if(rank != 0) return;
	if(system("cd src && scripts/bc2ast 1> ../tmp_bcs 2> /dev/null && cd .."))
	{
		fprintf(stderr,"AC Error: Was not able to generate bcs!\n");
		exit(EXIT_FAILURE);
	}
	const int different = system("diff tmp_bcs src/astaroth/DSL/local/boundconds.h");
	if(different)
	{
		if(system("mv tmp_bcs src/astaroth/DSL/local/boundconds.h"))
		{
			fprintf(stderr,"AC Error: Was not able move bcs to DSL/local!\n");
			exit(EXIT_FAILURE);
		}
		fprintf(stderr,"BCs different: recompiling!\n");
	}
	else
	{
		if(system("rm tmp_bcs"))
		{
			fprintf(stderr,"AC Error: Was not able to remove tmp_bcs!\n");
			exit(EXIT_FAILURE);
		}
	}
}
/***********************************************************************************************/
extern "C" void testBCs();     // forward declaration
/***********************************************************************************************/
extern "C" void initializeGPU(AcReal *farr, int comm_fint, double t, int nt_)  // MPI_Fint comm_fint
{
  //Setup configurations used for initializing and running the GPU code
  nt = nt_;
  comm_pencil = MPI_Comm_f2c(comm_fint);
  setupConfig(mesh.info);
/**
#if TRAINING
  #include "user_constants.h"
  //fprintf(stderr,"INDEX OF TAU.xx: %d\n",acGetTAU_Xx());
	{
		#include "user_constants.h"
		if (itauxx-1 != TAU.xx)
		{
			fprintf(stderr,"Mismatch of indeces for tauxx : %d,%d!!\n",itauxx,TAU.xx);
			exit(EXIT_FAILURE);
		}
		if (itauxy-1 != TAU.xy)
		{
			fprintf(stderr,"Mismatch of indeces for tauxy !!\n");
			exit(EXIT_FAILURE);
		}
		if (itauxz-1 != TAU.xz)
		{
			fprintf(stderr,"Mismatch of indeces for tauxz !!\n");
			exit(EXIT_FAILURE);
		}
		if (itauyy-1 != TAU.yy)
		{
			fprintf(stderr,"Mismatch of indeces for tauyy!!\n");
			exit(EXIT_FAILURE);
		}
		if (itauyz-1 != TAU.yz)
		{
			fprintf(stderr,"Mismatch of indeces for tauyz !!\n");
			exit(EXIT_FAILURE);
		}
		if (itauzz-1 != TAU.zz)
		{
			fprintf(stderr,"Mismatch of indeces for tauzz !!\n");
			exit(EXIT_FAILURE);
		}
	}
#endif
	**/

  if (rank==0 && ldebug) printf("memusage after pointer assign= %f MBytes\n", acMemUsage()/1024.);
#if AC_RUNTIME_COMPILATION
#include "cmake_options.h"

  //Not worth it to get this working inside the container
  const bool inside_container = ltraining;
  if(!inside_container) generate_bcs();
  MPI_Barrier(MPI_COMM_WORLD);
  acCompile(cmake_options,mesh.info);
  acLoadLibrary(rank == 0 ? stderr : NULL,mesh.info);
  acCheckDeviceAvailability();
  acLogFromRootProc(rank, "Done setupConfig && acCompile\n");
#else
  acLogFromRootProc(rank, "Done setupConfig\n");
#endif
  fflush(stdout);
  //TP: done after setupConfig and acCompile since we need maux_vtxbuf_index and acGetNumFields
  //TP: this is an ugly way to do this but works for now
  {
    const size_t z_offset  = (dimensionality == 2 && nzgrid == 1) ? NGHOST*mx*my : 0;
    for (int i = 0; i < mvar; ++i)
    {
      mesh.vertex_buffer[VertexBufferHandle(i)] = &farr[mw*i+ z_offset];
    }

    int n_aux_on_gpu = 0;
    for (int i = 0; i < mfarray; ++i)
    {
      if (maux_vtxbuf_index[i] != -1)
      {
	++n_aux_on_gpu;
        mesh.vertex_buffer[maux_vtxbuf_index[i]] = &farr[mw*i + z_offset];
      }
    }
    for (int i = 0; i < mfarray-mvar-maux; ++i)
    {
        mesh.vertex_buffer[mvar+n_aux_on_gpu+i] = &farr[mw*(mvar+maux+i) + z_offset];
    }
    //TP: for now for training we have all slots filled since we might want to read TAU components to the host for calculating validation error
    if (ltraining)
    {
    	for (int i = 0; i < acGetNumFields(); ++i)
    	{
	   if (mesh.vertex_buffer[i] == NULL)
	   {
    	    	mesh.vertex_buffer[i] = (AcReal*)malloc(sizeof(AcReal)*mw);
	   }
    	}
    }
  }
  checkConfig(mesh.info);
  if (rank==0 && ldebug) printf("memusage grid_init= %f MBytes\n", acMemUsage()/1024.);
  acGridInit(mesh);
  if (rank==0 && ldebug) printf("memusage after grid_init= %f MBytes\n", acMemUsage()/1024.);

  mesh.info = acGridDecomposeMeshInfo(mesh.info);
  //TP: important to do before autotuning
  acDeviceSetInput(acGridGetDevice(), AC_step_num,(PC_SUB_STEP_NUMBER)0);
  acDeviceSetInput(acGridGetDevice(), AC_dt,dt);
  acDeviceSetInput(acGridGetDevice(), AC_t,AcReal(t));
  acDeviceSetInput(acGridGetDevice(), AC_shear_delta_y,deltay);
		
  //if (ltest_bcs) testBCs();
  //TP: for autotuning
  afterSubStepGPU();
  autotune_all_integration_substeps();
  if (rank==0 && ldebug) printf("memusage before store config= %f MBytes\n", acMemUsage()/1024.);
  acStoreConfig(acDeviceGetLocalConfig(acGridGetDevice()), "PC-AC.conf");
  if (rank==0 && ldebug) printf("memusage after store config= %f MBytes\n", acMemUsage()/1024.);
  acGridSynchronizeStream(STREAM_ALL);
  if (rank==0 && ldebug) printf("memusage after store synchronize stream= %f MBytes\n", acMemUsage()/1024.);
  acLogFromRootProc(rank, "DONE initializeGPU\n");
  fflush(stdout);

  constexpr AcReal unit = 1.0;
  dt1_interface = unit/dt;
}
/***********************************************************************************************/
extern "C" void reloadConfig()
{
  setupConfig(mesh.info);
#if AC_RUNTIME_COMPILATION
  //save the current values of the vtxbufs since the device arrays are freed by acGridQuit
  copyFarray(mesh.vertex_buffer[0]);
  acGridQuit();
  const AcResult closed_res = acCloseLibrary();
  if (closed_res != AC_SUCCESS)
  {
	  if (rank == 0) fprintf(stderr,"Was not successful in closing Astaroth lib!\n");
	  exit(EXIT_FAILURE);
  }
#include "cmake_options.h"
  generate_bcs();
  MPI_Barrier(MPI_COMM_WORLD);
  acCompile(cmake_options,mesh.info);
  acLoadLibrary(rank == 0 ? stderr : NULL,mesh.info);
  acGridInit(mesh);
  acLogFromRootProc(rank, "Done setupConfig && acCompile\n");
  fflush(stdout);
  fflush(stderr);
  //TP: this is important that we don't overwrite the output buffer in middle of a timestep when the output buffer holds some meaning!
  autotune_all_integration_substeps();
  //TP: restore the vtxbuf values before quitting grid
  loadFarray();
#else
  acGridSynchronizeStream(STREAM_ALL);
  acDeviceUpdate(acGridGetDevice(), mesh.info);
  acGridSynchronizeStream(STREAM_ALL);
#endif
  acLogFromRootProc(rank, "DONE reloading on GPU\n");
  fflush(stdout);
  fflush(stderr);
}
/***********************************************************************************************/
extern "C" void updateInConfigArr(int index)
{
     if (mesh.info[(AcRealArrayParam)index] != nullptr)
        acDeviceLoadRealArray(acGridGetDevice(),STREAM_DEFAULT,mesh.info,static_cast<AcRealArrayParam>(index));
}
/***********************************************************************************************/
extern "C" void updateInConfigScal(int index, AcReal value)
{
     acDeviceLoadScalarUniform(acGridGetDevice(),STREAM_DEFAULT,static_cast<AcRealParam>(index),value);
     acDeviceSynchronizeStream(acGridGetDevice(),STREAM_DEFAULT);
}
/***********************************************************************************************/
extern "C" int updateInConfigArrName(char *name)
{
    int index = -1;
    for (int i=0; i<NUM_REAL_ARRAYS; i++){
       if (strcmp(get_array_info(static_cast<AcRealArrayParam>(i)).name,name)==0) index=i;
    }
    if (index>-1) updateInConfigArr(index);
    else
    {
       fprintf(stderr,"Astaroth WARNING: Did not entry named %s in config!!\n",name);
       fflush(stderr);
    }
    return index;
}
/***********************************************************************************************/
extern "C" int updateInConfigScalName(char *name, AcReal value)
{
    int index = -1;
    for (int i=0; i<NUM_REAL_PARAMS; i++){
       if (strcmp(realparam_names[i],name)==0) index=i;
    }
    if (index>-1) updateInConfigScal(index, value);
    else
    {
       fprintf(stderr,"Astaroth WARNING: Did not entry named %s in config!!\n",name);
       fflush(stderr);
    }
    return index;
}
/**********************************************************************************************/
extern "C" void finalizeGPU()
{
  // Deallocate everything on the GPUs and reset
  AcResult res = acGridQuit();
}
/***********************************************************************************************/
extern "C" void random_initial_condition()
{
  return;
  //acGridSynchronizeStream(STREAM_ALL);
  //AcMeshDims dims = acGetMeshDims(acGridGetLocalMeshInfo());
  //acGridLaunchKernel(STREAM_DEFAULT, randomize, dims.n0, dims.n1);
  //acGridSynchronizeStream(STREAM_ALL);
}
/***********************************************************************************************/
//TP: this is not written the most optimally since it needs two extra copies of the mesh where at least the tmp
//could be circumvented by temporarily using the output buffers on the GPU to store the f-array and load back from there
//but if we truly hit the mem limit for now the user can of course simply test the bcs with a smaller mesh and skip the test with a larger mesh
/***********************************************************************************************/
void sym_z(AcMesh mesh_in)
{
  const auto DEVICE_VTXBUF_IDX = [&](const int x, const int y, const int z)
  				{
					return acVertexBufferIdx(x,y,z,acGridGetLocalMeshInfo());
  				};
  AcMeshDims dims = acGetMeshDims(acGridGetLocalMeshInfo());
  for (size_t i = 0; i < dims.m1.x; i++)
  {
    for (size_t j = 0; j < dims.m1.y; j++)
    {
      for (size_t k = 0; k < dims.m1.z; k++)
      {
	if (
	   i >= NGHOST && i < dims.n1.x &&
	   j >= NGHOST && j < dims.n1.y &&
	   k >= NGHOST && k < dims.n1.z
	   ) continue;
	if (k >= NGHOST && k < dims.n1.z) continue;
        for (int ivar = 0; ivar < acGetNumFields(); ivar++)
        {
	  //BOT
	  if (k < NGHOST)
	  {
	  	const auto idx = DEVICE_VTXBUF_IDX(i,j,k);
	        const auto offset = NGHOST-k;
	  	const auto domain_z= NGHOST+offset;
	  	const auto domain_idx = DEVICE_VTXBUF_IDX(i,j,domain_z);
	        mesh_in.vertex_buffer[ivar][idx] = mesh_in.vertex_buffer[ivar][domain_idx];
	        //mesh_in.vertex_buffer[ivar][idx] = mesh.vertex_buffer[ivar][idx];
	  }
	  //TOP:
	  else
	  {
	  	const auto idx = DEVICE_VTXBUF_IDX(i,j,k);
	        const auto offset = k-mesh_in.info[AC_nlocal_max].z+1;
	  	const auto domain_z= mesh_in.info[AC_nlocal_max].z-offset;
	  	const auto domain_idx = DEVICE_VTXBUF_IDX(i,j,domain_z);
	        mesh_in.vertex_buffer[ivar][idx] = mesh_in.vertex_buffer[ivar][domain_idx];
	        //mesh_in.vertex_buffer[ivar][idx] = mesh.vertex_buffer[ivar][idx];
	  }
	}
      }
    }
  }
}
/***********************************************************************************************/
void check_sym_z(AcMesh mesh_in)
{
  const auto DEVICE_VTXBUF_IDX = [&](const int x, const int y, const int z)
  				{
					return acVertexBufferIdx(x,y,z,acGridGetLocalMeshInfo());
  				};
  if (rank == 1)
  {
  	//printf("HMM: %14e\n",mesh_in.vertex_buffer[0][DEVICE_VTXBUF_IDX(14,15,2)]);
  	//printf("HMM: %14e\n",mesh_in.vertex_buffer[0][DEVICE_VTXBUF_IDX(14,15,2)]);
  	//printf("HMM: %14e\n",mesh_in.vertex_buffer[0][DEVICE_VTXBUF_IDX(26,14,1)]);
  	//printf("HMM: %14e\n",mesh_in.vertex_buffer[0][DEVICE_VTXBUF_IDX(26,14,5)]);
  }
  AcMeshDims dims = acGetMeshDims(acGridGetLocalMeshInfo());
  for (size_t i = 0; i < dims.m1.x; i++)
  {
    for (size_t j = 0; j < dims.m1.y; j++)
    {
      for (size_t k = 0; k < dims.m1.z; k++)
      {
	if (
	   i >= NGHOST && i < dims.n1.x &&
	   j >= NGHOST && j < dims.n1.y &&
	   k >= NGHOST && k < dims.n1.z
	   ) continue;
	if (k >= NGHOST && k < dims.n1.z) continue;
        for (int ivar = 0; ivar < acGetNumFields(); ivar++)
        {
	  //BOT
	  if (k < NGHOST)
	  {
	  	const auto idx = DEVICE_VTXBUF_IDX(i,j,k);
	        const auto offset = NGHOST-k;
	  	const auto domain_z= NGHOST+offset;
	  	const auto domain_idx = DEVICE_VTXBUF_IDX(i,j,domain_z);
	        if (mesh_in.vertex_buffer[ivar][idx] !=  mesh_in.vertex_buffer[ivar][domain_idx])
	         	printf("WRONG\n");
	  }
	  //TOP:
	  else
	  {
	  	const auto idx = DEVICE_VTXBUF_IDX(i,j,k);
	        const auto offset = k-mesh_in.info[AC_nlocal_max].z+1;
	  	const auto domain_z= mesh_in.info[AC_nlocal_max].z-offset;
	  	const auto domain_idx = DEVICE_VTXBUF_IDX(i,j,domain_z);
	        mesh_in.vertex_buffer[ivar][idx] = mesh.vertex_buffer[ivar][domain_idx];
	        if (mesh_in.vertex_buffer[ivar][idx] !=  mesh_in.vertex_buffer[ivar][domain_idx])
	         	printf("WRONG\n");
	  }
	}
      }
    }
  }
}
/***********************************************************************************************/
void check_sym_x(const AcMesh mesh_in)
{
  const auto DEVICE_VTXBUF_IDX = [&](const int x, const int y, const int z)
  				{
					return acVertexBufferIdx(x,y,z,acGridGetLocalMeshInfo());
  				};
  if (rank == 1)
  {
  	//printf("HMM: %14e\n",mesh_in.vertex_buffer[0][DEVICE_VTXBUF_IDX(14,15,2)]);
  	//printf("HMM: %14e\n",mesh_in.vertex_buffer[0][DEVICE_VTXBUF_IDX(14,15,2)]);
  	//printf("HMM: %14e\n",mesh_in.vertex_buffer[0][DEVICE_VTXBUF_IDX(26,14,1)]);
  	//printf("HMM: %14e\n",mesh_in.vertex_buffer[0][DEVICE_VTXBUF_IDX(26,14,5)]);
  }
  AcMeshDims dims = acGetMeshDims(acGridGetLocalMeshInfo());
  for (size_t i = 0; i < dims.m1.x; i++)
  {
    for (size_t j = 0; j < dims.m1.y; j++)
    {
      for (size_t k = 0; k < dims.m1.z; k++)
      {
	if (
	   i >= NGHOST && i < dims.n1.x &&
	   j >= NGHOST && j < dims.n1.y &&
	   k >= NGHOST && k < dims.n1.z
	   ) continue;
	if (i >= NGHOST && i < dims.n1.x) continue;
        for (int ivar = UUY; ivar <= UUY; ivar++)
        {
	  //BOT
	  if (i < NGHOST)
	  {
	  	const auto idx = DEVICE_VTXBUF_IDX(i,j,k);
	 	const auto offset = NGHOST-i;
	  	const auto domain_x= NGHOST+offset;
	  	const auto domain_idx = DEVICE_VTXBUF_IDX(domain_x,j,k);
	 	if (mesh_in.vertex_buffer[ivar][idx] !=  mesh_in.vertex_buffer[ivar][domain_idx])
	 	{
	 		fprintf(stderr,"WRONG\n");
	 		exit(EXIT_FAILURE);
	 	}
	  }
	  //TOP:
	  else
	  {
	  	const auto idx = DEVICE_VTXBUF_IDX(i,j,k);
	 	const auto offset = i-mesh_in.info[AC_nlocal_max].x+1;
	  	const auto domain_x= mesh_in.info[AC_nlocal_max].x-offset;
	  	const auto domain_idx = DEVICE_VTXBUF_IDX(domain_x,j,k);
	 	if (mesh_in.vertex_buffer[ivar][idx] !=  mesh_in.vertex_buffer[ivar][domain_idx])
	 	{
	 		fprintf(stderr,"WRONG\n");
	 		exit(EXIT_FAILURE);
	 	}
	  }
	}
      }
    }
  }
}
/***********************************************************************************************/
void testBCs()
{
  //if (dimensionality != 3) return;
  // Set random seed for reproducibility
  srand(321654987);
  const auto DEVICE_VTXBUF_IDX = [&](const int x, const int y, const int z)
  				{
					return acVertexBufferIdx(x,y,z,acGridGetLocalMeshInfo());
  				};
  AcTaskGraph* rhs = acGetOptimizedDSLTaskGraph(AC_rhs);
  const bool all_periodic = acGridTaskGraphHasPeriodicBoundcondsX(rhs) && 
  			    acGridTaskGraphHasPeriodicBoundcondsY(rhs) && 
  			    acGridTaskGraphHasPeriodicBoundcondsZ(rhs); 
  if (all_periodic && !lshear) return;
  if (lshear) acLogFromRootProc(rank,"testBCS: deltay: %7e\n",deltay);
  auto bcs = acGetDSLTaskGraph(boundconds);

  AcMesh tmp_mesh_to_store;
  acHostMeshCopy(mesh, &tmp_mesh_to_store);
  
  acHostMeshRandomize(&mesh);

  AcMesh mesh_to_copy;
  acHostMeshCopy(mesh, &mesh_to_copy);

  acGridSynchronizeStream(STREAM_ALL);
  acDeviceLoadMesh(acGridGetDevice(), STREAM_DEFAULT, mesh);
  acGridSynchronizeStream(STREAM_ALL);

  int ivar1 = 1;
  int ivar2 = mvar;
  boundconds_x_c(mesh.vertex_buffer[0],&ivar1,&ivar2);
  boundconds_y_c(mesh.vertex_buffer[0],&ivar1,&ivar2);
  boundconds_z_c(mesh.vertex_buffer[0],&ivar1,&ivar2);

  acGridSynchronizeStream(STREAM_ALL);
  acGridExecuteTaskGraph(bcs,1);
  acGridSynchronizeStream(STREAM_ALL);

  acGridSynchronizeStream(STREAM_ALL);
  acDeviceStoreMesh(acGridGetDevice(), STREAM_DEFAULT, &mesh_to_copy);
  acGridSynchronizeStream(STREAM_ALL);

  //check_sym_x(mesh);
  //check_sym_x(mesh_to_copy);

  AcMeshDims dims = acGetMeshDims(acGridGetLocalMeshInfo());
  bool passed = true;
  AcReal max_abs_not_passed_val=-1.0;
  AcReal true_pair{};
  AcReal max_abs_relative_difference =-1.0;
  AcReal max_abs_value = -1.0;
  AcReal min_abs_value = 1.0;
  AcReal gpu_val_for_largest_diff{};
  AcReal true_val_for_largest_diff{};
  AcReal epsilon = (AcReal)pow(10.0,-12.0);
  int3 largest_diff_point{};
  //AcReal epsilon = 0.0;

  auto skip_x = nxgrid == 1 || acGridTaskGraphHasPeriodicBoundcondsX(rhs);
  auto skip_y = nygrid == 1 || acGridTaskGraphHasPeriodicBoundcondsY(rhs);
  auto skip_z = nzgrid == 1 || acGridTaskGraphHasPeriodicBoundcondsZ(rhs);

  auto start_x =  skip_x ? NGHOST : 0;
  auto start_y =  skip_y ? NGHOST : 0;
  auto start_z =  skip_z ? NGHOST : 0;

  auto end_x =  skip_x ? dims.n1.x : dims.m1.x;
  auto end_y =  skip_y ? dims.n1.y : dims.m1.y;
  auto end_z =  skip_z ? dims.n1.z : dims.m1.z;

  int num_of_points_where_different[NUM_VTXBUF_HANDLES]{};
  bool different_in[NUM_VTXBUF_HANDLES][27]{};
  //TP: stupid levels of safety
  memset(different_in,0,NUM_VTXBUF_HANDLES*27*sizeof(bool));
  int num_of_points = 0;

  for (size_t i = start_x; i < end_x; i++)
  {
    for (size_t j = start_y; j < end_y; j++)
    {
      for (size_t k = start_z; k < end_z; k++)
      {
	const int x_idx = i < NGHOST ? -1 : i >= dims.n1.x ? 1 : 0;
	const int y_idx = j < NGHOST ? -1 : j >= dims.n1.y ? 1 : 0;
	const int z_idx = k < NGHOST ? -1 : k >= dims.n1.z ? 1 : 0;

	const bool x_in = i >= NGHOST && i < dims.n1.x;
	const bool y_in = j >= NGHOST && j < dims.n1.y;
	const bool z_in = k >= NGHOST && k < dims.n1.z;
	if (
		x_in && y_in && z_in
	  ) continue;
	  //if (i < NGHOST | i > dims.n1.x  || j < NGHOST || j > dims.n1.y) continue;
	++num_of_points;
        for (int ivar = 0; ivar < mvar; ivar++)
        {

	  const int index = (x_idx+1) + 3*((y_idx+1) + 3*(z_idx+1));
          AcReal out_val = mesh_to_copy.vertex_buffer[ivar][DEVICE_VTXBUF_IDX(i, j, k)];
          AcReal true_val = mesh.vertex_buffer[ivar][DEVICE_VTXBUF_IDX(i, j, k)];
          AcReal abs_diff = fabs(out_val - true_val);
          if (fabs(true_val) > max_abs_value) max_abs_value = fabs(true_val);
          if (fabs(true_val) < min_abs_value) min_abs_value = fabs(true_val);
          if (fabs(abs_diff/true_val) > epsilon || (true_val == (AcReal)0.0 && fabs(out_val) > (AcReal)pow(0.1,13)) || (epsilon == (AcReal)0.0 && true_val != out_val))
          {
	    different_in[ivar][index] |= true;
            passed = false;
            num_of_points_where_different[ivar]++;
            if (max_abs_not_passed_val<abs(out_val)){
              max_abs_not_passed_val = abs(out_val);
              true_pair = true_val;
            }
            if (true_val != (AcReal)0.0){
              if (max_abs_relative_difference<(abs_diff/true_val)){
                max_abs_relative_difference=(abs_diff/true_val);
                gpu_val_for_largest_diff = out_val;
                true_val_for_largest_diff = true_val;
		largest_diff_point = (int3){(int)i,(int)j,(int)k};
              }
            }  
          }
          if (isnan(out_val))
          {
            acLogFromRootProc(rank,"nan before at %d,%d,%d,%s!\n!",i,j,k,acGetFieldName(Field(ivar)));
            acLogFromRootProc(rank,"%.7e\n",(double)out_val);
          }
        }
      }
    }
  }


  passed &= !has_nans(mesh);
  if (!passed)
  {
  	for (int ivar=0;ivar<mvar;ivar++)
        {
    		acLogFromRootProc(0,"ratio of values wrong for field: %s\t %f\n",acGetFieldName(Field(ivar)),(double)num_of_points_where_different[ivar]/num_of_points);
		for (int x = -1; x <= 1; ++x)
		{
			for (int y = -1; y <= 1; ++y)
			{
				for (int z = -1; z <= 1; ++z)
				{
	  				const int index = (x+1) + 3*((y+1) + 3*(z+1));
					if (different_in[ivar][index]) acLogFromRootProc(0,"different (%d,%d,%d)\n",x,y,z);
				}
			}
		}
	}
  	acLogFromRootProc(0,"max abs not passed val: %.7e\t%.7e\n",(double)max_abs_not_passed_val, (double)fabs(true_pair));
  	acLogFromRootProc(0,"max abs relative difference val: %.7e\n",(double)max_abs_relative_difference);
  	acLogFromRootProc(0,"Point where biggest rel diff: %d,%d,%d\n",largest_diff_point.x,largest_diff_point.y,largest_diff_point.z);
  	acLogFromRootProc(0,"largest difference: %.7e\t%.7e\n",(double)gpu_val_for_largest_diff, (double)true_val_for_largest_diff);
  	acLogFromRootProc(0,"abs range: %.7e-%7e\n",(double)min_abs_value,(double)max_abs_value);
    	acLogFromRootProc(0,"Did not pass BC test :(\n");
	fprintf(stderr,"Did not pass BC\n");

    	exit(EXIT_FAILURE);
  }
  fflush(stdout);

  acHostMeshDestroy(&mesh_to_copy);
  acHostMeshCopyVertexBuffers(tmp_mesh_to_store,mesh);
  acHostMeshDestroy(&tmp_mesh_to_store);
  acGridDestroyTaskGraph(bcs);
  acHostMeshDestroy(&tmp_mesh_to_store);
}
/***********************************************************************************************/
extern "C" void gpuSetDt(double t)
{
	acGridSynchronizeStream(STREAM_ALL);
 	acDeviceSetInput(acGridGetDevice(), AC_t,AcReal(t));
	beforeBoundaryGPU(false,0,t);
	if (!lcourant_dt)
	{
		fprintf(stderr,"gpuSetDt works only for Courant timestep!!\n");
		exit(EXIT_FAILURE);
	}
	//TP: not needed but for extra safety
  	acDeviceSetInput(acGridGetDevice(), AC_step_num, (PC_SUB_STEP_NUMBER) 0);
	const auto graph = acGetOptimizedDSLTaskGraph(AC_calculate_timestep);

	acGridExecuteTaskGraph(graph,1);
	acGridSynchronizeStream(STREAM_ALL);
	AcReal dt1_ = calc_dt1_courant(AcReal(t));
	set_dt(dt1_);
	dt1_interface = dt1_;
        acDeviceSwapBuffers(acGridGetDevice());
	loadFarray();
	//TP: not strictly needed but for extra safety
}
/***********************************************************************************************/
extern "C" void getGPUReducedVars(AcReal* dst)
{
        acDeviceSynchronizeStream(acGridGetDevice(),STREAM_DEFAULT);
#if LAXIONSU2BACK
	dst[0] = acDeviceGetOutput(acGridGetDevice(), AC_grand_sum);
	dst[1] = acDeviceGetOutput(acGridGetDevice(), AC_dgrant_sum);
	dst[2] = acDeviceGetOutput(acGridGetDevice(), AC_trdoteff2km_sum);
	dst[3] = acDeviceGetOutput(acGridGetDevice(), AC_trdoteff2m_sum);
	dst[4] = acDeviceGetOutput(acGridGetDevice(), AC_treff2km_sum);
	dst[5] = acDeviceGetOutput(acGridGetDevice(), AC_treff2m_sum);
	dst[6] = acDeviceGetOutput(acGridGetDevice(), AC_tldoteff2km_sum);
	dst[7] = acDeviceGetOutput(acGridGetDevice(), AC_tldoteff2m_sum);
	dst[8] = acDeviceGetOutput(acGridGetDevice(), AC_tleff2km_sum);
	dst[9] = acDeviceGetOutput(acGridGetDevice(), AC_tleff2m_sum);
#endif
#if LBACKREACT_INFL
	dst[0] = acDeviceGetOutput(acGridGetDevice(), AC_a2rhom_all__mod__backreact_infl);
	dst[1] = acDeviceGetOutput(acGridGetDevice(), AC_a2rhopm_all__mod__backreact_infl);
	dst[2] = acDeviceGetOutput(acGridGetDevice(), AC_a2rhophim_all__mod__backreact_infl);
	dst[3] = acDeviceGetOutput(acGridGetDevice(), AC_a2rhogphim_all__mod__backreact_infl);
	dst[4] = acDeviceGetOutput(acGridGetDevice(), AC_sige1m_all_nonaver__mod__backreact_infl);
	dst[5] = acDeviceGetOutput(acGridGetDevice(), AC_sigb1m_all_nonaver__mod__backreact_infl);
	dst[6] = acDeviceGetOutput(acGridGetDevice(), AC_ddotam_all__mod__backreact_infl);
	dst[7] = acDeviceGetOutput(acGridGetDevice(), AC_e2m_all__mod__backreact_infl);
	dst[8] = acDeviceGetOutput(acGridGetDevice(), AC_b2m_all__mod__backreact_infl);
#endif
#if LKLEIN_GORDON
	dst[0] = acDeviceGetOutput(acGridGetDevice(), AC_a2rhom_all__mod__klein_gordon);
	dst[1] = acDeviceGetOutput(acGridGetDevice(), AC_a2rhopm_all__mod__klein_gordon);
	dst[2] = acDeviceGetOutput(acGridGetDevice(), AC_a2rhophim_all__mod__klein_gordon);
	dst[3] = acDeviceGetOutput(acGridGetDevice(), AC_a2rhogphim_all__mod__klein_gordon);
	dst[4] = acDeviceGetOutput(acGridGetDevice(), AC_sige1m_all_nonaver__mod__klein_gordon);
	dst[5] = acDeviceGetOutput(acGridGetDevice(), AC_sigb1m_all_nonaver__mod__klein_gordon);
	dst[6] = acDeviceGetOutput(acGridGetDevice(), AC_ddotam_all__mod__klein_gordon);
	dst[7] = acDeviceGetOutput(acGridGetDevice(), AC_e2m_all__mod__klein_gordon);
	dst[8] = acDeviceGetOutput(acGridGetDevice(), AC_b2m_all__mod__klein_gordon);
#endif
}
/***********************************************************************************************/
